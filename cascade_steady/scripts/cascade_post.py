#!/usr/bin/env python3
"""
Per-case post-processing for the Chapter 2 cascade runs.

    python3 cascade_post.py CASE_DIR

Reads the blade surface pressure (surfaces function object, raw format) at every write time,
the blade force history, residuals, inlet state and y+. Compares with
  * NUMECA FINE/Turbo (reference/fineturbo_Tu*.csv, Tu 0.8% and 5.2%): FINE/Turbo node values are
    interpolated along the suction- and pressure-surface polylines to each OpenFOAM face centre;
    deviation = |p_OF - p_FT| / max(p_OF, p_FT) x 100 (the Chapter 2 metric);
  * the v1712 OpenFOAM results (reference/v1712_*.csv), face by face (version check);
  * Consigny and Richards' measured isentropic Mach number (reference/consigny_Mis.csv).
Writes CASE_DIR/post/: summary.json, surface_p.csv, surface_p.png, Mis.png, deviation.png, convergence.png
Transient runs (rhoPimpleFoam): the surface pressure is pMean, the forces are their time mean from avgStart
(with the rms and peak-to-peak of the fluctuation and its main frequency), and the loss coefficient comes from
the mean fields; the v1712 comparison then shows the transient mean against the steady v1712 result.
"""
import csv, glob, json, math, os, re, sys
import numpy as np
_trapz = getattr(np, "trapezoid", None) or np.trapz   # numpy 2 renamed trapz
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
REF = os.path.join(HERE, "..", "reference")
CHORD, P0_INLET = 0.08, 127500.0


def read_params(case):
    p = {}
    for k, v in re.findall(r"^\s*(\w+)\s+([^;\s]+)\s*;", open(os.path.join(case, "caseParameters")).read(), re.M):
        try:
            p[k] = float(v)
        except ValueError:
            p[k] = v
    return p


def numdirs(path):
    return sorted([d for d in os.listdir(path) if re.fullmatch(r"[0-9.eE+-]+", d)], key=float) if os.path.isdir(path) else []


def read_raw(path, field="p"):
    """raw surface file: '#' header (last one names the columns), then x y z values. Returns (xyz, values of
    the named field), or None if the file doesn't hold that field"""
    names, rows = None, []
    for line in open(path):
        if line.startswith("#"):
            t = line.lstrip("#").split()
            if t and t[0] == "x":
                names = t
            continue
        v = line.split()
        if v:
            rows.append([float(x) for x in v])
    if not rows:
        return None
    a = np.array(rows)
    if names and field in names:
        return a[:, :3], a[:, names.index(field)]
    if names is None and re.match(re.escape(field) + r"[_.]", os.path.basename(path)):
        return a[:, :3], a[:, -1]                                  # no column names: the file is named after the field
    return None


def surface_by_time(case, field="p"):
    """blade surface values of one field (p, or pMean for transient runs) at every write time that has it"""
    out = {}
    base = os.path.join(case, "postProcessing", "bladeSurface")
    for t in numdirs(base):
        for f in sorted(glob.glob(os.path.join(base, t, "**", "*.raw"), recursive=True)):
            r = read_raw(f, field)
            if r is not None:
                out[t] = r
                break
    return out


def read_dat(pattern):
    files = sorted(glob.glob(pattern))
    if not files:
        return None, None
    names, rows = None, []
    for line in open(files[-1]):
        if line.startswith("#"):
            t = line.lstrip("#").replace("\t", " ").split()
            if t and t[0] == "Time":
                names = t
            continue
        v = line.replace("(", " ").replace(")", " ").split()
        if v:
            rows.append(v)
    return names, rows


def profile_polylines():
    """blade profile from reference/blade_profile.csv (side, y, z; each side from the leading edge): the geometry for
    s/c when no reference solution is present. The third column (reference pressure) is empty."""
    f = os.path.join(REF, "blade_profile.csv")
    d = {"PS": [], "SS": []}
    for r in csv.DictReader(l for l in open(f) if not l.startswith("#")):
        d[r["side"]].append((float(r["y"]), float(r["z"]), float("nan")))
    return {k: np.array(v) for k, v in d.items()}


def polylines(tu):
    f = os.path.join(REF, f"fineturbo_Tu{tu}.csv")
    if not os.path.isfile(f):
        return None
    d = {"PS": [], "SS": []}
    for r in csv.DictReader(open(f)):
        d[r["side"]].append((float(r["y"]), float(r["z"]), float(r["p"])))
    for k in d:
        d[k] = np.array(d[k])
    # both polylines start at the leading edge: PS starts there; SS is stored trailing edge first
    le = d["PS"][0, :2]
    if np.hypot(*(d["SS"][0, :2] - le)) > np.hypot(*(d["SS"][-1, :2] - le)):
        d["SS"] = d["SS"][::-1]
    return d


def project(pt, line):
    """nearest point on a polyline: (distance, arc length from start, interpolated p)"""
    P = line[:, :2]
    seg = P[1:] - P[:-1]
    L = np.hypot(seg[:, 0], seg[:, 1])
    s0 = np.concatenate([[0], np.cumsum(L)])
    t = np.clip(((pt - P[:-1]) * seg).sum(1) / np.maximum(L ** 2, 1e-30), 0, 1)
    q = P[:-1] + seg * t[:, None]
    dist = np.hypot(*(q - pt).T)
    i = int(np.argmin(dist))
    return dist[i], s0[i] + t[i] * L[i], line[i, 2] + t[i] * (line[i + 1, 2] - line[i, 2])


def mach_is(p, p0=P0_INLET):
    return np.sqrt(np.maximum(5.0 * ((p0 / p) ** (0.4 / 1.4) - 1.0), 0.0))


def span_from_mesh(case):
    """thickness of the 2D slab (ASCII, binary or .gz points): along the axis of smallest extent, the
    distance between the median positions of the front and back points, so stray points don't count"""
    import gzip
    f = os.path.join(case, "constant", "polyMesh", "points")
    raw = open(f, "rb").read() if os.path.isfile(f) else gzip.open(f + ".gz", "rb").read()
    head = raw[:4096].decode("latin-1")
    m = re.search(rb"\n(\d+)\s*\n\(", raw)
    n, start = int(m.group(1)), m.end()
    if re.search(r"format\s+binary", head):
        a = np.frombuffer(raw[start:start + 24 * n], dtype="<f8").reshape(n, 3)
    else:
        a = np.array(re.findall(r"\(([-\d.eE+]+) ([-\d.eE+]+) ([-\d.eE+]+)\)", raw[start:].decode("latin-1"))[:n], float)
    k = int(np.argmin(a.max(0) - a.min(0)))
    x = a[:, k]
    mid = 0.5 * (x.min() + x.max())
    return float(np.median(x[x >= mid]) - np.median(x[x < mid]))


def total_pressure(p, U, T, R=8314.47 / 28.9, g=1.4):
    M2 = U * U / (g * R * T)
    return p * (1.0 + 0.5 * (g - 1.0) * M2) ** (g / (g - 1.0))


def checkmesh_stats(case):
    f = os.path.join(case, "log.checkMesh")
    if not os.path.isfile(f):
        return {}
    s = open(f).read()
    out = {}
    for key, pat in (("mesh_cells", r"cells:\s+(\d+)"), ("mesh_max_nonorth", r"non-orthogonality Max:\s*([\d.eE+-]+)"),
                     ("mesh_mean_nonorth", r"non-orthogonality Max:\s*[\d.eE+-]+\s+average:\s*([\d.eE+-]+)"),
                     ("mesh_max_skewness", r"Max skewness = ([\d.eE+-]+)"), ("mesh_max_aspect_ratio", r"Max aspect ratio = ([\d.eE+-]+)"),
                     ("mesh_failed_checks", r"Failed (\d+) mesh checks")):
        m = re.search(pat, s)
        if m:
            out[key] = float(m.group(1))
    out["mesh_ok"] = "Mesh OK" in s
    m = re.search(r"\*\*\*(.*)", s)
    if m:
        out["mesh_check_note"] = m.group(1).strip()             # AutoGrid5: edges not aligned with the 2D direction
    return out


def v1712_column(P):
    tu = {0.8: "A1", 3.0: "A2", 3.8: "A3", 5.2: "A4"}.get(round(P["Tu"], 1))
    if abs(P["beta"] - 30) > 1e-6 and P["solver"] == "rhoSimpleFoam" and P["inletType"] == "PIBC":
        return "v1712_angle_surface_p.csv", f"beta{int(round(P['beta']))}"
    if tu is None or abs(P["beta"] - 30) > 1e-6:
        return None, None
    kind = "inc" if P["solver"] == "simpleFoam" else P["inletType"]
    if kind == "VIBC" and tu != "A1":
        return None, None
    return "v1712_surface_p.csv", f"{tu}_{kind}"


def main():
    case = os.path.abspath(sys.argv[1])
    out = os.path.join(case, "post"); os.makedirs(out, exist_ok=True)
    P = read_params(case)
    mesh = str(P.get("meshType", "ag5"))
    res = dict(case=os.path.basename(case), solver=P["solver"], inlet=P["inletType"], Tu=P["Tu"], beta=P["beta"],
               mesh=mesh, Uin=P.get("Uin"), numerics=str(P.get("numerics", "v1712")), wall=str(P.get("wall", "lowRe")))
    res.update(checkmesh_stats(case))
    warn = []
    transient = P["solver"] == "rhoPimpleFoam"                    # time-accurate: time means from avgStart

    surf = surface_by_time(case, "pMean" if transient else "p")
    if not surf:
        sys.exit("no blade surface output found" + (" for pMean (did the run pass avgStart?)" if transient else ""))
    t_last = list(surf)[-1]
    xyz, p = surf[t_last]
    yz = xyz[:, 1:3]
    res["surface_times"] = list(surf)

    if mesh == "ag5":
        # reference face order (v1712 face centres): match each output face to it
        ref = list(csv.DictReader(open(os.path.join(REF, "v1712_surface_p.csv"))))
        ref_yz = np.array([[float(r["y"]), float(r["z"])] for r in ref])
        order = np.array([int(np.argmin(np.hypot(*(ref_yz - c).T))) for c in yz])
        if len(set(order)) != len(order):
            warn.append("face matching to the v1712 face order is not one-to-one")
        inv = np.argsort(order)                                   # output rows in v1712 face order
        p_by_time = {t: v[1][inv] for t, v in surf.items()}
        yz = yz[inv]
    else:                                                         # another mesh: no face-by-face version check
        p_by_time = {t: v[1] for t, v in surf.items()}
    p = p_by_time[t_last]

    # FINE/Turbo comparison (Tu 0.8 and 5.2 only)
    design = abs(P["beta"] - 30.0) < 1e-6                         # FINE/Turbo data exist for the design angle only
    lines = polylines({0.8: "0.8", 5.2: "5.2"}.get(round(P["Tu"], 1), "none")) if design else None
    side, s_c, p_ft = [], [], []
    geo = lines or polylines("0.8") or profile_polylines()        # geometry for s/c even without FT data
    for c in yz:
        dps, sps, pps = project(c, geo["PS"])
        dss, sss, pss = project(c, geo["SS"])
        if dps <= dss:
            side.append("PS"); s_c.append(sps / CHORD); p_ft.append(pps)
        else:
            side.append("SS"); s_c.append(sss / CHORD); p_ft.append(pss)
    side, s_c, p_ft = np.array(side), np.array(s_c), np.array(p_ft)
    res["faces"] = int(len(p)); res["faces_SS"] = int((side == "SS").sum()); res["faces_PS"] = int((side == "PS").sum())
    dev = None
    if lines:
        dev = np.abs(p - p_ft) / np.maximum(p, p_ft) * 100.0
        res.update(dev_mean=float(dev.mean()), dev_max_SS=float(dev[side == "SS"].max()), dev_max_PS=float(dev[side == "PS"].max()),
                   n_over_1pct=int((dev > 1.0).sum()),
                   dev_max_PS_at_s_c=float(s_c[side == "PS"][np.argmax(dev[side == "PS"])]))
        ps_len = s_c[side == "PS"].max()
        te = (side == "PS") & (s_c > 0.95 * ps_len)
        res["dev_max_excluding_PS_trailing_edge"] = float(dev[~te].max())

    # stations and peak isentropic Mach number (mesh study and mesh comparison)
    Mis = mach_is(p)
    for sd in ("SS", "PS"):
        m = side == sd
        o = np.argsort(s_c[m]); ss, pp = s_c[m][o], p[m][o]
        for f in (0.25, 0.5, 0.75):
            res[f"p_{sd}_{int(f * 100)}"] = float(np.interp(f * ss.max(), ss, pp))
    m = side == "SS"
    i = int(np.argmax(Mis[m]))
    res["Mis_peak_SS"], res["Mis_peak_SS_at_s_c"] = float(Mis[m][i]), float(s_c[m][i])

    # version check against v1712 (face by face), at 1000 iterations and at the last write
    fname, col = v1712_column(P) if mesh == "ag5" else (None, None)
    if fname:
        ref_v = np.array([float(r[col]) for r in csv.DictReader(open(os.path.join(REF, fname)))])
        for t in (k for k in p_by_time if abs(float(k) - 1000) < 1e-6 or k == t_last):
            d = np.abs(p_by_time[t] - ref_v) / np.maximum(p_by_time[t], ref_v) * 100.0
            res[f"v1712_dev_mean_at_{t}"] = float(d.mean()); res[f"v1712_dev_max_at_{t}"] = float(d.max())
        res["v1712_column"] = col

    # convergence: blade force change over the last 500 iterations; final residuals
    names, rows = read_dat(os.path.join(case, "postProcessing", "bladeForce", "*", "force.dat"))
    if rows:
        F = np.array([[float(v) for v in r[:4]] for r in rows])
        Fm = np.hypot(F[:, 2], F[:, 3])
        res["blade_force_N"] = [float(F[-1, 2]), float(F[-1, 3])]
        if transient:
            # time mean over the averaging window, the size of the fluctuation about it, and its main
            # frequency (trailing-edge shedding, if the flow sheds)
            m = F[:, 0] >= float(P["avgStart"]) - 1e-12
            if m.sum() < 20:
                warn.append("fewer than 20 force samples after avgStart: the means use the whole run")
                m = np.ones(len(F), bool)
            # the time step can vary (Courant-limited), so means are time-weighted and the spectrum is taken
            # on a uniform grid at the median step
            t = F[m, 0]
            tmean = lambda y: float(_trapz(y, t) / (t[-1] - t[0])) if t[-1] > t[0] else float(y.mean())
            Fy_avg, Fz_avg = tmean(F[m, 2]), tmean(F[m, 3])
            Fbar = math.hypot(Fy_avg, Fz_avg)
            res["force_change_last_500_pct"] = None
            res["averaging_window_s"] = [float(t[0]), float(t[-1])]
            res["time_step_s"] = [float(np.diff(F[:, 0]).min()), float(np.median(np.diff(t))), float(np.diff(t).max())]
            res["force_rms_pct"] = float(math.sqrt(max(tmean((Fm[m] - tmean(Fm[m])) ** 2), 0.0)) / Fbar * 100)
            res["force_peak_to_peak_pct"] = float((Fm[m].max() - Fm[m].min()) / Fbar * 100)
            dt = float(np.median(np.diff(t)))
            tu = np.arange(t[0], t[-1], dt)
            fy = np.interp(tu, t, F[m, 2]); fy -= fy.mean()
            spec = np.abs(np.fft.rfft(fy * np.hanning(len(fy))))
            freq = np.fft.rfftfreq(len(fy), dt)
            k = 1 + int(np.argmax(spec[1:]))
            res["force_peak_freq_Hz"] = float(freq[k])
        else:
            n = min(500, len(Fm) - 1)
            res["force_change_last_500_pct"] = float((Fm[-n:].max() - Fm[-n:].min()) / Fm[-1] * 100) if n > 0 else None
            # a steady run that settles into a small oscillation (e.g. trailing-edge shedding on a coarse mesh)
            # is reported by its mean over the last 500 iterations; for a converged run this equals the last value
            Fy_avg, Fz_avg = (float(F[-n:, 2].mean()), float(F[-n:, 3].mean())) if n > 0 else (float(F[-1, 2]), float(F[-1, 3]))
        try:
            span = span_from_mesh(case)
            res["span_m"] = span
            res["force_per_span_N_m"] = [Fy_avg / span, Fz_avg / span]
            res["force_per_span_last_N_m"] = [float(F[-1, 2]) / span, float(F[-1, 3]) / span]
        except Exception as e:                                    # points unreadable: forces stay per slice
            warn.append(f"span not read from the mesh ({e}); forces are per slice only")
        if res["force_change_last_500_pct"] is not None and res["force_change_last_500_pct"] > 0.01:
            warn.append(f"blade force still changing: {res['force_change_last_500_pct']:.3f}% over the last 500 iterations")
        plt.figure(figsize=(6, 3)); plt.plot(F[:, 0], Fm, "k-", lw=0.8)
        if transient:
            plt.axvline(float(P["avgStart"]), color="r", lw=0.8, ls="--")
        plt.xlabel("time (s)" if transient else "iteration"); plt.ylabel("|F| on the blade"); plt.title(res["case"], fontsize=9)
        plt.tight_layout(); plt.savefig(os.path.join(out, "convergence.png"), dpi=150); plt.close()
    names, rows = read_dat(os.path.join(case, "postProcessing", "residuals", "*", "solverInfo.dat"))
    if rows and names:
        last = rows[-1]
        for i, nme in enumerate(names):
            if nme.endswith("_initial") and i < len(last):
                try:
                    res[f"residual_{nme}"] = float(last[i])
                except ValueError:
                    pass

    # inlet state (area averages on row_1_INFLOW)
    names, rows = read_dat(os.path.join(case, "postProcessing", "inletState", "*", "surfaceFieldValue.dat"))
    if rows:
        v = [float(x) for x in rows[-1]]
        res["inlet_U"] = v[1:4]
        res["inlet_Umag"] = float(np.linalg.norm(v[1:4]))
        res["inlet_angle_deg"] = float(math.degrees(math.atan2(v[2], v[3])))
        # guard: the compressible pressure inlet should give the v1712 inlet speed at the same angle
        # (84.0 m/s at the design 30 deg); more than 5% away means the inlet has broken down
        if P["solver"] == "rhoSimpleFoam" and P["inletType"] == "PIBC" and os.path.isfile(os.path.join(REF, "v1712_angle_inlet_velocity.csv")):
            ref_u = {float(r["beta"]): float(r["Umag"])
                     for r in csv.DictReader(open(os.path.join(REF, "v1712_angle_inlet_velocity.csv")))}
            u0 = ref_u.get(float(P["beta"]))
            if u0 and abs(res["inlet_Umag"] - u0) > 0.05 * u0:
                warn.append(f"pressure-inlet velocity {res['inlet_Umag']:.1f} m/s against {u0:.1f} m/s in v1712 at "
                            f"beta {float(P['beta']):g} deg: the solution is not physical, check the inlet")
        if len(v) > 4:
            res["inlet_p"] = v[4]
        if len(v) > 5:
            res["inlet_T"] = v[5]
    # total-pressure loss from the mass-weighted inlet and outlet states (compressible runs):
    # omega = (p0_in - p0_out) / (p0_in - p_out), p0 isentropic from the mass-averaged p, |U|, T
    # (transient runs: the same from the time-mean fields UMean, pMean, TMean)
    st = {}
    pair = ("inletMean", "outletMean") if transient else ("inletMass", "outletMass")
    for nm in pair:
        names, rows = read_dat(os.path.join(case, "postProcessing", nm, "*", "surfaceFieldValue.dat"))
        if rows and len(rows[-1]) >= 6:
            v = [float(x) for x in rows[-1]]
            st[nm] = (np.linalg.norm(v[1:4]), v[4], v[5])
    if len(st) == 2:
        (Ui, pin, Ti), (Uo, pout, To) = st[pair[0]], st[pair[1]]
        p0i, p0o = total_pressure(pin, Ui, Ti), total_pressure(pout, Uo, To)
        res.update(p0_in_mass=float(p0i), p0_out_mass=float(p0o), p_out_mass=float(pout),
                   loss_coefficient=float((p0i - p0o) / (p0i - pout)),
                   exit_Mis_mass=float(mach_is(pout, p0i)))
    names, rows = read_dat(os.path.join(case, "postProcessing", "massFlow", "*", "surfaceFieldValue.dat"))
    if rows and "span_m" in res:
        res["mass_flow_per_span_kg_s_m"] = abs(float(rows[-1][1])) / res["span_m"]

    names, rows = read_dat(os.path.join(case, "postProcessing", "yPlus", "*", "yPlus.dat"))
    if rows:
        r = rows[-1]
        res["yPlus_min"], res["yPlus_max"], res["yPlus_mean"] = float(r[2]), float(r[3]), float(r[4])

    # tables and plots
    with open(os.path.join(out, "surface_p.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["face", "y", "z", "side", "s_over_c", "p", "p_fineturbo", "deviation_pct", "M_is"])
        for i in range(len(p)):
            w.writerow([i, f"{yz[i,0]:.8g}", f"{yz[i,1]:.8g}", side[i], f"{s_c[i]:.5f}", f"{p[i]:.6g}",
                        f"{p_ft[i]:.6g}" if lines else "", f"{dev[i]:.4f}" if dev is not None else "", f"{Mis[i]:.5f}"])
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    for sd, mk in (("SS", "-"), ("PS", "--")):
        m = side == sd; o = np.argsort(s_c[m])
        ax[0].plot(s_c[m][o], p[m][o] / 1e3, "k" + mk, lw=1, label=f"OpenFOAM v2606 {sd}")
        ax[1].plot(s_c[m][o], Mis[m][o], "k" + mk, lw=1, label=f"OpenFOAM v2606 {sd}")
        if lines:
            s_ft = np.concatenate([[0], np.cumsum(np.hypot(*np.diff(lines[sd][:, :2], axis=0).T))]) / CHORD
            ax[0].plot(s_ft, lines[sd][:, 2] / 1e3, "r" + mk, lw=0.8, label=f"FINE/Turbo {sd}")
            ax[1].plot(s_ft, mach_is(lines[sd][:, 2]), "r" + mk, lw=0.8, label=f"FINE/Turbo {sd}")
    exp = list(csv.DictReader(open(os.path.join(REF, "consigny_Mis.csv"))))
    for sd, mk in (("SS", "o"), ("PS", "s")):
        e = [(float(r["s_over_c"]), float(r["M_is"])) for r in exp if r["side"] == sd]
        ax[1].plot(*zip(*e), mk, mfc="none", color="b", label=f"Consigny and Richards {sd}")
    ax[0].set_xlabel("s/c"); ax[0].set_ylabel("p (kPa)"); ax[1].set_xlabel("s/c"); ax[1].set_ylabel("M_is")
    for a in ax:
        a.legend(fontsize=7, frameon=False)
    fig.suptitle(f"{res['case']}: {P['solver']} {P['inletType']}, Tu {P['Tu']}%, beta {P['beta']:g} deg", fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(out, "surface_p.png"), dpi=200); plt.close(fig)
    if dev is not None:
        plt.figure(figsize=(6, 3.2))
        for sd, mk in (("SS", "o"), ("PS", "s")):
            m = side == sd
            plt.plot(s_c[m], dev[m], mk, ms=3, mfc="none", label=sd)
        plt.axhline(1.0, color="r", lw=0.8, ls="--")
        plt.xlabel("s/c"); plt.ylabel("deviation from FINE/Turbo (%)"); plt.legend(frameon=False)
        plt.tight_layout(); plt.savefig(os.path.join(out, "deviation.png"), dpi=200); plt.close()

    res["warnings"] = warn
    json.dump(res, open(os.path.join(out, "summary.json"), "w"), indent=1)
    for k, v in res.items():
        if k not in ("warnings", "surface_times"):
            print(f"{k:38s} {v:.5g}" if isinstance(v, float) else f"{k:38s} {v}")
    for w_ in warn:
        print("WARNING:", w_)


if __name__ == "__main__":
    main()
