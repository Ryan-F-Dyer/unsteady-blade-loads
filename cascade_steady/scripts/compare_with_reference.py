#!/usr/bin/env python3
"""
Compare the result of this run (post/summary.json) with the run reported in the dissertation
(reference/as_run_summary.json).

    python3 scripts/compare_with_reference.py
"""
import json, os, sys

here = os.path.dirname(os.path.abspath(__file__))
new_f = os.path.join(here, "..", "post", "summary.json")
ref_f = os.path.join(here, "..", "reference", "as_run_summary.json")
if not os.path.isfile(new_f):
    sys.exit("no post/summary.json: run scripts/cascade_post.py first")
new, ref = json.load(open(new_f)), json.load(open(ref_f))

rows = [("blade force per unit span, pitchwise (N/m)", "force_per_span_N_m", 0),
        ("blade force per unit span, axial (N/m)", "force_per_span_N_m", 1),
        ("mass flow per unit span (kg/s/m)", "mass_flow_per_span_kg_s_m", None),
        ("inlet total pressure, mass average (Pa)", "p0_in_mass", None),
        ("inlet static pressure (Pa)", "inlet_p", None),
        ("total-pressure loss coefficient", "loss_coefficient", None),
        ("exit isentropic Mach number", "exit_Mis_mass", None),
        ("peak isentropic Mach number, suction side", "Mis_peak_SS", None),
        ("pressure at half the suction side (Pa)", "p_SS_50", None),
        ("pressure at half the pressure side (Pa)", "p_PS_50", None),
        ("largest y+ on the blade", "yPlus_max", None)]


def get(d, key, i):
    v = d.get(key)
    if v is None:
        return None
    return v[i] if i is not None else v


print(f"{'quantity':48s} {'this run':>12s} {'dissertation':>12s} {'difference':>11s}")
for label, key, i in rows:
    a, b = get(new, key, i), get(ref, key, i)
    if a is None or b is None:
        print(f"{label:48s} {'-' if a is None else format(a, '12.5g'):>12s} {'-' if b is None else format(b, '12.5g'):>12s}")
        continue
    print(f"{label:48s} {a:12.5g} {b:12.5g} {100.0 * (a - b) / b:+10.3f}%")
print("change of the blade force over the last 500 iterations: "
      f"{new.get('force_change_last_500_pct', float('nan')):.2g}% (dissertation run: {ref.get('force_change_last_500_pct', float('nan')):.2g}%)")
if new.get("warnings"):
    print("warnings of the post-processing:", "; ".join(new["warnings"]))
