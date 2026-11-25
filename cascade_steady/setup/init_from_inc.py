#!/usr/bin/env python3
"""
Initialise the compressible case from the incompressible run (the v1712 workflow).

    python3 init_from_inc.py INC_CASE COMP_CASE

Copies, from the latest time of INC_CASE into COMP_CASE/0:
  - the internal fields of U, p, nuTilda and nut
  - the boundary values of U and p on the inlet and outlet (row_1_INFLOW, row_1_OUTFLOW)
keeping the compressible boundary-condition types. The incompressible p carries the compressible
pressure numbers (rho = 1 convention, as in v1712), so it is copied as it is.

Why the boundary values: v1712 copied the whole incompressible U and p files, so the compressible
run started with inlet values consistent with the interior. Copying only the internal field left
the template's inlet p of 98 384 Pa against about 124 000 Pa in the first cells. That jump drove the
first momentum and pressure solves (first p residual 0.33 against 0.038 in v1712) and the
pressure-inlet runs settled on a non-physical state with 1.6 times the mass flow (2 Nov 2026).
"""
import os, re, sys

inc, comp = (os.path.abspath(x) for x in sys.argv[1:3])
times = [d for d in os.listdir(inc) if re.fullmatch(r"[0-9.eE+-]+", d) and float(d) > 0]
if not times:
    sys.exit(f"no solution times in {inc}")
latest = max(times, key=float)


def entry_end(text, i):
    """index just past the ';' that ends the entry starting at i (bracket depth tracked)"""
    depth, j = 0, i
    while j < len(text):
        c = text[j]
        if c == "(":
            depth += 1
        elif c == ")":
            depth -= 1
        elif c == ";" and depth == 0:
            return j + 1
        j += 1
    raise ValueError("unterminated entry")


def internal_field(text):
    i = text.index("internalField")
    return i, entry_end(text, i)


def patch_block(text, patch):
    """(start, end) of the {...} dictionary of a named patch inside boundaryField, or None"""
    b = text.index("boundaryField")
    m = re.search(r"\n[ \t]*" + re.escape(patch) + r"[ \t]*\n[ \t]*\{", text[b:])
    if not m:
        return None
    i = b + m.end() - 1
    depth, j = 0, i
    while j < len(text):
        if text[j] == "{":
            depth += 1
        elif text[j] == "}":
            depth -= 1
            if depth == 0:
                return i, j + 1
        j += 1
    raise ValueError(f"unbalanced braces in patch {patch}")


def value_entry(text, block):
    """(start, end) of the 'value ...;' entry inside a patch block, or None"""
    s, e = block
    m = re.search(r"\n[ \t]*value[ \t\n]", text[s:e])
    if not m:
        return None
    i = s + m.start() + 1
    i += len(text[i:]) - len(text[i:].lstrip())              # to the keyword itself
    return i, entry_end(text, i)


for f in ("U", "p", "nuTilda", "nut"):
    src = open(os.path.join(inc, latest, f)).read()
    dst_path = os.path.join(comp, "0", f)
    dst = open(dst_path).read()
    s0, s1 = internal_field(src)
    d0, d1 = internal_field(dst)
    dst = dst[:d0] + src[s0:s1] + dst[d1:]
    copied = ["internal field"]
    if f in ("U", "p"):
        for patch in ("row_1_INFLOW", "row_1_OUTFLOW"):
            sb, db = patch_block(src, patch), patch_block(dst, patch)
            if not (sb and db):
                continue
            sv, dv = value_entry(src, sb), value_entry(dst, db)
            if not (sv and dv):
                continue
            dst = dst[:dv[0]] + src[sv[0]:sv[1]] + dst[dv[1]:]
            copied.append(f"{patch} value")
    if dst.count("{") != dst.count("}"):
        sys.exit(f"{f}: unbalanced braces after copying - not written")
    open(dst_path, "w").write(dst)
    print(f"{f}: {', '.join(copied)} copied from {os.path.basename(inc)}/{latest}")
