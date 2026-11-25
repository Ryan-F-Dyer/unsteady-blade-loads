# cascade_steady

Steady compressible flow through the turbine cascade of Consigny and Richards (1982), on the medium level of the
snappyHexMesh study (Chapter 2 of the dissertation). This is the mesh and the solution on which the wake method is built.

## Set-up

| | |
|---|---|
| Solver | rhoSimpleFoam (SIMPLEC), started from an incompressible simpleFoam solution (`init/`) |
| Turbulence model | Spalart-Allmaras, wall resolved (y+ below 1) |
| Mesh | `shm5u_L2`: 48 705 cells, one cell thick (span 1 mm); 24 wall layers from 1.14 µm at a ratio of 1.24 |
| Inlet | Velocity 83.9544 m/s at 30° to the axial direction, total temperature 419 K, turbulence intensity 0.8% |
| Outlet | Static pressure 98 384 Pa |
| Blade | Isothermal wall at 296 K |
| Gas | Air, perfect gas |

The inlet velocity is the one that gives the measured inlet total pressure of 127.5 kPa on this mesh. Chord 80 mm,
pitch 53.6 mm. Axes: x span, y pitchwise, z axial.

## Run

    ./Allrun

About ten minutes on one core. `checkMesh` reports three highly skew faces on the periodic boundaries and
"Failed 1 mesh checks" for that reason; this is the mesh as used.

`Allrun` ends with a table of this run against the run reported in the dissertation.

## Results of the run reported in the dissertation

| Quantity | Value |
|---|---|
| Blade force per unit span, pitchwise / axial | 1 075.3 N/m / 1 310.8 N/m |
| Mass flow per unit span | 4.036 kg/s per metre |
| Inlet total pressure (mass average) | 127.50 kPa |
| Total-pressure loss coefficient | 0.0443 |
| Exit isentropic Mach number | 0.620 |

That run continued a converged solution at 84.0 m/s for 2 000 iterations. The case here runs 3 000 iterations from
the incompressible start. The case files are otherwise identical to those of that run.

## Files

| Path | Content |
|---|---|
| `0/`, `constant/`, `system/`, `caseParameters` | The compressible case. The dictionaries read their values from `caseParameters` |
| `init/` | The incompressible start case |
| `setup/setup_case.py`, `setup/template_rho`, `setup/template_inc` | The script and templates that wrote the case files. `setup.json` records its options |
| `setup/init_from_inc.py` | Copies the incompressible solution into `0/`, keeping the compressible boundary conditions |
| `setup/mesh/shm5u_L2/polyMesh.tar.gz` | The mesh |
| `scripts/cascade_post.py` | Blade force, surface pressure, isentropic Mach number, loss coefficient: `post/summary.json`, `post/surface_p.csv`, plots |
| `scripts/compare_with_reference.py` | This run against `reference/as_run_summary.json` |
| `reference/as_run_summary.json`, `reference/as_run_surface_p.csv` | The run reported in the dissertation |
| `reference/consigny_Mis.csv` | Isentropic Mach number measured by Consigny and Richards (1982) |
| `reference/blade_profile.csv` | The blade profile as it is in this mesh, for the surface distance s/c |

To write the case files again:

    cd setup
    python3 setup_case.py ../new_case --solver rhoSimpleFoam --inlet VIBC --Tu 0.8 --beta 30 --iterations 3000 \
        --init-iterations 1000 --numerics hc_u95 --wall lowRe --schemes ho_unlim --U 83.9544 --mesh shm5u_L2

## Reference

Consigny, H. and Richards, B.E. (1982). Short-duration measurements of heat-transfer rate to a gas turbine rotor
blade. *Journal of Engineering for Power*, 104(3).
