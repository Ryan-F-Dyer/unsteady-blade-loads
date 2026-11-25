# unsteady-blade-loads

A numerical method to predict the unsteady aerodynamic loads that the wakes of an upstream blade row impose on a
turbine blade, built on OpenFOAM v2606 and Python.

The repository holds the cases of the method as developed in the MEng dissertation
*A numerical method to predict unsteady aerodynamic loads on axial-flow turbomachinery blades*
(Ryan Dyer, University of Pretoria, 2026; supervisor Prof. P.S. Heyns).

## Cases

| Folder | Case | Dissertation | Status |
|---|---|---|---|
| `cascade_steady` | Steady compressible flow through the turbine cascade of Consigny and Richards (1982) on the selected snappyHexMesh mesh, with rhoSimpleFoam and the Spalart-Allmaras model | Chapter 2 | included |
| `cylinder_transient` | Vortex shedding from a circular cylinder with rhoPimpleFoam, on the selected mesh | Chapter 3 | being added |
| `wake_design_point` | The wake method at the design operating point: a wake prescribed at the inlet, non-reflecting inlet and outlet, force and pressure harmonics for a structural analysis | Chapters 4 and 5 | being added |

Each folder is an OpenFOAM case that runs by itself.

## Requirements

- OpenFOAM v2606 (www.openfoam.com). The cases were written for and run with this version.
- Python 3 with `numpy` and `matplotlib`, for the post-processing.
- Linux, or Windows with WSL.

## Running a case

    source /usr/lib/openfoam/openfoam2606/etc/bashrc     # or wherever OpenFOAM v2606 is installed
    cd cascade_steady
    ./Allrun
    ./Allclean        # removes everything Allrun wrote

The README in each folder gives the set-up, the run time and the results to expect.

## Layout of a case folder

| Path | Content |
|---|---|
| `0/`, `constant/`, `system/` | The OpenFOAM case as it was run |
| `setup/` | The script and templates that wrote the case files, and the mesh, packed |
| `scripts/` | Post-processing |
| `reference/` | Results of the run reported in the dissertation, to compare a new run against, and the published data used for comparison |

## Licence

MIT, see `LICENSE`.

OpenFOAM itself is distributed under the GNU General Public License. This offering is not approved or endorsed by
OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM®
and OpenCFD® trade marks.
