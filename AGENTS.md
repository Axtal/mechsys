# MechSys — instructions for future agents

This file is loaded automatically by DSH for work in this repository. It is
deliberately short and operational. Read
[`MECHSYS_HANDBOOK.md`](MECHSYS_HANDBOOK.md) before making substantive changes;
the handbook is the durable knowledge base for users and agents.

## First checks

1. Run `pwd`, `git status`, and inspect the relevant file before editing.
2. This checkout is `mechsys` at `/home/sergio/mechsys`; `mechsys` is a symlink
   to `lib/`.
3. **Do not assume `HEAD` is the optimized/correct baseline.** The DEM and
   CUDA optimization work is committed up to `0d9d696`; `HEAD` also contains
   the MPM fork merge `74407d4`. The pre-merge tip is branch
   `backup/pre-mario-merge`. Never run `git reset`, `git checkout -- .`, or
   delete generated/untracked files without explicit permission.
4. Read the current optimization reports before changing DEM, CUDA, LBM-DEM or
   MPM performance code:
   - [`CPU_DEM_OPTIMIZATION_CHANGES.md`](CPU_DEM_OPTIMIZATION_CHANGES.md)
   - [`CUDA_OPTIMIZATION_CHANGES.md`](CUDA_OPTIMIZATION_CHANGES.md)
   - [`REDUCED_TEST_CU_03_PERFORMANCE.md`](REDUCED_TEST_CU_03_PERFORMANCE.md)
   - [`MPM_OPTIMIZATION_CHANGES.md`](MPM_OPTIMIZATION_CHANGES.md)
5. **GPU capability is not guaranteed in an agent shell.** The `nvidia_*` kernel
   modules can be loaded while the agent sandbox exposes no `/dev/nvidia*`
   device nodes; then `nvidia-smi` reports "couldn't communicate with the
   NVIDIA driver" and CUDA initialization fails even though the host GPU is
   healthy. Run `ls /dev/nvidia*` before planning any CUDA case, and do not
   report a GPU build or run you could not execute. When the nodes are absent,
   restrict yourself to CPU targets (128 cores are available, no scheduler) and
   ask the user to run GPU cases outside the sandbox.

## Memory maintenance (mandatory after every prompt)

This checkout is a long-running engineering notebook as well as a codebase.
At the end of **every user prompt** in this repository, before the final
response, perform this checkpoint:

1. Ask: *did this prompt teach a new durable fact about MechSys?* A durable
   fact is architecture, a verified invariant, a command that works, a measured
   benchmark, a demonstrated failure mode, or a corrected workflow. Transient
   logs, one-off `/tmp` paths, and guesses are not durable facts.
2. If there is new knowledge, update the right file:
   - `AGENTS.md` — short operational rules that every future session must see;
   - `MECHSYS_HANDBOOK.md` — architecture, build/run procedures, exactness
     protocol, traps, and benchmark context;
   - `MEMORY.md` — distilled experience and the short verification checklist;
   - the relevant optimization report (`CPU_DEM_...`, `CUDA_...`,
     `REDUCED_TEST_CU_03_...`) — designs, measurements, and accuracy results;
   - `SKILL.md` and `.agents/skills/mechsys/SKILL.md` — task-triggered
     guidance. Keep these two copies in sync except for their mirror note.
3. Prefer editing an existing section over appending duplicate prose. State
   the exact command, input, thread count, and whether a number is wall time or
   solve time.
4. If nothing durable was learned, make no memory edit. Do not manufacture
   activity or copy transient output into the repository.
5. Run `git diff --check` after memory edits. Do not commit or push unless the
   user explicitly asks.

The final response should briefly mention which memory file changed, or say
that no memory update was needed.

## `tdem/ttt` triaxial runs

`ptype = state` restart is **broken** in the current binary: the saved
interaction list still references the six wall particles, and
`DelParticles(-2..-7)` + `GenBoundingBox` recreate them at the same indices with
stale contact records, giving a NaN at particle index `NP-6`.  Regenerate from
scratch instead of restarting, and never shorten `T0` below `2*Time` (a negative
`tspan` reverses the stress ramp).  Details in
[`MECHSYS_HANDBOOK.md`](MECHSYS_HANDBOOK.md), section 7.

## Build and test

The root CMake/Makefiles are already configured. Build from the test directory:

```bash
cd /home/sergio/mechsys/tdem
make -j8 test_periodic          # CPU Voronoi reference
make -j8 test_cu_01 test_cu_02 test_cu_03 test_cu_04

cd /home/sergio/mechsys/tlbmdem
make -j8 tlbmdem01 tlbmdem02
make -j8 tlbmdem_cu_01 tlbmdem_cu_02 tlbmdem_cu_03 tlbmdem_cu_04

cd /home/sergio/mechsys/tmpm
make -j8 tmpm01 tmpm02 tmpm03 tmpm04 tmpm05
make -j8 tmpm_cu_01 tmpm_cu_02 tmpm_cu_03 tmpm_cu_04 tmpm_cu_05

cd /home/sergio/mechsys/tlbmmpm
make -j8 tlbmmpm01 tlbmmpm02      # also compiles the MPM headers
make -j8 tlbmmpm_cu_01 tlbmmpm_cu_02
```

The `tmpm` cases are large (up to 216,000 particles and 580,000 steps). For a
smoke run, copy the source outside the checkout and reduce the grid and `Tf`;
`tmpm04` must keep its physical `tf = 200 s` to reach steady state. Measured
smoke timings and the full benchmark sizes are in `MECHSYS_HANDBOOK.md`,
section 9.

Typical runs:

```bash
# CPU reduced Voronoi benchmark: FileKey argument; .inp suffix is appended.
cd /home/sergio/mechsys/tdem
./test_periodic /home/sergio/mechsys/benchmarks/cpu_voronoi/tper_small 8

# CUDA test: Nproc, Verlet distance (Alpha), threads per block.
./test_cu_03 1 0.05 128
```

`tlbmdem01`/`tlbmdem02` are large. For a quick integration smoke test, reduce
their grid/Tf in a temporary copy; do not modify the checked-in test just to
make it fast.

## FLBM / GPU fluid cases

- FLBM lives in `lib/flbm/Domain.h`; the CUDA kernels are in
  `lib/flbm/lbm.cuh`. `D2Q9` with `iVec3_t(nx,ny,1)` is the 2-D GPU setup.
- The GPU streaming kernel `cudaFusedStreamMacro` zeroes `BForce` every step.
  A constant body force must be rewritten on the device (`dom.pBForce`) in the
  `Setup` callback before the collision; setting the host array once is not
  enough.
- `dom.Nthread` (default 256) sets the CUDA block size. The `Nproc` argument of
  `Domain::Solve` only affects host-side packing/transfers on the GPU path.
- `tclbm01` (2400x300, 40k steps) is flat for `Nthread=32..256`; 32 vs 256
  was 45.11 vs 45.72 s over three runs each (Welch p=0.11). Do not quote its
  small differences as a tuning gain.
- Velocity convention: in this tree **both** active paths report the raw
  second moment, `u = sum_k f_k c_k / rho` – CPU `StreamSC` and the GPU fused
  stream kernel. The `+0.5*dt*BForce` velocity lives in `StreamEFS_SRT`, which
  `Solve` does not call (verified 2026-11-17: plane Poiseuille, CPU vs GPU
  agree to 1.4e-16). So do not apply a blind "half-force" correction: measure
  it. In a forced channel with bounce-back walls the raw profile equals the
  analytical parabola plus a constant slip, see the handbook 5.1.
- Device population addressing uses `FLBM::FIDX(Ncells,Nneigh,il,k,ic)` from
  `lib/flbm/lbm.cuh`; do not add local copies of the direction-major index
  formula in test kernels.
- The device 3x3 alternative solver is `Mat3::SolAlt` in
  `lib/linalg/matvec.h`; do not add another `lbmmpm_sol_alt`-style copy in a
  CUDA kernel.
- LBMMPM CUDA (`lib/lbmmpm/lbmmpm.cuh`, `tlbmmpm/tlbmmpm_cu_01/02.cu`) is
  D3Q15-only and keeps `MPMDOM.OneStepCPI` on the CPU.  Build it with
  `--fmad=false` (the CMake targets do); without it the facet-edge `Gamma`
  classification differs from the CPU.  MechSys host/CPU code is double.
  CUDA precision is user-selectable with
  `USE_GPU_DOUBLE`; `mechsyscc` and CMake enable it by default.  `Solve`
  imprints once after
  `UpLoadDevice` before the first output; without it `frame 0000` has no
  initial body `Gamma`.  Compare HDF5 datasets, not whole files: `frame 0000`
  payloads are exact but HDF5 timestamps differ, while later frames are not
  bitwise.  The GPU `CD` diagnostic is not bitwise or quantitatively
  validated; `Re`/velocity and t01 deflection are the usable checks.
  `tlbmmpm_cu_01` is the elastic beam; `tlbmmpm_cu_02` is the settled cube
  (the only case with `CD(formula)` in its `.res`).  Its MPM mesh is
  24 tetrahedra/14 corners; body force is `0.25`.  CMAME 2021 reference:
  P. Zhang, S. Sun, Y. Chen, W. Cui, S.A. Galindo-Torres, "Coupled material
  point Lattice Boltzmann method for modeling fluid-structure interactions
  with large deformations," Comput. Methods Appl. Mech. Eng. 385 (2021)
  114040, doi:10.1016/j.cma.2021.114040.  The paper's built-in Eq. (38) is
  `Cd=128*Re^-0.8/5.5`; its cube validation uses a six-tet body and
  `b=2.5e-5`, so stock 24-tet/`b=0.25` is not directly paper-equivalent.
  A manual 6-tet/8-corner
  cube (12 BFac, exact volume) at `b=0.843` reaches
  `Re=100.03`, `CD=0.3026` vs Eq. (38) `0.5844` (-48.2%); at matched Re it
  is worse than the 24-tet body, so the mesh alone is not the fix.
  Velocity inlet (Zou-He, x=0) + zero-gradient outlet (x=Nx-1) was tried:
  at Re=100 six-tet CD=0.3201 (-45.2%); mini sweep Re=5..100 gives
  -26.2..-45.2%, no match, so no full sweep.  The current `nu=0.0020833`
  is ~1000x the CMAME 2021 low-viscosity value `2.0803e-6`;
  `Tau=0.500125` is unstable with the current SRT/coupling.  LES is
  implemented in `cudaCollideMPM` (`Sc=0.17` default, active; Sc=0 vs
  0.17 changes 1k CD 0.2682 vs 0.3117), but does not stabilize the low
  viscosity; all velocity-BC results
  already had LES on.  Original `tlbmmpm_cu_02` stock 24-tet body sweep:
  b=0.25 -> Re 37.99, b=0.80 -> 95.54, b=0.84 -> 98.91; b=0.84 is the
  best Re~100 point (CD=0.3551 vs Eq.38 0.5897).  The facet
  `-2 rho W Cs C` term is now always full; the former empirical
  `ForceCorrection` member and its 0.6 default were removed.  Historical
  FC=0.6 results (CD=0.5858 at Re99.2, 1.2803 at Re38) are no longer
  reproducible without re-adding that scaling.  The CMAME 2021 low-viscosity
  run (`nu=2.0803e-6`,
  `bforce=2.5e-5`)
  is unstable at low Tau: full run `Re=10.553`, `CD=-9170` vs Eq. (38)
  `3.533`; 24-tet also negative; Sc=2 gives 161.5.  Full cube `.res` run:
  `CD` 41.1 -> 0.913, formula 22.25 -> 1.268, final -28.0%.  A longer
  streamwise box helps modestly (at matched `Re=30`, error -25.5% for
  241x61x61 vs -24.1% for 361x61x61 and -23.2% for 481x61x61), but a wider
  transverse box hurts badly (-39.2%/-41.5% for 81/91 at `Re=30`).  Do not
  treat a uniformly larger box as a fix.  The legacy `ImprintLattice`
  differs by one assignment that overwrites the solid-neighbour population;
  it was reactivated in the CUDA facet kernel by user request, then removed
  again on 2026-11-08.  With the copy it was a GPU data race and full-force
  cube drag was `CD=1.88`/`5.17` at `b=0.84`/`0.25` (historical FC=0.6:
  `2.07`/`5.49`); without it the fresh stock `b=0.25` run is `Re=37.9866`,
  `CD=0.912101` vs Eq. (38) `1.268068` (-28.07%, Mach about 0.076), so do not
  use the copy for validation.  The legacy CPU
  `tlbmmpm02` at `b=0.25`/`0.84` with 64 threads gives `CD=4.878`/`1.965`
  vs Eq. (38) `1.250`/`0.581` (+290%/+238%), confirming the legacy cube
  itself was not accurate and that the GPU race is not the only cause.  The
  same copy is what matches the 10 mm MDPI flap deflection at Re=9
  (`Appl. Sci. 10(18):6177`), while no-copy is about 7x too stiff; the
  paper defines `Re=v*w/nu` with a uniform inlet `v` and flap width `w`.
  The MDPI model is **3D**: its abstract says three-dimensional transient
  simulations and Section 2.2.1 gives the tank as (150w,20w,12.5w); its
  bending lines are 2D projections, not a 2D simulation.  The old t01
  parabola used `umax=1.5*Re*nu/Ly` and `dVel[j]`, so it varied
  along z (length), not width.  t01 now uses the paper's uniform inlet
  `v=Re*nu/Ly` at Re=9; its no-copy deflection is 8.78 mm versus 4.99 mm for
  the old parabolic mean-Re=9 case (the parabola was nearly zero at the
  clamp).  The current slab mesh is `ndiv=64`: 3586 tets, 1624 surface
  triangles, minimum surface edge 1.937 mm, so every triangle remains larger
  than the 1.5 mm LBM cell.  Full refined run: mean `dx=9.314 mm`, wall
  1170 s.  `ndiv=80` gives min 1.646 mm (1.10 dx); do not refine past about
  `ndiv=80` without checking the generated facets.  `tlbmmpm_cu_01` takes an
  optional `argv[2]` domain size (default 0.3 m, `n=200`), `argv[3]` inlet
  velocity (default `v=1.815 m/s`), `argv[5]` slab thickness (default
  0.01 m), `argv[6]` `ndiv` and `argv[7]` MPM grid `beta`; the body/report
  centre is `xc=0.5*n*dx`.  Measured `t>=0.75 s` means: 0.2 m -> 8.839 mm
  (std 0.048), 0.3 m -> 9.314 mm (std 0.168), 0.5 m -> 7.867 mm
  (std 0.768).  Larger is not better; 0.3 m is the best compromise.  Current
  default is `Ldom=0.3`, `v=1.815 m/s` (`Re=44.27`), giving a 39.6 mm
  centreline tip; the response is sensitive near v=1.81--1.83, so do not
  overclaim a precise target velocity.  Final centreline plot
  `/home/sergio/DATA/t01_L030_v1815_uniform/t01_v1815_middle_line.png` and
  100-frame VisIt animation
  `/home/sergio/DATA/t01_L030_v1815_uniform/t01_v1815_fluid.gif`.  Do not
  tune a global `ForceCorrection` to satisfy both, and treat the copy as an
  unvalidated ghost-population closure.  It is phase-correct for
  post-streaming data, but the current CUDA copy sits inside the facet loop
  after `newF`; a consistent version needs a separate post-streaming
  `Fghost` pass used by both bounce-back and force.  VisIt animation of the
  refined (`ndiv=64`) run is
  `/home/sergio/DATA/t01_ndiv64_uniform/t01_ndiv64_fluid.gif`.  GPU body-force sweep: `b=2.5e-5`
  gives `Re=0.005`, `CD=4657` vs Eq. (38) 1505; `b=0.25` gives `Re=38`,
  0.913 vs 1.268.  At matched Re, CD depends on bforce, so instantaneous
  constant-force points are not a quasi-steady validation.  The
  2026-11-06 pass fuses Rho/Vel into `cudaCollideMPM` (no `cudaVelDen`
  launch), precomputes `1/Cs`/`1/rho`, and uses the fixed-point limiter
  without intermediate `Ftemp` stores; a single-scan limiter is not
  equivalent.  Measured `Solve`: full t02 window 4.993 -> 4.062 s, full t01
  window 10.664 -> 8.938 s.  Later-frame differences are within the GPU's own
  atomic run-to-run spread; `frame 0000` is exact.
- One-off cases should live outside the checkout and be compiled with
  `/home/sergio/mechsys/mechsyscc case.cu`; do not add them to `tflbm/` unless
  they belong in the library test suite.
- Reference case: `/home/sergio/tflbm10_porous/tclbm10.cu` (random 90%-porosity
  medium, +x body force). Measured on the RTX A6000: 512x256, 20,000 steps,
  `Nthread=256` -> 5.37 s solve / 5.64 s wall, physical `K = 1.08` (raw moment
  1.00).
- Permeability sweep: `/home/sergio/tflbm10_porous/permeability_sweep/`
  (`phi=0.50...0.85`, 512x512, 3 seeds). The fluid phase is random site
  percolation: `K=0` below `phi_c ~ 0.593`; percolating values are
  `K(0.60)=0.0422`, `K(0.70)=0.1163`, `K(0.85)=0.5441`. A `1024x1024` check
  changes `K` by <1% at `phi=0.60/0.65`.
- Random porous media / permeability: `/home/sergio/tflbm_porous/`
  (`porous.cu`, `run_sweep.sh`, `validate_slabs.sh`). With the raw second
  moment it reproduces `K=phi*H^2/12` exactly on periodic straight channels,
  so no `+0.5*nu*phi` shift belongs in `K`. Averaged `K` has an accuracy floor
  `|K|<~0.03` from D2Q9 diagonal corner leakage + staircase walls.
- Plane-Poiseuille validation: `/home/sergio/tflbm_poiseuille/`
  (`poiseuille.cu`, `analyze.py`, `plot.py`, `README.md`; the older
  `tflbm10_porous/validation/` tree no longer exists). Walls at `y=0.5` and
  `y=H+0.5` (half-way bounce-back: `H` = number of fluid rows, verified with
  `nsolid=1/3/6`), uniform `+x` force, `Dom.Sc=0` to keep LES off. The profile
  is exactly parabolic with `nu_fit = nu` and a constant slip
  `u_raw = a/(2 nu) s(H-s) - c(tau) a`, `c = 1/(8 nu) - 6 nu + 1/2` (fits 23
  runs to 4 digits, independent of `H`). At `tau=1` this makes
  `K_raw = <u>nu/fx = H^2/12` to `1e-11` for `H=8..128` **with no correction**;
  pointwise L2 error `7.13e-3/1.78e-3/4.46e-4/1.11e-4/2.79e-5` for
  `H=8/16/32/64/128` (second order, `u_slip/u_max = (1-48nu^2+4nu)/H^2`), and
  `|dK/K| <= 1.1e-3` for `0.56 <= tau <= 1.1` at H=32 (4x smaller at H=64).
  Scale `fx` with `H^2` to keep `u_max ~ 0.01`. The old half-force-based
  numbers (error `1.56e-2` at H=8 …) were produced by a tree whose reported
  velocity differed; re-measure before applying `u+0.5*dt*fx` to a porous `K`.
- Advection-diffusion point source: `/home/sergio/tflbm10_porous/solute_source/`.
  Two-layer `AdvectionDiffusion` (`nu={1/6,D}`): layer 0 is the fluid, layer 1
  the passive concentration. Set `BForce[1]=0.5*fx/Tau[1]` so the AD collision
  advects at `u_phys=u_raw+0.5*dt*fx/rho`; a source must add to both the
  populations and `Rho[1]`. The entropic limiter in `lib/flbm/lbm.cuh` is
  guarded in this working tree; without the guard, negative solute populations
  produce NaN in the zero-concentration far field at 512x512. `mass_fluid`
  (pore sum) is about 14% below `mass_total` (includes bounce-back populations
  on static solids), while `mass_total` matches `S(t-t_on)`. Production
  512x512, phi=0.70, seed=12345, fx=1e-2, D=0.02, S=1e-4, t=20k..60k: 67.4 s
  wall, `mass_total=4.0`, `Cmax=1.09e-3`; the plume wraps between 55k and 60k.
- Obstacle channel case: `/home/sergio/tflbm_triangle/` (`tclbm_obstacle.cu`,
  `--shape=rect|tri`, `Re=u_mean*D/nu`, 640x201, `D=40`). Rows 0 and `ny-1`
  are solid, so the walls are at `y=0.5` and `y=ny-1.5`, the channel width is
  `L=ny-2` and the symmetry mirror is `j <-> ny-1-j`; using `L=ny-1` makes the
  whole solution asymmetric. At Re=40 both wakes are steady and symmetric:
  `Lr/D=1.88` for the 60x40 rectangle and `2.88` for the triangular box
  (grid-converged 2.86, Mach-independent); no shedding up to Re=100 at 20%
  blockage.

## VisIt visualization

**Always use the project fluid colormap, `turbo`, for plots and animations**
(`colorTableName = "turbo"` in VisIt; `matplotlib.colormaps["turbo"]` in
Python).  Do not substitute `hot`, `bwr`, `jet`, or a custom blue-white-red
map unless the user explicitly asks for a different map.

VisIt 3.3.2 is installed at `/usr/local/visit/bin/visit`; `~/.profile` adds
`/usr/local/visit/bin` to `PATH`.  Use it headlessly:

```bash
/usr/local/visit/bin/visit -cli -nowin -s check.py data.xmf
```

A verified minimal script is:

```python
OpenDatabase(sys.argv[1])
md = GetMetaData(sys.argv[1])
print(md.GetFileFormat(), md.GetNumStates())
print([md.GetMeshes(i).name for i in range(md.GetNumMeshes())])
print([md.GetScalars(i).name for i in range(md.GetNumScalars())])
AddPlot("Mesh", "DEM_Faces")
DrawPlots()
SaveWindow()
```

VisIt 3.3.2 accepts MechSys' original one-dimensional HDF5 XDMF datasets;
do not change the dataset rank just to satisfy VTK's XDMF reader.  For a
temporal collection, keep every child `Grid Name` equal (for example
`DEM_Faces`) or VisIt exposes the first child's name as the mesh variable.
For pseudocolor, use the full variable name for a standalone file
(`DEM_Faces/Velocity`) and the short name for a collection (`Velocity`);
the Threshold variable must likewise be `DEM_Faces/Tag` versus `Tag`.
`SaveWindow` appends a frame number to the requested filename
(`visit.png0000.png`).

For the t01 LBM/MPM slab animation, create `lbm.visit`/`mpm.visit`
playlists, plot `speed=magnitude(Velocity_0)` with `turbo` on a `Slice`
plane (`y=0.15`, normal `+y`), overlay the deformed `MPM_mesh`, and call
`quit()` after the save loop.  A `Slice` is ~3 s/frame at 900x900; a `Clip`
of the same 200^3 grid is much slower.  Current target run:
`/home/sergio/DATA/t01_L030_v1815_uniform/t01_v1815_fluid.gif`; fixed speed
range 0--2.0 m/s for v=1.815 m/s.  The direct `position.res` plot is
`/home/sergio/DATA/t01_L030_v1815_uniform/t01_v1815_position_file.png`
(`plot_position_res.py`); the file has normalized `(x-xc)/Lz` and `z/Lz`
columns, 731 rows including the clamped z<0 section.  The `Delta z` vs
`Delta x` version is
`/home/sergio/DATA/t01_L030_v1815_uniform/t01_v1815_position_delta.png`
(`plot_position_delta.py`); initial coordinates come from the matching
frame-0 particle ids.  The requested `z` vs `Delta x` plot is
`/home/sergio/DATA/t01_L030_v1815_uniform/t01_v1815_position_z_dx.png`
(`plot_position_z_dx.py`).  For a half-thickness slab (`Lx=0.005 m`) keep the
MPM x padding at `PadX=Lz`, not `10*Lx`; the narrower box crashes with
`nhigh out of bounds`.  Completed half run:
`/home/sergio/DATA/t01_L030_v1815_Lx005_padLz/`, mean tip 83.34 mm, final
83.11 mm, wall 1081.20 s.  Refining the thickness grid
(`ndiv=80,beta=0.75` or `beta=0.5`) is unstable, so that run is
under-resolved/qualitative.  Both middle planes are overlaid in
`/home/sergio/DATA/t01_L030_v1815_uniform/t01_v1815_middle_planes_10_5mm.png`
(`plot_both_middle_planes.py`), using `z` vs `x-xc` for both.  Its
100-frame animation is
`/home/sergio/DATA/t01_L030_v1815_Lx005_padLz_video/t01_Lx005_fluid.gif`
(poster `t01_Lx005_fluid_poster.png`); the video rerun took 1137.84 s and
VisIt rendering 137.70 s.  Raw `*.h5`/`*.xmf` frames for both t01 runs
were deleted; rerun with `RenderVideo=1` to regenerate fields.

## Exactness rules (critical)

- Compare against a **clean worktree of `HEAD`**, not against the modified
  working tree.
- For CPU DEM, exactness is checked at `Nproc=1` by comparing **HDF5 datasets**,
  not whole-file MD5. HDF5 container metadata is not deterministic even
  baseline-to-baseline. Expected result is `neq []`.
- CPU MPM follows the same rule: build the reduced driver against a copy of the
  pre-change `lib/mpm` headers and compare HDF5 datasets at `Nproc=1`
  (bitwise).  At `Nproc>1` the node sums can be reordered; a healthy long run
  differs by <=1e-13 relative in `PVelocity`/`Sigma`.
- At `Nproc>1`, CPU force accumulation order is thread-order dependent. Small
  last-digit differences are expected, not automatically a regression.
- CUDA `test_cu_03/04` use `atomicAdd`; the original baseline is itself
  non-deterministic. Do not claim bitwise equality for those GPU runs.
- Run the real test window, not only a short smoke test, before making a final
  accuracy claim. State the window and thread count.

Useful comparison pattern:

```python
import h5py, numpy as np
b = h5py.File("baseline.h5", "r")
o = h5py.File("optimized.h5", "r")
for k in b:
    if isinstance(b[k], h5py.Dataset) and not np.array_equal(b[k][()], o[k][()]):
        print("DIFF", k)
```

## Hot-path invariants

These are easy to break accidentally:

1. **Preserve candidate-list order.** Contact forces are accumulated in list
   order. Sorting `L` changes floating-point results even if contacts are
   identical. Sorted distances, if used, must remain separate from `L`.
2. **Keep the baseline `Distance` arithmetic.** Do not replace full distance
   calculations with cached normals or squared-distance variants in
   `_update_contacts`. `Distance2` is only an early rejection in the force loop.
3. **The CPU Verlet cutoff uses `4.0 * ContactDispGlobal()`.** A factor of two
   was observed to skip real contacts (`Ld` just outside the bound that became
   a contact after motion). `ContactDispGlobal()` accumulates the
   `MaxDisplacement` values reported each step; do not lower the factor
   without a new safe derivation and an exact `Nproc=1` comparison.
4. **Do not sort the contact list or change original `(i,j)` push order.**
5. **CUDA movement:** keep the separate translate/rotate fallback for custom
   movement and for cases where the fused movement kernel is not valid.
6. **CUDA thread blocks:** `Nthread=512` fails for the current `test_cu_03`
   register usage. Valid tested range is up to 384; large test optimum is
   around 128, reduced test optimum around 32.
7. **CRLF:** `tdem/test_cu_04.cu` has CRLF line endings. Do not rewrite it with
   an editor/tool that normalizes line endings; it creates a huge meaningless
   diff.
8. **Shared headers affect both CPU and CUDA.** After changing `lib/dem/*`,
   build at least one CPU target, one CUDA target, and the LBM-DEM targets when
   relevant.
9. **No physics shortcuts.** Broadphase rejections must be conservative; if in
   doubt, remove the rejection rather than risk a missed contact.
10. **Clean warnings.** The only expected third-party warning is Blitz
    `tinymat.h(60): warning #191-D`. Do not add new warnings in MechSys code.

11. **MPM hot path.** `Node` is `alignas(64)` with
    `Mass/Mn/Fn/valid/lck` in the first cache line; `VNodes` is built by first
    touch inside `ParticleToNode`; `Node::Reset` deliberately leaves `Vn`/`Sn`
    alone; and `FillShapeCache3D` must keep the same multiplication order as
    `GIMP3D`/`BSpline*3D`.  Changing any of these breaks either the speedup or
    the `Nproc=1` bitwise exactness.  See
    [`MPM_OPTIMIZATION_CHANGES.md`](MPM_OPTIMIZATION_CHANGES.md).

12. **Rotation caches are derived state.** `Particle::wrot` and `Qc` are
    refreshed by `CalcProps`, `InitializeVelocity`, `Rotate`,
    `Domain::UpdateContacts`, and after a user `Setup` callback. The CUDA
    upload copies these fields; `DnloadParticle` copies them back. A custom
    GPU kernel that changes `Q`/`w` must refresh `DPar[].wrot`/`Qc`. If host
    code assigns `Q` or `w` directly and then calls an interaction force
    without going through the normal paths, call `Particle::UpdateRotCache()`
    first.

13. **CUDA MPM hot path.** Hot node fields are SoA arrays
    (`pMass/pMn/pFn/pVn/pValid`); the scatter builds the touched-node list
    with `atomicExch` + `atomicAdd`, the previous list is reset first, and the
    counter is ping-ponged. Do not remove `pValidNumberPrev` as "redundant":
    a single-counter trial stayed deterministic but moved the non-CPDI
    `tmpm_cu_01` trajectory (final `.res` difference `8.74e-3 m`) and gave no
    speedup; CPDI stayed exact. Reset/solve launch over `Nnodes` and exit early on
    the device count. Keep the final `cudaDeviceSynchronize()` in `Solve`.
    `ParticleCU` now holds 64 nodes per particle (cubic B-spline limit);
    quadratic/GIMP normally use 27. GPU sums use atomics, so never promise
    bitwise equality. The `BodyMesh!=NULL` CPDI branch is implemented with
    `cudaCornerToNode`, `cudaParticleToNodeCPDI`, `cudaMoveCorners` and
    `cudaNodeToParticleCPDI`; keep the launch order matching `OneStepCPI`.
    LBMMPM still runs `MPMDOM.OneStepCPI` on the CPU.  The checked-in
    `tmpm01`/`tmpm_cu_01` pair takes `[Nproc] [cpdi|mpm]` (default non-CPDI).
    CPDI uses `ndiv=125`, `beta=1`, and loads only boundary corners in the tip
    slice; non-CPDI keeps `ndiv=50` and loads only `ForcePar`.  The CUDA version keeps `ForceCor` as a device index list and
    sets corner `h` with `SetCornerForce`.  The 2026-11-06 CPDI pass caches
    `An=Fn/Mass` (`cudaComputeAn`) and `GradvecN=Gradvec/(6*Vt)` in
    `ParticleCU`, keeps particle/corner state in registers, and drops the
    redundant valid-node rebuild in `cudaParticleToNodeCPDI` because the corner
    scatter already covers every node.  It cut `ndiv=125` `Solve` from
    17m23s to 4m03s; CPDI `.res` is unchanged to `1e-16`, while HDF5 fields
    shift by round-off.  The non-CPDI scatter evaluates separable 1-D GIMP/
    B-spline factors (`MPMShapeFunction1DCU`) and its gather also uses
    `An=Fn/Mass`; `Solve` fell from 1m54.8s to 21.03s, with the final `.res`
    matching the previous GPU baseline to `1e-13 m`.  HDF5 payloads are not
    promised bitwise after these caches; the observed `PVelocity` shift is
    about `7e-8` relative.  For this rectangular beam, the CPDI tip load
    selects corners by `x` alone because the last cross-section is an external
    face; do not generalize that shortcut to a body mesh that can have
    interior vertices in the last slice.  See
    [`CUDA_MPM_OPTIMIZATION_CHANGES.md`](CUDA_MPM_OPTIMIZATION_CHANGES.md).

## Performance rules

- Measure wall time, not only printed per-step solver time. Small tests are
  dominated by construction/output overhead.
- Report baseline and optimized on the same machine, input, and thread count.
- Separate solve time from total wall time.
- For CPU Voronoi, use the reduced input in
  `benchmarks/cpu_voronoi/tper_small.inp` for exactness; use
  `tper_med.inp` for a larger timing check.
- For CUDA, re-benchmark after struct/layout changes; do not reuse stale
  numbers. See `REDUCED_TEST_CU_03_PERFORMANCE.md` for the reduced sweep.
- Recommend a thread count based on the actual problem size; the large and
  reduced CUDA workloads have different optima.

## Preferred working style

- Keep documentation in the existing root markdown files (`MECHSYS_HANDBOOK.md`,
  `MEMORY.md`, `CPU_DEM_OPTIMIZATION_CHANGES.md`, ...). Do not add a `docs/`
  tree or other new top-level structure; MechSys upstream has none.
- Reply in English unless the user asks otherwise.
- Give explicit commands, file paths, and measured numbers.
- Distinguish “compiles”, “runs”, “matches baseline”, and “is faster”.
- If a verification was not done, say so.
- The per-prompt memory checkpoint above is mandatory; do not treat it as
  optional cleanup.
