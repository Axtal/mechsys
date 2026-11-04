# MechSys Handbook

Durable working knowledge for users and coding agents working on the MechSys
checkout in this workspace. This is a living document: update it when a hard
lesson, benchmark result, or verification procedure changes.

> Automatic agent instructions live in [`AGENTS.md`](AGENTS.md). This handbook
> is the detailed reference those instructions point to.

## 1. What MechSys is

MechSys is an open C++ programming library for mechanical-system simulation.
The parts exercised in this workspace are:

- DEM (discrete element method) for granular and particle systems
- LBM (lattice Boltzmann method) for fluids
- LBM-DEM coupling for fluid-particle interaction
- mesh, geometry, and supporting numerical utilities

The library is mostly header-based/template-based. Test programs are small
standalone drivers under `tdem/`, `tlbmdem/`, `tlbm/`, `tflbm/`, and similar
directories. There is no single “MechSys executable”; each test is the entry
point for its feature.

Repository root:

```text
/home/sergio/mechsys
├── lib/                    # actual headers, also reachable as ./mechsys -> lib/
├── tdem/                   # DEM tests, including the CUDA tests
├── tlbmdem/                # LBM-DEM tests, CPU and CUDA
├── benchmarks/             # reduced benchmark inputs and measurement data
├── CMakeLists.txt          # configured build
├── Makefile                # configured top-level Makefile
├── CUDA_OPTIMIZATION_CHANGES.md
├── REDUCED_TEST_CU_03_PERFORMANCE.md
├── CPU_DEM_OPTIMIZATION_CHANGES.md
├── AGENTS.md
├── MEMORY.md
└── MECHSYS_HANDBOOK.md
```

`mechsys -> lib/` is intentional: code includes headers as
`<mechsys/dem/...>` or `<mechsys/...>` while the real files live in `lib/`.

## 2. Build environment in this workspace

Known working configuration:

| Item | Value |
|---|---|
| Compiler | `/usr/bin/c++` (GCC 15.2.0) |
| OpenMP | enabled (`-fopenmp`) |
| CUDA | 12.4, `nvcc` V12.4.131 |
| GPU | NVIDIA RTX A6000, 46 GB |
| CUDA architecture | `sm_86` / `compute_86` |
| CUDA real type | `double` (`-DUSE_GPU_DOUBLE`) |
| HDF5 | `/home/sergio/pkg/hdf5-1.8.15-patch1` |
| Blitz++ | `/home/sergio/pkg/blitz-0.9` |
| Voro++ | `/home/sergio/pkg/voro++-0.4.5` |
| gmsh | `/home/sergio/pkg/gmsh-4.15.2` |
| igraph | `/home/sergio/pkg/igraph-0.8.2` |

The root is already CMake-configured. Do not reconfigure casually; the
generated `Makefile`, `CMakeCache.txt`, and `CMakeFiles/` are part of the
working state and are normally ignored by git.

### Build commands

```bash
cd /home/sergio/mechsys/tdem
make -j8 test_periodic
make -j8 test_cu_01 test_cu_02 test_cu_03 test_cu_04

cd /home/sergio/mechsys/tlbmdem
make -j8 tlbmdem01 tlbmdem02
make -j8 tlbmdem_cu_01 tlbmdem_cu_02 tlbmdem_cu_03 tlbmdem_cu_04
```

A single test can be built by name. The configured Makefiles also expose the
usual CMake targets.

### Running

`test_periodic` takes a **FileKey prefix**, not a complete filename; it appends
`.inp` and output suffixes:

```bash
cd /home/sergio/mechsys/tdem
./test_periodic /home/sergio/mechsys/benchmarks/cpu_voronoi/tper_small 8
```

`test_cu_03` takes `Nproc`, `Alpha` (Verlet distance), and `Nthread`:

```bash
cd /home/sergio/mechsys/tdem
./test_cu_03 1 0.05 128
```

`tlbmdem01` and `tlbmdem02` are very large. They are useful build/integration
checks, but a full run can exceed several minutes. For a smoke test, make a
temporary reduced copy (smaller grid and `Tf`) rather than changing the
checked-in test.

## 3. DEM architecture

### 3.1 Main types and files

| File | Responsibility |
|---|---|
| `lib/dem/domain.h` | `Domain`, particle container, interaction lists, `Solve`, output |
| `lib/dem/particle.h` | `Particle`: position, velocity, quaternion, vertices/edges/faces, movement/rotation |
| `lib/dem/interacton.h` | `CInteracton`, feature-list construction, contact force/friction |
| `lib/dem/distance.h` | feature distance functions, broadphase `Overlap` tests |
| `lib/dem/face.h`, `edge.h`, etc. | geometry primitives and cached geometry |
| `lib/dem/dem.cuh`, CUDA headers | device data structures and kernels |
| `lib/lbmdem/Domain.h` | CPU/GPU LBM-DEM coupling |

### 3.2 DEM solve cycle

The conceptual cycle in `Domain::Solve` is:

1. reset particle forces/torques (including fixed `Ff`/`Tf`);
2. evaluate all interactions in the current list;
3. accumulate force/torque per particle;
4. translate and (if enabled) rotate particles;
5. compute maximum vertex displacement since the last list rebuild;
6. if displacement exceeds `Alpha`, rebuild the interaction/Verlet lists;
7. advance `Time` and `iter` once.

The order matters. In particular, contact lists are rebuilt from the state
after movement, and `MaxDisplacement` compares current vertices against
`Vertso`, which is reset by `ResetDisplacements()` during
`Domain::UpdateContacts()`.

### 3.3 Feature pairs and Verlet lists

For non-spherical particles, a contact pair is evaluated through feature
lists:

- `Lee`: edge-edge
- `Lvf`: vertex-face
- `Lfv`: face-vertex
- plus torus/cylinder lists where relevant

`CInteracton::UpdateContacts(alpha)` builds candidate lists using a broadphase
and a distance threshold of roughly `Rsum + 2*alpha`. Forces are then computed
from the candidate lists.

`FrictionMap_t` is an `std::unordered_map<size_t, Vec3_t>` storing tangential
displacement history per feature pair. Friction history must survive contact
list rebuilds while the pair remains present.

### 3.4 CPU Voronoi hot path (current optimized state)

The CPU optimization report is
[`CPU_DEM_OPTIMIZATION_CHANGES.md`](CPU_DEM_OPTIMIZATION_CHANGES.md). The
important implementation facts:

- `Particle` stores `wrot = Rotation(w,Q)` and `Qc = Conjugate(Q)`. They are
  refreshed in `Particle::UpdateRotCache()` from `CalcProps`,
  `InitializeVelocity`, `Rotate`, `Domain::UpdateContacts`, and after a user
  `Setup` callback. Interactions
  read `P1->wrot`, `P2->wrot`, `P1->Qc`, `P2->Qc`; there are no
  per-interaction cache fields. Domain-side refreshes use
  `UpdateRotCacheIfChanged()` so setup callbacks that do not touch `Q`/`w`
  pay only a comparison; `Rotate`/`CalcProps` refresh unconditionally.
- `Particle::Rotate` fuses the old/new orientation passes and refreshes the
  caches after installing the new `Q`, mirroring the GPU kernels.
- `Overlap` has safe AABB/plane broadphase with a `1e-12` margin.
- `Distance2` is used only as an early rejection in the force loop. Actual
  contact geometry uses the baseline `Distance` implementation.
- Each feature list stores its build-time distance (`Leed`, `Lvfd`, etc.).
- The force loop keeps original list order and skips a candidate when
  `Ld[k] > Rsum + 4.0*ContactDispGlobal()`.
- `ContactDispGlobal()` adds the `maxdis` value reported each step since the
  last rebuild and is reset in `Domain::UpdateContacts()`. `MaxDisplacement`
  itself is measured from the last `ResetDisplacements()`, so this is a
  conservative accumulated value, not a mathematically tight feature bound.
- `_update_contacts` uses full `Distance`, original nested-loop order, and does
  not sort the contact list.
- CPU forces use per-thread `MtData::F`/`T` arrays and a final reduction; the
  per-interaction particle locks were removed.
- The interaction loop uses OpenMP `schedule(guided)`.

**The exactness-critical parts are:** original list order, full `Distance` in
contact rebuild, and the conservative cutoff factor. Do not “clean up” these
without re-verifying at `Nproc=1`.

### 3.5 CPU benchmark baseline

Reduced input: 480 Voronoi particles, `L=10`, `n=10`, `T0=Tf=0.01`,
`dt=5e-5`, 400 total steps. Stored at
`benchmarks/cpu_voronoi/tper_small.inp`.

Measured on this machine:

| Threads | Baseline wall | Optimized wall | Speedup |
|---:|---:|---:|---:|
| 1 | 13.94 s | 9.27 s | 33.5% wall |
| 8 | 3.20 s | 2.55 s | 20.3% wall |

Solve-only totals: `Nproc=1` 12.73 s → 8.06 s (36.7%); `Nproc=8`
1.99 s → 1.36 s (32.0%).

Medium check: 20³ seeds, about 4,000 cells, `Nproc=8`. In the latest matched
run, total wall was 26.06 s → 21.64 s (about 17% faster) and the solve phase
was 16.40 s → 12.00 s (about 27% faster). An earlier run under different
machine load gave 32.55 s → 25.71 s (about 21% wall / 30% solve), so treat
absolute medium timings as run-dependent.

A separate A/B run isolated the particle-level rotation-cache change on the
medium input: two-solve total 12.13 s → 12.02 s (about 0.9%), second solve
6.62 s → 6.48 s (about 2%). The reduced 480-particle case is noise-limited.

At `Nproc=1` every final HDF5 dataset is bit-identical to the `HEAD` baseline
for the reduced run. At `Nproc>1`, small differences are expected because the
original accumulation order is thread-dependent.

### 3.6 VisIt visualization

VisIt 3.3.2 is available at `/usr/local/visit/bin/visit`; `~/.profile`
adds its `bin` directory to `PATH`.  Use the CLI headlessly when no GUI is
needed:

```bash
/usr/local/visit/bin/visit -cli -nowin -s visit_check.py data.xmf
```

The verified Python API for a MechSys XDMF file is:

```python
OpenDatabase(f)
md = GetMetaData(f)
print(md.GetFileFormat(), md.GetNumStates())
meshes = [md.GetMeshes(i).name for i in range(md.GetNumMeshes())]
scalars = [md.GetScalars(i).name for i in range(md.GetNumScalars())]
AddPlot("Mesh", "DEM_Faces")
DrawPlots()
SaveWindow()
```

MechSys writes `Verts`/`FaceCon` as one-dimensional HDF5 datasets even
though the XDMF DataItems advertise `(N,3)`.  **VisIt 3.3.2 accepts this
MechSys convention.**  Only the XMF path must be the HDF basename; do not
reshape the datasets solely for VTK's stricter reader.

A single `WriteXDMF` call is one state.  For animation, write per-frame
pairs and combine them into a temporal collection whose child grids all have
`Grid Name="DEM_Faces"`.  VisIt then reports the mesh as `DEM_Faces`; if the
children use different names, the first name becomes the mesh variable.

Plotting pseudocolor with transparent walls:

```python
scalars = [md.GetScalars(i).name for i in range(md.GetNumScalars())]
pvar = next(s for s in scalars if s.endswith("Velocity"))
tagvar = next(s for s in scalars if s.endswith("Tag") and "Faces" in s)
AddPlot("Pseudocolor", pvar)          # DEM_Faces/Velocity standalone
AddOperator("Threshold")
t = ThresholdAttributes()
t.SetDefaultVarName(tagvar)           # DEM_Faces/Tag standalone, Tag collection
t.SetListedVarNames(tagvar)
t.SetLowerBounds(0.5); t.SetUpperBounds(1.0e9)
SetOperatorOptions(t)

AddPlot("Mesh", "DEM_Faces")
m = MeshAttributes(); m.SetMeshColor((0.55,0.55,0.55,255))
m.SetMeshColorSource(0); m.SetOpacity(0.18)
SetPlotOptions(m)
AddOperator("Threshold")
w = ThresholdAttributes(); w.SetDefaultVarName(tagvar)
w.SetListedVarNames(tagvar); w.SetLowerBounds(-1.0e9); w.SetUpperBounds(-0.5)
SetOperatorOptions(w)
DrawPlots()
```

MPM uses a different grid name and attribute set.  Each `WriteXDMF` frame is
a standalone single-state file with grid `MPM_points` (`Polyvertex` topology)
and node attributes `Tag`, `Velocity`, `Pressure`, `Sigma`, `Strain`.  A
verified point-cloud render is:

```python
OpenDatabase("granular_0019.xmf")
DefineScalarExpression("speed", "magnitude(Velocity)")
AddPlot("Pseudocolor", "speed")
pc = PseudocolorAttributes()
pc.colorTableName = "hot"
pc.limitsMode = 0; pc.min = 0.0; pc.max = 1.2   # fixed scale across frames
pc.pointSizePixels = 3
SetPlotOptions(pc)
DrawPlots()
ResetView()
v = View3DAttributes()
v.viewNormal = (0.75, 0.35, 0.56); v.viewUp = (0.0, 1.0, 0.0)
v.focus = (0.13, 0.045, 0.10); v.parallelScale = 0.20; v.perspective = 0
SetView3D(v)
sw = SaveWindowAttributes()
sw.fileName = "granular_0019"
sw.format = sw.PNG          # not the string "png"
sw.width = 900; sw.height = 900
SetSaveWindowAttributes(sw)
SaveWindow()
```

For a **fixed, pressure-proportional color scale**, use the plot's min/max
flags; `limitsMode` itself only accepts `0` (`OriginalData`) or `1`
(`ActualData`), not a user-specified mode:

```python
AddPlot("Pseudocolor", "Pressure")
pc = PseudocolorAttributes()
pc.colorTableName = "turbo"         # blue low -> red high
pc.limitsMode = 0
pc.minFlag = 1; pc.maxFlag = 1
pc.min = -30000.0; pc.max = 110000.0
SetPlotOptions(pc)
```

A 0.5 s water-column run (`tmpm_cu_03`, 134,820 points, 40 frames) took
8 min 43 s and produced `Pressure` from 11-774 Pa initially, through a
+105 kPa / -27 kPa transient at 0.25 s, to -24/+81 kPa at 0.49 s.  The
fixed -30/+110 kPa scale keeps the color mapping linear and comparable
between frames.

The periodic square-duct case (`tmpm_cu_04`) is different: with its
`EOS_LINEAR` EOS, initially uniform density and x-periodic domain, the
`Pressure` dataset is **identically zero in every frame**.  A pressure
colormap is therefore meaningless for this case; visualize the axial
velocity instead (`DefineScalarExpression("Vx","Velocity[0]")`).  A full
`Tf = 200 s` run (216,000 points, 50 frames) took 22 min 57 s solve and
reached `Vx_max = 7.68731e-4 m/s` versus the analytical `7.36712e-4 m/s`
(+4.35%).  The duct is low-Reynolds (`Re ~ 0.74`), so the profile develops
smoothly with no pressure signal in this formulation.

In VisIt 3.3.2 `SaveWindowAttributes` has no `resolution` member, and
`outputDirectory` is not honored from a CLI script: the PNG is written in the
process launch directory (or use `os.chdir`).  `SaveWindow` still appends
`0000` to the name.

`SaveWindow` appends `0000`, `0001`, ... to the requested name.  For a time
series, call `SetTimeSliderState(i)` and `DrawPlots()` before each save.
VisIt CLI may return a nonzero status after a successful save if the script
leaves a warning; check the saved files and the log, not only the exit code.

### 3.7 Periodic boundaries with finite-size particles

`DEM::Domain` uses `Xmin/Xmax`, `Ymin/Ymax`, `Zmin/Zmax`; `Solve()` builds
`Per`, `ResetDisplacements()` wraps free centres, and `CInteracton` applies
`Per` only when both particles report `IsFree()`.  The reference setup is
`tdem/test_periodic.cpp`:

- the periodic domain is several particle lengths across, independent of the
  particle size;
- fixed plates are created with `AddPlane`, then constrained with
  `FixVeloc()` and **not** `FixFree`;
- the plates are wider than the periodic cell: `Lx*2` by
  `Ly + 4*maxd`, where `maxd=dom.MaxDim()`.

For the avalanche driver's rounded boxes, setting the periodic length equal to
the 6 cm slot width (`periodic_Ly=0.06 m` for 0.0456 m boxes) repeatedly
aborted with a free-box/free-box overlap at `y≈±0.03 m`.  Enlarging only the
fixed planes, or setting `FixFree` on them, did not fix it.  A separate
`periodic_Ly=0.20 m` (4.4 box lengths), with the same nominal 6 cm insertion
width, non-free planes, and plane span `periodic_Ly+4*L`, completed 800 events
cleanly (`./ricepile rice_narrow_lx3_periodic_mu025.inp 32`, final simulated
time 414.207 s, exit 0).  With burn=100 and fixed `xmin=5`, the wider periodic
cell gives `S_path` alpha=2.302+/-0.123 (n_tail=429, KS p=6.8e-4) for the
800-event seed 12345; that burn starts during pile growth (`h=0.095 m` at
event 100), and selecting only `h>=0.30 m` lowers `S_path` alpha to
2.144+/-0.174 (184 events).  Pooling a second 402-event seed 67890 gives
alpha=2.384+/-0.118 (n_tail=532, KS p=0.0149), or `N_moved`
alpha=2.743+/-0.231 (n_tail=218, KS p=9.2e-10).  Final pile height is 0.365 m
for seed 12345 and 0.227 m for seed 67890, versus 0.663 m for the walled
mu=0.25 run.  There are fewer large avalanches: seed 12345 has S>=50: 11 vs 24
and N>=50: 1 vs 19; the pooled sample has S>=50: 12 and N>=50: 1.  The periodic exponents are
therefore not closer to the rice value 1.12-1.21, and the poor KS p-values mean
the power-law form is not a good description of either periodic tail.  Do not
blame linked-cell candidate generation: at `1.5Lx=0.06845 m` the run aborts at
event 10 whether `Domain::Beta` is 0.5, 1, 2, or 4, and whether
`Domain::Alpha` is 0.00225, 0.005, 0.01, or 0.02; `Beta=20` with `Alpha=0.05`
puts all boxes in one cell and still aborts at `2Lx` on a boundary-straddling
pair; `3Lx` survives 300 events but aborts near event 835; its first 832 events give
`S_path` alpha=2.094+/-0.095 (n_tail=512, KS p=0.115).  The remaining issue is
the finite-box periodic contact/force geometry.
Use the periodic length in the lost-particle test as well; otherwise boxes
outside the nominal slot are wrongly parked.  Changing the periodic length changes the
effective pile width and therefore the avalanche statistics.

### 3.8 Save/load and continuation

`Domain::Save`/`Domain::Load` are a complete DEM checkpoint: positions,
velocities, both Verlet references `xb` and `wb`, orientations, constraint
flags (`vxf..wzf`, `FixFree`, `Eroded`, `Closed`, `Cluster`), properties,
particles, interaction history and per-tag friction overrides.

`Load` sets `Initialized=true` and intentionally does **not** reinitialize the
loaded particles: it preserves the saved Verlet references `xb` and `wb`, which is what
makes an in-memory continuation and a load continuation bit-identical.
`FreePar`, `NoFreePar` and `Particle::Index` are derived state: `Solve()`
clears and rebuilds them after `Initialize()`, so they are intentionally empty
right after `Load`.  Code that inspects or extends a loaded domain before the
first `Solve()` must rebuild the partition itself.

Two consequences:

- A particle added **after** a load is not initialized by the next `Solve()`
  (the domain is already marked initialized).  Call
  `p->InitializeVelocity(dt)` after adding or repositioning it, or the
  uninitialized `xb` makes the first Verlet step invalid.
- `Load` appends to the current particle array.  Load into an empty domain
  unless appending is intended.

Regression command:

```bash
cd tdem
make -j8 test_saveload
./test_saveload 1 300 300
```

The test compares positions, velocities, `xb`, `wb`, orientations, all
constraint flags and sphere friction history.  After the continuation
`Solve()` it also compares `FreePar`/`NoFreePar` and `Particle::Index`, which
`Solve()` owns.
Expected result after 300 saved steps and 300 more steps:

```text
max|dxb|=0.000e+00 max|dw|=0.000e+00 max|dwb|=0.000e+00   flags=MATCH   BIT IDENTICAL
after 300 more      BIT IDENTICAL
free partition      ld: 3 free / 1 fixed   dom: 3 free / 1 fixed   MATCH
RESULT: PASS - the reloaded domain is bit-for-bit the same
```

A freshly loaded domain is also stable before any insertion: loading an
802-particle wall/box state and solving 0.20 s gave
`z=[0.00996,0.34906]`, no NaNs, no particles below the base.

## 4. CUDA DEM

The CUDA optimization report is
[`CUDA_OPTIMIZATION_CHANGES.md`](CUDA_OPTIMIZATION_CHANGES.md), with the
reduced thread/Verlet study in
[`REDUCED_TEST_CU_03_PERFORMANCE.md`](REDUCED_TEST_CU_03_PERFORMANCE.md).

### 4.1 Data structures

- `DynParticleCU` stores the dynamic particle state, including cached rotation
  vectors and cached conjugate quaternions (`wrot`, `Qc`) rather than separate
  device arrays.
- The host prepares `Particle::wrot`/`Qc` in `CalcProps`/`InitializeVelocity`;
  `UploadParticle` copies them to the device, and `DnloadParticle` copies the
  device values back. The upload no longer calls `Rotation`/`Conjugate`.
  `UpLoadDevice` still refreshes the host caches if `Q`/`w` changed, and the
  device movement kernels refresh the device caches after rotation.
- The layout remains AoS (array of structs). AoS was retained after comparison;
  keep it unless a measured, verified SoA migration is requested.
- Device interaction lists are separate from host `Interactons`.

### 4.2 Kernels and dispatch

The main CUDA flow in `Domain::Solve` is:

1. `pReset`
2. `pForceVV`
3. `pForceEE`
4. `pForceVF`
5. `pForceFV`
6. `pMove` (fused) or `pTranslate` + `pRotate` fallback
7. `MaxD`
8. optional `UpdateContactsDevice`

Host-selectable function pointers allow custom force laws and movement. The
default dispatch uses the specialized kernels. The fallback path must remain
valid for custom movement; an earlier attempt to always fuse movement broke
`test_cu_04`.

`MaxD` writes one partial maximum per block; the host reduces only the block
entries, not the whole per-vertex array.

### 4.3 CUDA performance facts

For the large `test_cu_03` workload (77,003 particles), representative
before/after measurements were about 99.06 ms/step → 84.51 ms/step, roughly
14.7%. These values predate the final struct/movement refactor: re-benchmark
before quoting them as current.

Thread-count behavior for that large case:

- `Nthread=64`: about 83–84 ms/step
- `Nthread=128`: best, about 82.8 ms/step
- `192`: about 83–84 ms/step
- `224–384`: about 85–86 ms/step
- `512`: launch failure (too many resources)

For the reduced 1.6k-particle test, the optimum is much smaller: about 32
threads per block, and a Verlet distance around `0.01–0.02` rather than `0.05`.
This is a different workload; do not transfer the large-case optimum blindly.

### 4.4 CUDA accuracy caveat

The original CUDA `test_cu_03/04` use `atomicAdd`, so the baseline itself is
non-deterministic run-to-run. Bitwise equality with a single baseline run is
not a meaningful requirement. Use energy trends, tolerance-based comparisons,
or deterministic reduced configurations, and state exactly what was compared.

## 5. LBM, LBM-DEM, and clocks

`lib/lbmdem/Domain.h` couples `LBMDOM` and `DEMDOM`. The historical problem was
that LBM and DEM updated iteration/time counters in different places, so
`iter`/`Time` could disagree, especially in LBMDEM.

Current rule: the enclosing solve loop is the single authority. Advance
`Time` and `iter` once per step, and publish the same values to CUDA through
`SetClock`. Do not add independent increments in the DEM or LBM substeps.

For CUDA, `SetClock<<<1,1>>>(...)` writes the authoritative clock after the
step. CPU paths use the same loop variables.

LBM-DEM tests are large. Build them and use a reduced temporary copy for smoke
tests. A reduced 41×21×21 LBM grid with `Tf=100` completed 6,250 iterations in
verification.

### 5.1 GPU FLBM porous-medium case

A standalone 2D D2Q9 GPU case for body-force-driven flow through a random
porous medium lives outside the checkout, so it is not a CMake target:

```text
/home/sergio/tflbm10_porous/tclbm10.cu   (plus Makefile, README, results)
```

Build and run it with the configured standalone compiler:

```bash
cd /home/sergio/tflbm10_porous
/home/sergio/mechsys/mechsyscc -o tclbm10 tclbm10.cu
./tclbm10 [nx] [ny] [porosity] [fx] [nsteps] [seed] [Nthread]
./tclbm10 512 256 0.90 1.0e-4 20000 12345 256
```

It flags `round((1-porosity)*nx*ny)` random nodes as solid and drives the fluid
with a constant `+x` body force. The GPU streaming pass resets `BForce` to zero,
so the constant force is written back on the device by the `Setup` callback
before every collision; it cannot be set only once on the host. The case prints
the raw second-moment velocity and a `K = nu*<ux>/fx` estimate, with the
`0.5*dt*fx` half-force variant used by the CPU path for comparison.

Measured on the RTX A6000 (double precision, 512×256×1, 20,000 steps,
`Nthread=256`): solver 5.37 s, total wall 5.64 s, realised porosity 0.900002,
steady `<ux>_fluid = 6.677e-4`, `<ux>_Darcy = 6.009e-4`, `K = 1.08` lattice
units (raw second moment 1.00). The raw GPU moment carries a constant
`-0.5*nu*phi` bias in `K`; the physical value uses
`u_physical = u_raw + 0.5*dt*fx*phi`. Two identical 128×64 runs produced
identical `.res` time series (no atomics in this path), and `K` changes by less
than 0.02% over `fx = 1e-5 … 1e-3`, confirming the Darcy regime.

A permeability sweep over `phi = 0.50 … 0.85` (512×512, 3 seeds, `fx=1e-3`,
20,000 steps) is in `/home/sergio/tflbm10_porous/permeability_sweep/`; every
case has its own folder with `.res`, `.h5`, and `.xmf`. Because the solid phase
is random site percolation, the fluid does not span the domain below
`phi_c ≈ 0.593`, where `K = 0`; the percolating values are `K(0.60)=0.0422`,
`K(0.65)=0.0709`, `K(0.70)=0.1163`, `K(0.75)=0.1889`, `K(0.80)=0.3121`, and
`K(0.85)=0.5441`. A 1024×1024 check changes `K` by less than 1% at
`phi=0.60/0.65`.

The half-force convention is independently validated by a plane-Poiseuille
case in `/home/sergio/tflbm10_porous/validation/`. With walls at `y=0.5` and
`y=H+0.5`, the exact permeability is `H^2/12`. The corrected `K` converges to
that value with approximately second-order accuracy: relative error `1.56e-2`
at `H=8`, `3.91e-3` at `H=16`, `9.77e-4` at `H=32`, `1.99e-4` at `H=64`, and
`1.63e-5` at `H=128`. At `H=32`, varying the force over a 16x range changes `K`
by only `2e-12` relative. Scale the force with `H^2` to keep `u_max ~ 0.01`; a
high-Mach run (`H=128`, `fx=1e-4`, `u_max ~ 1.23`) is not valid.


### 5.2 GPU FLBM advection-diffusion point source

A two-layer D2Q9 GPU case for a passive solute in the same random porous medium
is in `/home/sergio/tflbm10_porous/solute_source/`:

```text
solute_source.cu   Makefile   plot_solute.py   README.md
production_phi_0.70/
```

The domain uses `Array<double> nu = {1/6, D}` and
`FLBM::Domain(..., AdvectionDiffusion)`: layer 0 is the fluid (density and
velocity), layer 1 is the solute concentration. The medium is the same random
site-percolation medium as `permeability_sweep/phi_0.70/seed_12345`.

Two details are non-obvious:

1. **Solute advection velocity.** The AD collision forms
   `u_adv = u_raw + dt*Tau[1]/rho*BForce[1]`, while the physical fluid velocity
   is `u_phys = u_raw + 0.5*dt*fx/rho`. Set
   `BForce[1] = 0.5*fx/Tau[1]` on every fluid node in the device `Setup`
   callback. With `BForce[1] = 0`, the solute advects with the raw GPU moment,
   slower by `0.5*fx/rho`; at `fx=1e-2` this is comparable to the Darcy
   velocity.
2. **Source injection.** A source must add `dC` to both the layer-1
   populations and `Rho[1]`; the AD collision reads the previous macroscopic
   concentration before colliding. The production case uses a 9x9 Gaussian
   stencil (`sigma=2`) over fluid nodes, normalised so the total source is
   `S` per step. A single-node source is much sharper and can drive the
   entropic limiter negative.

The entropic limiter in `lib/flbm/lbm.cuh` (`cudaCollideAD`) was guarded in
this working tree. Without the guard, a negative population coinciding with
`F == Feq` makes `tau*F/(F-Feq)` equal `-inf`, then `alpha*(F-Feq)` becomes
`-inf*0 = NaN`; at 512x512 the zero-concentration far field eventually turns
NaN. The guard only applies the limiter ratio for `denom > 0` and `temp > 0`,
clamps negative layer-1 populations when the limiter cannot fix them, and
clamps any remaining negative layer-1 population after the two iterations.
The fluid layer is not post-clamped. Total solute mass is conserved after the
fix.

Production run: 512x512, `phi=0.70`, seed 12345, `fx=1e-2`, `D=0.02`,
`S=1e-4` per step, source on at `t=20,000`, 60,000 steps, `Nthread=256`.
Wall time 67.37 s (solve 66.67 s) on the RTX A6000. The source is placed on a
fast path near `(56,203)` with local `ux=0.0614`. At `t=60,000`,
`mass_total=4.0000`, `mass_fluid=3.4332`, `Cmax=1.089e-3`, centroid
`(219.1,198.3)`, and `var_x=9579`, `var_y=616`. The plume is channelised and
reaches the periodic right boundary between 55k and 60k; the 55k snapshot is
the last unwrapped field. `mass_fluid` (pore sum) is about 14% below
`mass_total` because the static bounce-back populations on solid nodes retain
part of the scalar; the total is the conserved quantity and matches
`S(t-t_on)`.

Outputs are `solute_source.res` (time, mass_fluid, mass_total, Cmax, centroid,
variance, source velocity) and `solute_XXXX.xmf/.h5` snapshots. HDF5 datasets
are `Density` (fluid), `Concentration` (solute), `Velocity_0`, `Velocity_1`,
and `Gamma`. `Domain::WriteXDMF` writes `Nx/Ny/Nz` inside the layer loop, so
`Nl=2` emits harmless "name already exists" HDF5 diagnostics; the case
silences them with `H5Eset_auto(H5E_DEFAULT, NULL, NULL)`.

Build and reproduce:

```bash
cd /home/sergio/tflbm10_porous/solute_source
make
make run          # ~67 s on RTX A6000
make plot
```


## 6. Accuracy verification procedure

### 6.1 Recreate a clean baseline

Use a separate git worktree so the modified working tree is untouched:

```bash
cd /home/sergio/mechsys
git worktree add /tmp/mechsys_base 94f64b8   # or the chosen baseline commit
ln -s lib /tmp/mechsys_base/mechsys
```

If the project Makefile in the worktree cannot resolve includes, compile the
baseline test manually with the same definitions and include paths as the
configured target. The working binary from the current tree is not a baseline.

### 6.2 Compare datasets, not files

HDF5 container metadata can change between two identical baseline runs, so
whole-file MD5 is not a valid correctness test. Open both files and compare
each dataset. XMF output can be checked separately (historically the XMF
checksums matched even when container metadata differed).

For CPU DEM force changes:

1. run both binaries with `Nproc=1`;
2. compare all final-state datasets;
3. expect no differences (`neq []`);
4. also compare printed energies;
5. then run `Nproc=8` for timing and note that only timing, not exactness, is
   expected there.

### 6.3 Keep output names separate

Tests write fixed `FileKey` names. Copy or rename outputs between runs so a
later run does not overwrite the baseline. Remove stale `_a`, `_b`, `_bf_`, and
`.res` files before a controlled comparison.

## 7. Known traps and hard-won lessons

These caused real failures during optimization work.

### Contact/cutoff

- With a cutoff bound of `2*ContactDispGlobal()`, a real contact was skipped:
  the build-time distance was inside the `2*alpha` list but outside the
  tightened force cutoff after motion. The observed ratio was about 3.4×.
  The current bound is `4×`, verified exact on the reduced CPU test.
- An invalid edge-edge lower bound reduced `cee` contacts and changed energy.
  It was removed.
- A planar face-bounds broadphase rejected valid periodic contacts and changed
  energy from `0.0397229` to `0.0416998`. Face bounds were reverted; only the
  conservative plane test remains.
- Sorting candidate lists by distance changes floating-point accumulation.
  Keep `L` in original order; use separate sorted arrays only if you do not
  change the force order.
- Replacing the full `Distance` with cached-normal/squared variants produced
  tiny but real final-state differences. Do not do that for contacts.

### CUDA

- `Nthread=512` fails for `test_cu_03` register usage.
- Fusing movement unconditionally broke `test_cu_04`; retain the separate
  translate/rotate fallback.
- Changing `MaxD` signature or clock plumbing broke `tlbmdem`; update
  `lib/lbmdem/Domain.h` when shared APIs change.
- A CUDA cache-initialization launch outside `#ifdef USE_CUDA` broke CPU
  `tlbmdem01/02`. Guard device-only code correctly.
- `test_cu_04.cu` uses CRLF; preserve it.

### Build and verification

- Building a baseline in a git worktree can fail on missing include paths;
  use the configured flags manually or create the `mechsys -> lib` symlink.
- The full `test_cu_03` run is long; a temporary shorter-`Tf` copy is useful,
  but restore the source afterward.
- The full `tlbmdem01` is large and can exceed a 300 s timeout; a reduced copy
  is the practical smoke test.
- `Array::Find()` returns `-1` on failure; assigning it to `size_t` is a bug.
  Check the signed result before using it as an index.
- `std::is_same`-based feature counters can be wrong because feature
  containers may hold pointers (`Array<Edge*>`, etc.), not values. Prefer
  explicit type traits or per-container logic.
- Warning policy: the only expected third-party warning is Blitz
  `tinymat.h(60): warning #191-D`. Keep MechSys sources warning-free.

### `tdem/ttt` triaxial test: restart does not work

`ttt <key> <ncores>` reads `<key>.inp` from the current directory.  With
`ptype = state` it restarts from `state.hdf5`.  **That path is broken in the
current binary** (`HEAD` `4a567cc`): the run aborts at the first step of the
hold, printing `Index = 49300`, producing NaN, and dying with
`terminate called after throwing an instance of 'Fatal*'` (exit 134).

Cause, from `ttt.cpp` around line 640:

- `Domain::Load` restores the interaction list, and that list still contains
  the six wall (platen) particles at indices `NP-6 .. NP-1`.
- `DelParticles(-2..-7)` then deletes those walls.
- `GenBoundingBox(-2, ...)` recreates six walls *at the same indices*, so the
  restored contact records now point at freshly built particles while carrying
  stale `Kn/Kt/Gn/Gt/Mu`, and the force comes out NaN.

It fails even for a state written by the same binary (`gen_b.hdf5`, `Time=100`,
walls tagged `-2..-7`) and for the untouched original input, so it is not a
consequence of editing the `.inp`.

Other `ttt` facts:

- `Domain::Save(key)` writes `Time`, `Xmin/Xmax`, ... and is the nominally
  restartable format.  Older `*_b.hdf5` files (e.g. the Oct-15
  `fine/ttt_b.hdf5`) lack `Time`/`Xmin` and also restart to NaN -- do not use
  them as restart states.
- `T0` is the total compression time: stage a ramps `0 -> T0/2`, stage b holds
  `T0/2 -> T0`.  `dat.tspan = T0/2 - dom.Time` is evaluated *before* stage a,
  so **never shorten `T0` below `2*Time` on a restart**: a negative `tspan`
  drives the stress servo backwards.
- `dtOut` is global to both stages.  For a fine cadence on the hold only, run
  from a fresh packing, delete the stage-a `_a_*` frames as they appear, and
  stop the process once the wanted number of `_b_*` frames exist.
- `WriteBF` writes `<key>_b_bf_NNNN.h5` (datasets `Branch`, `Normal`,
  `Position`, `Rolling`, `Tangential`, `ID1`, `ID2`; about 9.2 MB per frame at
  `NP=49306`).  There are no velocity fields in those files.
- The stored cadence `dtOut=0.05` is `0.90` of the contact vibration period
  `T = 2*pi*sqrt(mu_red/Kn) = 0.0557` (`m = rho*(4/3)*pi*R^3 = 1.5708`,
  `mu_red = 0.7854`, `Kn = 1e4`), so increments read at `dt <~ 0.05` are
  aliased.  Measured on the 49300-sphere packings: the local exponent of
  `<df^2>` rises to `1.8` and `K2` collapses below `T`; see
  `/home/sergio/ttt_compare/README.md`.

## 8. Current working-tree state

As of 2026-11-03:

- `HEAD` is `74407d4`, the merge of `maliaowave/mechsys` (`ad01691`) into
  `0d9d696`. The pre-merge tip is preserved as branch
  `backup/pre-mario-merge`; the fork is available as remote-tracking branch
  `mario/main` (`git fetch ~/mario/mechsys main:refs/remotes/mario/main`).
- The DEM/DEM-CUDA/FLBM optimization work described below is committed in the
  commits up to `0d9d696`.  Uncommitted in the tracked working tree: the CPU
  MPM optimization (`lib/mpm/Node.h`, `lib/mpm/shape_functions.h`,
  `lib/mpm/Domain.h`), the new `MPM_OPTIMIZATION_CHANGES.md`, and the memory
  documents (`AGENTS.md`, `MECHSYS_HANDBOOK.md`, `MEMORY.md`, `SKILL.md`,
  `.agents/skills/mechsys/SKILL.md`).  Do not commit or push them unless asked.
- Untracked build products: `tflbm/tclbm06..09`,
  `tlbmdem/tlbmdem0[12]`, `tlbmdem/tlbmdem_cu_0[1-4]`,
  `tlbmmpm/tlbmmpm01`, `tlbmmpm/tlbmmpm02`, `tmpm/tmpm0[45]`,
  `tmpm/tmpm_cu_0[1-5]`, `MPM_OPTIMIZATION_CHANGES.md`,
  `MECHSYS_CHANGES.txt` and the other MPM binaries from the build. They are
  regenerable and must not be deleted without asking.

## 9. MPM and tmpm

The Material Point Method library is `lib/mpm/`; its benchmark drivers are
`tmpm/tmpmNN.cpp` (CPU) and `tmpm/tmpm_cu_NN.cu` (CUDA). `A_WITH_TMPM=ON` in
this checkout, so all targets come from the root build system:

```bash
cd /home/sergio/mechsys/tmpm
make -j8 tmpm01 tmpm02 tmpm03 tmpm04 tmpm05
make -j8 tmpm_cu_01 tmpm_cu_02 tmpm_cu_03 tmpm_cu_04 tmpm_cu_05
```

`tlbmmpm01/02` include `<mechsys/lbmmpm/Domain.h>`, which pulls in the MPM
headers, so build them after changing `lib/mpm/*`:

```bash
cd /home/sergio/mechsys/tlbmmpm && make -j8 tlbmmpm01 tlbmmpm02
```

### 9.1 Library entry points

| File | Responsibility |
|---|---|
| `lib/mpm/Domain.h` | background grid, particle-to-node/node-to-particle, `Solve`, body mesh |
| `lib/mpm/Particle.h` | material-point state, constitutive models, `CalcProps` |
| `lib/mpm/shape_functions.h` | GIMP, quadratic and cubic B-splines, support radius and node range |
| `lib/mpm/eos.h` | Tait/Cole (`EOS_TAIT`) and linear (`EOS_LINEAR`) acoustic EOS |
| `lib/mpm/fluid_filter.h` | vNRL artificial viscosity for weakly-compressible water |
| `lib/mpm/mpm.cuh` | CUDA kernels used by `tmpm_cu_*` |

- `dom.ShapeFunction` defaults to `SHAPE_GIMP`; the newer cases set
  `MPM::SHAPE_BSPLINE_QUADRATIC`, and `SHAPE_BSPLINE_CUBIC` also exists.
- `Particle::Constitutive == 4` is the weakly-compressible fluid (`Mu`, `Cs`,
  `EOSType`, `EOSexp`, `AVLin`, `AVQuad`, `FluidDx`).
- `Constitutive == 5` is mu(I) rheology. Soil models are 1 = Drucker-Prager,
  2 = Bui Drucker-Prager, 3 = hypoplastic, driven by `K`, `G`, `Friction`,
  `Dilation`, `Cohesion`.
- `dom.Gn` is the nodal dissipative (hourglass-damping) constant, not gravity.
  Body acceleration is set in the case `Setup` callback.

### 9.2 Checked-in benchmark sizes

| Case | Physics | Particles | Grid nodes | Time steps |
|---|---|---|---|---|
| `tmpm01`/`tmpm_cu_01` | 3-D column, fixed base, end force | 5,000 | 125,000 | ~86,000 |
| `tmpm02`/`tmpm_cu_02` | dry granular column, Drucker-Prager | ~26,000 | ~106,000 | ~580,000 |
| `tmpm03`/`tmpm_cu_03` | closed-box water column, Tait EOS + vNRL | ~135,000 | ~358,000 | ~259,000 |
| `tmpm04`/`tmpm_cu_04` | periodic square-duct flow vs analytic | 216,000 | ~39,000 | ~42,400 |
| `tmpm05`/`tmpm_cu_05` | granular column, mu(I) rheology | ~26,000 | ~106,000 | ~500,000 |

All of these are long. For a smoke test, copy the source outside the checkout,
reduce `ndiv` / `num_particles_x` / `particles_per_diameter` and, except for
`tmpm04`, `Tf`; compile the copy with
`/home/sergio/mechsys/mechsyscc copy.cpp` (or `.cu`). `tmpm04` needs the
physical `tf = 200 s` to reach steady state (viscous diffusion time is about
25 s at `h = 5 cm`), so reduce only `particles_per_diameter` there.

Measured on 2026-11-03 (RTX A6000, 128-core host; wall time includes startup
and output):

| Smoke reduction | CPU wall | GPU wall |
|---|---|---|
| `tmpm01`: `ndiv 50->10`, `Tf 5->0.02` | 0.04 s | 1.49 s |
| `tmpm02`: `num_particles_x 40->10`, `Tf 1.5->0.02` | 0.87 s | 3.72 s |
| `tmpm03`: `particles_per_diameter 70->20`, `Tf 5->0.02` | 0.39 s | 7.34 s |
| `tmpm04`: `particles_per_diameter 60->10`, `tf 200->1` | 0.07 s | 1.78 s |
| `tmpm05`: `num_particles_x 40->10`, `Tf 1.5->0.01` | 0.43 s | 2.35 s |

All ten exit 0 and write 36-100 HDF5 frames with finite `Position`,
`PVelocity`, `Pressure`, `Sigma`, and `Strain`.

`tmpm04` at `particles_per_diameter = 20` (8,000 particles, 14,200 steps,
`tf = 200`) ran in 22.9 s CPU wall. `Vx_max` converged to `8.57e-4 m/s` against
the case's own analytic `7.37e-4 m/s` (`Rerror = 16.3%`); that resolution has
only five cells across the pipe radius, so this is a convergence check, not an
accuracy claim. The checked-in 60-ppd case (216,000 particles, 42,400 steps,
96 OpenMP threads, `tf = 200`) completed in 1932 s (32 min 12 s) wall before
the CPU optimization and 1676 s (27 min 56 s) after it (Section 9.3), with
`Vx_max = 7.68731e-4 m/s`, i.e. `Rerror = 4.35%` from `t = 100 s` to the end.
That is the checked-in accuracy of this benchmark; it does not converge
further with time.

At the tested reduced resolution (20 ppd, 3,160 particles, 294 steps),
`tmpm_cu_03` and `tmpm03` produced bitwise-identical HDF5 datasets before the
CUDA optimization.  After it the GPU atomics can reorder sums; do not
generalize bitwise identity to other GPU MPM cases without checking.

`tmpm_cu_02.cu` was repaired after the fork merge: its `dtOut` declaration had
been commented out while the variable was still printed, and `Solve` used
`Tf/100` instead. It now mirrors `tmpm02.cpp`
(`double dtOut = std::min(Tf/100.0,5000.0*dt)`, passed to `Solve`).

**`tmpm_cu_02` is not kinetically equilibrated at its checked-in `Tf = 1.5 s`.**
A full run (25,280 points, 50 frames, `Nthread=256`) took 30 min 37 s in
`Solve` (1837 s wall). The deposit shape stabilizes by about 0.6 s
(`x in [-0.126, 0.314] m`, `z` the same, `y_max = 0.075 m`) and changes by
less than a millimetre to 1.47 s. However, a thin basal layer keeps sliding:
at 1.47 s `vmax = 1.009 m/s`, `vrms = 0.104 m/s`, median speed `0.017 m/s`.
The fixed ground nodes plus `Gn = 0.001` do not bring the material to rest.
Do not present the 1.5 s state as static equilibrium; it is a steady deposit
shape with persistent basal flow. A true static state needs basal friction,
stronger damping, or a much longer run.

### 9.3 CPU performance changes

The grid transfer in `ParticleToNode`/`NodeToParticle` was optimized without
changing the arithmetic; full design and measurements are in
[`MPM_OPTIMIZATION_CHANGES.md`](MPM_OPTIMIZATION_CHANGES.md).

- `Node` is `alignas(64)` and keeps `Mass/Mn/Fn/valid/lck` in one cache line.
- `Node::Reset` no longer clears `Vn` or `Sn` (both are overwritten/unused).
- `VNodes` is built during the scatter (first touch per thread), not by a scan
  of all `Nnodes`.
- Per-particle 1-D shape weights are cached and reused for the 27/64 nodes.
- The NaN guard is outside the node loop; single-thread runs skip the locks.
- The node array is constructed in parallel so pages are NUMA-local.

Measured speedups (`Solve` wall time unless stated; 96 threads on the
128-thread host):

| Case | Before | After | Speedup |
|---|---:|---:|---:|
| `tmpm03` reduced, 20 ppd, 998 steps | 2.03 s | 0.89 s | 2.29x |
| `tmpm03` reduced, 50 ppd, 150 steps | 0.79 s | 0.40 s | 1.98x |
| `tmpm01` full, 86,391 steps, 100 frames (wall) | 92.42 s | 33.77 s | 2.74x |
| `tmpm04` full, 42,437 steps, 100 frames (wall) | 1932.08 s | 1676.27 s | 1.15x |

`Nproc=1` output is dataset-bit-identical for GIMP (`tmpm01`),
Drucker-Prager (`tmpm02`), Tait water (`tmpm03`) and mu(I) (`tmpm05`).  At
`Nproc>1` only the floating-point reduction order can change; full `tmpm01`
at 96 threads differs only in `PVelocity` by `2.8e-14` (7.6e-11 relative),
and `tmpm04` frames differ by at most `~1e-14` relative in
`PVelocity`/`Sigma`.

`tmpm04` is the dense limit: 216,000 particles on 38,880 nodes (~5.5 per
cell), so each node is updated ~150 times per step and the scatter is bound
by cache-line ownership.  On a short 5-frame `tmpm04` run, 8/32/96 threads
take 98.1/60.6/71.3 s after the optimization; for such dense grids fewer
threads than `0.75*MAXTHREADS` can be faster.  An eight-group partial-grid
prototype reached 60.8 s at 96 threads but regressed sparse cases, so it was
not kept.

### 9.4 CUDA MPM performance changes

The GPU path was optimized after the CPU work; full design, phase profile and
numbers are in
[`CUDA_MPM_OPTIMIZATION_CHANGES.md`](CUDA_MPM_OPTIMIZATION_CHANGES.md).
Files changed: `lib/mpm/mpm.cuh`, the `USE_CUDA` parts of
`lib/mpm/Domain.h`, and `lib/mpm/Particle.h`.

- Hot node fields (`Mass`, `Mn`, `Fn`, `Vn`, `valid`) are separate device
  arrays (SoA), so neighbouring nodes read/write consecutive addresses.  The
  cold `NodeCU` array still holds `Vnf` and the fixed-velocity flags.
- The scatter builds the compact list of touched nodes itself
  (`atomicExch` + `atomicAdd`).  The old per-step `thrust::fill`,
  `thrust::copy_if`, `thrust::distance` and count `cudaMemcpy` are gone, along
  with `bAllIndexes`/`bValidMask`.  The previous list is reset first; the
  device count is ping-ponged.
- Reset and solve are launched over `Nnodes` blocks and exit early on the
  device count, so the host never needs the valid-node count.
- `ParticleCU` now stores 64 nodes per particle (`4^3`, enough for cubic
  B-spline and for quadratic/GIMP at the usual particle size), down from 125;
  `sizeof(ParticleCU)` drops from 6024 B to 3584 B.
- A final `cudaDeviceSynchronize()` makes `Solve` timing cover queued work.

Reduced four-frame benchmarks (`Nthread=256`, RTX A6000) improved
`tmpm_cu_01/03/04/05` by 4/8/9/10%; `tmpm_cu_02` regressed by about 1%
(short steps, frequent particle deletion/re-upload).  These are `Solve`
times; see the report for wall times and the phase split.

The GPU scatter uses floating-point atomics, so bitwise identity with the
pre-change binary is not guaranteed.  On the reduced cases the final-frame
differences were at most `~2e-5` absolute in `Pressure`/`Sigma` (about
`1e-8` relative) and `~1e-12` in `PVelocity`; `tmpm_cu_02` was bitwise
identical.  At the previously tested reduced resolution, `tmpm_cu_03` and
`tmpm03` were bitwise identical before this CUDA optimization; do not expect
that after it.

Caveat: the 64-entry arrays cover up to four nodes per dimension.  A GIMP
particle roughly larger than 1.5 cells would need five nodes per dimension
(125 entries); use a smaller particle size or widen the arrays in
`Particle.h` for that configuration.

## 10. Documentation map

| Document | Use |
|---|---|
| `AGENTS.md` | automatic, concise agent rules |
| `MEMORY.md` | distilled engineering experience and a short verification checklist |
| `SKILL.md` | visible, downloadable copy of the project skill |
| `.agents/skills/mechsys/SKILL.md` | project-local DSH skill; discovered automatically when this checkout is the project root |
| `MECHSYS_HANDBOOK.md` | this durable reference |
| `CPU_DEM_OPTIMIZATION_CHANGES.md` | CPU Voronoi changes, exactness, benchmarks |
| `CUDA_OPTIMIZATION_CHANGES.md` | CUDA kernels, dispatch, clocks, accuracy, performance |
| `REDUCED_TEST_CU_03_PERFORMANCE.md` | reduced CUDA thread/Verlet sweeps |
| `MPM_OPTIMIZATION_CHANGES.md` | CPU MPM changes, exactness, benchmarks, dense-grid limit |
| `CUDA_MPM_OPTIMIZATION_CHANGES.md` | CUDA MPM SoA/valid-list changes, phase profile, GPU benchmarks and 64-node limit |
| `benchmarks/cpu_voronoi/` | reproducible CPU Voronoi inputs |

## 11. Maintaining this knowledge

Maintenance is mandatory, not opportunistic. At the end of every user prompt
in this repository, evaluate whether the prompt produced a durable fact. If it
did, update the owning document before replying:

1. **Invariant or short operational rule** -> `AGENTS.md`, which is always
   loaded. Keep it terse; this file is the last place for a long explanation.
2. **Architecture, command, exactness procedure, or trap** -> this handbook.
3. **Distilled lesson or verification pattern** -> `MEMORY.md`.
4. **Design or benchmark result** -> the relevant
   `CPU_DEM_OPTIMIZATION_CHANGES.md`,
   `CUDA_OPTIMIZATION_CHANGES.md`, or
   `REDUCED_TEST_CU_03_PERFORMANCE.md`.
5. **Task-triggered guidance** -> both `SKILL.md` and
   `.agents/skills/mechsys/SKILL.md`, kept in sync.
6. **Known trap** -> the “Known traps” section here, plus a one-line rule in
   `AGENTS.md` if it is severe enough to affect every session.

Rules for the update:

- Prefer an existing section over a duplicate paragraph.
- Record durable facts only: architecture, verified invariants, working
  commands, measured benchmarks, demonstrated failure modes. Do not promote
  transient `/tmp` paths, one-off logs, guesses, or unverified numbers.
- For a benchmark, record the exact input, thread count, time window, machine,
  and whether the value is wall or solve time.
- If nothing durable was learned, make no edit; do not manufacture activity.
- Run `git diff --check`; do not commit or push unless explicitly asked.
- Say in the final response which memory file changed, or that no update was
  needed.
