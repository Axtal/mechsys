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

### 3.2.1 Integrator and energy tracks (position Verlet / velocity Verlet)

`DEM::Domain` carries two runtime switches; both are plain members that must be
set **before** `Solve`, and neither is part of a saved state:

| member | values | effect |
| --- | --- | --- |
| `UseVelocityVerlet` | `false` (default) / `true` | position Verlet / velocity Verlet |
| `EnergyBalance` (set with `SetEnergyBalance(bool)`) | `false` (default) / `true` | baseline contact kernels / kernels that also record `Epot`, `dEvis`, `dEfric`, `Fndpot`, `Ftdpot`, `Nc`, `NormalOverlap` |

`Solve` splits one iteration into three stages chosen once by
`SelectStages()`: a pre-force stage (`NULL` for position Verlet; half kick,
drift, rotational update and a list-validity check on the post-drift geometry
for velocity Verlet), the force pass, and a post-force stage (velocity
finalisation for velocity Verlet, translate/rotate plus the list check for
position Verlet). The loop body itself tests nothing:

```c++
if (pPreForce != NULL) (this->*pPreForce)();
(this->*pCalcForce)();
(this->*pPostForce)();
```

The three calls are indirect *host* calls, one each per iteration, measured at
**5.3 ns/iteration** against 1.8 ns for direct calls (a throwaway `-O3`
program: two loops of 2e8 iterations, one calling a member function directly and
one through a member pointer, so it measures call overhead only; it is not worth
keeping in the tree). The cheapest step
in the suite is `tdem/test_cu_01` at 73 us, so the overhead is 0.007 % there and
below 1e-4 % on coupled runs. `LBMDEM::Domain::Solve` drives the same stages
(`DEMDOM.SelectStages()` before its loop) and exposes
`SetDEMIntegrator(bool)` / `SetDEMEnergyBalance(bool)`; its CPU coupling loop
remains position-Verlet only and refuses velocity Verlet with a `Fatal`.

The default contact kernels are the pre-existing ones **verbatim**: an
unmodified position-Verlet run must stay bit-identical to the pre-change
library. Verified on this tree: `tdem/test_cu_04` reports
`0 / 0.00914244 / 0.00244635` with and without the new options, and the coupled
Stokes case (`stokes_main.cu 1 256 0 0`) writes a `cur_stokes_force.res`
byte-identical to the pre-change build. The energy-track kernels therefore only
add recording; they do not change the contact law.

**Trap:** in `tdem/test_cu_04` the phase-2 wall prescription (`MoveTopPar`) is
integrator-specific. Position Verlet keeps the original
`Translate(Verts,Vertso,Par[ip],DPar[ip],trans)` call, because its step recovers
the prescribed displacement from the translated position; velocity Verlet
instead prescribes `v.z = Dz/dt` and `xb.z = x.z - Dz`, because the drift reads
`v_half*dt`. Applying both mechanisms moves the wall twice. The choice is the
`vv` field of `UdCu`, not a per-step test inside the kernel.

**Position Verlet vs velocity Verlet, measured.** `tdem/test_cu_04` with a
local reduced copy of `tpercu.inp` (sphere packing, `R=0.05`, `dt=2e-5`, 34
parameters, named `cu04.inp` because the driver reads `<key>.inp`; `tpercu.inp`
itself is the large Voronoi case) as `test_cu_04 cu04 1 4.0
0.0 <vv> 0` at `dt = 2e-5, 1e-5, 5e-6, 2.5e-6` (the input's `dt` line, 501+1001
steps at the stock value). Reported observable: the kinematic energy, i.e. the
only energy `CalcEnergy` fills here (`Potential energy` prints 0 because gravity
enters as `Ff`, not as a potential).

| phase | PV, dt=2e-5 | PV, dt/8 | PV dt->0 | VV, dt=2e-5 | VV, dt/8 | VV dt->0 |
| --- | --- | --- | --- | --- | --- | --- |
| end of compression `T0` | 0.00914244 | 0.00904751 | 0.00904955 (+1.03 % at stock dt) | 0.00905745 | 0.00905192 | 0.00904988 (+0.08 %) |
| end of `T0+Tf` | 0.00244635 | 0.00195421 | 0.00194917 (+25.5 % at stock dt) | 0.0117339 | 0.0117279 | 0.0117313 (dt-flat) |

PV converges at an observed order of 2.2-2.6 and needs about `dt/4` to be within
0.1 % of the limit. VV is already dt-flat at the stock step (spread <=0.15 %, the
noise floor of a chaotic packing) and its `T0` limit matches PV's to 0.004 %.
The 6x difference in the second row is *not* an integrator error: the residual
kinetic energy of a static jammed packing is a chaotic rattle and the two schemes
settle into different microscopic attractors (the damping is read from the
half-step velocity in VV and from the step-start velocity in PV). The
macroscopic packing agrees at every dt: height `1.004566 vs 1.004565` at `T0`
and `0.894047 vs 0.894053` at the end (2e-6 relative), coordination
`Z=2.250 vs 2.228`, mean normal overlap `4.33e-4 vs 4.38e-4`, max overlap
`3.49e-3` in both.

Cross-checks with an answer: `tdem/test_cu_01` (sliding block, 1e6 steps) gives
`Observed L = 51.0019` against `Theoretical L = 51.0204` (-0.036 %) for **both**
integrators, final kinematic energy `0.0200495` vs `0.0200358`; `tdem/test_cu_03`
(1000-particle sheared packing, 1e5 steps) gives a steady-state total kinematic
energy of `32105.9` (PV) vs `32154.3` (VV), 0.15 % apart, at 799 vs 804 steps
between list updates.

**The CPU path has the same switch.** `SelectStages` chooses the stages
independently of `USE_CUDA`: on the CPU, `CalcForcesVV` installs
`v_half`/`w_half`, accumulates directly into `Particle::F/T` under per-particle
locks and keeps the `CalcForceSphere` fast path; `PreForceVV` does
`UpdateVelocityHalf`, `UpdateAngularVelocityHalf`, `TranslateVelVerlet` and
`RotateVelVerlet` plus the post-drift list check; `PostForceVV` does
`UpdateVelocityFull`, `UpdateAngularVelocityFull` and the post-move check.
`test_01`, `test_02` and `test_periodic` now take `[integrator] [energy]`
(`test_periodic`: `<key> [Nproc] [dthoot] [integrator] [energy]`) beside
`test_saveload`. Measured on CPU:

- `test_01 0 0` vs `test_01 1 0` (cube sliding on a plane, 1e6 steps,
  `OMP_NUM_THREADS` default): `Observed L = 51.0012` vs `51.0013` against
  `Theoretical L = 51.0204` (-0.038 % both), solve 2.674 vs 2.602 s.
- `test_periodic tper 4 4.0 <vv> <eb>` on a reduced 6^3 sphere packing
  (`T0=0.01`, `Tf=0.02`, `dt=5e-5`, hertz): `PV 0.120176/1.05864`,
  `VV 0.120678/1.13789`. With the energy tracks on, the CPU VV run is
  **bit-identical** to plain VV (unlike the GPU, where the extra `atomicAdd`s
  reorder the force sums): `Evis=2.666e-4` and `Efric=2.329e-6` after `T0`,
  `Evis=7.463e-4` and `Efric=1.328e-4` after `Tf`, 87146 sphere contacts.
- `test_02 <vv> <eb>` (cube rolling on a turntable, `Tf=30`, `dt=1e-5`,
  3e6 steps) prints the analytic orbit `xt`/`yt`; the deviation of the integrated
  position from it is the cleanest CPU accuracy probe: final
  `1.175e-2` (PV) vs `1.795e-3` (VV), mean `5.82e-3` vs `9.41e-4`, i.e. **VV is
  6.5x closer to the analytic orbit at the same step**. Total energy
  `269131` vs `269136` (2e-5 relative).
- CPU cost, same packing, per step: `Nproc=1` PV 14.758 ms vs VV 14.735 ms
  (equal); `Nproc=4` PV 6.956 ms vs VV 4.322 ms (**VV 38 % faster**). The
  reason is structural: the legacy CPU PV force pass zeroes and reduces
  per-thread `MTD[i].F/T` arrays (`O(Nproc*N)` traffic per step plus a
  reduction), while the VV force pass accumulates in place. The two paths
  therefore have different summation orders; neither is bitwise reproducible
  at `Nproc>1`.
- The CPU coupling loop in `LBMDEM::Domain` still refuses velocity Verlet
  (`Fatal: LBM-DEM velocity Verlet requires the CUDA path`), so
  `tlbmdem01`/`tlbmdem02` are position-Verlet only.

**Cost.** Within noise on real workloads, because VV trades the fused move for a
half-kick and a pre-force list check:

| workload | steps | PV | VV |
| --- | --- | --- | --- |
| `test_cu_03`, 1000 particles | 1e5 | 1.8492 ms/step | 1.8606 ms/step (+0.6 %) |
| `test_cu_04` phase 2, `dt=2.5e-6` | 8000 | 99.11 us/step | 98.23 us/step (-0.9 %) |
| coupled Stokes 64^3 | 1251 | 3.349 ms/step | 3.325 ms/step (-0.7 %) |
| `test_cu_01`, single body | 1e6 | 72.66 s | 75.36 s (+3.7 %) |

The single-body case is latency bound (73 us/step), so VV's two extra kernel
launches cost 2.7 us there; do not quote that as a throughput loss.

**Test switches.** Every DEM and LBM-DEM test takes the same two trailing
arguments, `<integrator> <energy>` (`0`/`1` each; both default to `0`). No test
prints them itself: `Domain::Solve` reports both in its banner, immediately
after the Verlet distance, so the library -- not each driver -- says what ran:

```
  Verlet distance                           =  0.05
  Velocity Verlet integrator                =  on
  Energy balance tracking                   =  on
```

The coupled `LBMDEM::Domain::Solve` prints the same two lines for its DEM
sub-domain. Keep the reporting in `Solve`: a driver that echoes the flags is
duplicating the library and goes stale when the defaults change.

```bash
cd /home/sergio/mechsys/tdem
./test_01 1 0                   # CPU: <integrator> <energy>
./test_02 1 0
./test_periodic tper 4 4.0 1 0  # <key> <Nproc> <dthoot> <integrator> <energy>
./test_saveload 1 100 100 1 1   # <Nproc> <steps before save> <steps after load>
./test_cu_01 1 10.0 1 0         # <Nproc> <Tf> <integrator> <energy>
./test_cu_02 1 10.0 1 0
./test_cu_03 1 0.05 128 10.0 1 1  # <Nproc> <Alpha> <Nthread> <Tf>
./test_cu_04 tpercu 1 4.0 0.0 1 1  # <key> <Nproc> <dthoot> <qin>; reads <key>.inp
cd ../tlbmdem
./tlbmdem01 1 1 0                  # <Nproc>; the CPU coupling refuses VV
./tlbmdem_cu_01 1 1 1
./tlbmdem_cu_02 1 256 1 1          # <Nproc> <nthread>
./tlbmdem_cu_03 1 256 1 1
./tlbmdem_cu_04 1 256 1 1
./tlbmdem02 1 256 1 1
```

**The coupled solver has the switch on the CUDA path only.** `LBMDEM::Domain`
forwards it to its DEM sub-domain (`SetDEMIntegrator`/`SetDEMEnergyBalance`) and
`LBMDEM::Domain::Solve` drives the same three stages, putting the fluid
imprint at the top of the step; the CPU coupling loop does not drive the
stages and refuses velocity Verlet
(`Fatal: LBM-DEM velocity Verlet requires the CUDA path`). All six tlbmdem
programs accept the two trailing arguments. Verified with reduced copies (grid
and `Tf` cut, same physics parameters, `RenderVideo` off), PV vs VV:

| case | reduced setup | result |
| --- | --- | --- |
| `tlbmdem_cu_01` (sphere in a channel, drag) | `Tf=200`, 301x61x61 | `_force.res` **identical in every column over 200 rows**: `Re=33.3418`, `CD=2.24298`, `CDsphere=2.00553` |
| `tlbmdem_cu_02` (Stokes terminal velocity) | 100^3, `Tf=100` | `\|Flbm\| 2.42949` vs `2.44483`, `\|v\| 0.0429552` vs `0.0429540` (**2.5e-5 relative** in the terminal velocity) |
| `tlbmdem_cu_03` (settling sphere) | `Tf=100` | `_force.res` **identical**, `Ft=0.162860` in both |
| `tlbmdem_cu_04` (gravity settling) | `Tf=2e3` | both integrators run to completion, 18.65 s wall each |
| `tlbmdem01`/`tlbmdem02` (CPU coupling) | — | `Fatal: LBM-DEM velocity Verlet requires the CUDA path` |

**HDF5 `name already exists` on multi-layer output (fixed).** `FLBM::Domain`
has three HDF5 writers -- `WriteXDMF`, `WriteXDMF_DEM` and the restart `Save`
-- and all three created the per-file grid datasets (`Nx`/`Ny`/`Nz`, and
`Ncell` in `Save`) *inside* their `for (j=0;j<Nl;j++)` layer loop. With
`Nl>1` every extra layer re-created them, so HDF5 dumped a `H5Dcreate2 ...
name already exists` trace per frame per layer and kept going: the files were
complete, the errors were failed duplicate creates. It shows up in any two-layer
case (`tlbmdem_cu_04` has `nu(2)`, layers 0/1 with different densities;
`tclbm03`; the `tlbmdem_cu_02` layers), with either integrator, and is not a
velocity-Verlet problem. The three writers now guard those datasets with
`if (j==0)`, so one set is written per file. Verified: the reduced two-layer
coupled case went from 600 error lines to 0 over 200 frames, the real
`tlbmdem_cu_04 24 256 1 1` produces 0, and `tclbm03` (which exercises `Save` as
well) writes `F_0/F_1/Ncell/Nx/Ny/Nz/Solid_0/...` with 0 errors.

Output volume to expect: at 100^3 a `_lbm_%04d.h5` frame is ~34 MB, so
`tlbmdem_cu_04`'s 400 frames are >10 GB of LBM frames alone, beside the `_dem_`
and `_dem_bf_` sets (150 s of running produced 11 GB).

The switch is a plain domain member in every case (`dom.UseVelocityVerlet`,
`dom.SetEnergyBalance(eb)`; `dom.SetDEMIntegrator(vv)`/
`dom.SetDEMEnergyBalance(eb)` for a coupled domain), so any other driver sets
the same two fields before `Solve`. `EnergyBalance` only adds recording to the
contact kernels; a position-Verlet run with it off must stay bit-identical to
the pre-change library, and a CPU run with it on is bit-identical to the same
run with it off because the host path has no atomics (the GPU path does, so its
force sums reorder and the trajectories drift apart).

`tdem/test_saveload <Nproc> <K> <M> <integrator> <energy>` runs one
combination per invocation (three `Solve` calls: the pre-save run, the
continuation, and the reloaded copy), so check all four with
`for k in 0 1; do for m in 0 1; do ./test_saveload 1 100 100 $k $m; done; done`
-- each must print `RESULT: PASS`; with the energy tracks on it also compares every contact's
`Epot`/`dEvis`/`dEfric`/`Fndpot`/`Ftdpot`/`Nc`/`NormalOverlap` bitwise. Note
that `PairtoCInt` is an `unordered_map`, so its iteration order (and therefore
the last bit of a domain total summed from it) can differ after a restart; the
test sums the totals in sorted-key order for the comparison.

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
and node attributes `Tag`, `Velocity`, `Pressure`, `Sigma`, `Strain`.  Use the
project **fluid colormap, `turbo`** (`colorTableName = "turbo"`) for every
fluid, solid, and animation plot unless the user explicitly requests another
map.  A verified point-cloud render is:

```python
OpenDatabase("granular_0019.xmf")
DefineScalarExpression("speed", "magnitude(Velocity)")
AddPlot("Pseudocolor", "speed")
pc = PseudocolorAttributes()
pc.colorTableName = "turbo"         # project fluid colormap
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
pc.colorTableName = "turbo"         # project fluid colormap, blue low -> red high
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

**Note (2026-11-17):** `/home/sergio/tflbm10_porous/` no longer exists on this
machine (nor `permeability_sweep/`, `validation/`, `solute_source/`); the
numbers below are kept as measurements, but there is no source tree to rerun.
The reproducible plane-Poiseuille validation now lives in
`/home/sergio/tflbm_poiseuille/` (section 5.1.1).

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
units (raw second moment 1.00). Two identical 128×64 runs produced identical
`.res` time series (no atomics in this path), and `K` changes by less than
0.02% over `fx = 1e-5 … 1e-3`, confirming the Darcy regime.

The velocity convention for those numbers: **both** active paths report the raw
second moment `u = sum_k f_k c_k / rho` (`StreamSC` on the CPU, the fused
stream kernel on the GPU). The `+0.5*dt*BForce` term only exists in
`StreamEFS_SRT`, which `Solve` does not call, so the "CPU path adds the
half-force velocity" statement of earlier revisions was wrong (verified
2026-11-17 against the CPU and GPU Poiseuille profiles, which agree to
`1.4e-16`). Any half-force correction must therefore be *measured*, not assumed:
in a forced channel the raw profile is the analytical parabola minus a constant
slip, see 5.1.1. The porous `K` values above were derived with the old
"`u_raw + 0.5*dt*fx*phi`" reading and are about `0.5*nu*phi` high in `K`; treat
them as relative/ordering data, not as absolute validated permeabilities.

A permeability sweep over `phi = 0.50 … 0.85` (512×512, 3 seeds, `fx=1e-3`,
20,000 steps) is in `/home/sergio/tflbm10_porous/permeability_sweep/`; every
case has its own folder with `.res`, `.h5`, and `.xmf`. Because the solid phase
is random site percolation, the fluid does not span the domain below
`phi_c ≈ 0.593`, where `K = 0`; the percolating values are `K(0.60)=0.0422`,
`K(0.65)=0.0709`, `K(0.70)=0.1163`, `K(0.75)=0.1889`, `K(0.80)=0.3121`, and
`K(0.85)=0.5441`. A 1024×1024 check changes `K` by less than 1% at
`phi=0.60/0.65`.

### 5.1.1 Reproducible plane-Poiseuille validation

`/home/sergio/tflbm_poiseuille/` (`poiseuille.cu`, `analyze.py`, `summary.py`,
`plot.py`, `profile_table.py`, `README.md`, `results.md`, four PNGs). Same
recipe as the porous case: D2Q9, `nx x (H+2*nsolid)`, `IsSolid` walls, constant
`+x` force rewritten into `BForce` (device `pBForce`) by `Setup` every step,
`Dom.Sc = 0` so the Smagorinsky term does not change `tau`, `fx = 8 nu u_max/H^2`
with `u_max ~ 0.01`, `Tf = 40 H^2` (the viscous time is `H^2/nu`).

```bash
cd /home/sergio/tflbm_poiseuille
/home/sergio/mechsys/mechsyscc poiseuille.cu -o poiseuille_gpu
./poiseuille_gpu 32 8 0.01 0.16666666666666666 256 40 1
python3 profile_table.py data/H032_nu0.166667
```

Walls are at `y = 0.5` and `y = H+0.5` (half-way bounce-back, `H` = number of
fluid rows, not `H+1`; `nsolid=1/3/6` gives bit-identical profiles), and the
exact permeability is `H^2/12`. Measured on the RTX A6000 (`Nthread=256`,
double, `nx=8`), `nu=1/6` (`tau=1`):

| `H` | steps | L2 rel. | max rel. | `nu_fit` | `K_raw` | `K/H^2*12 - 1` | wall |
|--:|--:|--:|--:|--:|--:|--:|--:|
| 8 | 2,560 | 7.131e-3 | 5.291e-3 | 0.1666667 | 5.333333 | -1.3e-14 | 0.31 s |
| 16 | 10,240 | 1.783e-3 | 1.307e-3 | 0.1666667 | 21.333333 | -5.0e-15 | 0.45 s |
| 32 | 40,960 | 4.457e-4 | 3.258e-4 | 0.1666667 | 85.333333 | +4.5e-13 | 1.80 s |
| 64 | 163,840 | 1.114e-4 | 8.140e-5 | 0.1666667 | 341.333333 | +3.1e-12 | 6.3 s |
| 128 | 655,360 | 2.786e-5 | 2.035e-5 | 0.1666667 | 1365.333333 | +1.4e-11 | 24.6 s |

The steady profile is **exactly parabolic in the cell index** (quadratic fit
residual `1e-13`) with the input viscosity and a *uniform* offset, i.e. a
spurious wall slip:

```text
u_raw(y) = a/(2 nu) (y-y_w)(y_t-y) - c(tau) a ,   a = fx/rho,
c(tau)   = 1/(8 nu) - 6 nu + 1/2          (dx = dt = 1, raw second moment)
```

`c` fits 23 independent runs to four digits, is independent of `H` (checked
`H=32` vs `H=64`), linear in `fx`, independent of the initial condition (rest,
analytical parabola, plug flow give the same fixed point to `1e-17`), and
identical on CPU and GPU (`1.4e-16`, round-off from `NonEq/tau` vs
`NonEq*itau`); repeated GPU runs are bit-identical. It vanishes at
`nu = 0.1919` (`tau = 1.0757`), which is why the pointwise error is smallest
near `tau = 1.1`. Consequences:

- the pointwise error is *only* that constant, so it is second order in the
  grid: `u_slip/u_max = (1-48 nu^2+4 nu)/H^2`;
- at `tau=1` the cell-averaged bulk velocity gives `K_raw = H^2/12` to `1e-11`
  with **no** correction (the slip happens to cancel the `dx^2/12` sampling
  bias of the cell average), and `|dK/K| <= 1.1e-3` for `0.56 <= tau <= 1.1`
  at `H=32`, four times smaller at `H=64`;
- the naive half-force velocity `u_raw + 0.5 a dt` is *not* the right
  correction: it makes the `tau=1` bulk error `+0.098%` instead of `0` and only
  halves the pointwise error at `tau=0.65`.

Controls: doubling `fx` doubles the profile exactly; `nx=32` instead of `8` is
identical; `nsolid=1/3/6` identical; CPU `StreamSC` and the GPU kernel agree to
round-off. Scale the force with `H^2` to keep `u_max ~ 0.01`; high-Mach runs
are invalid (the old `H=128`, `fx=1e-4`, `u_max ~ 1.23` case).


### 5.1.2 Random porous media: permeability sweep (validated)

`/home/sergio/tflbm_porous/` (`porous.cu`, `run_sweep.sh`,
`run_size_study.sh`, `validate_slabs.sh`, `analyze_porous.py`,
`plot_porous.py`, `README.md`, figures, `permeability.csv`). Fully periodic
`N x N` box, `round((1-phi)N^2)` cells tagged solid from a reproducible
`xorshift64*` stream seeded with `seed`, constant `+x` force re-written into
`BForce`/`pBForce` by `Setup` every step, `Sc=0`, `fx=5e-5`, `Tf=20000`.

```bash
cd /home/sergio/tflbm_porous
/home/sergio/mechsys/mechsyscc porous.cu -o porous_gpu
./porous_gpu 256 0.70 12345 5e-5 0.16666666666666666 256 20000
N=512 ./run_sweep.sh          # phi=0.50..0.90, 3 seeds, ~9.5 s per 512^2 run
./validate_slabs.sh           # exact-channel validation
python3 analyze_porous.py --csv && python3 plot_porous.py
```

`K = nu <u>_Darcy/a`, `<u>_Darcy = (1/N^2) sum_all u_x`, `u_x` = raw second
moment. **Convention settled:** periodic straight channels (solid row every
`SLABS` rows, gap `H=SLABS-1`) give `K=phi*H^2/12` exactly with the raw moment
(rel. error `1e-13` for `H=3/7/15`), while `K+0.5*nu*phi` is high by exactly
`0.5*nu*phi`. So the raw moment is the physical velocity and the older
`u+0.5*dt*fx` numbers (including the `tflbm10_porous` sweep: `K(0.90)=1.08`,
`K(0.70)=0.1163`, `K(0.60)=0.0422`) are that variant; subtract `0.5*nu*phi`
to compare (`K(0.90)=1.00`, `K(0.70)=0.0581`, `K(0.60)=-0.0077`).

Measured (`N=512`, 3 seeds, mean `+-` half-range): `phi=0.50 -2.93e-2`,
`0.55 -2.21e-2`, `0.60 -7.65e-3`, `0.65 +1.74e-2`, `0.70 +5.81e-2`,
`0.75 +1.27e-1`, `0.80 +2.47e-1`, `0.85 +4.76e-1`, `0.90 +9.99e-1`. The
seed spread is `<1%`. `K` is independent of `N=128..1024` at `phi=0.60/0.90`
(the staircase error is below the scatter).

**Accuracy floor:** below `phi_c=0.5927` the 4-connected fluid is disconnected
(largest cluster `<2%`) but the D2Q9 diagonal links still transport mass
through corner contacts (8-connected: one cluster, 98%), so the averaged `K` is
`-0.02..-0.03` instead of `0`, independent of `fx`, `nu` and `N`. Treat
`|K| <~ 0.03` as unresolved; `phi >= 0.70` is reliable. Do not fit a
near-threshold critical exponent with this data. For a physically cleaner
medium use a 4-connected solid mask (no diagonal-only contacts) or a finer
grid, and normalise by the percolating cluster only.

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


### 5.3 GPU FLBM: channel flow past an obstacle (triangle or rectangular box)

A standalone case for flow past a single obstacle in a plane channel at a small
Reynolds number is in `/home/sergio/tflbm_triangle/` (`tclbm_obstacle.cu`,
`animate.py`, `README.md`, the frames and the `.res`).  One binary runs both
shapes:

```bash
cd /home/sergio/tflbm_triangle
/home/sergio/mechsys/mechsyscc -o tclbm_obstacle tclbm_obstacle.cu
# rectangular box: 84 frames over the startup transient
./tclbm_obstacle --shape=rect --Lx=60 --Tf=25000 --dtOut=300 --render=1 --prefix=tclbm_rect
# triangle, converged field
./tclbm_obstacle --shape=tri --Tf=60000.5 --dtOut=60000 --render=1 --prefix=steady
python3 animate.py gif tclbm_rect --shape=rect --Lx=60 --dtout=300 --out tclbm_rect.gif
```

Walls are rows 0 and `ny-1`; the inlet is a Zou-He parabolic velocity profile
and the outlet a zero-gradient velocity condition, both adapted from
`tflbm/tclbm01.cu`.  `u_mean = (2/3) u_max` and `Re = u_mean*D/nu`, with `D`
the cross-stream extent of the obstacle (20.1% blockage at `ny=201`, `D=40`).
`--D=0` removes the obstacle and gives a control run.  The masks use closed
bounds, so `D=40` is 41 rows thick and `Lx=60` is 61 columns long; `ny` is odd
so both shapes are exactly symmetric about row `j=100`.

The easy trap is the wall geometry: with rows 0 and `ny-1` solid, the half-way
bounce-back walls sit at `y=0.5` and `y=ny-1.5`, so the fluid rows are
`1..ny-2`, the channel width is `L=ny-2` and the symmetry mirror is
`j <-> ny-1-j`.  Using `L=ny-1` shifts the inlet profile by one row and makes
the whole solution asymmetric even though grid and obstacle are symmetric; it
showed up as a frozen ~2% asymmetry in `ux` that persisted downstream.  The
no-obstacle control run is the quickest check: with the correct profile the
parabolic solution is preserved with `max|uy| = 6e-6` (1e-4 of `u_max`, the
symmetric residual of the Zou-He inlet; centreline `uy ~ 1e-16`) and the mirror
asymmetry stays at `1e-13`.

Measured at 640x201, `Nthread=256`, `Re=40` (bulk velocity 0.04, `nu=0.04`,
`Tau=0.62`):

- rectangular box 60x40: 80,000 steps in 19.7 s wall, 2501 solid nodes,
  realized `Re=40.00`; steady and exactly symmetric wake (centreline probe
  `uy ~ 1e-16`, asymmetry `1.3e-13 u_max`), `Lr/D = 1.88` with a peak reverse
  velocity of `0.15 u_mean`.  `Lr/D` overshoots to 2.05 at `t ~ 8000` and
  relaxes by `t ~ 20 000`; the side faces stay attached (`ux > 0` at `j=121`),
  separation is at the trailing edge.  Doubling the resolution (`D=80`,
  `Lx=120`, 1280x401, `nu=0.08`) gives `Lr/D = 1.81`, 3.7% below the coarse
  value, and `Lr` is 74/75 cells at `u_max = 0.03/0.10`, so this shape is
  grid- and Mach-converged too.
- triangular box (apex upstream, base 40): 60,000 steps in 15.0 s wall, 726
  solid nodes, steady symmetric wake with `Lr/D = 2.88` and a peak reverse
  velocity of `0.36 u_mean`.  Grid convergence (D=80, 1280x401, `nu=0.08`)
  gives `Lr/D = 2.86`, and `Lr` is 115 cells at `u_max = 0.03, 0.06, 0.10`
  (same `Re`, `nu` scaled with `u_mean`), so the scalar is grid- and
  Mach-converged.  `Lr/D` is 1.65 at Re=20, 3.75 at Re=60, 4.4 at Re=80 and
  4.8 at Re=100 (20% blockage), and 4.4 at Re=40 with 10% blockage (`ny=401`).
- neither obstacle sheds in this confined configuration: the symmetric steady
  wake persists to `Re=100` for the triangle, so a Karman street cannot be
  produced without unconfining the domain.
- the generalised binary reproduces the earlier triangle-only binary bit for
  bit (all HDF5 datasets of frame 79 at `t=19 750`), so `--shape` did not
  change the numerics.

`animate.py gif` reads the `Velocity_0`/`Gamma` frames (`x` fastest, three
velocity components per node) and writes a two-panel speed/vorticity GIF with a
fixed axes box (no layout jitter) and a figure height derived from the panel
aspect, so `--xmax` crops stay equal-aspect.  There is no
`ffmpeg`/`convert`/`gifsicle` on this machine; assemble GIFs with
`PIL.Image.save(..., save_all=True)`.


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
- The single-body DEM checks are paired by number and are *not* interchangeable
  with `test_periodic`: `tdem/test_01` (CPU sliding block, Tf=10, prints
  `Observed L`) with `tdem/test_cu_01` (GPU, same physics), and
  `tdem/test_02` (CPU turning table, Tf=30, `dt=1e-5` -> 3e6 steps) with
  `tdem/test_cu_02` (GPU). The "02" pair also writes a `<FileKey>_trayectory.res`
  trajectory (`Time x y xt yt`, `xt`/`yt` being the ideal rolling orbit) via the
  `Report` callback, which is the cheapest way to check a rolling-contact change
  without comparing HDF5 frames.
- `test_01`/`test_02`/`test_periodic` pass no `Nproc` to `Solve`, so they inherit
  the OpenMP default and run multithreaded (`test_01` shows `CPU 18.8s` against
  a `2.57s` wall). That reorders force accumulation; for an exactness comparison
  run them with `OMP_NUM_THREADS=1` and check that the printed `CPU` time equals
  the wall time before trusting a bitwise result.
- `Array::Find()` returns `-1` on failure; assigning it to `size_t` is a bug.
  Check the signed result before using it as an index.
- `std::is_same`-based feature counters can be wrong because feature
  containers may hold pointers (`Array<Edge*>`, etc.), not values. Prefer
  explicit type traits or per-container logic.
- Warning policy: the only expected third-party warning is Blitz
  `tinymat.h(60): warning #191-D`. Keep MechSys sources warning-free.
- **The `benchmarks/cpu_voronoi/` directory does not exist in this checkout**
  (referenced by section 3.5 and by `SKILL.md`). To reproduce the reduced
  exactness run, derive it from `tdem/tper.inp`: set `ptype = voronoi`,
  `Lx=Ly=Lz=10`, `nx=ny=nz=10`, `T0=Tf=0.01` (400 steps at `dt=5e-5`). The
  unmodified `tper.inp` is `sphereboxhcp` with `Tf=100`, so it exercises the
  **sphere** contact path, not the polyhedral one.
- `test_periodic` writes both `<key>_a.hdf5`/`<key>_b.hdf5` (full state,
  one dataset per particle field) and `<key>_a.h5` / `<key>_a_0000.h5`
  (19 flat visualisation arrays). The `.hdf5` files carry ~15k–430k datasets;
  a sphere run is far larger than a Voronoi run.
- `Interactons/Fdvv`, `Interactons/Fdr` and per-particle `Diam` are written
  from **uninitialised buffers** in `Domain::Save` (`new double[3*NI]` filled
  only for sphere pairs, written unconditionally). They differ between two
  identical baseline runs (`Fdr` 3e-3, `Diam` denormal). Exclude them from any
  HDF5 dataset comparison; they are pre-existing, not a regression signal.
- `PAngacceleration` at frame 0 came from an uninitialised `Particle::wa` in
  the baseline (reproducible garbage); it is deterministic and zero when `wa`
  is initialised.

### `MostlySpheres` segfaults on sphere-sphere contacts

`Domain::ResetContacts` builds the pair list with

```cpp
if (!MostlySpheres) PairtoCInt[hash] = new CInteractonSphere(...);
...
if (PairtoCInt[hash]->UpdateContacts(Alpha,Per,iter)) ...
```

With `MostlySpheres = true` the key is never inserted, `operator[]` returns a
default-constructed **NULL** pointer, and `->UpdateContacts` segfaults. It is
reproducible with two touching spheres at `Nproc=1`, in both the current
`HEAD` and in PR #2, for contact law 0; `gdb` stops in
`DEM::Domain::ResetContacts()` under `GOMP_parallel`. The linear
`Domain::CalcForceSphere` fast path is therefore unreachable for contacting
spheres, and `ttt <key> <ncores> 1` cannot work as written. PR #2 fixes only
the `ContactLaw==1` case (it adds `|| ContactLaw==1` to the creation test);
the linear case still crashes. Treat `MostlySpheres` as broken until
`PairtoCInt[hash]` is null-checked, and prefer `MostlySpheres = false` +
`CInteractonSphere` (as `PID_v2`, `TX_PID_periodic` and the new tests do).

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
- Post-processing facts worth not rediscovering: the increment flatness
  `m4/m2^2 = K4/(K2^2 dt) = kappa4+3` of the force magnitude diverges as
  `dt^-1` (`4.5 -> 79.8`, ordered; `4.0 -> 75.0`, random) between `dt=0.25` and
  `dt=0.005`, which is the signature of contact birth/death jumps, not of a
  diffusion.  So no `dt -> 0` Fokker--Planck limit exists and Pawula's theorem
  yields only a deterministic transport equation there.
- `/home/sergio/ttt_compare/pawula.py` exceeds the 600 s bash cap when it
  processes both packings in one pass (it reads ~9 GB of `bf` frames per run).
  It takes a run selector (`python3 pawula.py ord` / `rnd`) and merges into
  `pdfform/pawula.json`.  `fine_lag.py` handles both in one pass.

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
  `tlbmmpm/tlbmmpm01`, `tlbmmpm/tlbmmpm02`, `tlbmmpm/tlbmmpm_cu_0[12]`,
  `tmpm/tmpm0[45]`, `tmpm/tmpm_cu_0[1-5]`, `MPM_OPTIMIZATION_CHANGES.md`,
  `MECHSYS_CHANGES.txt` and the other MPM binaries from the build. They are
  regenerable and must not be deleted without asking.
- New untracked LBMMPM CUDA sources are `lib/lbmmpm/lbmmpm.cuh` and
  `tlbmmpm/tlbmmpm_cu_01/02.cu`; keep them with the tracked changes until a
  commit is requested.

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
| `tmpm01`/`tmpm_cu_01` (non-CPDI default) | GIMP beam, fixed base, end force | 5,000 | 125,000 | ~86,000 |
| `tmpm01`/`tmpm_cu_01` (`cpdi`) | CPDI tetrahedral beam, fixed base, end force | 44,032 | 238,328 | ~216,000 |
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
| `tmpm01`/`tmpm_cu_01` CPDI `ndiv=50`, `Tf 5->0.5`, 10 frames | 51.8 s | 23.8 s |
| `tmpm02`: `num_particles_x 40->10`, `Tf 1.5->0.02` | 0.87 s | 3.72 s |
| `tmpm03`: `particles_per_diameter 70->20`, `Tf 5->0.02` | 0.39 s | 7.34 s |
| `tmpm04`: `particles_per_diameter 60->10`, `tf 200->1` | 0.07 s | 1.78 s |
| `tmpm05`: `num_particles_x 40->10`, `Tf 1.5->0.01` | 0.43 s | 2.35 s |

All ten exit 0 and write 36-100 HDF5 frames with finite `Position`,
`PVelocity`, `Pressure`, `Sigma`, and `Strain`.

The CPDI smoke row above is the historical `ndiv=50` mesh, reducing only
`Tf`; setup is about 17 s on both paths.  The checked-in CPDI case is now
`ndiv=125`; the smoke is kept only as the CPU/GPU parity check.  `Solve` was
33.88 s CPU versus 5.44 s GPU for the 8,640-step window.  That parity check
used `Fb=0.1`; the current checked-in tip force is `Fb=5e-2`.  With `Nproc=1`, all 11 HDF5 frames had identical
dataset payloads except one `PVelocity` value at frame 7 (`4.4e-16` absolute);
the GPU sums use atomics, so this is the expected class of difference.

Halving the tip force (`Fb 0.10 -> 0.05`, both sources) halves the absolute
deflection, `2.297e-3 -> 1.148e-3 m`, but not the relative error: the
full-`Tf=5` GPU run gives `1.14806452e-3 m` versus analytical
`9.75102376e-4 m`, ratio `1.1774` (the `Fb=0.1` ratio was `1.1780`).  The
`y/Fb` profiles collapse, so force scaling is not an accuracy fix.  The
baseline and half-force result files are
[`tmpm01_Comparison_F01.res`](tmpm/tmpm01_Comparison_F01.res) and
[`tmpm01_Comparison_F005.res`](tmpm/tmpm01_Comparison_F005.res).

CPDI body-mesh refinement with the current boundary-corner load (GPU,
`Tf=5`, `Fb=0.05`, `beta=1`, `Nproc=1`).  `ResizeDomainMesh(...,1.0)` sets the
background spacing to the largest body-cell dimension, so the underlying grid
is refined with the body mesh:

| `ndiv` | particles | grid nodes | `Dx` (m) | Tip (m) | MPM/analytic tip | MAE (m) |
|---:|---:|---:|---:|---:|---:|---:|
| 50 | 3,315 | 17,576 | 7.551e-3 | 8.0449e-4 | 0.8250 | 6.842e-5 |
| 75 | 10,045 | 50,653 | 5.266e-3 | 8.7912e-4 | 0.8878 | 4.554e-5 |
| 100 | 23,180 | 125,000 | 3.962e-3 | 9.4365e-4 | 0.9486 | 2.065e-5 |
| 125 | 44,032 | 238,328 | 3.204e-3 | 9.7016e-4 | 0.9694 | 1.217e-5 |

Refinement moves the tip toward the analytical value and reduces MAE by about
5.6x from `ndiv=50` to 125.  The checked-in CPDI case is `ndiv=125`;
non-CPDI remains `ndiv=50`.  The pre-optimization GPU `Solve` was about
17m23s at 125; the 2026-11-06 hot-path pass reduced it to 4m02.85s (see
below).  An earlier `Tf=2` refinement test (ratios
1.144/1.292/1.440) used the old combined particle+corner load; it does not
describe the current CPDI load path.

### 9.3 CPDI CUDA hot-path optimization (2026-11-06)

The CPDI algorithm and the six-launch order are unchanged.  The full
`ndiv=125` case (44,032 particles, 238,328 nodes) dropped from
**17m23.33s to 4m02.85s** (`4.30x`); the reduced `ndiv=50` case dropped from
**47.53s to 26.74s** (`1.78x`).

| CPDI case | particles | grid nodes | before `Solve` | after `Solve` | speed-up |
|---|---:|---:|---:|---:|---:|
| `ndiv=50` | 3,315 | 17,576 | 47.53 s | 26.74 s | 1.78x |
| `ndiv=125` | 44,032 | 238,328 | 1043.33 s | 242.85 s | 4.30x |

The deflection is unchanged: tip `9.70156508e-4 m`, analytical
`1.00081114e-3 m`, ratio `0.9693702`, MAE `1.2168888e-5 m`.  The optimized
`.res` differs by at most `1e-16 m`; `ndiv=50` is bitwise identical.  HDF5
fields shift by round-off (about `2e-8` relative in `Sigma`/`PVelocity`), not
by a changed formula.

What changed in `lib/mpm/mpm.cuh`, `lib/mpm/Domain.h` and
`lib/mpm/Particle.h`:

- `cudaComputeAn` precomputes `An=Fn/Mass` once per step; the gather no longer
  divides for every contribution.
- `ParticleCU::GradvecN[4]` caches `Gradvec/(6*Vt)` from `CalcVolCPI` and
  `UploadParticle`; scatter/gather use it directly.
- Particle `v/L/x`, corner `x`, and scatter `v/m` stay in registers.
- CPDI particle scatter no longer marks nodes: the corner scatter already
  visits every node it can touch, so the redundant `atomicExch`/`atomicAdd`
  and `SelIndexes` write are gone.
- Corner scatter skips momentum/force atomics for exactly zero corner forces
  and skips node atomics for zero shape values; the `Nodes`/`Sfunc` slots are
  still filled.
- The non-CPDI path has a separate pass, described in Section 9.4.

`beta` is the ratio `Dx/dmax`: below 1 is a finer background grid, above 1 is
coarser.  At the current `ndiv=125`, `Fb=0.05`, full `Tf=5`:

| `beta` | grid nodes | `Dx` (m) | tip/analytic | MAE (m) | GPU `Solve` |
|---:|---:|---:|---:|---:|---:|
| 1.00 | 238,328 | 3.204e-3 | 0.9694 | 1.217e-5 | 17m23s |
| 1.25 | 117,649 | 4.005e-3 | 1.1652 | 8.878e-5 | 17m15s |
| 1.50 | 68,921 | 4.806e-3 | 1.1392 | 7.653e-5 | 17m17s |

`beta=1.0` is clearly best; coarsening overshoots.  The earlier positive
impression of `beta=1.5` came from the old combined particle+corner load and
does not transfer to the boundary-corner path.  Those beta timings are
pre-optimization; at `ndiv=125` the optimized `Solve` is about 4m03s, so use
beta for accuracy, not as a time knob.

`tmpm01`/`tmpm_cu_01` now take an optional path argument:

```bash
./tmpm01 [Nproc] [cpdi|mpm]
./tmpm_cu_01 [Nproc] [cpdi|mpm]
```

The default is `mpm` (standard particle-grid, `ndiv=50`); pass `cpdi` to
enable the refined CPDI mesh (`ndiv=125`, `beta=1`).  In CPDI, the end load is
applied **only to boundary corners in the tip slice**.  The last cross-section
of the rectangular beam is an external face, so selecting corners by `x` alone
is sufficient; no `BFac`/`BoundaryMesh` mask is needed.  `ForcePar` is empty.  Non-CPDI applies the
load only to `ForcePar`, with `ForceCor` empty.

Final frame at `Fb=0.05`, full `Tf=5`:

| Path / load | Particles | Grid nodes | Tip (m) | Tip/analytic | MAE (m) |
|---|---:|---:|---:|---:|---:|
| Non-CPDI, particles, `ndiv=50` | 5,000 | 125,000 | 9.3456e-4 | 0.9623 | 1.539e-5 |
| CPDI, boundary corners, `ndiv=50` (44) | 3,315 | 17,576 | 8.0449e-4 | 0.8250 | 6.842e-5 |
| CPDI, boundary corners, `ndiv=75` (99) | 10,045 | 50,653 | 8.7912e-4 | 0.8878 | 4.554e-5 |
| CPDI, boundary corners, `ndiv=100` (144) | 23,180 | 125,000 | 9.4365e-4 | 0.9486 | 2.065e-5 |
| CPDI, boundary corners, `ndiv=125` (229, checked-in CPDI) | 44,032 | 238,328 | 9.7016e-4 | 0.9694 | 1.217e-5 |
| CPDI, old particles + all tip corners, `ndiv=50` | 3,315 | 17,576 | 1.1481e-3 | 1.1774 | 6.191e-5 |

The non-CPDI path remains closest to the analytical tip, but CPDI converges
toward it as the body mesh is refined.  The checked-in CPDI path uses
`ndiv=125`, `beta=1`, with no particle load.  The optimized non-CPDI final
GPU `.res` matches the previous GPU baseline to `1e-13 m`; CPDI matches to
`1e-16`.  CPU/GPU HDF5 payloads are not promised bitwise after the
separable-shape/`An` caches; the observed `PVelocity` shift is at the `1e-8`
relative level.  This is an accuracy comparison, not a general GPU speed
claim.

`tmpm04` at `particles_per_diameter = 20` (8,000 particles, 14,200 steps,
`tf = 200`) ran in 22.9 s CPU wall. `Vx_max` converged to `8.57e-4 m/s` against
the case's own analytic `7.37e-4 m/s` (`Rerror = 16.3%`); that resolution has
only five cells across the pipe radius, so this is a convergence check, not an
accuracy claim. The checked-in 60-ppd case (216,000 particles, 42,400 steps,
96 OpenMP threads, `tf = 200`) completed in 1932 s (32 min 12 s) wall before
the CPU optimization and 1676 s (27 min 56 s) after it (Section 9.5), with
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

### 9.4 Non-CPDI CUDA pass (2026-11-06)

The default GIMP beam (5,000 particles, 125,000 nodes, 86,393 steps) is now
**5.46x** faster: `Solve` fell from **1m54.82s to 21.03s**.  The final `.res`
matches the previous GPU baseline to `1e-13 m`: tip `9.34560567e-4 m`,
analytical `9.71136179e-4 m`, ratio `0.9623373`, MAE `1.5385435e-5 m`.
HDF5 fields shift by round-off (about `7.4e-8` relative in `PVelocity`,
`2.0e-9` in `Sigma`).

What changed:

- GIMP and the B-splines are separable.  `cudaParticleToNode` evaluates the
  1-D shape factors once per axis (`MPMShapeFunction1DCU`) and combines them,
  instead of evaluating the full 3-D shape for every node.  For a 3x3x3
  support that is 9 1-D evaluations instead of 81.  The direct 3-D path is
  retained as a fallback when an axis exceeds 8 nodes.
- The non-CPDI gather consumes the same precomputed `An=Fn/Mass` array as
  CPDI, so it no longer divides at every particle-node contribution.
- The gather keeps `v`, `x` and `L` in registers; the scatter keeps `v/m` in
  registers.

The CPDI kernels and their launch order are unchanged by this pass.

### 9.5 CPU performance changes

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
| pre-CPDI `tmpm01` full (5,000 particles), 86,391 steps, 100 frames (wall) | 92.42 s | 33.77 s | 2.74x |
| `tmpm04` full, 42,437 steps, 100 frames (wall) | 1932.08 s | 1676.27 s | 1.15x |

`Nproc=1` output is dataset-bit-identical for pre-CPDI GIMP (`tmpm01`),
Drucker-Prager (`tmpm02`), Tait water (`tmpm03`) and mu(I) (`tmpm05`).  At
`Nproc>1` only the floating-point reduction order can change; full pre-CPDI
`tmpm01` at 96 threads differs only in `PVelocity` by `2.8e-14` (7.6e-11
relative), and `tmpm04` frames differ by at most `~1e-14` relative in
`PVelocity`/`Sigma`.

`tmpm04` is the dense limit: 216,000 particles on 38,880 nodes (~5.5 per
cell), so each node is updated ~150 times per step and the scatter is bound
by cache-line ownership.  On a short 5-frame `tmpm04` run, 8/32/96 threads
take 98.1/60.6/71.3 s after the optimization; for such dense grids fewer
threads than `0.75*MAXTHREADS` can be faster.  An eight-group partial-grid
prototype reached 60.8 s at 96 threads but regressed sparse cases, so it was
not kept.

### 9.6 CUDA MPM performance changes

The GPU path was optimized after the CPU work; full design, phase profile and
numbers are in
[`CUDA_MPM_OPTIMIZATION_CHANGES.md`](CUDA_MPM_OPTIMIZATION_CHANGES.md).
Files changed: `lib/mpm/mpm.cuh`, the `USE_CUDA` parts of
`lib/mpm/Domain.h`, `lib/mpm/Particle.h`, `lib/mpm/Corner.h` and
`lib/mpm/shape_functions.h`.

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
  B-spline and for quadratic/GIMP at the usual particle size), down from 125.
  The grid-path size is 3584 B; the CPDI fields (`Corners`, `Gradvec`, `Vt`)
  add 128 B, so the current structure is 3712 B.
- `MPM::Domain::Solve` now has a CUDA CPDI/body-mesh branch.  It launches
  `cudaCornerToNode`, `cudaParticleToNodeCPDI`, `cudaSolveNode`,
  `cudaMoveCorners` and `cudaNodeToParticleCPDI` in the same order as
  `OneStepCPI`, and downloads corners plus particles at output frames.  The
  checked-in `tmpm_cu_01` mirrors `tmpm01` (`AddRectangularBeamMesh`,
  `ResizeDomainMesh(...,1.0)`) and applies the tip load with a device
  `ForceCor` list and `SetCornerForce`.  This is separate from the LBMMPM
  coupling in §9.5, which still uses the CPU `MPMDOM.OneStepCPI`.
- A final `cudaDeviceSynchronize()` makes `Solve` timing cover queued work.

Reduced four-frame pre-CPDI benchmarks (`Nthread=256`, RTX A6000) improved
`tmpm_cu_01/03/04/05` by 4/8/9/10%; `tmpm_cu_02` regressed by about 1%
(short steps, frequent particle deletion/re-upload).  `tmpm_cu_01` is now the
CPDI beam, so its old four-frame timing is historical.  These are `Solve`
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

### 9.7 CUDA LBMMPM

The LBMMPM CUDA path mirrors the CPU coupling in `lib/lbmmpm/Domain.h`; the
kernels are in the new `lib/lbmmpm/lbmmpm.cuh`.  It currently supports
**D3Q15 only** (`LBMDOM.Nneigh == 15` is checked in `Solve`).  In this
two-domain coupling the MPM material-point update (`MPMDOM.OneStepCPI`)
remains on the CPU; only the LBM and the immersed-boundary coupling run on
the device.  Standalone `MPM::Domain::Solve` has a separate CUDA CPDI path
(§9.4).  The detailed data map, upload
schedule and launch order are in `CUDA_LBMMPM_STRUCTURE.txt`.

Host-side work in `LBMMPM::Domain`:

- `UpLoadDevice` uploads `F/Rho/Vel/BForce/IsSolid`, the body-cell and facet
  vertex ids, corner positions/velocities/forces, and builds two compact work
  lists: `(body cell, candidate lattice cell)` and `(facet, candidate cell)`.
  The candidate cells come from the same AABB bounds as the CPU loops.
- Immediately after `UpLoadDevice`, `Solve` calls `UpdateDeviceCorners` and
  `ImprintLatticeCUDA` once before the first output.  Without that extra
  imprint, `frame 0000` has `Gamma=0` in the initial body region while the CPU
  has already imprinted, so the first HDF5 output differs.
- `UpdateDeviceCorners` re-uploads moving corner `x/v` and clears `bCornerH`.
- `ImprintLatticeCUDA` zeroes `Gamma` with `cudaMemsetAsync` and launches, in
  order: `cudaImprintCells`, `cudaRefill`, `cudaMarkFacets` and
  `cudaImprintFacets`.  The separate marking pass finalises `Gamma` before the
  forcing pass reads the solid neighbour; without it a single mixed kernel
  races marking against the `Gamma[ns]` read.
- `DownloadCornerForces` copies the atomic corner forces back into
  `MPMDOM.Corners[*]->h`.
- Do not restore the legacy `ImprintLattice` assignment
  `F[nxs][nys][nzs][k] = F[nx][ny][nz][k]`.  The current CPU and CUDA paths
  both use the original solid-neighbour population; the assignment was
  disabled in `e1ad057` and raises the current CPU cube `CD` from 6.83 to
  31.5 at `Re~5` (formula 6.41).
- The per-step LBM sequence is `cudaCollideMPM<15>`, pointer swap,
  `cudaStreamMPM` (periodic), pointer swap.  `cudaCollideMPM` reconstructs
  `Rho`/`Vel` itself, so there is no separate `cudaVelDen` launch.
- Device `F` is direction-major.  Use `FLBM::FIDX(Ncells,Nneigh,il,k,ic)` from
  `lib/flbm/lbm.cuh` for device addressing; the boundary kernels compute the
  cell index `n` once and pass it to `FIDX`.  Do not add local copies of the
  index formula.

`cudaImprintFacets` intentionally keeps the CPU strict tests
(`r>0, s>0, t>0, s+t<1`, normal dot `>0`).  The CPU free `SolAlt` and the
device `Mat3::SolAlt` method in `lib/linalg/matvec.h` use the same arithmetic;
on facet edges the CPU returns exact zeros and the device inverse leaves
~1e-16 noise, which changes which surface cells are marked.  The CMake targets therefore compile the CUDA
sources with **`--fmad=false`**; with FMA contraction the reduced and full t02
cases disagree on 145 and 176 `Gamma` cells respectively, and with the flag
`Gamma` is bitwise identical in those tests.  For one-off builds use

```bash
NVCC_PREPEND_FLAGS=--fmad=false /home/sergio/mechsys/mechsyscc case.cu
```

Build the two checked-in CUDA drivers from the root-configured tree:

```bash
cd /home/sergio/mechsys/tlbmmpm
make -j8 tlbmmpm_cu_01 tlbmmpm_cu_02
```

`tlbmmpm_cu_01.cu` is the elastic-beam case with device inlet/outlet/top
boundary kernels; `tlbmmpm_cu_02.cu` is the settled-cube case with a device
`BForce = Rho*fluidforce` write in `Setup`.  Both use the same full grids and
`Tf` as `tlbmmpm01/02`, so they are production-length runs; reduce a copy for a
smoke test.  The CPU examples remain unchanged.

Measured on the RTX A6000 (`Nthread=256`, 8 CPU threads for the CPU side,
CUDA 12.4, double precision):

| Case | Grid / steps | CPU solve (wall) | GPU solve (wall) | Speedup (solve) |
|---|---|---:|---:|---:|
| `t02` reduced, `Tf=100`, `dt=0.2` | 61x21x21, 500 | 0.937 s (1.02 s) | 0.617 s (0.81 s) | 1.5x |
| `t02` full grid, window `Tf=200` | 241x61x61, 1000 | 39.62 s (39.93 s) | 5.80 s (6.16 s) | 6.8x |
| `t01` reduced beam, `Tf=0.02`, `ndiv=40` | 40x40x40, 277 | 1.285 s (1.42 s) | 0.806 s (1.26 s) | 1.59x |
| `t01` full-grid window, `Tf=0.002`, `ndiv=40` | 200x200x200, 28 | 14.69 s (16.62 s) | 6.35 s (8.63 s) | 2.31x |

The 2026-11-06 hot-path pass keeps this order and tunes the kernels.  On the
same RTX A6000 with `Nthread=256`, `-maxrregcount=128`, `--fmad=false`:

| Case | Grid / steps | Before `Solve` | After `Solve` | Speedup |
|---|---|---:|---:|---:|
| `t02` full window, `Tf=200` | 241x61x61, 1001 | 4.993 s | 4.062 s | 1.23x |
| `t01` full window, `Tf=0.01`, `ndiv=40` | 200x200x200, 139 | 10.664 s | 8.938 s | 1.19x |

The pass fuses the velocity/density reconstruction into `cudaCollideMPM`,
precomputes `1/Cs` and `1/rho`, removes the intermediate `Ftemp` stores from
the entropic fixed-point loop, replaces stream `%` wrap with conditional
add/subtract, and uses `cudaMemsetAsync` for `Gamma`.  A single-scan
`alpha=1` limiter is faster but is **not** equivalent: a negative-noneq
direction can become negative only after alpha is reduced, so the fixed-point
loop must be retained.  Register caps 64/72/80/96/160/192, `Nthread=128/512`
and a shared-memory `NonEq` variant gave no repeatable gain.

Later frames are not bitwise.  Repeated baseline GPU runs already differ by
about `7.2e-4` relative in `Density_0` and `2.2e-3` in `Velocity_0`; the
optimized-vs-baseline difference is the same order, while `frame 0000` is
byte-identical.  The reduced t01 final beam deflection differs by `5.6e-6 m`
(`5.6e-5` relative).  Do not treat small later-frame LBM differences as an
optimization regression without a baseline-vs-baseline rerun.

The reduced t02 GPU run is 7.5x faster than the same case at `Nproc=1`
(4.693 s solve).  On the reduced grid `Nthread=128` was 0.70 s and `Nthread=256`
0.62 s in the final build; the difference is smaller than the run-to-run
spread.  The full-grid window is 1000 of the checked-in 100,000 steps; it is a
per-step benchmark, not the full case.

Accuracy in those windows:

- Reduced t02 at `t=50`: `Re` 26.021 (CPU) vs 25.789 (GPU), -0.9%; `CD` 8.369
  vs 7.66, -8.5%.
- Full t02 at `t=180.2`: `Re` 94.720 vs 94.679, -0.04%; `CD` 0.6727 vs 0.5109,
  -24%.  At `t=20.2` the `CD` gap is already -8.5%.
- Reduced t01 at `t=0.015`: `Re` 12.202 vs 12.257, +0.45%; beam deflection
  `6.526e-3` vs `6.259e-3`, -4.1%; `CD` 0.6741 vs 0.6355, -5.7%.
- Full-grid t01 at `t=0.001013` (28 steps): the printed `Re` and both
  deflections were identical (`Re=0.330572695`, deflection `0`); this is a
  short startup window before measurable bending, not a full `Tf=1.5` run.
- `Gamma` was bitwise identical at the compared frames in both t02 windows.
- At `frame 0000` all common LBM and MPM dataset payload bytes are identical
  after the initial imprint fix.  The HDF5 files are still not byte-identical:
  the raw files differ by 12 bytes (LBM) or 22 bytes (MPM) of HDF5 metadata,
  mostly object timestamps; compare datasets rather than whole-file MD5.
- After one step the payloads are not byte-identical.  Reduced t02 `frame 0001`
  has `Gamma` exact but `Density_0`/`Velocity_0` differ (8706/28862 bytes).
  Full t02 `frame 0001` differs in 2.95M/13.4M density/velocity bytes; t01
  reduced `frame 0001` has MPM `Sigma/Strain/PVelocity/Position/Vertex`
  differences.  The initial `frame 0000` exactness therefore does not imply a
  bitwise trajectory.

The velocity fields and `Re` agree closely, but the integrated corner force
(`CD`) is not a validated GPU quantity: the CPU OpenMP facet loop and the GPU
atomic corner sums update `F`/`h` in different orders, and this diagnostic is
sensitive to local cancellation.  `atomicAdd` also makes repeated GPU runs
non-bitwise (reduced `CD` varied by ~0.2% between runs).  Do not present the
GPU `CD` as matching the CPU; the t01 deflection is the better mechanics check.

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
| `CUDA_LBMMPM_STRUCTURE.txt` | CUDA LBMMPM host/device data map, upload schedule, launch order, kernel responsibilities and extension checklist |
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
