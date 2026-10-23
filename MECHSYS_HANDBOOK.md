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

## 8. Current working-tree state

As of the last update of this handbook:

- `HEAD` is `94f64b8`.
- Modified files:
  - `lib/dem/distance.h`
  - `lib/dem/domain.h`
  - `lib/dem/interacton.h`
  - `lib/dem/particle.h`
- Added reports:
  - `CPU_DEM_OPTIMIZATION_CHANGES.md`
  - `CUDA_OPTIMIZATION_CHANGES.md`
  - `REDUCED_TEST_CU_03_PERFORMANCE.md`
- Added benchmark inputs:
  - `benchmarks/cpu_voronoi/tper_small.inp`
  - `benchmarks/cpu_voronoi/tper_med.inp`

These changes are uncommitted. Do not treat them as `HEAD`, and do not delete
them while experimenting.

## 9. Documentation map

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
| `benchmarks/cpu_voronoi/` | reproducible CPU Voronoi inputs |

## 10. Maintaining this knowledge

When you learn something that would save a future agent from repeating a
failure:

1. add it to the “Known traps” section;
2. if it is an invariant, add it to `AGENTS.md`;
3. if it is a benchmark, record the exact input, thread count, window, and
   whether it is wall or solve time;
4. do not record transient `/tmp` paths or one-off logs as authoritative;
5. keep this handbook concise enough to be read in one sitting.
