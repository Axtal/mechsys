# MechSys Engineering Memory

This is a distilled record of the practical experience gained while optimizing
the DEM, CUDA DEM, and LBM-DEM parts of MechSys. It is written for people who
download the repository and need to understand the code quickly, and for
coding agents that want more than a source-file tour.

For operational agent rules use [`AGENTS.md`](AGENTS.md). DeepSeek Harness
also discovers the project-local skill at
[`.agents/skills/mechsys/SKILL.md`](.agents/skills/mechsys/SKILL.md); a visible
copy is [`SKILL.md`](SKILL.md). The skill loads this material on a
MechSys-related task instead of keeping all of it in the always-on context.
For the full architecture, commands, and verification
procedure use [`MECHSYS_HANDBOOK.md`](MECHSYS_HANDBOOK.md). The detailed
optimization reports are:

- [`CPU_DEM_OPTIMIZATION_CHANGES.md`](CPU_DEM_OPTIMIZATION_CHANGES.md)
- [`CUDA_OPTIMIZATION_CHANGES.md`](CUDA_OPTIMIZATION_CHANGES.md)
- [`REDUCED_TEST_CU_03_PERFORMANCE.md`](REDUCED_TEST_CU_03_PERFORMANCE.md)

This file is intentionally shorter than the handbook. It records the things
that are easy to get wrong and the shortest path to a trustworthy result.

## 1. Orientation

- The real headers are in `lib/`; `mechsys` is a symlink to `lib/`. Include
  headers as `<mechsys/dem/...>`.
- There is no single MechSys executable. The tests under `tdem/`, `tlbmdem/`,
  `tlbm/`, `tflbm/`, etc. are the entry points.
- The root is already CMake-configured in this checkout. Build individual
  targets from the test directory:
  ```bash
  cd tdem
  make -j8 test_periodic
  make -j8 test_cu_01 test_cu_02 test_cu_03 test_cu_04

  cd ../tlbmdem
  make -j8 tlbmdem01 tlbmdem02
  make -j8 tlbmdem_cu_01 tlbmdem_cu_02 tlbmdem_cu_03 tlbmdem_cu_04
  ```
- The optimization work described here was developed as uncommitted changes
  in one checkout. Do not assume a given `HEAD` contains it; check `git status`
  and the optimization reports before comparing against a baseline.
- Read an existing file before replacing it. In particular, `tdem/test_cu_04.cu`
  uses CRLF line endings; do not rewrite it with a tool that normalizes them.

## 1a. VisIt

- Binary: `/usr/local/visit/bin/visit` (VisIt 3.3.2), added to `PATH` by
  `~/.profile`.
- Headless: `/usr/local/visit/bin/visit -cli -nowin -s script.py file.xmf`.
- VisIt accepts MechSys' one-dimensional HDF5 datasets in XDMF; leave that
  convention alone.  The XMF must reference the HDF by basename.
- Temporal collections: give every child grid `Name="DEM_Faces"` so the
  mesh variable is `DEM_Faces`.  Use `SetTimeSliderState(i)` to animate.
- Pseudocolor: `AddPlot("Pseudocolor","DEM_Faces/Velocity")`; threshold with
  `DEM_Faces/Tag` for a standalone file.  For a collection use
  `Velocity` / `Tag`.  Add a second `Mesh` plot of the walls with
  `MeshAttributes.SetOpacity(0.18)` and a negative-tag threshold.
- VisIt 3.3.2 CLI attribute names: `ThresholdAttributes` takes
  `lowerBounds`/`upperBounds` **tuples** (not `lowerBound`/`upperBound`), and
  the 3-D view class is `View3DAttributes` (not `View3DAtts`).
- `SaveWindow` appends a frame number to the requested filename.
- A nonzero CLI exit after a successful save can be a warning; inspect the
  log and PNGs.
- A positive-`Tag` threshold alone includes parked/recycled boxes: the
  avalanche driver keeps them at `x = -10 - 0.1*tag`, `z = -0.4` with positive
  tags.  They are storage, not grains in the container.  For a clean render,
  add a `Box` operator cropped to the slot (or exclude `z < -0.2`); do not
  diagnose a wall leak from a positive-tag image.  VisIt 3.3.2 exposes
  `BoxAttributes.SetMinx/SetMaxx/...`; applying it after the tag threshold
  works on standalone frames.  On a temporal collection, VisIt 3.3.2 produced
  intermittently blank frames with that operator; render each frame separately
  for an active-only animation.

## 2. The DEM solve cycle

The order in `DEM::Domain::Solve` is important:

1. reset particle forces and torques (including fixed `Ff`/`Tf`);
2. evaluate all current interactions;
3. accumulate force/torque per particle;
4. translate and, if enabled, rotate particles;
5. compute `maxdis` from the current vertices and `Vertso`;
6. rebuild interaction/Verlet lists when `maxdis > Alpha`;
7. advance `Time` and `iter` once.

`MaxDisplacement()` is measured from the last `ResetDisplacements()`, which is
called while rebuilding contacts. This is why the CPU cutoff uses accumulated
`ContactDispGlobal()` rather than a single-step bound.

### 2a. Domain-wide energy and bounds are not active-only

`Domain::CalcEnergy` loops over **all** `Particle` objects.  A box that has
already left the chute and is free-falling below the base still contributes
translational and rotational kinetic energy, even though it is no longer part
of the pile.  In the avalanche driver this produced `Ekin=3.47 J` while every
active box had `vmax=0.0034 m/s`, causing a false “pile still moving” stop.
For quiescence or avalanche criteria, filter with the same `IsActive()`/
parked logic as the driver and compute the kinetic sum over active boxes only.

A per-smallest-box energy threshold is also too strict for a large pile:
hundreds of boxes at a few mm/s can exceed one box's energy.  Use
`0.5*M_active*v_settle^2` for the post-check, plus an explicit `vmax` bound.
Do **not** pass that mass-scaled value as `Domain::Solve`'s `minEkin`: the
solver can then stop while one fast box is still moving and the post-check
will reject the event.  Keep the solver's strict smallest-box threshold.

`Domain::BoundingBox` likewise scans free particles across the whole domain;
inspect and clamp its result before using it as a pile-top estimate.

### 2b. Reset a recycled particle's orientation

A recycled DEM box can keep the orientation it had when it tumbled out.  In a
narrow slot, a long axis across the gap makes the oriented `MaxY()` larger
than the clear width, so the insertion code stops with “recycled box ... is
too wide for the slot.”  This is a geometric dead end, not an avalanche
criterion.  Before reusing a parked particle, rotate its vertices back to the
body frame and reset `Q`/`Qc` (see `AlignParkedBox` in the avalanche driver).
Do not reuse the old quaternion blindly.  Resetting only when the box is
parked is not sufficient: a 417-event rerun still found a parked box with a
non-identity quaternion at reuse.  Reset again immediately before measuring
the oriented width at insertion.

### 2c. Periodic boundaries in DEM

`DEM::Domain` supports periodic boxes through `Xmin/Xmax`, `Ymin/Ymax`,
`Zmin/Zmax`.  `Domain::Solve` computes `Per` from those bounds;
`Domain::ResetDisplacements` wraps free particles, and the linked-cell
neighbour search uses the same periodic flags.  `CInteracton` applies `Per`
only when both particles are free, so fixed walls are never periodic.  To make
the narrow `y` direction periodic, remove the two side planes, set
`Ymin=-0.5*cell_gap`, `Ymax=0.5*cell_gap` before `Solve`, and remove the
side-wall term from the insertion-width and lost-particle tests.  The fixed
base/feed walls keep their normal contact.

## 3. CPU Voronoi hot path: what must stay true

The exactness-critical rules are:

1. **Do not sort candidate lists.** Contact forces are accumulated in list
   order; sorting `L` changes the last bits even when the contacts are the same.
2. **Use the original `Distance` arithmetic for real contacts.** `Distance2`
   may be used only as an early rejection in the force loop. Replacing full
   contact distances with cached normals or squared-distance variants changes
   results.
3. **Keep the conservative Verlet cutoff:** `Rsum + 4.0*ContactDispGlobal()`.
   A factor of two was observed to skip a real contact that became a contact
   after motion. Do not lower the factor without a new derivation and an exact
   `Nproc=1` comparison.
4. **Do not change the original `(i,j)` contact push order.**
5. **Broadphase rejections must be conservative.** If a bound is not clearly
   safe, remove it rather than risk a missed contact.

The current optimizations that are compatible with those rules:

- per-thread force/torque accumulators plus an OpenMP reduction, instead of a
  particle lock per interaction;
- build-time feature distances stored beside each candidate list, used for
  rejection before the full distance calculation;
- safe AABB/plane broadphase in `Overlap`;
- fused old/new orientation passes in `Particle::Rotate`;
- particle-level rotation caches (`wrot = Rotation(w,Q)`, `Qc = Conjugate(Q)`);
- one friction-map lookup per interaction with a reference to the stored
  tangential displacement.

The primary CPU reference input is
`benchmarks/cpu_voronoi/tper_small.inp`; `tper_med.inp` is the larger timing
check.

## 4. Rotation caches: CPU and GPU

`Particle` owns:

```cpp
Vec3_t       wrot;        // Rotation(w,Q)
Quaternion_t Qc;          // Conjugate(Q)
Vec3_t       wrot_source; // w used for the cache
Quaternion_t Qc_source;   // Q used for the cache
```

The cache is refreshed by:

- `Particle::CalcProps()`;
- `Particle::InitializeVelocity()`;
- `Particle::Rotate(double dt)`, immediately after installing the new `Q`;
- `Domain::UpdateContacts()` and after a user `Setup` callback, using
  `UpdateRotCacheIfChanged()`.

The interaction force path reads `P1->wrot`, `P2->wrot`, `P1->Qc`, and
`P2->Qc`. There are no `Wrot1/Wrot2/Qc1/Qc2` interaction fields.

On CUDA:

- `UploadParticle` **copies** the host caches into `DynParticleCU`; it does not
  recompute them.
- `DnloadParticle` copies the device caches back and updates the source
  snapshots.
- `Domain::UpLoadDevice` calls `UpdateRotCacheIfChanged()` before packing, so a
  direct upload after a host edit cannot send stale values.
- The device `Rotate`/`MoveParticles` kernels refresh the device caches after
  changing `Q`/`w`.

If you write a custom GPU movement/rotation kernel, you own the device cache:
after changing `Q`/`w`, execute the equivalent of

```cpp
Rotation (DPar[ic].w,DPar[ic].Q,DPar[ic].wrot);
Conjugate(DPar[ic].Q,DPar[ic].Qc);
```

If host code assigns `Q` or `w` directly and then calls a force routine outside
the normal solve paths, call `Particle::UpdateRotCache()` first.

## 5. CUDA DEM: design and limits

- `DynParticleCU` stays AoS. The particle, vertex, and interaction arrays are
  AoS. A full SoA rewrite was considered and not adopted because the force
  kernels consume several fields of the same particle and the vertex records
  are compact `real3` values.
- The default force dispatch uses the cached `CalcForceEEPre`, `CalcForceVFPre`,
  and `CalcForceFVPre` kernels. `pForceVV` remains on its normal path. Users can
  restore the original kernels or supply their own through the function
  pointers.
- `pMove = MoveParticles` is the default combined translation/rotation path.
  `pMove = NULL` keeps the separate `pTranslate` + `pRotate` fallback. Keep that
  fallback: making fused movement unconditional broke `test_cu_04` for custom
  movement.
- `MaxD` writes one partial maximum per block; the host reduces only those
  block entries.
- The outer solve loop owns `Time`/`iter`. CUDA receives the authoritative
  values through `SetClock` once per step. Do not add independent clock
  increments inside DEM or LBM substeps.
- `Nthread=512` fails for the current large `test_cu_03` register usage. The
  tested valid range is up to 384; the large workload optimum is around 128,
  while the reduced workload optimum is around 32.

CUDA `test_cu_03` and `test_cu_04` use `atomicAdd`, so the baseline itself is
not bitwise deterministic. Do not claim bitwise equality for those GPU runs.
Use energy/trend checks, a deterministic reduced configuration, or a host-side
inspection of the copied fields.

## 6. LBM-DEM and clocks

`lib/lbmdem/Domain.h` couples `LBMDOM` and `DEMDOM`. The historical failure was
that LBM and DEM advanced `iter`/`Time` independently, so they could disagree.

Current rule: the enclosing solve loop is the single authority. Advance `Time`
and `iter` once per outer step and publish the same values to CUDA through
`SetClock`. Do not increment them in the DEM or LBM substep.

The LBM-DEM tests are large. Build them as integration checks; use a temporary
reduced grid/time copy for a smoke run rather than modifying the checked-in
test just to make it fast.

For GPU FLBM, note two separate conventions:

- the streaming kernel resets `BForce` every step, so a constant body force must
  be rewritten on the device in the `Setup` callback;
- the raw GPU second moment carries a `-0.5*nu*phi` bias in permeability;
  `u_physical = u_raw + 0.5*dt*fx*phi`.

## 7. How to prove a CPU change is exact

1. Recreate a clean baseline from the chosen commit in a separate worktree.
2. Run baseline and modified programs with `Nproc=1`.
3. Compare **HDF5 datasets**, not whole-file MD5. Container metadata is not
   deterministic even between two baseline runs.
4. Expect `neq []` and identical printed energies.
5. Then run the relevant `Nproc=8` case for timing; small thread-order
   differences are expected there because the original force accumulation is
   thread-dependent.

Useful pattern:

```python
import h5py, numpy as np
b = h5py.File("baseline.h5", "r")
o = h5py.File("optimized.h5", "r")
for k in b:
    if isinstance(b[k], h5py.Dataset) and not np.array_equal(b[k][()], o[k][()]):
        print("DIFF", k)
```

Do not compare the modified working tree against itself. Also keep output names
separate (`_a`, `_b`, `_bf_`) so a later run does not overwrite the baseline.

## 8. Hard-won traps

- **HDF5 whole-file hashes are useless for exactness.** Compare datasets.
- **`Array::Find()` returns `-1` on failure.** Assigning it to `size_t` is a
  real bug; check the signed result first.
- **Feature containers often hold pointers** (`Array<Edge*>`, not values).
  `std::is_same`-based type traits can silently misclassify them.
- **Planar face bounds** once rejected valid periodic contacts and changed
  energy from `0.0397229` to `0.0416998`. They were reverted; keep only the
  conservative plane test.
- **An invalid edge-edge lower bound** reduced the edge-edge contact count and
  changed energy. It was removed.
- **`TheFileKey == NULL`** must not be passed to `String::Printf("%s", ...)`;
  the solve path handles it.
- **Changing `MaxD` or the clock API affects LBM-DEM.** Update
  `lib/lbmdem/Domain.h` when shared APIs change.
- **Guard device-only code with `#ifdef USE_CUDA`.** A cache-init launch outside
  the guard broke CPU `tlbmdem01/02`.
- **Preserve `tdem/test_cu_04.cu` CRLF endings.**
- **Do not delete uncommitted optimization changes or generated benchmark
  inputs** while experimenting.
- **An overhanging grain at the open end is not a wall leak.**  The avalanche
  driver uses `exit_margin=0.010 m`, so a box can remain active until its
  center is 1 cm past `cell_Lx`.  Check the side walls (`|y|`, `z`) separately
  before calling the run faulty.

## 9. Benchmarks: what the numbers mean

The reduced CPU Voronoi case is 480 particles, `L=10`, `n=10`,
`T0=Tf=0.01`, `dt=5e-5`, 400 total steps. It is useful for exactness and for
relative CPU measurements, but construction and output dominate at high thread
counts.

Representative results from this optimization work:

| Case | Result |
|---|---:|
| Reduced CPU, `Nproc=1`, solve | 12.73 s -> 8.06 s |
| Reduced CPU, `Nproc=8`, solve | 1.99 s -> 1.36 s |
| Medium CPU, 4k cells, `Nproc=8`, wall | 26.06 s -> 21.64 s |
| Medium CPU, solve phase | 16.40 s -> 12.00 s |
| Large CUDA `test_cu_03` | ~99.06 -> ~84.51 ms/step |

The CUDA figure predates the final struct/dispatch refactor; re-benchmark it
before quoting it as a current number. The medium absolute times vary with
machine load. The particle-level rotation
cache follow-up was small in wall time (about 1% on the medium two-solve total,
about 2% on the second solve); its main value is removing per-interaction work
and matching the GPU cache lifecycle.

For the reduced CUDA test, the thread optimum is about 32 and the useful Verlet
distance is around `0.01-0.02`, not the large-case `0.05`/128. Always state the
workload, thread count, Verlet distance, and whether the number is wall time or
solve time.

## 10. Short verification checklist

After changing `lib/dem/*`:

```bash
cd tdem
make -j8 test_periodic test_cu_01 test_cu_02 test_cu_03 test_cu_04 test_dynamics

cd ../tlbmdem
make -j8 tlbmdem01 tlbmdem02 tlbmdem_cu_01

# CPU exactness: clean HEAD baseline, Nproc=1, compare HDF5 datasets
cd ../tdem
./test_periodic /path/to/input_prefix 1

# For a CUDA smoke test, copy test_cu_03.cu to a temporary file, reduce Tf
# (for example 10.0 -> 1.0e-4), build the copy, and run it. The checked-in
# full test is long; do not treat a shortened run as a full accuracy run.
```

Then state separately:

- what compiled;
- what ran;
- what matched the baseline;
- what was faster, with the measured input and thread count.

That distinction prevents most misleading optimization reports.

## 11. Memory maintenance duty

This repository is intentionally a long-running notebook. At the end of every
user prompt, before the final response, check whether the prompt produced a
durable fact. If so, update the owning file:

| Kind of knowledge | File |
|---|---|
| short rule every session must see | `AGENTS.md` |
| architecture, command, exactness, trap | `MECHSYS_HANDBOOK.md` |
| distilled experience, checklist | `MEMORY.md` |
| design or benchmark result | the relevant optimization report |
| task-triggered guidance | `SKILL.md` and `.agents/skills/mechsys/SKILL.md` |

The two skill copies must stay in sync. Record only durable, verified facts;
exact commands, inputs, thread counts, and time windows belong with any
measurement. Do not add transient `/tmp` logs or guesses. If nothing durable
was learned, make no edit and say so in the final response. Run
`git diff --check`; never commit or push unless explicitly asked.
