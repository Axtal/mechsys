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
