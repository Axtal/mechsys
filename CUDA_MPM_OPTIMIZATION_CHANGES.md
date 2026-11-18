# CUDA MPM optimisation changes

Date: 2026-11-04
Tree: `/home/sergio/mechsys` (uncommitted working-tree changes)
GPU: NVIDIA RTX A6000, CUDA 12.4, driver as installed on this host

This note covers the GPU (CUDA) MPM path only: `lib/mpm/mpm.cuh`,
`lib/mpm/Domain.h` (the `USE_CUDA` parts) and `lib/mpm/Particle.h`.
The CPU MPM work is documented separately in `MPM_OPTIMIZATION_CHANGES.md`.

## Result

All timings are the `Elapsed time` printed by `Domain::Solve` (seconds,
`Nthread=256`, four HDF5 frames per run).  The wall time includes process
start-up and output.  The reduced sources are the checked-in
`tmpm/tmpm_cu_0*.cu` cases with a smaller final time, so the relative
change can be measured in minutes rather than hours.

| case | baseline Solve | new Solve | speed-up | baseline wall | new wall |
|---|---:|---:|---:|---:|---:|
| tmpm_cu_01 pre-CPDI (5 k particles, 125 k nodes, 345 steps) | 0.611 | 0.583 | 1.05x | 0.98 | 0.94 |
| tmpm_cu_02 (25 k active, 43.7 k nodes, 7 735 steps) | 24.545 | 24.859 | 0.99x | 24.69 | 24.98 |
| tmpm_cu_03 (135 k active, 358 k nodes, 1 035 steps) | 25.048 | 23.372 | 1.07x | 25.31 | 23.65 |
| tmpm_cu_04 (216 k particles, 38.9 k nodes, 849 steps) | 34.360 | 31.408 | 1.09x | 34.62 | 31.67 |
| tmpm_cu_05 (25 k active, 43.7 k nodes, 587 steps) | 2.850 | 2.601 | 1.10x | 2.98 | 2.78 |

`tmpm_cu_02` is the only regression (about 1%).  It has many short steps
and frequent particle deletion/re-upload; the other cases gain 4-10%.
The 5 k-particle pre-CPDI case is dominated by per-step launch and output
cost, so its 4% is close to the noise floor.  The checked-in case has since
been switched to the CPDI tetrahedral beam, so these `tmpm_cu_01` timings are
historical.

## 2026-11-06 CPDI hot-path pass

This pass keeps the CPDI algorithm and launch order unchanged and removes
work that was being repeated per particle-node contribution.  The checked-in
CPDI case (`ndiv=125`, 44,032 particles, 238,328 nodes, 9,151 corners, 229
boundary tip corners) improved from **17m23.33s to 4m02.85s** (4.30x).  The
reduced `ndiv=50` case improved from **47.53s to 26.74s** (1.78x):

| CPDI case | particles | grid nodes | before Solve | after Solve | speed-up |
|---|---:|---:|---:|---:|---:|
| `ndiv=50` | 3,315 | 17,576 | 47.53 s | 26.74 s | 1.78x |
| `ndiv=125` (checked-in) | 44,032 | 238,328 | 1043.33 s | 242.85 s | 4.30x |

The final CPDI deflection is unchanged: `ndiv=125` tip `9.70156508e-4 m`,
analytical `1.00081114e-3 m`, ratio `0.9693702`, MAE `1.2168888e-5 m`.
The optimized `.res` differs from the pre-optimization one by at most
`1e-16 m`; `ndiv=50` is bitwise identical.  HDF5 fields are not bitwise
identical because the normalized-gradient cache changes the last few bits:
at `ndiv=125` the largest relative differences are about `2.3e-8` in
`Sigma`/`MSigma`, `6.9e-8` in `PVelocity`, and `1.9e-8` in
`Strain`/`MStrain`.  These are round-off-level and do not change the final
deflection.

Changes:

- `cudaComputeAn` precomputes `An[n] = Fn[n]/Mass[n]` once per step.  The
  gathers previously divided by `Mass` for every particle-node contribution.
- `ParticleCU::GradvecN[4]` stores `Gradvec/(6*Vt)`; it is filled by
  `CalcVolCPI` and `UploadParticle`.  Scatter and gather use it directly
  instead of doing a division for each of up to 32 contributions.
- Particle `v`, `L`, `x`, corner `x`, and scatter `v/m` are kept in local
  registers inside the kernels instead of being read/written through the
  global `ParticleCU`/`CornerCU` pointer each iteration.
- CPDI particle-to-node no longer rebuilds the touched-node list: the preceding
  `cudaCornerToNode` already visits every node the particle scatter can touch.
  The redundant `atomicExch`/`atomicAdd` pair and `SelIndexes` write were
  removed from `cudaParticleToNodeCPDI`.
- `cudaCornerToNode` skips the six momentum/force atomics when the corner
  force is exactly zero, and skips all node atomics for zero shape values.
  The `Nodes`/`Sfunc` slots are still filled for the particle kernels.
- The non-CPDI path has a separate pass, described below.

## 2026-11-06 non-CPDI pass

The default `tmpm_cu_01` path (5,000 GIMP particles, 125,000 nodes, 86,393
steps) improved from **1m54.82s to 21.03s** (**5.46x**):

| non-CPDI case | particles | grid nodes | before `Solve` | after `Solve` | speed-up |
|---|---:|---:|---:|---:|---:|
| `ndiv=50`, GIMP | 5,000 | 125,000 | 114.82 s | 21.03 s | 5.46x |

Profiling showed the old particle scatter at about `108s` of the run; the
separable shape evaluation cut that phase to about `17s`.  The final `.res`
matches the pre-optimization GPU result to at most `1e-13 m` (atomic ordering
prevents a general bitwise promise):
tip `9.34560567e-4 m`, analytical `9.71136179e-4 m`, ratio `0.9623373`, MAE
`1.5385435e-5 m`.  HDF5 fields are not bitwise identical because the
precomputed nodal acceleration changes accumulation order: the largest relative
differences are about `7.4e-8` in `PVelocity`, `2.0e-9` in `Sigma`, and
`1.4e-9` in `Strain`.

Changes:

- GIMP/B-spline are separable.  `cudaParticleToNode` now evaluates the 1-D
  factors once per axis (`MPMShapeFunction1DCU`, normally 3-5 values) and
  combines them, instead of calling `MPMShapeFunction3DCU` 3 times per node.
  For the standard 3x3x3 support this is 9 1-D evaluations instead of 81.
  A direct 3-D fallback remains for support wider than 8 nodes per axis.
- `cudaComputeAn` now also feeds the non-CPDI gather, so `Fn/Mass` is computed
  once per touched node rather than once per particle-node contribution.
- The non-CPDI gather keeps `v`, `x` and `L` in registers and writes the
  particle back once; the scatter keeps `v/m` in registers.
- The corner/CPDI kernels and launch order are unchanged by this pass.

The phase profile that motivated the pass (reduced `ndiv=50`, before this
pass) was: reset `1.04s`, corner scatter `3.83s`, particle scatter `17.11s`,
nodal solve `0.85s`, corner motion `1.77s`, gather `30.11s`.  The gather was
the largest cost and was division-bound; caching `An` cut it from `30.1s` to
`18.0s`, and normalized gradients cut the particle scatter/gather pair to
about `10.1s + 11.1s` before the final `GradvecN` cache.

## How the time was measured

`nsys` profiling is unavailable on this host (`importer binary ... not
found`) and `ncu` is blocked by `ERR_NVGPUCTRPERM` (no permission to read
GPU performance counters).  The phase split was therefore measured with
host `std::chrono` timers around each launch plus
`cudaDeviceSynchronize()`, compiled only in a profiling copy
(`-DMPM_GPU_PROFILE`); nothing profiling-related is committed.

Profiling copy of the reduced `tmpm_cu_03` (seconds over the whole run):

| phase | before | after |
|---|---:|---:|
| reset previous valid nodes | 0.204 | 0.174 |
| particle -> node scatter | 9.262 | 8.492 |
| mask fill + compaction + count copy | 0.038 | removed |
| nodal solve | 0.033 | 0.026 |
| node -> particle gather | 10.140 | 9.927 |
| sum of GPU phases | 19.677 | 18.619 |

The remainders are the HDF5 output frames (`DnLoadDevice` copies the
whole `ParticleCU` array), the `Setup` callback and process overhead.
That is why the reduced runs show a smaller relative gain than the GPU
phases alone.

## Changes

1. **Hot node data in structure-of-arrays form.**
   `NodeCU` is a large cold struct (208 B, including the unused `Sn`).
   Scatter and gather now use separate device arrays for the hot fields:
   `Mass`, `Mn`, `Fn`, `Vn` and the `valid` flag.  A warp reading
   neighbouring nodes now reads consecutive addresses instead of a
   208-byte stride.  The cold `NodeCU` array is still used for the
   fixed-velocity flags and `Vnf`; its hot fields are no longer read by
   the kernels.

2. **Fused selection-list construction.**
   The old sequence was, every step:
   `cudaResetNode` over the previous list, `thrust::fill` over all
   `Nnodes`, scatter, `thrust::copy_if` over all `Nnodes`,
   `thrust::distance` (host synchronisation), `cudaMemcpy` of the count,
   then the solve and gather.
   The scatter now appends each newly touched node itself
   (`atomicExch` on the per-node `valid` flag plus `atomicAdd` on a device
   counter).  The previous list is reset at the start of the next step and
   the counter is ping-ponged (`pValidNumber` / `pValidNumberPrev`), so
   there is no mask fill, no compaction and no per-step host copy.
   `bAllIndexes`, `bValidMask` and `pValidMask` are gone.

3. **Reset and solve over the previous list.**
   `cudaResetNode` resets only the nodes touched in the previous step
   (`O(valid)` rather than `O(Nnodes)`).  Both reset and solve are
   launched with `Nnodes` blocks and exit early on the device counter, so
   the host never has to know the valid count.  A final
   `cudaDeviceSynchronize()` at the end of the time loop makes the
   reported `Solve` time cover all queued work (and prevents the process
   from appearing to finish early).

4. **Smaller `ParticleCU`.**
   The per-particle node list/shape arrays were reduced from 125 entries
   to 64 entries (`4^3`, the maximum for cubic B-spline; quadratic and
   GIMP use at most 4 nodes per dimension for the usual particle size).
   `sizeof(ParticleCU)` drops from 6024 B to 3584 B, which reduces the
   per-frame host/device copy and improves cache behaviour in scatter and
   gather.  This is the change with the largest effect.  The CPDI fields
   below add 128 B, so the current size is 3712 B; the optimised grid path
   still benefits from the 64-entry reduction.

5. **CPDI/body-mesh path on the GPU.**
   Under `USE_CUDA`, `Domain::Solve` previously had an empty
   `if (BodyMesh!=NULL)` branch.  It now runs the CPDI update on the device
   in the same sequence as `Domain::OneStepCPI`:

   | stage | kernel |
   |---|---|
   | reset previously touched nodes | `cudaResetNode` |
   | corner -> node scatter | `cudaCornerToNode` |
   | particle -> node scatter | `cudaParticleToNodeCPDI` |
   | nodal solve | `cudaSolveNode` |
   | move corners | `cudaMoveCorners` |
   | node -> particle gather | `cudaNodeToParticleCPDI` |

   `ParticleCU` now also carries `Corners[4]`, `Gradvec[4]` and `Vt`.
   `ParticleCU::CalcVolCPI` mirrors `Particle::CalcVolCPI` using the four
   device corner positions.  A new `CornerCU` mirrors `MPM::Corner` and is
   uploaded as `bCornersCU`; `UpLoadDevice` copies its position,
   previous position, velocity, fixed velocity, hydraulic force and fixed
   flags.  Corner `Nodes` and `Sfunc` are computed on the device every step.
   `DnLoadDevice` copies particles and corners back to the host at output
   frames and calls `UpdateMesh()` so the body mesh bounds/normals match the
   downloaded corner positions.  The checked-in `tmpm_cu_01.cu` mirrors
   `tmpm01.cpp`: the standard non-CPDI path is the default; pass `cpdi` to select
   `AddRectangularBeamMesh` + `ResizeDomainMesh(...,1.0)`.  In CPDI
   the tip of the rectangular beam is an external face, so the host selects
   corners by their `x` coordinate alone (the last mesh slice), fixes base
   corners, stores only the tip corner indexes in a device `ForceCor` vector,
   and applies the load with `SetCornerForce`; `ForcePar` is empty.  No
   `BFac` boundary mask or setup-time `BoundaryMesh()` call is needed.  This
   shortcut is valid because the beam tip is a single external cross-section;
   a general body mesh with possible interior vertices in the last slice
   would still need a boundary/edge test.  The non-CPDI path
   applies `SetForce` to `ForcePar` and leaves `ForceCor` empty.

   The zero-shape case is handled explicitly: when a corner lies exactly on
   a grid node, some of its eight shape values are exactly zero.  Both the
   CPU and CUDA CPDI scatters/gathers now skip those contributions instead
   of dividing by a zero nodal mass.  Nonzero contributions are unchanged.

   This GPU CPDI path is for standalone `MPM::Domain::Solve`.  The LBMMPM
   coupling still calls the CPU `MPMDOM.OneStepCPI` and has its own device
   corner arrays in `lib/lbmmpm`; see `CUDA_LBMMPM_STRUCTURE.txt`.

## Correctness

The reduced cases were run with the same inputs and the four output HDF5
frames were compared field by field against the pre-change baseline.  The
GPU scatter uses floating-point atomics, so bitwise identity is not
guaranteed in general.  Observed differences:

- `tmpm_cu_02`: all frames bitwise identical.
- pre-CPDI `tmpm_cu_01`: all frames bitwise identical on the final run (an
  earlier run had a 9.1e-13 `Strain` difference, i.e. run-to-run atomic
  noise).
- `tmpm_cu_03`: `Pressure`/`Sigma` up to 3.1e-5 absolute (about 1e-8
  relative), `PVelocity` up to 4.7e-10.
- `tmpm_cu_04`: `Sigma` up to 8e-18 absolute, `PVelocity` up to 5e-20.
- `tmpm_cu_05`: `Pressure`/`Sigma` up to 9e-10 absolute, `PVelocity`
  up to 2e-12.

These are at the level expected from different atomic summation order;
no systematic drift was observed.  All five full-size `tmpm/tmpm_cu_*.cu`
sources compile with the modified headers.

### CPDI exactness check

A reduced CPDI case was built in both CPU and CUDA form with
`AddRectangularBeamMesh` (44 particles, 22 corners) on a 0.237 m grid, then
run for 20 steps with five output frames.  All HDF5 dataset payloads
(`Position`, `PVelocity`, `Sigma`, `Strain`, `Pressure`, `Con`, `Vertex`,
`PTag`, `MTag`, `MSigma`, `MStrain`) were bitwise identical in every frame.
The case is small, so it checks the algorithm/data transfer, not GPU
speed-up.

The checked-in `tmpm01`/`tmpm_cu_01` pair supports an optional CPDI path
(`cpdi` argument; `AddRectangularBeamMesh`, 44,032 particles, 238,328 grid
nodes, 9,151 corners at `ndiv=125`) and defaults to non-CPDI.
The CPU/GPU parity smoke used `Fb=0.1`; the checked-in tip force is now
`Fb=5e-2`.  With `Tf` reduced from 5 to 0.5 (8,640 steps, 10 output frames,
`Nproc=1`), all HDF5 datasets matched between CPU and CUDA except one
`PVelocity` value at frame 7 (`4.4e-16` absolute).  `Solve` took 33.88 s on
the CPU and 5.44 s on the GPU in that smoke window; startup/output dominate
the wall time.  Full `Tf=5` half-force GPU run: tip `1.14806452e-3 m`,
analytical `9.75102376e-4 m`, ratio `1.1774`; halving force halves deflection
but does not change the relative error.

CPDI body-mesh refinement with the current boundary-corner load (full
`Tf=5`, `Fb=0.05`, `beta=1`, GPU) gives:

| `ndiv` | particles | grid nodes | `Dx` (m) | tip/analytic | MAE (m) |
|---:|---:|---:|---:|---:|---:|
| 50 | 3,315 | 17,576 | 7.551e-3 | 0.8250 | 6.842e-5 |
| 75 | 10,045 | 50,653 | 5.266e-3 | 0.8878 | 4.554e-5 |
| 100 | 23,180 | 125,000 | 3.962e-3 | 0.9486 | 2.065e-5 |
| 125 | 44,032 | 238,328 | 3.204e-3 | 0.9694 | 1.217e-5 |

Refinement moves the tip toward the analytical value and cuts MAE about 5.6x
from `ndiv=50` to 125.  The checked-in CPDI case uses `ndiv=125`; non-CPDI
remains `ndiv=50`.  GPU `Solve` is about 17m23s at 125; body-mesh work
dominates.  An earlier `Tf=2` test (ratios 1.144/1.292/1.440) used the old
combined particle+corner load and does not describe this path.

`beta` in `ResizeDomainMesh` is `Dx/dmax`: below 1 is finer, above 1 is
coarser.  At the current `ndiv=125`, `Fb=0.05`, full `Tf=5`:

| `beta` | grid nodes | `Dx` (m) | tip/analytic | MAE (m) | GPU `Solve` |
|---:|---:|---:|---:|---:|---:|
| 1.00 | 238,328 | 3.204e-3 | 0.9694 | 1.217e-5 | 17m23s |
| 1.25 | 117,649 | 4.005e-3 | 1.1652 | 8.878e-5 | 17m15s |
| 1.50 | 68,921 | 4.806e-3 | 1.1392 | 7.653e-5 | 17m17s |

`beta=1.0` is clearly best; coarsening overshoots.  The old `beta=1.5`
observation came from the old combined particle+corner load and does not
transfer to this path.  Runtime is body-dominated, so beta does not buy time
here.

`tmpm01`/`tmpm_cu_01` now accept `[Nproc] [cpdi|mpm]`; the default is non-CPDI.
CPDI uses `ndiv=125` (`229` boundary tip corners) and leaves `ForcePar` empty.
Non-CPDI keeps `ndiv=50`, selects `ForcePar`, and leaves `ForceCor` empty.

At `Fb=0.05`, full `Tf=5`:

| Path / load | Particles | Grid nodes | Tip (m) | Tip/analytic | MAE (m) |
|---|---:|---:|---:|---:|---:|
| Non-CPDI, particles, `ndiv=50` | 5,000 | 125,000 | 9.3456e-4 | 0.9623 | 1.539e-5 |
| CPDI, boundary corners, `ndiv=50` | 3,315 | 17,576 | 8.0449e-4 | 0.8250 | 6.842e-5 |
| CPDI, boundary corners, `ndiv=75` | 10,045 | 50,653 | 8.7912e-4 | 0.8878 | 4.554e-5 |
| CPDI, boundary corners, `ndiv=100` | 23,180 | 125,000 | 9.4365e-4 | 0.9486 | 2.065e-5 |
| CPDI, boundary corners, `ndiv=125` (checked-in CPDI) | 44,032 | 238,328 | 9.7016e-4 | 0.9694 | 1.217e-5 |
| CPDI, old particles + all tip corners, `ndiv=50` | 3,315 | 17,576 | 1.1481e-3 | 1.1774 | 6.191e-5 |

The optimized non-CPDI final GPU `.res` matches the previous GPU baseline to
`1e-13 m`.  CPU/GPU HDF5 payloads are no longer promised bitwise because
the separable-shape and `An` caches change the last bits; the observed
differences are at the `1e-8` relative level in `PVelocity`.

## 2026-11-06 `pValidNumberPrev` check

`pValidNumberPrev` looks redundant: reset only needs the count of the list
built in the previous step, and that count is still in `pValidNumber` until
the memset.  I tested the single-counter version by resetting with
`pValidNumber`, zeroing it, scattering into it, and deleting the swap.  It is
**not a safe drop-in removal**:

- CPDI `ndiv=50` stayed bitwise identical.
- Non-CPDI `tmpm_cu_01` was deterministic across reruns but diverged from the
  baseline `.res` from output row 80 onward; the final maximum difference was
  `8.74e-3 m`.
- Replacing `cudaMemsetAsync` with synchronous `cudaMemset`, adding an explicit
  `cudaDeviceSynchronize()` after reset, and restoring a dummy second
  allocation all reproduced the same divergence, so this is not a memset
  ordering or pointer-layout effect.
- Removing the second counter did not improve the 21.00 s non-CPDI run.

Keep the two-counter ping-pong.  It costs one `cudaMalloc` and a host pointer
swap, and it defines which list the next sparse reset clears.  Any future
single-buffer redesign must be validated against the full non-CPDI case, not
only CPDI.

## Caveats

- `ParticleCU` now holds at most 64 nodes per particle.  Cubic B-spline
  needs exactly 64; quadratic and GIMP at the standard particle size use
  27.  A GIMP particle significantly larger than a cell (roughly
  `|Dx| > 1.5 * grid Dx`) could need 5 nodes per dimension (125).  Such a
  configuration must either use a smaller particle size or widen the
  arrays back to 125 (at the cost of the 64-entry gain).
- The reported numbers are for four output frames.  Full production runs
  with many frames pay the `DnLoadDevice` copy of `ParticleCU`, which is
  now smaller but still proportional to the particle count.
- `tmpm_cu_02` is about 1% slower; the other measured cases improve.

## What was tried and did not help

- packed `NodeHot` arrays (AoS or 64-byte aligned) for gather: slower;
- an early `Fn/Mass` precompute attempt without the other CPDI changes was
  neutral/slower; the current separate `cudaComputeAn` pays for itself now
  that the gathers no longer divide per contribution and the normalized
  gradients/local caches are in place;
- scalar/separable shape-function cache in scatter: neutral;
- splitting gather and constitutive update into two kernels: neutral
  (after the update work was actually included);
- `-maxrregcount` 64/80/96/112 (registers 128/118 in the original
  kernels): no gain or slower;
- block sizes 128/512 instead of 256: within noise (512 slightly better
  on some cases, worse on the 5 k-particle case).

## Reproducing

Reduced sources and the baseline binaries used here live in
`/tmp/gpu_mpm` (not part of the repository).  To rebuild a full case:

```sh
cd /home/sergio/mechsys
./mechsyscc tmpm/tmpm_cu_04.cu -o /tmp/tmpm_cu_04
/tmp/tmpm_cu_04 96
```

The first argument is the OpenMP thread count used by the host parts; the
CUDA path uses `Domain::Nthread` (256 by default) for the GPU blocks.
