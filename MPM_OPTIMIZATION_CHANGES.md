# MPM CPU optimization

This document describes the CPU optimizations made to the material point
method in `lib/mpm/`, using the `tmpm` benchmarks as reference.  Only the
grid-based CPU transfer (`Domain::ParticleToNode` / `Domain::NodeToParticle`)
is changed; the CUDA kernels, the CPDI/body-mesh path and the constitutive
models are left untouched.

The `Nproc=1` results are dataset-bit-identical to the pre-optimization code
for every case checked (GIMP `tmpm01`, Drucker-Prager `tmpm02`, Tait water
`tmpm03` and mu(I) `tmpm05`).  At `Nproc>1` the momentum sums may be evaluated
in a different order, so last-bit differences are possible; the largest
measured deviation is reported in Section 3.

## 1. Reference environment and workloads

- Host: 2x Intel Xeon Gold 6430 (64 physical cores, 128 logical CPUs, 2 NUMA
  sockets), 62 GB `/tmp` tmpfs.
- Compiler: the normal repository build (`g++ -O3 -fopenmp`, GCC 15).
- All timings below are wall-clock `Domain::Solve` times (the stopwatch
  printed by `Solve`), except the full `tmpm01`/`tmpm04` rows, which are
  `/usr/bin/time` wall times including HDF5 output.

| Workload | Particles | Nodes | Steps | Notes |
|---|---:|---:|---:|---|
| W1 column collapse (`tmpm03`, 20 ppd) | 10,740 | 39,304 | 998 | no output |
| W2 column collapse (`tmpm03`, 50 ppd) | 49,400 | 140,608 | 150 | no output |
| W3 `tmpm01` full (GIMP) | 5,000 | 125,000 | 86,391 | 100 HDF5 frames |
| W4 `tmpm04` full (Tait water, duct) | 216,000 | 38,880 | 42,437 | 100 HDF5 frames |

W1 and W2 are reduced copies of the checked-in case built with
`mechsyscc` in a scratch directory; the only edits are
`particles_per_diameter` and `Tf`, and the call
`Solve(Tf,dt,Tf,NULL,NULL,NULL,false,Nproc)` to remove the output from the
timing.  W3 and W4 are the checked-in programs,
so their numbers include the usual report output.

### Profile of the optimized code

Timers around `ParticleToNode` and `NodeToParticle` (the latter includes the
node update, the particle gather and `CalcVol`/`CalcStress`/`Reset`), W1 at
1/8/32/96 threads, seconds of `Solve`:

| Threads | `ParticleToNode` | `NodeToParticle` | Sum | `Solve` |
|---:|---:|---:|---:|---:|
| 1 | 3.19 | 2.43 | 5.62 | 5.62 |
| 8 | 1.73 | 0.46 | 2.19 | 2.19 |
| 32 | 0.72 | 0.18 | 0.90 | 0.90 |
| 96 | 0.69 | 0.24 | 0.93 | 0.93 |

The scatter is 74-80% of the step at 32-96 threads and stops scaling between
32 and 96 threads (0.72 s -> 0.69 s on W1; 0.350 s -> 0.320 s on W2).  The
gather is 20-26% and also flattens (W1 0.18 s at 32 vs 0.24 s at 96; W2
0.114 s vs 0.079 s).  This is why the changes concentrate on the scatter and
why dense grids remain limited.

## 2. Changes

| Area | Change |
|---|---|
| Node layout | `Node` is now `alignas(64)` with `Mass`, `Mn`, `Fn`, `valid` and the OpenMP lock in the first 64-byte line.  A scatter update touches one cache line instead of four, and adjacent nodes no longer share a line |
| Node reset | `Reset` no longer clears `Vn` (always overwritten by `UpdateNode`) or `Sn` (never read by the CPU transfer) |
| Valid-node list | `VNodes` is built during the scatter from the first thread that touches each node, instead of a serial scan of all `Nnodes` after the transfer |
| Shape functions | Per-particle separable 1-D weights (`FillShapeCache3D`) replace a full 3-D GIMP/B-spline evaluation per node.  The multiplication order is unchanged, so the cached path is bitwise identical |
| NaN check | Moved out of the inner particle-node loop (it depended only on the particle, not on the node) |
| Single-thread path | The node locks are skipped when `Nproc==1`, where they are pure overhead |
| NUMA first touch | The node array is allocated with `posix_memalign` and constructed in an OpenMP loop, so each page is first touched by the thread that will later update it |

## 3. Exactness

Reduced runs with HDF5 output were compared dataset by dataset against the
pre-optimization library built from the same sources.

| Case | Model | Threads | Result |
|---|---|---:|---|
| W1 (`tmpm03`, 20 ppd) | Tait/Cole water | 1, 8, 96 | bitwise identical |
| `tmpm01` reduced | GIMP elastic | 1 | bitwise identical |
| `tmpm02` reduced | Drucker-Prager | 1 | bitwise identical |
| `tmpm05` reduced | mu(I) rheology | 1 | bitwise identical |
| W3 full `tmpm01` | GIMP elastic | 96 | `Position`, `Pressure`, `Sigma`, `Strain`, `PTag` bitwise; `PVelocity` maxabs `2.8e-14` on `3.7e-4` (7.6e-11 relative) |

The `Nproc=1` identity is expected: the cached shape values use the same
arithmetic, the node sums are accumulated in particle order, and `VNodes`
order does not affect any per-node computation.  At `Nproc>1` only the order
of the floating-point reductions can change.

The CUDA path is unchanged: `mpm.cuh`, `NodeCU` and the `tmpm_cu_*` drivers
were not edited.  As a check that the host `Node` layout/reset change does not
disturb it, the reduced `tmpm_cu_01` (`Tf = 0.02`, 70 frames) was rebuilt with
the new headers and run on the RTX A6000: 0.41 s wall, 0.078 s `Solve`, all
finite.  Its final HDF5 frame is bitwise identical to the pre-change GPU run
(`PTag`, `PVelocity`, `Position`, `Pressure`, `Sigma`, `Strain`).  That is one
reduced case; the GPU kernels use atomics in general, so this is not a general
bitwise claim for GPU MPM.

## 4. Benchmarks

W1 (10,740 particles, 998 steps), seconds of `Solve`:

| Threads | Before | After | Speedup |
|---:|---:|---:|---:|
| 1 | 11.03 | 5.58 | 1.98x |
| 8 | 4.09 | 2.16 | 1.90x |
| 32 | 2.37 | 0.87 | 2.73x |
| 96 | 2.03 | 0.89 | 2.29x |

W2 (49,400 particles, 150 steps), seconds of `Solve`:

| Threads | Before | After | Speedup |
|---:|---:|---:|---:|
| 1 | 8.86 | 4.95 | 1.79x |
| 8 | 1.98 | 1.24 | 1.60x |
| 32 | 1.09 | 0.45 | 2.42x |
| 96 | 0.79 | 0.40 | 1.98x |

W3 full `tmpm01` (5,000 particles, 86,391 steps, 100 frames, 96 threads):

| Before | After | Speedup |
|---:|---:|---:|
| 92.42 s | 33.77 s | 2.74x |

W4 full `tmpm04` (216,000 particles, 42,437 steps, 100 frames, 96 threads):

| Before | After | Speedup |
|---:|---:|---:|
| 1932.08 s | 1676.27 s | 1.15x |

W4 is the worst case: its grid has only 38,880 nodes for 216,000 particles
(about 5.5 particles per cell), so every node receives ~150 locked updates
per step and the scatter is dominated by cache-line ownership transfers
between threads.  The serial and medium-thread gains are still present, but
Amdahl's law limits the 96-thread gain.  A short 5-frame W4 run (2,120 steps)
takes 75.20 s before and 71.32 s after; that run shows the compute part of the
gain without the long-run NUMA effects.

## 5. What did not help

| Attempt | Measurement |
|---|---|
| Caching `Fn/Mass` per node in an extra `An` field after `FixNode` | no change at 1, 8, 96 threads; reverted |
| Removing the node locks at `Nproc>1` (diagnostic, incorrect) | slower: the race lets several threads RMW the same line repeatedly, while the lock makes one thread own the line for all seven writes |
| `numactl --interleave=all` after adding the parallel node allocation | no further gain; the allocation already distributes the pages |
| Eight per-group partial node grids for the dense scatter (prototype) | short W4 at 96 threads 71.3 s -> 60.8 s, but sparse W2 at 32 threads 0.45 s -> 0.60 s; with a branch that keeps the old path below 64 threads the net W4 gain was only 5-7%, so it was not kept |
| More than eight partial groups | 4 groups 66.4 s, 8 groups 60.8 s, 16 groups 70.5 s on short W4 |

The dense-limit experiment also showed that fewer threads can be faster when
the grid is much denser than 1 particle per cell: on short W4, 8/32/96
threads took 98.1/60.6/71.3 s after the optimization (96 threads pay for
cache-line ownership transfers on the 38,880 shared nodes).

## 6. Reproduce

Build the library tests normally:

```
make -j8 tmpm01 tmpm02 tmpm03 tmpm04 tmpm05
```

The reduced drivers used for the tables are `tmpm03.cpp` with
`particles_per_diameter = 20` (W1) or `50` (W2), `Tf = 0.045` (W1) or
`0.00406` (W2), and

```
dom.Solve(Tf,dt,Tf,NULL,NULL,NULL,false,Nproc);
```

The full cases are the checked-in `tmpm/tmpm01` and `tmpm/tmpm04`
executables.  Keep `tf=200` for `tmpm04`: reducing it changes the physical
window (viscous diffusion takes ~25 s), while reducing
`particles_per_diameter` changes the spatial resolution.
