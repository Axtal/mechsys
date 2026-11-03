# Performance Study: Reduced `test_cu_03`

## 1. What was measured

The user-reduced `tdem/test_cu_03.cu` was benchmarked on the available GPU.
The reduced test generates:

| Quantity | Value |
|---|---:|
| Particles | 1,601 |
| Vertices | 37,726 |
| Edges | 56,587 |
| Faces | 22,062 |
| Vertex-vertex interactions | 0 |
| Edge-edge interactions | 173k–1.39M, depending on Verlet distance |
| Vertex-face interactions | 66k–440k, depending on Verlet distance |
| Face-vertex interactions | 66k–439k, depending on Verlet distance |

Hardware and build:

- NVIDIA RTX A6000
- CUDA 12.4, `sm_86`
- `real = double` (`-DUSE_GPU_DOUBLE`)
- DEM CUDA path with the cached `wrot/Qc` and `*Pre` force kernels

The full test uses `Tf = 10` and `dt = 1e-5`, i.e. 1,000,000 steps. That is
still a long run for interactive tuning, so the performance study used
shortened, otherwise identical runs:

- `Tf = 0.1`: 10,000 steps
- `Tf = 0.2`: 20,000 steps
- `dt = 1e-5`
- A temporary benchmark copy added `Tf` as a command-line argument and
  measured the wall time around `dom.Solve`.
- The setup/initialization cost was measured separately with `Tf = 0` under
  the same parameters and subtracted. "ms/step" means:

```text
(total Solve time - Tf=0 setup time) / number of steps
```

This removes particle generation, initial contact construction, and device
upload from the per-step number. It still includes the final output performed
by `Solve`; that overhead is small relative to 10k–20k steps.

One anomaly is worth recording: `Nthread = 512` is not usable in this build.
The run failed at launch with `cudaErrorInvalidDevice: invalid device ordinal`
after the kernel configuration was accepted; 512-thread launches had also
previously exceeded the kernel's launch resources. The sweep therefore stops
at 384 threads.

## 2. Thread-count sweep

Fixed parameters:

```text
Alpha (Verlet distance) = 0.01
Tf                       = 0.1 (10,000 steps)
Nthread                  = 32 ... 384
```

| Nthread | Total (s) | Setup (s) | ms/step |
|---:|---:|---:|---:|
| **32** | **14.197** | 0.776 | **1.342** |
| 64 | 14.261 | 0.768 | 1.349 |
| 96 | 14.574 | 0.767 | 1.381 |
| 128 | 14.955 | 0.674 | 1.428 |
| 160 | 15.441 | 0.662 | 1.478 |
| 192 | 15.948 | 0.664 | 1.528 |
| 224 | 15.869 | 0.682 | 1.519 |
| 256 | 16.532 | 0.705 | 1.583 |
| 320 | 17.363 | 0.710 | 1.665 |
| 384 | 17.919 | 0.699 | 1.722 |

Repeated runs of the three smallest configurations at the same settings:

| Nthread | Run 1 (s) | Run 2 (s) | Run 3 (s) | Median (s) |
|---:|---:|---:|---:|---:|
| 32 | 14.170 | 14.106 | 14.157 | 14.157 |
| 64 | 14.236 | 14.318 | 14.288 | 14.288 |
| 128 | 14.834 | 14.870 | 14.902 | 14.870 |

![Thread-count sweep](benchmarks/test_cu_03_reduced/nthread_sweep.png)

### Interpretation

- **32 threads is the best block size for this reduced problem.** It is about
  0.5–1% faster than 64 and about 5–6% faster than 128.
- Compared with the previous default of 256, 32 threads is about **15% faster**
  per step.
- The curve is essentially monotonic: larger blocks are progressively worse.
  The small 192/224 inversion is within run-to-run noise. There is no useful
  "sweet spot" at 128 or 256 for this particle count.
- The absolute differences are real but modest around 32/64; 32 and 64 are
  effectively tied within run-to-run noise.
- This is a **size-dependent** result. The larger 77k-particle version had its
  optimum around 128 threads. After reducing the particle count, the kernel
  grid became much smaller, and 32-thread blocks provide more blocks and
  better SM occupancy for these scattered interaction kernels.
- Larger blocks also pay more tail/quantization cost when the number of
  feature interactions is only in the hundreds of thousands.

## 3. Verlet-distance sweep

Fixed parameters:

```text
Nthread = 32
Tf      = 0.1 (10,000 steps) and 0.2 (20,000 steps)
Alpha   = 0.0025 ... 0.2
```

"Updates" is the approximate number of contact-list rebuilds in the run,
derived from the average number of iterations between updates. A missing
average means no rebuild occurred during that interval.

### 10,000-step runs (Tf = 0.1)

| Alpha | Total (s) | Setup (s) | ms/step | Avg. steps/update | Approx. updates | EE interactions | VF interactions |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.0025 | 16.876 | 0.686 | 1.619 | 453 | ~22 | 172,834 | 66,242 |
| 0.005 | 14.597 | 0.659 | 1.394 | 907 | ~11 | 233,729 | 88,548 |
| **0.01** | **14.084** | 0.703 | **1.338** | 1,902 | ~5 | 300,686 | 111,758 |
| 0.02 | 14.523 | 0.767 | 1.376 | 4,254 | ~2 | 371,472 | 134,752 |
| 0.05 | 17.282 | 0.917 | 1.636 | 9,513 | ~1 | 510,707 | 177,229 |
| 0.1 | 22.054 | 1.078 | 2.098 | 0 | 0 | 745,743 | 248,495 |
| 0.2 | 36.017 | 1.651 | 3.437 | 0 | 0 | 1,385,622 | 440,196 |

### 20,000-step runs (Tf = 0.2)

| Alpha | Total (s) | Setup (s) | ms/step | Avg. steps/update | Approx. updates | EE interactions | VF interactions |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.005 | 43.461 | 0.659 | 2.140 | 316 | ~63 | 214,031 | 81,291 |
| 0.01 | 36.084 | 0.703 | 1.769 | 637 | ~31 | 289,363 | 107,741 |
| **0.02** | **32.918** | 0.767 | **1.608** | 1,301 | ~15 | 365,220 | 132,675 |
| 0.05 | 35.417 | 0.917 | 1.725 | 3,294 | ~6 | 506,649 | 175,969 |

![Verlet-distance sweep](benchmarks/test_cu_03_reduced/alpha_sweep.png)

### Interpretation

A Verlet distance controls two competing costs:

1. **Too small:** the contact list is rebuilt very often. At `Alpha = 0.0025`
   there is a rebuild about every 453 steps, and the run is slower than
   `Alpha = 0.01` even though its contact list is smaller.
2. **Too large:** the contact list contains many more candidate pairs. At
   `Alpha = 0.2` the edge-edge list is about 4.6x larger than at `Alpha = 0.01`
   (1.386M vs 301k), and the step is about 2.6x slower.

For the shorter window, the optimum is around **0.01**. The broad optimum is
approximately **0.005–0.02**; all three are within a few percent of each other.

For the longer window, the optimum shifts to **0.02**:

- 0.02: 1.608 ms/step
- 0.01: 1.769 ms/step (+10%)
- 0.05: 1.725 ms/step (+7%)
- 0.005: 2.140 ms/step (+33%)

The shift happens because the particles move and settle differently as the
simulation progresses. The average rebuild interval at `Alpha = 0.01` drops
from 1,902 steps in the first 0.1 s to 637 steps over 0.2 s. In the longer run,
`Alpha = 0.005` is forced into about 63 rebuilds, which dominates its smaller
contact lists.

The practical conclusion is therefore **not** that 0.01 is a universal
constant. It is that:

- the default-like `Alpha = 0.05` is too large for this reduced system;
- `Alpha = 0.01–0.02` is the useful tuning range;
- 0.02 is the safer choice for longer runs, while 0.01 is best for the first
  part of the collapse;
- the optimum should be rechecked for the actual simulated duration.

## 4. Interaction between Nthread and Alpha

The thread-count preference does not reverse when the contact list becomes
larger. With `Tf = 0.1`:

| Alpha | 32 threads (s) | 64 threads (s) | 128 threads (s) | 256 threads (s) |
|---:|---:|---:|---:|---:|
| 0.05 | 17.430 | 17.489 | 18.127 | 20.011 |
| 0.10 | 22.192 | 22.345 | 22.911 | 24.893 |

With the longer `Tf = 0.2` run at `Alpha = 0.02`:

| Nthread | Total (s) | ms/step |
|---:|---:|---:|
| 32 | 32.944 | 1.609 |
| 64 | 33.180 | 1.621 |
| 128 | 34.357 | 1.680 |
| 256 | 37.757 | 1.850 |

So **32 remains best**, 64 is within about 1%, and 128/256 lose progressively
more as the work per kernel increases.

## 5. Practical conclusions

For the reduced `test_cu_03` configuration:

1. **Use 32 threads per block** (or 64 if a slightly larger block is preferred
   for consistency with other tests). Do not keep the old 256-thread default
   for this particle count: it costs about 15% per step in the Alpha = 0.01
   case and about 13% in the longer Alpha = 0.02 case.
2. **Use a Verlet distance around 0.01–0.02**, not 0.05. For the short window
   0.01 is best; for the longer window 0.02 is best. `0.0025` is too small
   because contact rebuilding dominates, and `0.1–0.2` is far too large
   because the candidate interaction count explodes.
3. The best thread count is **problem-size dependent**. The larger 77k test
   preferred 128 threads; the reduced 1.6k test prefers 32. Re-run a short
   sweep whenever the particle count or contact density changes materially.
4. The full `Tf = 10` run is still long. At roughly 1.3–1.6 ms/step it is on
   the order of 20–27 minutes for 1e6 steps, depending on the chosen Alpha and
   the evolving contact state. The short benchmarks are suitable for tuning,
   but the final choice should be checked over the real simulation window.
5. **Do not interpret the reduced-case optimum as a general DEM optimum.**
   The old large-case conclusions (128 threads best, about 14.7% overall gain
   from the force-kernel optimizations) and these reduced-case conclusions
   answer different workloads.

## 6. Suggested benchmark settings

For quick tuning of the reduced test, make a temporary copy with `Tf=0.1`
(the checked-in source still runs `Tf=10`). The copy accepts the same
`Nproc Alpha Nthread` arguments plus optional adaptive controls
`flag window passes warm-up`:

```bash
# fixed 32 threads
./test_cu_03_short 1 0.01 32

# adaptive, starting from 256; sweep 30-step windows, 2 passes, 100 warm-up
./test_cu_03_short 1 0.01 256 1 30 2 100
```

The reduced copy differed from the checked-in `tdem/test_cu_03.cu` only in
`Tf` and the output file key.

## 7. Raw data

The CSV files next to this report contain the plotted data:

- `benchmarks/test_cu_03_reduced/nthread_sweep.csv`
- `benchmarks/test_cu_03_reduced/alpha_sweep.csv`

All runs are single measurements except the 32/64/128 repeatability set shown
in Section 2. The main sources of uncertainty are GPU clock variation and the
fact that the evolving particle configuration changes the contact list even
when the numerical parameters are fixed. Differences below about 2% should not
be treated as significant without more repetitions.


## 8. Adaptive controllers (rolled back)

`AdaptiveNthread` and `AdaptiveAlpha` were both implemented experimentally and
then rolled back. No adaptive field or controller code remains in DEM, FLBM, or
LBMDEM. The fixed-`Nthread` measurements in this report are still valid; the
adaptive tables from the earlier draft are not part of the current code.

For any future Verlet-distance controller: increasing `Alpha` requires an
immediate contact rebuild, decreasing it is safe only while the accumulated
displacement is below the new value, and the linked-cell `Beta` cap must be
respected.
