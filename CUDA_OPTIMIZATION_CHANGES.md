# CUDA DEM / LBM-DEM Optimization and Cleanup

This document describes the CUDA changes made to the DEM and LBM-DEM
implementations in this working tree, including the optimizations, the
resulting API changes, the clock-consistency refactor, and the associated
build and accuracy checks.

The changes are relative to the original CUDA implementation in the
repository. The main modified files are:

- `lib/dem/dem.cuh`
- `lib/dem/domain.h`
- `lib/dem/particle.h`
- `lib/dem/distance.h`
- `lib/lbmdem/Domain.h`
- `lib/flbm/Domain.h`
- `lib/flbm/lbm.cuh`
- `lib/mesh/mesh.h`
- `lib/util/array.h`
- `tdem/test_cu_03.cu`, `tdem/test_cu_04.cu`
- `tlbmdem/tlbmdem_cu_02.cu`, `tlbmdem/tlbmdem_cu_04.cu`

## 1. Reference Environment and Workload

The performance measurements in this document were made on:

- GPU: NVIDIA RTX A6000, 47.69 GB
- CUDA: 12.4
- Architecture: `sm_86`
- Registers: 65,536 per SM
- Data type: `real = double`, compiled with `-DUSE_GPU_DOUBLE`

The main reference test is `tdem/test_cu_03`, which uses 77,003 particles and:

| Quantity | Value |
|---|---:|
| Vertices | 1,815,920 |
| Edges | 2,723,878 |
| Faces | 1,061,963 |
| Vertex-vertex contacts | 448,164 |
| Edge-edge feature interactions | 26,872,324 |
| Vertex-face feature interactions | 9,305,675 |
| Face-vertex feature interactions | 9,303,317 |

The test is dominated by the feature interaction kernels, not by particle
array access.

## 2. Cached Rotation Quantities

### 2.1 New fields in `DynParticleCU`

`lib/dem/particle.h` now stores two additional per-particle quantities:

```cpp
struct DynParticleCU
{
    real3 x, xb, v, w;
    real3 wrot;   // Rotation(w,Q)
    real3 wa, wb;
    real3 F, Flbm;
    real4 Q;
    real4 Qc;     // Conjugate(Q)
};
```

The fields are prepared on the host by `Particle::CalcProps()`,
`Particle::InitializeVelocity()`, and `Particle::UpdateRotCache()`.
`UploadParticle` now copies them with the rest of the particle state instead
of evaluating the quaternion functions again:

```cpp
DPc.wrot.x = Par.wrot(0);
DPc.wrot.y = Par.wrot(1);
DPc.wrot.z = Par.wrot(2);
DPc.Qc.w   = Par.Qc(0);
DPc.Qc.x   = Par.Qc(1);
DPc.Qc.y   = Par.Qc(2);
DPc.Qc.z   = Par.Qc(3);
```

`DnloadParticle` performs the inverse copy into `Par.wrot`/`Par.Qc` and
updates the host source snapshots. `UpLoadDevice` also calls
`UpdateRotCacheIfChanged()` before packing, so a direct upload after a host
edit cannot send stale caches. This removes the need for separate device
arrays for the rotated angular velocity and the conjugate quaternion; the
first force evaluation has valid cache values.

The low-level `UploadParticle` function itself is a pure copy. Callers that
use it directly, outside `Domain::UpLoadDevice`, must prepare the host cache
first.

### 2.2 Cache refresh during rotation

The cache is refreshed after the new quaternion is computed, in both default
rotation paths:

- `Rotate`
- `MoveParticles`

immediately after:

```cpp
DPar[ic].Q = Qnew;
```

with:

```cpp
Rotation (DPar[ic].w,DPar[ic].Q,DPar[ic].wrot);
Conjugate(DPar[ic].Q,DPar[ic].Qc);
```

This is the same arithmetic and the same inputs that were previously used by
`MaxD` after the movement kernel. Moving it into the rotation path makes the
ownership explicit: **`wrot` and `Qc` are refreshed exactly when `w` or `Q`
changes**.

The initial host-side values come from `InitializeVelocity`/`CalcProps` and are
copied by `UploadParticle`; there is no separate device-side
cache-initialization kernel. `DnloadParticle` copies the final device caches
back so a later upload remains consistent.

### 2.3 `MaxD` simplification

`MaxD` no longer has a `cacheOnly` argument and no longer updates `wrot/Qc` or
the clock. Its signature is now:

```cpp
__global__ void MaxD(real3 const * Verts, real3 const * Vertso,
                     real * maxd, dem_aux const * demaux);
```

It performs only the maximum-displacement reduction. The previous
`MaxD(..., true)` startup call was removed.

## 3. Specialized DEM Force Kernels

### 3.1 Motivation

The original edge-edge, vertex-face, and face-vertex force kernels recomputed
the particle rotation quantities inside the interaction loop:

```cpp
Rotation (DPar[i].w,DPar[i].Q,rotated_w);
Conjugate(DPar[i].Q,q_conjugate);
```

These quantities depend only on the particle, not on the interaction. Since a
particle participates in many feature interactions, this repeated the same
quaternion work many times per force evaluation.

### 3.2 New kernels

Three specialized kernels were added in `lib/dem/dem.cuh`:

- `CalcForceEEPre`
- `CalcForceVFPre`
- `CalcForceFVPre`

They read the cached values:

```cpp
real3 t1 = DPar[i1].wrot;
real3 t2 = DPar[i2].wrot;
...
Rotation(Tt,DPar[i1].Qc,T);
```

The arithmetic of the force law itself is unchanged. Only the source of the
rotation quantities changed: from per-interaction recomputation to the cached
per-particle values.

### 3.3 Default dispatch

`DEM::Domain` still exposes the force function pointers:

```cpp
ForceVV_ptr_t pForceVV;
ForceEE_ptr_t pForceEE;
ForceVF_ptr_t pForceVF;
ForceFV_ptr_t pForceFV;
```

The default assignments in the constructor are now:

```cpp
pForceEE = CalcForceEEPre;
pForceVF = CalcForceVFPre;
pForceFV = CalcForceFVPre;
```

Users can still assign the original kernels or their own kernels after
construction. The new kernels are the default because they are the optimized
path for the standard force laws.

`CalcForceVV` remains on its existing path; the vertex-vertex kernel was not
changed to use the cache in this work.

## 4. Combined Movement Kernel

### 4.1 `MoveParticles`

A new combined movement kernel was added:

```cpp
__global__ void MoveParticles(real3 * Verts, ParticleCU const * Par,
                              DynParticleCU * DPar, dem_aux const * demaux,
                              bool doRotate, void * extraparams);
```

It contains the translation arithmetic of `Translate` followed by the
rotation arithmetic of `Rotate`, selected by `doRotate`. The arithmetic is
line-for-line equivalent to the original two kernels.

The corresponding function pointer is:

```cpp
typedef void (*Move_ptr_t)(real3 *, ParticleCU const *, DynParticleCU *,
                           dem_aux const *, bool, void *);
```

and `DEM::Domain` has:

```cpp
Move_ptr_t pMove;   // NULL means use pTranslate and pRotate separately
```

The constructor assigns:

```cpp
pMove = MoveParticles;
```

### 4.2 Fallback behavior

`Domain::Solve` uses:

```cpp
if (pMove)
    pMove<<<...>>>(..., RotPar, pExtraParams);
else
{
    pTranslate<<<...>>>(...);
    if (RotPar) pRotate<<<...>>>(...);
}
```

Setting `pMove = NULL` restores the old two-kernel path. This preserves the
ability to use a custom translation and rotation independently, or a custom
combined movement law.

### 4.3 Cache ownership for custom movement

Because `MaxD` no longer refreshes `wrot/Qc`, a custom `pMove`, `pRotate`, or
other custom kernel that changes `Q`/`w` is responsible for refreshing:

```cpp
Rotation (DPar[ic].w,DPar[ic].Q,DPar[ic].wrot);
Conjugate(DPar[ic].Q,DPar[ic].Qc);
```

The default `MoveParticles` and `Rotate` do this automatically. If rotation is
disabled (`doRotate == false`), `Q`/`w` do not change, so the existing cache
remains valid.

## 5. Maximum-Displacement Reduction

`MaxD` previously wrote one displacement value per vertex. The host then ran
`thrust::max_element` over the entire per-vertex array every step.

The current kernel reduces within each block and writes:

```cpp
maxd[blockIdx.x] = block_max;
```

The host scans only the actual number of launched blocks:

```cpp
size_t nmaxblocks = demaux.nverts/Nthread+1;
thrust::max_element(bMaxDCU.begin(), bMaxDCU.begin()+nmaxblocks);
```

The maximum of a set is exact and independent of reduction order, so this does
not change the Verlet/contact-update condition. It substantially reduces the
number of values written and read by the host each step.

## 6. Clock Consistency

### 6.1 Problem

`dem_aux` and `lbm_aux` both contain:

```cpp
size_t iter;
real   Time;
```

Before this refactor they were advanced in different places:

- `dem_aux` was incremented inside `DEM::MaxD`, and sometimes again on the
  host in `DEM::Domain::Solve` / `LBMDEM::Domain::Solve`.
- `lbm_aux` was incremented inside `cudaStream2`, `cudaStreamPF2`,
  `cudaStreamPFI2`, `cudaStreamSW2`, and `cudaTick`, depending on the solver.
- The host `Time` was incremented separately in the solver loop.
- In `LBMDEM`, `DEMDOM.Time` and `LBMDOM.Time` were not advanced at all.

This produced multiple clocks that could diverge, especially in the coupled
LBM-DEM solver.

### 6.2 Single authority

The host is now the single owner of the simulation clock:

- The host advances `Time` and `iter` once per outer step.
- The host writes the same values into the host copies of `dem_aux` and
  `lbm_aux`.
- A small kernel sets the device copies to those absolute values.

New device helpers:

```cpp
__global__ void DEM::SetClock(dem_aux * demaux, size_t iter, real Time);
__global__ void FLBM::cudaSetClock(lbm_aux * lbmaux, size_t iter, real Time);
```

They **set** the fields; they do not increment them.

All kernel-side clock increments were removed:

- `MaxD` no longer advances `dem_aux`.
- `cudaStream2`, `cudaStreamPF2`, `cudaStreamPFI2`, and `cudaStreamSW2` no
  longer advance `lbm_aux`.
- `cudaTick` and `LaunchTick` were removed.

### 6.3 Standalone DEM

At the end of each DEM step:

```cpp
Time += Dt;
iter++;
demaux.iter = iter;
demaux.Time = Time;
SetClock<<<1,1>>>(pdemaux,demaux.iter,demaux.Time);
```

`UpLoadDevice` also initializes:

```cpp
demaux.iter = iter;
demaux.Time = Time;
```

so a restarted domain starts from the correct host values.

### 6.4 Standalone LBM

`FLBM::Domain::Solve` keeps a local `cuda_iter` and, after each step:

```cpp
Time += dt;
cudaSetClock<<<1,1>>>(plbmaux,++cuda_iter,Time);
```

### 6.5 LBMDEM

At the end of each coupled outer step, all four views are written from the
same `(iter, Time)` pair:

```cpp
Time += dt;
DEMDOM.Time = Time;
LBMDOM.Time = Time;
DEMDOM.iter++;
iter++;
DEMDOM.demaux.iter = iter;
DEMDOM.demaux.Time = Time;
DEM::SetClock<<<1,1>>>(DEMDOM.pdemaux, iter, Time);
FLBM::cudaSetClock<<<1,1>>>(LBMDOM.plbmaux, iter, Time);
```

A temporary two-step diagnostic run confirmed:

```text
host=0      DEMhost=0      DEMdev=0      LBMhost=0      LBMdev=0
host=5.008  DEMhost=5.008  DEMdev=5.008  LBMhost=5.008  LBMdev=5.008  iter=313
host=10     DEMhost=10     DEMdev=10     LBMhost=10     LBMdev=10     iter=625
```

The diagnostic print was removed after verification.

## 7. Distance-Function Optimizations

Two changes were made in `lib/dem/distance.h`.

### 7.1 `DistanceEE`

The four candidate edge-edge distances are now compared using squared
lengths:

```cpp
real l1 = dotreal3(s1,s1);
real l2 = dotreal3(s2,s2);
real l3 = dotreal3(s3,s3);
real l4 = dotreal3(s4,s4);
```

Square root is monotonic, so the selected candidate is the same, but four
`sqrt` operations per feature are avoided.

### 7.2 `DistanceVF`

The vertex-face edge loop no longer computes the first edge twice or relies on
the first `norm()` result as the initial minimum. It initializes:

```cpp
double ld = 1.0e300;
Xf = make_real3(0.0,0.0,0.0);
S  = make_real3(0.0,0.0,0.0);
```

and compares squared edge distances:

```cpp
double lt = dotreal3(st,st);
if (lt < ld) { ... }
```

This removes the duplicated first-edge `DistanceVE` call and the extra square
roots while preserving the selected edge.

## 8. Compilation-Warning Cleanup

The following warnings were fixed:

- `lib/dem/domain.h`
  - HDF5 `hsize_t` was printed with `%zu`; now `%llu` with an explicit cast.
  - `CriticalDt()` could fall through without returning; now returns `0.0`.
- `lib/mesh/mesh.h`
  - `InsideCell(...)` could fall through without returning; now returns
    `false`.
  - Removed the unused `vtot`.
- `lib/util/array.h`
  - Removed an always-false `IdxToRemove<0` check for a `size_t`.
  - Fixed `DelVal`: `Find()` returns `long` with `-1` as the not-found
    sentinel, but the result was stored in `size_t`, turning `-1` into a huge
    index. It now remains `long` and is checked before conversion.
- `tdem/test_cu_03.cu`
  - Removed unused local `Nthread` and `g`.
- `tdem/test_cu_04.cu`
  - Removed unused `load`.
- `tlbmdem/tlbmdem_cu_02.cu`, `tlbmdem/tlbmdem_cu_04.cu`
  - Removed unused `dat` references.

The only remaining warning in a full `tdem` build is from the external Blitz++
header:

```text
/home/sergio/pkg/blitz-0.9/blitz/tinymat.h(60):
warning #191-D: type qualifier is meaningless on cast type
```

It is third-party code and was not modified.

## 9. Data Layout: AoS Retained

No structure-of-arrays rewrite was made.

The particle, vertex, and interaction arrays remain array-of-structures:

- `ParticleCU` and `DynParticleCU` are AoS.
- `Verts` is an array of `real3`, i.e. AoS with a 24-byte stride.
- `DynInteractonCU`, `ComInteractonCU`, and `InteractonCU` are AoS.
- The edge/face index arrays are flat scalar arrays, not SoA vector arrays.

For the force kernels, the relevant access pattern is:

- `DPar[i1]` / `DPar[i2]` use scattered particle indices.
- The kernels consume several fields of the same particle.
- `Verts[iv]` accesses use a compact 24-byte record and consume x/y/z
  together.

A full SoA conversion would increase load instruction count and lose useful
spatial locality in these paths. The adopted change is the smaller,
targeted one: cache the two derived rotation quantities in the existing AoS
particle record. A hot/cold split of `DynParticleCU` remains a possible future
experiment but is not part of this change.

## 10. Accuracy

### 10.1 Original accuracy observations

The following checks were made during the optimization work:

| Test | Result |
|---|---|
| `test_cu_01` | Bit-identical fields and vertices |
| `test_cu_02` | 1 ULP in one `w`, `wb`, and `Q`; positions, vertices, forces, `wa`, and energies identical |
| `test_cu_03` | `x`, `xb`, `v`, `Q`, and vertices identical; `w`, `wb`, `wa`, `F`, `Erot` last-bit differences; energy `0.155521` |
| `test_cu_04` | Energy `0.589492`; last-bit differences in `F`/`wa` and one vertex |

`test_cu_03` and `test_cu_04` are **not bitwise deterministic in the baseline**
because the force kernels use `atomicAdd`. Two baseline runs already differ in
the last bits. Exact cross-version equality is therefore not meaningful for
those tests; the relevant check is that the differences remain at the
last-bit level and do not change the physical result.

### 10.2 Checks after the final cache/clock changes

- `test_cu_01` completed with energy `0.0200495`.
- `test_cu_02` completed with energy `269130`.
- A one-step polyhedral `test_cu_03` run completed without NaN/overlap errors,
  with energy `9.25173e-05`.
- A host-side `UploadParticle`/`DnloadParticle` round-trip check confirmed that
  the copied `wrot`/`Qc` values are bit-identical to direct
  `Rotation`/`Conjugate` calls.
- A reduced 1,600-particle `test_cu_03` case ran 10 GPU steps with
  `Nthread=32` after the copy change and finished with finite energy
  `0.000691404` (not a bitwise baseline comparison; the GPU force path uses
  `atomicAdd`).
- Full `tlbmdem_cu_03` completed all 62,501 iterations normally in about
  3m40s.

The full `test_cu_03` (10 simulated time units, 1e6 steps) was not rerun after
the final cache move because of its runtime; the one-step polyhedral run was
used to exercise the cached `Pre` kernels. The final cache move preserves the
same arithmetic and inputs as the previous `MaxD` cache update.

## 11. Performance

Representative measurements during the optimization were:

| Version | Time per step |
|---|---:|
| Original CUDA | 99.06 ms |
| Optimized CUDA | 84.51 ms |

This is approximately a **14.7% speedup** for the `test_cu_03` workload.

A thread-count sweep gave:

| `Nthread` | Time per step |
|---:|---:|
| 64 | ~83.0–83.7 ms |
| 128 | **82.79 ms** |
| 192 | ~83.0–83.7 ms |
| 224–384 | ~85.1–86.0 ms |
| 512 | launch failure: too many resources requested |

The recommended default is therefore **128 threads** for this workload, with
the range 64–192 being close. `Nthread=256` remains a conservative default for
general use. The 512-thread configuration is not valid for the current
register usage.

These numbers were collected before the final struct-field and movement
refactor; the final code should be re-benchmarked with the same workload before
using the values as a final performance claim.

## 12. Build and Run

From the repository root:

```bash
cd tdem
make -j8 test_cu_01 test_cu_02 test_cu_03 test_cu_04

cd ../tlbmdem
make -j8 tlbmdem_cu_01 tlbmdem_cu_02 tlbmdem_cu_03 tlbmdem_cu_04
make -j8 tlbmdem01 tlbmdem02
```

The CUDA targets use the project CMake configuration, including the
`sm_86`/`compute_86` code generation and `-DUSE_GPU_DOUBLE`.

Typical runs:

```bash
cd tdem
./test_cu_01
./test_cu_02

# test_cu_03: Nproc, Alpha, Nthread
./test_cu_03 1 0.05 128

# test_cu_01..04 are long; use a shorter Tf while developing
```

For LBM-DEM:

```bash
cd tlbmdem
./tlbmdem_cu_03 1
```

`tlbmdem_cu_01` and `tlbmdem_cu_02` are substantially longer; reduced-`Tf`
runs were used for smoke testing during this work.

## 13. API and Customization Notes

The following function pointers remain available after construction:

```cpp
pForceVV, pForceEE, pForceVF, pForceFV
pTranslate, pRotate, pMove, pReset
```

- `pForceEE/VF/FV` default to the cached `*Pre` kernels.
- Assigning the original `CalcForceEE/VF/FV` restores the previous force
  path.
- `pMove = MoveParticles` is the default combined movement path.
- `pMove = NULL` uses `pTranslate` and `pRotate` separately.
- Custom movement/rotation kernels that change `Q`/`w` must refresh
  `DPar[].wrot` and `DPar[].Qc` themselves; otherwise the next force
  evaluation will use stale cached values.

The `dem_aux`/`lbm_aux` clocks are set by the host after each outer step. A
custom kernel reading `demaux->Time`/`iter` sees the current start-of-step
values; the clock is advanced once at the end of the step.

## 14. File-by-File Index

| File | Changes |
|---|---|
| `lib/dem/dem.cuh` | `Move_ptr_t`; `CalcForceEEPre/VFPre/FVPre`; `MoveParticles`; cache refresh in `Rotate`/`MoveParticles`; `SetClock`; simplified `MaxD` |
| `lib/dem/domain.h` | `pMove`; new default dispatch; combined movement with fallback; `MaxD` call/signature; host-owned clocks; HDF5 and `CriticalDt` warning fixes; opt-in adaptive-`Nthread` startup sweep |
| `lib/dem/particle.h` | `wrot`/`Qc` fields; initialization in `CalcProps`/`InitializeVelocity`; copies in `UploadParticle`/`DnloadParticle`; fused CPU rotation |
| `lib/dem/distance.h` | Squared-distance comparisons in `DistanceEE`; `DistanceVF` edge-loop cleanup |
| `lib/lbmdem/Domain.h` | `pMove` path; `MaxD` block scan; removal of duplicate clock updates; unified DEM/LBM clocks |
| `lib/flbm/Domain.h` | Host-driven `cuda_iter`; `cudaSetClock` after each step |
| `lib/flbm/lbm.cuh` | `cudaSetClock` replaces `cudaTick`; removed clock writes from stream kernels |
| `lib/mesh/mesh.h` | Missing `return false`; removed unused variable |
| `lib/util/array.h` | Unsigned comparison cleanup; `DelVal` sentinel fix |
| `tdem/test_cu_03.cu` | Optional adaptive-`Nthread` arguments; warning cleanup |
| `tdem/test_cu_04.cu` | Warning cleanup (CRLF preserved) |
| `tlbmdem/tlbmdem_cu_02.cu`, `tlbmdem/tlbmdem_cu_04.cu` | Unused-variable warning cleanup |

## 15. Known Limitations

1. Full `test_cu_03` and `test_cu_04` runs are long; not every full run was
   repeated after the final cache relocation. The one-step polyhedral check
   and the existing accuracy results should be repeated if a strict bitwise
   audit is required.
2. `test_cu_03`/`test_cu_04` cannot be made bitwise reproducible against the
   baseline because the baseline itself is nondeterministic due to atomics.
3. The thread-count values are workload-specific. Re-run the sweep for a
   different particle count or shape distribution.
4. The remaining Blitz++ warning is external and was left untouched.
5. Adaptive `Nthread` changes block size during the startup sweep, so it
   changes the atomic accumulation order and can diverge from a fixed-`Nthread`
   trajectory. It is a performance mode, not a reproducibility mode.

## 16. Adaptive controllers (rolled back)

### 16.1 Adaptive `Nthread`

`AdaptiveNthread` was implemented for `DEM::Domain`, `FLBM::Domain` and
`LBMDEM::Domain`, benchmarked, and then rolled back. No `AdaptiveNthread`
field, CLI option, sweep, or occupancy probe remains; `Nthread` is a fixed
launch parameter again.

The startup sweep was not worth the complexity or the extra steps. On the
reduced DEM test it selected 32, but the LBM cases were effectively flat:

- `tclbm01` (2400x300 D2Q9, 40,000 steps, no file output): `Nthread=32..256`
  all took 45.1-45.7 s solve time. Three runs each at 32 and 256 gave 45.11 vs
  45.72 s; Welch p=0.11 and a 95% CI of [-1.49,+0.26] s, so the apparent 1.3%
  difference is noise.
- Reduced `tlbmdem_cu_02` (200^3, `Tf=20`, 1,250 steps, no file output):
  fixed 32-256 all took 33.62-33.74 s wall; adaptive selected 96 at 33.86 s.

The candidate windows were also sequential and therefore exposed to ordinary
GPU load and clock drift. Use a fixed `Nthread` and benchmark the actual
geometry/workload if the default is suspect.

### 16.2 Verlet-distance adaptation: not implemented

An experimental `AdaptiveAlpha` controller was implemented and benchmarked,
then rolled back on request. No `AdaptiveAlpha` field or controller code
remains. The safety rules for any future attempt are:

- **Increasing `Alpha` requires an immediate contact rebuild.** Existing
  `Interactons`/linked-cell lists were built with the old, smaller margin. If
  `Alpha` is increased without rebuilding, the next rebuild is delayed and
  pairs between the old and new margins are never considered.
- **Decreasing `Alpha` is safe** as long as the new value is at least the
  current `maxdis` between rebuilds; changing it at a rebuild is simplest
  because the displacement counter resets.
- **Linked cells impose a cap.** If `Alpha > 2*Beta*MaxDmax`, update `Beta`
  and the cell dimensions before rebuilding, as at the start of `Solve`;
  otherwise the +/-1-cell neighborhood can miss pairs.
- Keep the existing trigger `maxdis > Alpha`. Adaptive tuning must never
  replace this safety check.

Adaptive `Alpha` changes which candidate contacts exist and therefore changes
the trajectory at the floating-point level. It cannot be claimed bitwise
reproducible. Validate energy and contact counts against a fixed-parameter
run, not with a bitwise GPU test.
