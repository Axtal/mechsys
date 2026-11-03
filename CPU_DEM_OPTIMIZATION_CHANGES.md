# CPU DEM Optimization for Voronoi Polyhedra

This document describes the CPU-side DEM optimizations made in
`lib/dem/`, using `tdem/test_periodic.cpp` as the reference. The work targets
the Voronoi-polyhedron path, but the hot loops are shared by the other CPU DEM
particle types.

The changes preserve the baseline CPU result exactly for the verification run.
The CUDA-side work is described separately in
[`CUDA_OPTIMIZATION_CHANGES.md`](CUDA_OPTIMIZATION_CHANGES.md).


## What is in the current change set

| Area | Change |
|---|---|
| Rotation state | `Particle` now owns `wrot = Rotation(w,Q)` and `Qc = Conjugate(Q)`; refreshed in `CalcProps`, `InitializeVelocity`, `Rotate`, `UpdateContacts`, and after a `Setup` callback. CUDA upload/download copies them |
| Interaction state | Removed the per-interaction `Wrot1/Wrot2/Qc1/Qc2` fields; `CalcForce` reads the particle fields |
| Rotation geometry | One fused vertex pass instead of the old-orientation/new-orientation double pass |
| Feature broadphase | Safe AABB/plane rejection in `Overlap`; `Distance2` early rejection in the force loop |
| Verlet pruning | Build-time feature distances plus a conservative `4.0*ContactDispGlobal()` cutoff |
| Threading | Per-thread force/torque accumulation and an OpenMP reduction instead of particle locks |
| Friction | One `find`/`emplace` lookup and a reference to the stored tangential displacement |

The exactness run described in Section 3 was repeated after the particle-level
rotation-cache change; it remains dataset-bit-identical at `Nproc=1`.

## 1. Reference Environment and Workload

- CPU: the same machine used for the other benchmarks in this workspace
- Compiler: the repository's normal `g++ -O3 -fopenmp` build
- Reference program: `tdem/test_periodic`
- Reference input: a reduced Voronoi input derived from `tdem/tpercu.inp`

The reduced input uses:

| Parameter | Value |
|---|---:|
| Particle type | `voronoi` |
| Particles | 480 |
| `Lx`, `Ly`, `Lz` | 10 |
| `nx`, `ny`, `nz` | 10 |
| Fraction | 0.5 |
| `dt` | 5e-5 |
| `T0` | 0.01 |
| `Tf` | 0.01 |
| Verlet | 0.05 |

`Solve(T0)` and `Solve(T0+Tf)` each advance 200 steps, so the run has 400
integration steps. The final state and energy are written to the usual HDF5
and XMF output files. The reduced and medium inputs are stored under
`benchmarks/cpu_voronoi/`.

The important point is that the reduced test still exercises the real Voronoi
feature structure: edge-edge (`Lee`), vertex-face (`Lvf`), and face-vertex
(`Lfv`) candidate lists, contact-list rebuilding, friction history, and rigid
body rotation. On the full `tpercu.inp` the same loops dominate, but with
roughly 32k Voronoi cells rather than 480.

## 2. Main Optimizations

### 2.1 Particle-level rotation caches

`CInteracton` previously evaluated `Rotation(w,Q)` and `Conjugate(Q)` for both
particles repeatedly for every feature pair. These quantities depend only on
one particle, not on the interaction.

`Particle` now stores:

```cpp
Vec3_t       wrot;   // Rotation(w,Q)
Quaternion_t Qc;     // Conjugate(Q)
```

and refreshes them in `Particle::UpdateRotCache()`:

```cpp
Rotation  (w,Q,wrot);
Conjugate (Q,Qc);
```

The cache is maintained at these points:

- `Particle::CalcProps()` after the orientation is calculated;
- `Particle::InitializeVelocity()`, the normal start of a Verlet run;
- `Particle::Rotate(double dt)` after `Q` is updated, exactly mirroring the
  GPU `Rotate`/`MoveParticles` kernels;
- `Domain::UpdateContacts()`, so the first force evaluation and direct
  `UpdateContacts`/`CalcForce` use current values;
- after a user `Setup` callback, which is allowed to change `Q` or `w`.

The domain-side refreshes use `UpdateRotCacheIfChanged()`, which compares the
cached source `w`/`Q` and skips the quaternion work when nothing changed.
`Rotate`, `CalcProps`, and `InitializeVelocity` refresh unconditionally.

The CUDA transfer path now treats the caches as part of the particle state:
`UploadParticle` copies `Par.wrot`/`Par.Qc` to the device instead of calling
`Rotation`/`Conjugate` again, and `DnloadParticle` copies the device values
back and updates the source snapshots. `UpLoadDevice` still calls
`UpdateRotCacheIfChanged()` before packing, so a direct upload after a host
edit cannot send stale caches. The device remains responsible for refreshing
the caches after `Q`/`w` change during movement.

`CInteracton::CalcForce` now simply reads `P1->wrot`, `P2->wrot`, `P1->Qc`, and
`P2->Qc`. The per-interaction `Wrot1/Wrot2/Qc1/Qc2` fields are gone. The work
is now done once per particle per step (and at initialization/contact
rebuild), instead of at the start of every interaction.

### 2.2 Fused rigid-body rotation

The original `Particle::Rotate(double dt)` performed two complete passes over
the vertices, edges, faces, and cylinders:

1. remove the old orientation,
2. install the new orientation.

The intermediate face normals and edge data are not used. The two vertex
passes are now fused, and the edge/face/cylinder update happens once:

```cpp
Vec3_t xt = *Verts[i]-x;
Vec3_t y; Rotation(xt,temp,y); y += x;
Vec3_t xt2 = y - x;
Vec3_t z; Rotation(xt2,Qnew,z);
*Verts[i] = z + x;
```

This halves the vertex and topology work in the rotation step while producing
the same final geometry, quaternion, and rotational energy.

### 2.3 Safe feature broadphase

`lib/dem/distance.h` now rejects impossible feature pairs before evaluating
the full distance:

- vertex-edge: axis-aligned bounding-box test on the segment, then the
  existing enclosing-sphere test;
- edge-edge: axis-aligned separating-axis test on the two segments, then the
  enclosing-sphere test;
- vertex-face: signed plane-distance test, using the maintained `Face::Nor`,
  then the existing face-centroid sphere test.

The broadphase tests include a small `1e-12` safety margin. They are only
allowed to reject pairs that are provably outside the contact range; the exact
`Distance` function still determines the closest points and the force.

`Distance2` helpers were added for the hot force path. They compute the same
squared feature distance without constructing `Xi`/`Xf`/`S`, and are used only
as a cheap rejection:

```cpp
if (Distance2((*A[i]),(*B[j]),Per) >= Rsum*Rsum) continue;
```

If the pair can be a contact, the original `Distance` is still called to form
the contact geometry and to keep the baseline arithmetic.

### 2.4 Verlet lists with build-time distances

Each `CInteracton` stores, alongside every Verlet candidate list, the feature
distance measured when the list was built:

```cpp
Array<double> Leed, Lvfd, Lfvd, Lvtd, Ltvd, Lvcd, Lcvd;
```

The lists keep their original `(i,j)` insertion order. That is essential:
changing the order of candidate features changes the floating-point
accumulation order of the contact force, even when the contact set is
identical.

During the force loop, a candidate whose build-time distance is already
outside the current conservative range is skipped:

```cpp
double cutoff = Rsum + 4.0*ContactDispGlobal();
...
if (Ld.Size()==L.Size() && Ld[k] > cutoff) continue;
```

`ContactDispGlobal()` adds `maxdis` every step, where `maxdis` is the
maximum vertex displacement measured from the last `ResetDisplacements()`;
it is reset by `Domain::UpdateContacts()`. The factor of
four is deliberately conservative: it leaves margin for relative motion of
the two features and for the fact that `MaxDisplacement` is measured from the
last contact rebuild rather than independently for every feature. A factor of
two was found insufficient on the reduced Voronoi test; with four the
candidate reduction is still large and the baseline state is reproduced
exactly.

When the list is rebuilt, `_update_contacts` uses the exact `Distance`
function and preserves the original nested-loop order. It does **not** sort the
contact list, because sorting would change the force summation order.

### 2.5 Lock-free force accumulation

The CPU force loop previously used one `omp_lock_t` per particle:

```cpp
omp_set_lock  (&Interactons[i]->P1->lck);
Interactons[i]->P1->F += Interactons[i]->F1;
...
```

Each OpenMP worker now accumulates into a private thread-local force/torque
array in `MtData`:

```cpp
size_t tid = omp_get_thread_num();
MTD[tid].F[Interactons[i]->I1] += Interactons[i]->F1;
MTD[tid].T[Interactons[i]->I1] += Interactons[i]->T1;
MTD[tid].F[Interactons[i]->I2] += Interactons[i]->F2;
MTD[tid].T[Interactons[i]->I2] += Interactons[i]->T2;
```

After the parallel loop the arrays are reduced into `Particles[ip]->F/T`.
For `Nproc=1` this is bit-identical to the original sequential accumulation;
the lock-free path removes the per-interaction particle locks for
multi-threaded runs.

The interaction loop uses `schedule(guided)` because the cost of an
interaction is now data-dependent: a pair may have many feature candidates,
but the distance cutoff skips most of them.

### 2.6 Single-lookup friction history

The friction history is an `unordered_map` keyed by the feature pair. The
contact loop previously performed `count()` followed by `operator[]` for every
contact. It now performs one lookup and keeps a reference to the stored
tangential displacement:

```cpp
auto fit = FMap.find(p);
if (fit==FMap.end()) fit = FMap.emplace(p,OrthoSys::O).first;
Vec3_t & fd = fit->second;
...
Ftnet += Kt*fd;
```

The arithmetic and the stored friction state are unchanged.

## 3. Correctness Verification

The reference baseline was built from the repository `HEAD` in a separate
worktree, using the same headers, flags, libraries, and input. Both programs
were run with `Nproc=1`, and all HDF5 datasets in the final state were compared
numerically.

Result: **all final-state datasets are bit-identical** to the baseline for the
reduced Voronoi run, including:

- `Position`, `Velocity`, `Quaternion`, `Verts`
- `PVelocity`, `PAngVelocity`, `PAngacceleration`, `PUForce`, `PEkin`
- `AngVelocity`, `Velocity`, `FaceCon`
- energy output (`0.0397229` and `5.98108`)

The latest particle-level cache implementation was verified separately after
the change: the complete final-state dataset comparison again returned no
differences (`neq []`). The test programs rebuilt after the change were
`test_periodic`, `test_cu_01`, `test_dynamics`, `tlbmdem01`, and `tlbmdem02`.

The HDF5 container MD5 can differ between two runs because of file metadata,
so the comparison is done dataset by dataset, not by whole-file hash.

The exact-order contact list, the original `Distance` arithmetic, and the
conservative displacement cutoff are the three details needed for this
bit-level agreement.

## 4. Benchmarks

Each time includes program startup, particle construction, two solves, and the
normal output writing. Timings are grouped by input and by the run that
produced them; the medium case was repeated to keep the baseline and current
binary under the same machine conditions.

### 4.1 Baseline versus optimized

| Threads | Baseline wall | Optimized wall | Wall speedup |
|---:|---:|---:|---:|
| 1 | 13.94 s | 9.27 s | 1.50x (33.5% faster) |
| 8 | 3.20 s | 2.55 s | 1.25x (20.3% faster) |

The simulation phase itself improves more than the total wall time because
the reduced run has a fixed construction/output overhead:

| Threads | Baseline solve time | Optimized solve time | Solve speedup |
|---:|---:|---:|---:|
| 1 | 12.73 s | 8.06 s | 1.58x (36.7% faster) |
| 8 | 1.99 s | 1.36 s | 1.47x (32.0% faster) |

### 4.2 Particle-level rotation-cache A/B

The move from the earlier per-interaction caches to the current particle-level
caches was measured against the binary built immediately before that edit.
Both runs used the same input, thread count, and build configuration.

| Case | Before | After | Change |
|---|---:|---:|---:|
| Reduced, `Nproc=1`, two-solve total | 8.03 s | 8.02 s | noise |
| Medium, `Nproc=8`, first solve | 5.51 s | 5.54 s | noise |
| Medium, `Nproc=8`, second solve | 6.62 s | 6.48 s | about 2% faster |
| Medium, `Nproc=8`, two-solve total | 12.13 s | 12.02 s | about 0.9% faster |
| Medium, `Nproc=8`, total wall | 21.83 s | 21.66 s | about 0.8% faster |

The reduced 480-particle case is noise-limited. The main benefit of this
follow-up is structural: the interaction path no longer carries or refreshes
four derived values per pair, and the cache lifecycle now matches the GPU
`DynParticleCU::wrot`/`Qc` design.

### 4.3 Medium Voronoi check

To make sure the reduced result was not an artifact of the 480-particle test,
a second input was made from `tpercu.inp` with `Lx=Ly=Lz=20`, `nx=ny=nz=20`,
`fraction=0.5` (about 4,000 Voronoi cells), `T0=Tf=0.01`, and
`dt=5e-5`. The same two-solve sequence was run with `Nproc=8`:

| Quantity | Baseline | Optimized |
|---|---:|---:|
| Construction/initial output | 3.02 s | 3.01 s |
| First 200-step solve | 8.22 s | 5.52 s |
| Second 200-step solve | 8.18 s | 6.47 s |
| Total wall | 26.06 s | 21.64 s |

That is about 17% less total wall time and about 27% less solve time on the
larger Voronoi case in this run. An earlier run on the same input was slower
for both binaries (32.55 s → 25.71 s / about 21% wall and 30% solve); the
machine load varies, so the relative A/B comparison in Section 4.2 is the
more controlled measurement of the cache follow-up. The `Nproc=8` energies
differ in the last printed digits because the original CPU force accumulation
is already order-dependent across threads; the exact `Nproc=1` comparison
above is the correctness check.

### 4.4 Optimized thread scaling

The optimized reduced run gives the following solve times:

| Threads | First solve | Second solve | Sum |
|---:|---:|---:|---:|
| 1 | 3.57 s | 4.49 s | 8.06 s |
| 2 | 1.98 s | 2.46 s | 4.44 s |
| 4 | 1.18 s | 1.42 s | 2.60 s |
| 8 | 0.618 s | 0.737 s | 1.36 s |

The reduced test is dominated by the unavoidable particle/element
construction and output work at high thread counts. The force loop is where
the optimization is concentrated; on a larger Voronoi run the force-loop
share is larger, so the benefit is expected to be larger than the 8-thread
total-wall result above.

## 5. Relevant Files

- `lib/dem/interacton.h`: particle-cache reads, friction lookup,
  Verlet-distance arrays, cutoff, contact rebuild
- `lib/dem/particle.h`: `wrot`/`Qc` fields, `UpdateRotCache()`,
  `UpdateRotCacheIfChanged()`, fused rotation, CUDA upload/download copies
- `lib/dem/domain.h`: cache refresh on contact rebuild/setup, per-thread
  force/torque accumulation, guided schedule, accumulated contact displacement
- `lib/dem/distance.h`: safe broadphase, squared-distance helpers
- `tdem/test_periodic.cpp`: reference test (not modified)

## 6. Build and Integration Checks

- `make test_periodic -j8` in `tdem/` succeeds without new warnings from the
  modified DEM code.
- `test_cu_01` and `test_dynamics` also rebuild after the shared
  `Particle`/`Interacton` changes.
- The reduced `test_periodic` run completes normally and prints the expected
  energies, with `neq []` against the `HEAD` baseline at `Nproc=1`.
- A host-side `UploadParticle`/`DnloadParticle` check confirmed that the copied
  `wrot`/`Qc` fields are bit-identical to calling `Rotation`/`Conjugate`
  directly, including a download-then-upload round trip.
- A reduced 1,600-particle `test_cu_03` case completed 10 GPU steps with
  `Nthread=32` and finite energy after the copy change. GPU output is not a
  bitwise baseline comparison because the force kernels use `atomicAdd`.
- The CPU LBM-DEM targets `tlbmdem01` and `tlbmdem02` and the CUDA LBM-DEM
  target `tlbmdem_cu_01` still build with the modified headers.
- The full `tlbmdem01` run (896,761 LBM cells) enters its solve phase without
  a crash but is too large to use as a quick smoke test on this machine; it
  was stopped by a 300 s timeout. A reduced copy of the same test
  (`41 x 21 x 21` LBM cells, `Tf=100`) completed **6250 iterations** with the
  modified DEM/LBM-DEM headers, confirming the CPU integration path still
  runs end to end. The reduced copy is a temporary smoke test, not part of the
  repository.

## 7. Limitations

- The optimization is benchmarked on the reduced 480-particle Voronoi input.
  The original `tpercu.inp` has roughly 32k cells and was not run to completion
  here; the hot loops are the same, but the absolute speedup will depend on
  the particle count and on how often the contact list is rebuilt.
- The displacement cutoff is a conservative heuristic. It is verified on the
  reduced run; if a simulation uses unusually large per-step rotations or a
  different contact margin, the factor can be increased or the cutoff can be
  disabled without changing the rest of the implementation.
- `Particle::wrot`/`Qc` are derived state. Normal `Solve` paths refresh them,
  but code that assigns `Q`/`w` directly and then calls a force routine without
  `UpdateContacts()` or `Rotate()` must call `Particle::UpdateRotCache()` first.
- A custom CUDA movement/rotation kernel that changes device `Q`/`w` must
  refresh `DPar[].wrot`/`DPar[].Qc` itself; the host copy is only an initial
  value and the built-in movement kernels maintain the device caches.
- The particle-level rotation-cache follow-up is a small timing change on these
  inputs. Its value is the eliminated per-interaction work and the closer
  correspondence with the GPU implementation, not a large wall-time win.
