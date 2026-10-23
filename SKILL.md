---
name: mechsys
description: Use when a task mentions MechSys, tdem, tlbmdem, DEM/Particle/Voronoi optimization, CUDA DEM tests, LBM-DEM, or the MechSys simulation library in this workspace.
whenToUse: Any request to inspect, build, benchmark, debug, optimize, or document this MechSys checkout, especially its CPU DEM, CUDA DEM, LBM-DEM, or FLBM paths.
---

# MechSys

This is the visible, downloadable copy of the MechSys skill. DSH discovers
`.agents/skills/mechsys/SKILL.md`; keep the two files in sync.

This is the project-local skill for the MechSys checkout. It is discovered
automatically from `.agents/skills/mechsys`; do not install a separate plugin
just to provide this static guidance.

The authoritative context is:

1. `AGENTS.md` — automatic, concise agent rules.
2. `MECHSYS_HANDBOOK.md` — architecture, build/run commands, exactness
   procedure, benchmarks, and known traps.
3. `MEMORY.md` — distilled experience and a short verification checklist.
4. The optimization reports:
   - `CPU_DEM_OPTIMIZATION_CHANGES.md`
   - `CUDA_OPTIMIZATION_CHANGES.md`
   - `REDUCED_TEST_CU_03_PERFORMANCE.md`

Read the relevant documents before making substantive changes, and update them
when a durable lesson or corrected benchmark is learned.

## Non-negotiable rules

- Do not assume a given checkout or `HEAD` contains the optimized state.
  Check `git status`, the reports, and compare against a clean baseline
  created from the chosen commit.
- CPU DEM exactness is checked at `Nproc=1` by comparing HDF5 **datasets**, not
  whole-file MD5. Expect `neq []` and matching printed energies.
- Preserve candidate-list order and the original `Distance` arithmetic in the
  contact path. Do not sort the lists.
- The CPU Verlet cutoff is `Rsum + 4.0*ContactDispGlobal()`. Do not lower it
  without a new safe derivation and an exact comparison.
- Broadphase bounds must be conservative. If in doubt, remove a
  rejection rather than risk a missed contact.
- CUDA `test_cu_03/04` use `atomicAdd`; the baseline is non-deterministic. Do
  not claim bitwise equality for those GPU runs. `Nthread=512` fails for the
  current large test; tested range is up to 384.
- Keep the separate `pTranslate` + `pRotate` fallback. Making fused movement
  unconditional broke `test_cu_04`.
- The outer solve loop owns `Time`/`iter`; publish them to CUDA with `SetClock`.
  Do not add independent clock increments in DEM or LBM substeps.
- After changing `lib/dem/*`, build at least one CPU target, one CUDA target,
  and the relevant LBM-DEM target.

## Rotation caches

`Particle` owns `wrot = Rotation(w,Q)` and `Qc = Conjugate(Q)`. They are
refreshed by `CalcProps`, `InitializeVelocity`, `Rotate`, `UpdateContacts`, and
after a user `Setup` callback.

- CUDA `UploadParticle` copies the host caches; it does not recompute them.
- CUDA `DnloadParticle` copies the device caches back.
- The built-in device movement kernels refresh the device caches after
  changing `Q`/`w`.
- A custom GPU movement/rotation kernel that changes `Q`/`w` must refresh
  `DPar[].wrot`/`DPar[].Qc` itself.
- Host code that assigns `Q`/`w` directly and bypasses the normal paths must
  call `Particle::UpdateRotCache()` before using the caches.

## Build and run

From the repository root:

```bash
cd tdem
make -j8 test_periodic
make -j8 test_cu_01 test_cu_02 test_cu_03 test_cu_04
./test_periodic ../benchmarks/cpu_voronoi/tper_small 8

cd ../tlbmdem
make -j8 tlbmdem01 tlbmdem02
make -j8 tlbmdem_cu_01
```

`test_periodic` takes a FileKey prefix (the `.inp` suffix is appended).
`test_cu_03` takes `Nproc`, Verlet distance, and threads per block:
`./test_cu_03 1 0.05 128`. Full LBM-DEM and CUDA tests are long; use a
temporary reduced copy for smoke tests and do not modify the checked-in test
just to shorten it.

For one-off cases, use the configured standalone compiler:
`./mechsyscc case.cu` from the repository root (or the equivalent
`mechsyscc` in the current installation).

## FLBM / GPU fluid cases

- `lib/flbm/Domain.h` is the fluid LBM; `lib/flbm/lbm.cuh` holds the CUDA
  kernels. Use `D2Q9` with `iVec3_t(nx,ny,1)` for a 2-D GPU run.
- The GPU streaming kernel resets `BForce` every step. A constant body force
  must be written back to `dom.pBForce` on the device in the `Setup` callback
  before collision; setting the host array once is not enough.
- `dom.Nthread` controls the CUDA block size. The `Nproc` argument of
  `Domain::Solve` affects host packing/transfers on the GPU path.
- CPU and GPU body-force conventions differ: CPU `StreamSC` adds
  `0.5*dt*BForce` to the reported velocity, while the GPU fused kernel returns
  the raw second moment. The raw GPU moment carries a `-0.5*nu*phi` bias in
  permeability; use `u_physical = u_raw + 0.5*dt*fx*phi`.
- Keep one-off FLBM cases outside the checkout; do not add them to `tflbm/`
  unless they belong in the library test suite.

## Reporting

Report measured wall times, separate solve time from total wall time, and
explicitly distinguish what compiled, what ran, what matched a baseline, and
what was faster. If a verification was not done, say so.
