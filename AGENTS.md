# MechSys — instructions for future agents

This file is loaded automatically by DSH for work in this repository. It is
deliberately short and operational. Read
[`MECHSYS_HANDBOOK.md`](MECHSYS_HANDBOOK.md) before making substantive changes;
the handbook is the durable knowledge base for users and agents.

## First checks

1. Run `pwd`, `git status`, and inspect the relevant file before editing.
2. This checkout is `mechsys` at `/home/sergio/mechsys`; `mechsys` is a symlink
   to `lib/`.
3. **Do not assume `HEAD` is the optimized/correct baseline.** The working tree
   currently contains uncommitted DEM optimization changes. Never run
   `git reset`, `git checkout -- .`, or delete generated/uncommitted files
   without explicit permission.
4. Read the current optimization reports before changing DEM, CUDA, or LBM-DEM
   performance code:
   - [`CPU_DEM_OPTIMIZATION_CHANGES.md`](CPU_DEM_OPTIMIZATION_CHANGES.md)
   - [`CUDA_OPTIMIZATION_CHANGES.md`](CUDA_OPTIMIZATION_CHANGES.md)
   - [`REDUCED_TEST_CU_03_PERFORMANCE.md`](REDUCED_TEST_CU_03_PERFORMANCE.md)

## Build and test

The root CMake/Makefiles are already configured. Build from the test directory:

```bash
cd /home/sergio/mechsys/tdem
make -j8 test_periodic          # CPU Voronoi reference
make -j8 test_cu_01 test_cu_02 test_cu_03 test_cu_04

cd /home/sergio/mechsys/tlbmdem
make -j8 tlbmdem01 tlbmdem02
make -j8 tlbmdem_cu_01 tlbmdem_cu_02 tlbmdem_cu_03 tlbmdem_cu_04
```

Typical runs:

```bash
# CPU reduced Voronoi benchmark: FileKey argument; .inp suffix is appended.
cd /home/sergio/mechsys/tdem
./test_periodic /home/sergio/mechsys/benchmarks/cpu_voronoi/tper_small 8

# CUDA test: Nproc, Verlet distance (Alpha), threads per block.
./test_cu_03 1 0.05 128
```

`tlbmdem01`/`tlbmdem02` are large. For a quick integration smoke test, reduce
their grid/Tf in a temporary copy; do not modify the checked-in test just to
make it fast.

## FLBM / GPU fluid cases

- FLBM lives in `lib/flbm/Domain.h`; the CUDA kernels are in
  `lib/flbm/lbm.cuh`. `D2Q9` with `iVec3_t(nx,ny,1)` is the 2-D GPU setup.
- The GPU streaming kernel `cudaFusedStreamMacro` zeroes `BForce` every step.
  A constant body force must be rewritten on the device (`dom.pBForce`) in the
  `Setup` callback before the collision; setting the host array once is not
  enough.
- `dom.Nthread` (default 256) sets the CUDA block size. The `Nproc` argument of
  `Domain::Solve` only affects host-side packing/transfers on the GPU path.
- CPU and GPU body-force conventions differ: CPU `StreamSC` adds
  `0.5*dt*BForce` to the reported velocity, while the GPU fused kernel returns
  the raw second moment. The raw GPU value carries a constant `-0.5*nu*phi`
  bias in `K`; use `u_physical = u_raw + 0.5*dt*fx*phi` for permeability and
  state which convention a number uses.
- One-off cases should live outside the checkout and be compiled with
  `/home/sergio/mechsys/mechsyscc case.cu`; do not add them to `tflbm/` unless
  they belong in the library test suite.
- Reference case: `/home/sergio/tflbm10_porous/tclbm10.cu` (random 90%-porosity
  medium, +x body force). Measured on the RTX A6000: 512x256, 20,000 steps,
  `Nthread=256` -> 5.37 s solve / 5.64 s wall, physical `K = 1.08` (raw moment
  1.00).
- Permeability sweep: `/home/sergio/tflbm10_porous/permeability_sweep/`
  (`phi=0.50...0.85`, 512x512, 3 seeds). The fluid phase is random site
  percolation: `K=0` below `phi_c ~ 0.593`; percolating values are
  `K(0.60)=0.0422`, `K(0.70)=0.1163`, `K(0.85)=0.5441`. A `1024x1024` check
  changes `K` by <1% at `phi=0.60/0.65`.
- Plane-Poiseuille validation: `/home/sergio/tflbm10_porous/validation/`
  (`poiseuille.cu`, `run_validation.py`). Walls at `y=0.5` and `y=H+0.5`,
  uniform `+x` force. The half-force corrected `K` converges to `H^2/12`:
  relative error `1.56e-2` (H=8), `3.91e-3` (H=16), `9.77e-4` (H=32),
  `1.99e-4` (H=64), `1.63e-5` (H=128); force linearity over 16x at H=32 is
  `2e-12`. Scale `fx` with `H^2` to keep `u_max ~ 0.01`; `H=128` at `fx=1e-4`
  is high-Mach and fails. This validates the half-force convention used above.
- Advection-diffusion point source: `/home/sergio/tflbm10_porous/solute_source/`.
  Two-layer `AdvectionDiffusion` (`nu={1/6,D}`): layer 0 is the fluid, layer 1
  the passive concentration. Set `BForce[1]=0.5*fx/Tau[1]` so the AD collision
  advects at `u_phys=u_raw+0.5*dt*fx/rho`; a source must add to both the
  populations and `Rho[1]`. The entropic limiter in `lib/flbm/lbm.cuh` is
  guarded in this working tree; without the guard, negative solute populations
  produce NaN in the zero-concentration far field at 512x512. `mass_fluid`
  (pore sum) is about 14% below `mass_total` (includes bounce-back populations
  on static solids), while `mass_total` matches `S(t-t_on)`. Production
  512x512, phi=0.70, seed=12345, fx=1e-2, D=0.02, S=1e-4, t=20k..60k: 67.4 s
  wall, `mass_total=4.0`, `Cmax=1.09e-3`; the plume wraps between 55k and 60k.

## Exactness rules (critical)

- Compare against a **clean worktree of `HEAD`**, not against the modified
  working tree.
- For CPU DEM, exactness is checked at `Nproc=1` by comparing **HDF5 datasets**,
  not whole-file MD5. HDF5 container metadata is not deterministic even
  baseline-to-baseline. Expected result is `neq []`.
- At `Nproc>1`, CPU force accumulation order is thread-order dependent. Small
  last-digit differences are expected, not automatically a regression.
- CUDA `test_cu_03/04` use `atomicAdd`; the original baseline is itself
  non-deterministic. Do not claim bitwise equality for those GPU runs.
- Run the real test window, not only a short smoke test, before making a final
  accuracy claim. State the window and thread count.

Useful comparison pattern:

```python
import h5py, numpy as np
b = h5py.File("baseline.h5", "r")
o = h5py.File("optimized.h5", "r")
for k in b:
    if isinstance(b[k], h5py.Dataset) and not np.array_equal(b[k][()], o[k][()]):
        print("DIFF", k)
```

## Hot-path invariants

These are easy to break accidentally:

1. **Preserve candidate-list order.** Contact forces are accumulated in list
   order. Sorting `L` changes floating-point results even if contacts are
   identical. Sorted distances, if used, must remain separate from `L`.
2. **Keep the baseline `Distance` arithmetic.** Do not replace full distance
   calculations with cached normals or squared-distance variants in
   `_update_contacts`. `Distance2` is only an early rejection in the force loop.
3. **The CPU Verlet cutoff uses `4.0 * ContactDispGlobal()`.** A factor of two
   was observed to skip real contacts (`Ld` just outside the bound that became
   a contact after motion). `ContactDispGlobal()` accumulates the
   `MaxDisplacement` values reported each step; do not lower the factor
   without a new safe derivation and an exact `Nproc=1` comparison.
4. **Do not sort the contact list or change original `(i,j)` push order.**
5. **CUDA movement:** keep the separate translate/rotate fallback for custom
   movement and for cases where the fused movement kernel is not valid.
6. **CUDA thread blocks:** `Nthread=512` fails for the current `test_cu_03`
   register usage. Valid tested range is up to 384; large test optimum is
   around 128, reduced test optimum around 32.
7. **CRLF:** `tdem/test_cu_04.cu` has CRLF line endings. Do not rewrite it with
   an editor/tool that normalizes line endings; it creates a huge meaningless
   diff.
8. **Shared headers affect both CPU and CUDA.** After changing `lib/dem/*`,
   build at least one CPU target, one CUDA target, and the LBM-DEM targets when
   relevant.
9. **No physics shortcuts.** Broadphase rejections must be conservative; if in
   doubt, remove the rejection rather than risk a missed contact.
10. **Clean warnings.** The only expected third-party warning is Blitz
    `tinymat.h(60): warning #191-D`. Do not add new warnings in MechSys code.

11. **Rotation caches are derived state.** `Particle::wrot` and `Qc` are
    refreshed by `CalcProps`, `InitializeVelocity`, `Rotate`,
    `Domain::UpdateContacts`, and after a user `Setup` callback. The CUDA
    upload copies these fields; `DnloadParticle` copies them back. A custom
    GPU kernel that changes `Q`/`w` must refresh `DPar[].wrot`/`Qc`. If host
    code assigns `Q` or `w` directly and then calls an interaction force
    without going through the normal paths, call `Particle::UpdateRotCache()`
    first.

## Performance rules

- Measure wall time, not only printed per-step solver time. Small tests are
  dominated by construction/output overhead.
- Report baseline and optimized on the same machine, input, and thread count.
- Separate solve time from total wall time.
- For CPU Voronoi, use the reduced input in
  `benchmarks/cpu_voronoi/tper_small.inp` for exactness; use
  `tper_med.inp` for a larger timing check.
- For CUDA, re-benchmark after struct/layout changes; do not reuse stale
  numbers. See `REDUCED_TEST_CU_03_PERFORMANCE.md` for the reduced sweep.
- Recommend a thread count based on the actual problem size; the large and
  reduced CUDA workloads have different optima.

## Preferred working style

- Reply in English unless the user asks otherwise.
- Give explicit commands, file paths, and measured numbers.
- Distinguish “compiles”, “runs”, “matches baseline”, and “is faster”.
- If a verification was not done, say so.
- When adding durable knowledge, update this file or
  `MECHSYS_HANDBOOK.md`, not an ephemeral note.
