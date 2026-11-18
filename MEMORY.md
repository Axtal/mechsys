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

## 0. Reviewing an incoming DEM pull request

The upstream PR workflow that has actually been used here:

```bash
git fetch origin 'refs/pull/*/head:refs/remotes/pr/*'   # PR refs are fetchable
git worktree add --detach /tmp/pr2_base origin/main      # baseline
git worktree add --detach /tmp/pr2_head pr/2             # the PR under review
ln -sfn lib /tmp/pr2_base/mechsys; ln -sfn lib /tmp/pr2_head/mechsys
cmake -S . -B . -DA_WITH_TDEM=ON -DA_MAKE_USE_GPU_DOUBLE=ON \
      -DA_USE_GMSH_STATIC=ON -DA_USE_HDF5=ON -DA_USE_OMP=ON -DA_USE_VTK=OFF
```

`git merge-tree --write-tree origin/main pr/N` dry-runs the merge without
touching a worktree. `gh` is **not** installed and the GitHub REST API is
rate-limited from this host; `git ls-remote origin` does work.

Verification done for PR #2 (2026-11, Marco Previtali) and reusable:

- CPU `test_periodic` at `Nproc=1`, Voronoi and `sphereboxhcp`, 400 and 2000
  steps, HDF5 datasets compared recursively with `h5py.visititems` (the files
  are group-per-particle, so a root-level loop silently compares nothing).
- CUDA `test_cu_01` (Tf=10) and `test_cu_02` (Tf=30) full windows, `test_cu_03`
  reduced, `test_cu_04` from `tests/dem_contact/inputs/test_cu_04.inp`. Build
  the PR's regressions with the group target `dem_contact_regression_binaries`
  (individual target names are `dem_test_<name>_cpu`). After 2026-11-13 the PR
  as reviewed adds no CMake option: `tests/dem_contact` is configured from the
  root file whenever `A_WITH_TDEM` is on (the default), the two
  `A_WITH_*` regression/coupling switches are gone, and the experimental
  coupled VV session (`lib/lbmdem/experimental/`, `tlbmdem_cu_benchmark`,
  `dem_test_coupling`) was deleted with them, so a stale build directory must
  be reconfigured and the registered test count is 10. The PR's own
  `scripts/verify_dem_contribution.py` then passes every stage.
- LBM-DEM: a reduced `tlbmdem01` copy (61x31x31, `Tf=100` -> 6250 steps)
  compiled with each tree's own `mechsyscc` is a cheap coupling check.

**Comparison hygiene.** Sort every HDF5 dataset diff by field name, and exclude
`Diam`, `Interactons/Fdvv`, `Interactons/Fdr` — they come from uninitialised
buffers and differ between two runs of the *same* binary. Both trees must be
run in their own output directory because `test_periodic` reuses a fixed
FileKey prefix. Do not write these runs to a small tmpfs: one sphere
`<key>_a.hdf5` has ~430k datasets and a full `test_cu_02` writes 300 files.

**PR #2 measured outcome** (Marco Previtali, three commits on `789a2de`,
"velocity Verlet, sphere Coulomb modes, contact output"). Mergeable, builds
clean (only the Blitz `tinymat.h` warning), and its eleven regressions pass
(88 CPU analytic contact cases, 32 fast-path cases, 5 CPU ctests).

- The default DEM numerics are preserved to floating-point reassociation only.
  `Nproc=1`, 2000 steps: Voronoi `x` 1.9e-14, `v` 6.0e-9, `w` 1.8e-7 relative;
  sphere `x` 2.4e-8, `v` 5.0e-5, `w` 1.5e-4 relative. Kinetic energy agrees to
  6 digits in both. CUDA `test_cu_01` `Observed L` is identical, `test_cu_02`
  final `Position` 1.6e-8 relative. LBM-DEM `tlbmdem01` (reduced, 6250 steps)
  is **bit-identical** including `_force.res`. There is no bitwise `neq []`
  for the CPU DEM paths, so the handbook's exactness gate is not met.
- **The whole sphere-path divergence is one change:** the new
  `Fdvv = 0` when `delta<=0` in `CInteractonSphere::CalcForce`. Rebuilding the
  PR with *only* that reset removed gives, over the same 2000 steps,
  `x` 1.7e-13, `v` 3.4e-10, `w` 6.0e-10, `Q` 3.3e-12 relative against the
  baseline - i.e. pure reassociation, the same magnitude as the Voronoi path.
  It also restores the saved friction-history row count exactly (baseline
  24467, PR 20266, PR-without-reset 24467). Everything else on the CPU sphere
  path (`F` reordering, `Epot`, `dEfric`, mode plumbing) is roundoff only.
- Two further default-path changes are *reporting* only (`Epot`/`dEfric` feed
  `CalcEnergy` and the `.res` file, never the dynamics): `Epot` is now
  accumulated, so `tper_sph` total energy moves 2.45 -> 2.55 (4%),
  `test_cu_02` grows by `Epot`, and CUDA `test_cu_03` (polyhedral, which had
  **no** `Epot` in `lib/dem/dem.cuh` at all in the baseline) reports
  `Epot=0.0739` instead of 0; `dEfric` becomes
  `friction_limit*|plastic_disp|` instead of `Kt*dot(Fdvv,vt)*dt`, but the
  `dem.cuh` polyhedral kernel keeps the old `dotreal3(ft,vt)*dt` form, so the
  diagnostic is now inconsistent between paths.
- **Do not read a dynamics change out of `test_cu_03` energies.** Two runs of
  the *same* baseline binary at `Nproc=1, Alpha=0.05, Nthread=128, Tf=10`
  (1,000,001 steps, 959-972 s wall) gave `Ekin` 0.0790183 and 0.0082049 - a
  10x run-to-run spread from `atomicAdd`. The PR's 0.00636631 is inside that
  spread. Only `Epot` moving off zero is conclusive there.
- `MostlySpheres` dispatch changed. `PairtoCInt` is now also populated when
  `ContactLaw==1`, and the linear `CalcForceSphere` fast path is gated on
  `ContactLaw==0`. Do **not** describe this as "hertz silently used a linear
  law" in the old code: with `MostlySpheres=true` a sphere-sphere pair leaves
  `PairtoCInt` empty, and `ResetContacts` then does
  `PairtoCInt[hash]->UpdateContacts(...)` on a NULL map entry, so the baseline
  **segfaults** for both laws (verified with a two-sphere probe, both trees,
  `gdb` frame `DEM::Domain::ResetContacts()`). The PR only makes
  `MostlySpheres`+hertz work; `MostlySpheres`+linear still segfaults in both.
  `CalcForceSphere`'s linear law is effectively dead code for contacting
  spheres.
- `test_cu_04`'s `MoveTopPar` is the one test-side physics change in the PR: it
  replaced the baseline `Translate(Verts,Vertso,Par[ip],DPar[ip],trans)` wall
  push with a `v.z=Dz/dt` + `xb.z=x.z-Dz` prescription so that velocity Verlet
  also receives `Dz`. The two mechanisms disagree (`5e-3` in `PVelocity`) and
  applying both moves the wall twice, so in this checkout the branch is on the
  integrator (`UdCu.vv`): PV keeps the baseline call and reproduces
  `0 / 0.00914244 / 0.00244635` bit-for-bit, VV uses the prescription. Without
  that branch the ported PV path silently differs from the baseline in phase 2
  (`0.012008` instead of `0.00244635`).
- The PR's `WriteBF` writes the rolling force into the existing `Rolling`
  dataset (the baseline put the normal force there) and widens which contacts
  are exported, so `_bf_` files are not backward compatible. `Rolling` is only
  emitted when `n_rl>0` (contacts where *both* particles have `Verts.Size()==1`),
  so a sphere on a quad plane never produces it; the four added datasets
  (`NormalOverlap`, `TangentialOverlap`, `DashpotForceNormal`,
  `DashpotForceTangential`) appear in every `_bf_` file.

**PR #2 on `test_01`/`test_02` (the canonical single-body DEM cases).** The
pairing is `tdem/test_01` = CPU sliding block, `tdem/test_cu_01` = GPU sliding
bodies, `tdem/test_02` = CPU turning table, `tdem/test_cu_02` = GPU turning
table; the two "02" variants write `<FileKey>_trayectory.res` via a `Report`
callback (columns `Time x y xt yt`, 300 rows at `dtOut=0.1`). Defaults already
match across trees (Tf=10, Tf=30, Nproc=1), so run the baseline with no
arguments and the PR as `test_cu_0N 1`.

- `test_02` CPU at `OMP_NUM_THREADS=1` is **bit-identical** between baseline and
  PR across all 300 `.res` rows and all 300 frames: `max|dx|=max|dy|=0` over
  3e6 steps, final `x=-3.726376620`, `y=-0.581608691`. The ball rolls in
  unbroken contact, so the PR's `delta<=0` history reset never fires.
- `test_01` CPU is roundoff only: final `Position` 1.8e-11 and `Velocity`
  3.7e-11 relative; `AngVelocity`/`PAngacceleration` look large in relative
  terms only because those components are ~1e-10.
- GPU `test_cu_02`: the **baseline** departs from the CPU reference by
  `max|dr|=1.37e-4` (285/300 rows), while the PR agrees with the CPU to
  `1.3e-8`. The PR improves GPU/CPU agreement on this case by four orders of
  magnitude; do not present this as the PR breaking GPU results.
- `PEkin` is the one CPU field that never matches: it is
  `float(P->Ekin + P->Erot)`, cached inside `Translate`/`Rotate`, and the PR
  refreshes those caches at a different point, so it differs on every frame
  after the first by <=9.8e-7 relative *while `Position`/`Velocity`/`Quaternion`
  are bit-identical*. Do not treat a `PEkin` diff as a dynamics regression.
- Measured at `Nproc=1`/`OMP_NUM_THREADS=1` (wall, single runs, sequential):
  CPU `test_01` 2.37 s base vs 2.58 s PR; CPU `test_02` 5.24 s vs 5.80 s;
  GPU `test_cu_01` 81.30 s vs 72.41 s; GPU `test_cu_02` 221.10 s vs 207.86 s.
  The CPU cases are two particles, so the 9-11% CPU delta is option plumbing,
  not a trend; the GPU numbers are single runs on a shared A6000.

**`test_cu_04` with `tper.inp`: the noisy-case protocol.** The correct command
is `./test_cu_04 <filekey> <Nproc>` (`tdem/test_cu_04 tper 1`): it reads
`<filekey>.inp` from the current directory and writes `<filekey>_b_walls.res`
with `Time sx ey`. Do **not** confuse it with `tdem/test_cu_03`, which takes
`Nproc Alpha Nthread` and writes no `.res` at all — `./test_cu_03 tper 4`
silently means `Nproc=atoi("tper")=0` and `Alpha=atof("4")=4.0`, i.e. a wrong
run that still starts. `tper.inp` packs **18,616 spheres** (the log's
`Num of particles`; `nx/ny/nz` are *not* the particle count, so do not guess the
size from the input file), with `L=5.2`, `R=0.1`, `dt=5e-5`, `T0=10`, `Tf=100`,
200001 + 2200000 steps; the
full run takes ~25 min per tree on the GPU at `Nproc=1`, so keep `tper.inp` and
every output in a scratch directory — `tdem/tper_b_walls.res` in the checkout is
an untracked **reference** left by an earlier full run and will be overwritten
by a run started from `tdem/`.

This case is chaotic, and the baseline is **not** reproducible run to run
(`cmp` of two baseline `.res` files differs; `atomicAdd` force accumulation).
Do not compare final values: the last `sx` point differs by 5% between two
baseline runs. Compare time averages over the last N points instead, and quote
the baseline's own repeat as the noise floor. Measured at `Nproc=1`, rms over
all 1000 rows: base-vs-base2 `sx` 0.0487 (max 0.181), base-vs-reference 0.0521,
base-vs-PR 0.0525 — i.e. the PR is indistinguishable from a re-run of the
baseline. Time-averaged `sx` over the last 100/300/500 points: base-PR
4.3e-4 / 6.4e-4 / 6.3e-3 against base-base2 1.5e-2 / 7.3e-4 / 1.3e-4, so
averaging collapses a 12% pointwise scatter to <=0.1% and the residual is at
the level of the baseline's own variability. `ey` has a small mean and larger
relative fluctuation; use the same treatment.

`test_cu_04` is also **~45-60% slower in the PR on this case**, and this is not
a build artifact (identical `CUDA_FLAGS`, `-O3 -DUSE_CUDA`, `sm_86`, equal step
counts). Reduced input (`T0=2`, `Tf=3`, 100k steps), alternating repeats:
baseline 77.01 / 76.93 s, PR 122.61 / 122.46 s; full run 1477 / 1498 s vs
2226 s. Note `USER` time rises with `WALL`; that is CUDA driver spin-wait, so it
is *not* evidence of host-side cost.

Decomposed with a stable 40k-step input (`T0=Tf=1`, ~8 s of packing included):
baseline 37.6-38.8 s, PR 55.0-55.7 s, i.e. +45%. Each number below is a
temporary rebuild (all reverted; both worktrees verified clean with
`git diff --stat`):

| variant | wall | share |
|---|---|---|
| baseline | 37.8 s | - |
| baseline with `ComInteractonCU` padded 64 -> 176 B | 42.3 s | ~9% |
| PR with the 8 new `Reset` writes disabled | 45.8 s | ~17% saved |
| PR with the diagnostic atomics off | 54.8 s | ~0% |

So: `ComInteractonCU` grew 64 -> 176 B and `DynInteractonCU` 136 -> 152 B
(hurts every strided access, ~9%); the `Reset` kernel now zeroes 8 extra fields
for **every candidate pair** each step, not just for active contacts (~17%);
the remaining ~18% is the force-kernel additions (writing `Ft=0`/`InContact`/
`PrevDelta` on the `delta<=0` path for every non-contacting candidate, plus
new branches). The per-contact atomics I first suspected are irrelevant.

All of it is **reporting-only state**, so it is removable: `Epot`, `dEvis`,
`dEfric`, `Fndpot` and `Ftdpot` can be derived at output time from `Fnnet`/
`Ftnet`/`NormalOverlap` (energies are only needed in `CalcEnergy`, which runs at
output frames), which shrinks the struct back to ~88 B and removes both the
extra `Reset` writes and the extra force-kernel stores. Keep `Nc` and
`NormalOverlap` (12 B) if the contact-count output needs them.

Profilers are unavailable on this box: `ncu` returns `ERR_NVGPUCTRPERM` (no
performance-counter permission), `nsys` produces a `.qdstrm` but the
`QdstrmImporter` is not installed, and `libcudart` is linked statically so an
`LD_PRELOAD` shim cannot intercept the CUDA runtime API. Isolate by rebuilding
variants and timing them.

**Fix that works: keep the baseline kernel bodies as the default and make the
PR versions opt-in.** The CUDA contact kernels are already dispatched through
function pointers (`pForceVV`, `pReset`, `pForceEE/VF/FV`) assigned in the
`Domain` constructor, so the PR should have *added* `*_EnergyBalance` kernels
instead of modifying the originals. Implemented on the `pr/2` worktree: original
bodies keep the plain names, the PR versions are `CalcForceVV_EnergyBalance`,
`CalcForceVV_Hertz_EnergyBalance`, `Reset_EnergyBalance`,
`CalcForceEEPre_EnergyBalance`, `CalcForceVFPre_EnergyBalance`,
`CalcForceFVPre_EnergyBalance`; `Domain::SetEnergyBalance(bool)` switches them
(default false), and the four energy-balance datasets plus their XDMF attributes
are written only when the flag is set. Call this feature "energy balance", never
"diagnostic". On the CPU path `EmitEnergyBalance()` returns true
unconditionally, because the CPU force
code fills those fields anyway and the PR's CPU output must not change. The
change is uncommitted on local branch `pr2-optional-contact-kernels` in the
permanent worktree `/home/sergio/mechsys-pr2` (branch point `pr/2` = `74ca119`),
because `/home/sergio/mechsys` is on `main` at `ca43286` and does not contain
the PR code at all. That worktree is plain-configured and `tdem/test_cu_04`
builds and runs there (0 compiler errors; 40.67 s on the 40k-step input with
`Epot = 0` and the baseline 7 datasets); its four modified files differ from the
timed copies only in trailing blanks and one relocated comment.

Same 40k-step input as above (`T0=Tf=1`, GPU, `Nproc=1`):

| build | wall |
|---|---|
| baseline | 37.8-38.8 s |
| PR as submitted | 55.0 s |
| PR + opt-in kernels, default | **42.1 s** |
| PR + `SetEnergyBalance(true)` | 54.8 s |

Full `tper.inp` window (2.2e6 steps): baseline 1477 / 1498 s, PR 2226 s,
default kernels 1661 s; the `_walls.res` time averages of the default agree
with the baseline at the level of a baseline repeat.

The default now reports `Epot = 0` and writes the baseline 7 datasets, so the
behaviour matches the baseline exactly, and the kernel cost is fully recovered.
The leftover ~4 s (11%) is the `sizeof(ComInteractonCU)` growth alone: padding
the *baseline* struct to 176 B gives 42.3 s, i.e. the same wall time. Removing
that needs the energy-balance fields moved to a separate optional array
(`ComInteractonDiagCU`) so the default struct stays at 64 B; a kernel switch
alone cannot recover it. Add the hook with care in `tdem/test_cu_04.cu`: it has
CRLF line endings, so edit it as bytes (`open(...,'rb')`/`'wb'`).

**Velocity Verlet (`UseVelocityVerlet`).** Every DEM and LBM-DEM test takes it
as a trailing `<integrator> <energy>` argument pair, and `Domain::Solve` prints
`Velocity Verlet integrator = on|off` and `Energy balance tracking = on|off`
right after the `Verlet distance` line (the coupled `LBMDEM::Domain::Solve`
prints the same for its DEM sub-domain). Report it from the library, never from
the individual drivers -- duplicated per-test echoes go stale.

Cost on `tper.inp`, GPU, `Nproc=1`, wall: 40k-step input 42.01 s default vs
42.99 s VV (3 interleaved repeats each, +2.3%; compression phase +5.0%, shearing
main phase +0.6%), full 2.2e6-step window 1661.2 s vs 1701.0 s (+2.4%). On the
2-particle `test_cu_02` (3e6 steps) it is 205.6/207.6 s vs 223.8 s (+8.8%),
because VV's extra per-step kernels cost about 6 us/step and there is almost no
force work to amortize them against.

Accuracy: the `tper` packing is chaotic, so use the averaged `.res`. VV's
whole-window mean `sx` is 1.141080 vs reference 1.145234 (-0.36%) and default
1.149644 (+0.39%); pointwise rms vs the reference is 0.0619 against a
baseline-repeat floor of 0.0487, and the last-100/300/500/1000 averages sit
1.8e-2/1.5e-2/1.2e-2/4.2e-3 below the reference (the default sits within
+/-5.7e-3). That is the size of the default integrator's own discretization
error, not a VV defect; resolving which is closer to the truth needs a
`tper` dt-refinement, which was not run.

The clean accuracy measurement is `test_cu_02`, whose `.res` carries the
analytic rolling trajectory: over 30 s at `dt=1e-5` the default integrator is
6.76e-3 rms off the analytic curve (max 1.17e-2) and VV is 1.05e-3 (max 1.84e-3),
i.e. **VV is 6.4x more accurate**, and repeated default runs are bit-identical,
so that ratio is not noise. Figure:
`/home/sergio/DATA/pr2_tper4/tt_integrator_error.png`.

**VV is rejected by the coupled LBM-DEM solver.** `LBMDEM::Domain::Solve` begins
with `if (DEMDOM.UseVelocityVerlet || DEMDOM.CheckpointReady) throw new
Fatal("LBM-DEM Solve uses legacy integration; use the opt-in experimental
benchmark coupling API for supported VV cases. Coupled restart is unsupported")`,
so neither `tlbmdem01/02` nor `tlbmdem_cu_01/02` can run with VV, and in
particular the Stokes comparison in their `_force.res`
(`|Flbm|`, `|v|`, `6*pi*nu*R*|v|`) is unavailable with VV. The only supported
route is `tlbmdem/tlbmdem_cu_benchmark.cu`, which is gated behind
`-DA_WITH_EXPERIMENTAL_VV_COUPLING=ON` (default OFF) and must be configured in;
usage `./tlbmdem_cu_benchmark vv|pv 01|02|03|tiny:FIXTURE STEPS CADENCE OUTDIR
--experimental` (the trailing `--experimental` token is mandatory, the output
directory must not exist). Fixture `02`/`moving` is the Stokes case: its
`report.csv` carries `Flbmx/y/z` and a `reference` column equal to
`6*pi*nu*R*|v|`.

Measured there (`Nproc=1`, D3Q15, 200^3 lattice, nu=1.6, R=1.8, dt=0.016,
driving force `m*0.1 = 2.4428`, so ideal Stokes terminal velocity 0.0450), 2000
steps: force balance `|Flbm|/Ff` reaches 0.9990 (PV) and 0.9995 (VV); drag at
matched velocity differs by 0.1-0.4% (VV lower); extrapolated terminal velocity
0.0402 (PV) vs 0.0409 (VV); effective `6*pi*mu*R` is 62.2 (PV) and 59.8 (VV)
against the ideal 54.3, i.e. both stay 10-15% above ideal Stokes because the
fixture sphere is only 3.6 cells across. So VV does not degrade the Stokes
comparison, and the absolute offset is spatial resolution, not integration.
Timing is unmeasurable there: the coupled 200^3 D3Q15 step is 65-77 ms, so VV's
~6 us/step is under 0.01%, and repeated walls scatter by about 6% (67.9/74.9,
65.6/66.6 and 150.9/142.1 s for pv/vv at 600/600/2000 steps).

**PR #2 inventory, for review or rebase.** Base is `789a2de` ("MPM by Mario
merged into mechsys"), head `74ca119`, three commits (`8f867a3` feature,
`4308abf` Coulomb projection, `74ca119` build/docs), 51 files, +8658/-397. The
branch is *not* rebased: current `main` (`ca43286`) carries 18 files of MPM
CUDA work the branch lacks, so `git diff main pr/2` shows `tlbmmpm_cu_0*.cu` and
the MPM headers as deletions; diff against the merge base to see the PR itself.
Build options: `A_WITH_DEM_CONTACT_REGRESSIONS` (OFF) adds `tests/dem_contact`
with ten ctests; adding `A_WITH_EXPERIMENTAL_VV_COUPLING=ON` makes eleven by
registering `dem_coupling_smoke`. `dem_test_history_cpu/gpu` are built but not
registered with ctest -- `tests/dem_contact/check_histories.py` runs them from
`scripts/verify_dem_contribution.py`, so they are covered, just not by ctest.

**PR #2's `lib/dem/checkpoint.h` (386 lines) is a VV restart, not a duplicate of
Save/Load.** `Domain::Save`/`Load` is the HDF5 output dump: per particle `x, xb,
v, w, wb, Q, V, I, m, Diam, Dmax, Rho/SR/Flags/Index, Tag, Ff, Tf` plus the
polyhedral geometry; per contact `Fdvv`, `Fdr` and the polyhedra maps EE/VF/FV
(as `key_EE/val_EE`, ...) with stored `Kn/Kt/Gn/Gt/Mu`; plus `Time`, `Iter` and
the bbox. That is FP64 and bit-exact (`tdem/test_saveload`). `Load` restores all of it, so the live
sphere and 3-map polyhedral friction history is already in Save/Load --
`tdem/test_saveload.cpp` (main, commit `bd2409a`) proves the round trip is
bit-for-bit and then continues both runs. I verified it is a full **VV** restart
too: copying that test with `UseVelocityVerlet=true` on both the original and
the reloaded domain gives max|dx|=dv=dxb=dw=dwb=dQ=0 after 100 further steps,
`Fdvv`/`Fdr` identical, PASS. What the checkpoint actually adds, corrected: the
new contact scalars `PrevDelta`/`InContact` (live only with
`SphereFirstContactCorrection`/`SphereTensileCutoff`, both default off;
`PrevDelta` also feeds the Hertz `Ft` upload on CUDA), the VV arrays
(`v_half`, `w_half`, `wa`, `A`, `Wdot`, `F`, `T`) which are reconstructible and
demonstrably exact without them, `vt_prev` (**never read** anywhere in the
library -- dead state), the VT/TV/VC/CV maps (torus/cylinder, and the checkpoint
refuses tori/cylinders, so unreachable), the candidate list by persistent id,
`NextPersistentId`, and the whole run configuration (`Dt`, `Alpha`, `Beta`, `Evis/Efric/Wext`, thermostat
range, `ContactLaw`, `SphereCoulombMode`, `RotPar`, tensile, activation,
`FricCoeff`). Format is boost property_tree text at `setprecision(17)`, with
`magic/version/precision/phase` and all-or-nothing `need()` refusals; it refuses
non-stock kernels, bonds, `MostlySpheres`, torus/cylinders and `UserData`, needs
a fresh empty Domain and the same `dt`, and only saves at a completed VV full
step. Fidelity is tested by `check_histories.py`: 32-step reference vs 16 +
resume for 16, matching within `1e-10 + 1e-9*|x|` on CPU and GPU. My review
position: the capability is needed, the third serialization format is not -- it
fits in `Save`/`Load`. `Save` is an explicit snapshot (`<FileKey>.hdf5`,
`H5F_ACC_TRUNC`), *not* the per-frame visualization writer: that is `WriteBF`,
which writes float32 `_bf_<idx>.h5` plus XDMF, and `Solve` no longer calls
`Save` at all. `Load` already restores 13 props per particle from a flat `(NP,13)`
array and tolerates missing datasets through `H5LTfind_dataset` (`/Time`,
`/Iter`), so additive datasets keep old files readable.

Folded in, in `/home/sergio/mechsys-pr2` (branch
`pr2-optional-contact-kernels`): `Save` also writes `/RestartVersion`, the run
scalars, `/RestartPersistentId`, `/RestartForce` (NP,12: `F,T,A,Wdot`),
`/Interactons/PrevDelta`, `/Interactons/InContact`, and `/ContactOrder`+
`/ContactData` (the active list with its 17 exported values); `Load` restores
them. There is no separate restart API: the PR's `lib/dem/checkpoint.h` and its
`SaveCheckpoint`/`LoadCheckpoint` wrappers are gone, and the helpers plus
`AssignPersistentIds`/`SaveRestartState`/`LoadRestartState` live at the end of
`lib/dem/domain.h` (one file format, nothing extra for a caller to know).  The
helpers are plain `inline` functions in `DEM` -- the helper namespace wrapper was
removed on request, and they are not class members either, since they are
stateless and only take the `hid_t`.
Two facts the
386-line version carried that a bare "80 lines" estimate misses, both learned
the hard way: (1) a report taken right after `Load` prints per-particle `F`/`T`
and the per-contact decomposition, so the last force evaluation **and** the
active contact list are part of the visible state -- without them
`check_histories.py` fails on its first row with `F=0`; (2) on CUDA
`DnLoadDevice` does not bring `T`/`F`/`A`/`Wdot` back to the host and
`UpLoadDevice` throws `VV device upload requires initialized acceleration
state`, so the writer mirrors `bParticlesCU`/`bDynParticlesCU`/`bA`/`bWdot`
first. `/EnergyBalance` must be stored too: it selects the kernel variant, and
without it a resumed run stops refreshing `Fndpot` (it stays frozen at the
restored value). No device history dump is needed: `UpLoadDevice` seeds
`Ft`/`Fr`/`PrevDelta`/`InContact` from the restored CPU contacts, so
`RestoreCheckpointDevice` and `CheckpointVV..` can be deleted. Verified after
all the deletions: `test_saveload` (PV) and the same source with
`UseVelocityVerlet = true` are bit-identical across save+load on CPU and CUDA,
an old 73 KB file still loads, and `tdem/test_cu_04` with
`UseVelocityVerlet=1 EnergyBalance=1` runs to 501+1501 iterations with
`Total energy = 0.00913382 / 0.0117137`.

**The opt-in must be explicit in every driver that wants the decomposition.**
`tests/dem_contact/*` was the PR's regression tree and is deleted; the drivers
that remain take the flag from the command line (`tdem/test_cu_04` has
`[integrator] [energy_balance]` as its 5th/6th optional arguments) or call
`d.SetEnergyBalance(true)`. Without it the baseline kernels never fill
`Fndpot`/`NormalOverlap`, so a reader of the HDF5 file sees
`object 'NormalOverlap' doesn't exist`.

**Coupled LBM-DEM:** velocity Verlet is available on the CUDA path only
(`SetDEMIntegrator`), and all six tlbmdem programs take `<integrator>
<energy>`. Reduced-copy checks: `tlbmdem_cu_01` and `tlbmdem_cu_03` produce
`_force.res` files identical between PV and VV, `tlbmdem_cu_02` agrees on the
terminal velocity to 2.5e-5, and the CPU coupling exits with
`Fatal: LBM-DEM velocity Verlet requires the CUDA path`. `tlbmdem_cu_04` turns `RenderVideo` on in its own source, so a full run writes
>10 GB; and any `Nl>1` FLBM output used to emit `H5Dcreate2 ... name already
exists` because `WriteXDMF`, `WriteXDMF_DEM` and `Save` all created
`Nx`/`Ny`/`Nz`/`Ncell` inside their layer loop -- fixed with `if (j==0)`, so
those traces are not evidence of a data problem (the frames were complete). Those coupled frames are large (a 100^3 field set is ~30 MB per
frame, 400 frames = 12 GB); smoke-test with a reduced copy, `Render=true` ->
`false` and a small output cadence, and check `df -h` first.

**CPU velocity Verlet exists and is worth using.** The stage selection is not
`USE_CUDA`-specific: `CalcForcesVV`/`PreForceVV`/`PostForceVV` have full host
bodies (`TranslateVelVerlet`, `RotateVelVerlet`, `UpdateVelocityHalf/Full`,
`UpdateAngularVelocityHalf/Full`). `test_saveload`, `test_01`, `test_02` and
`test_periodic` expose it. Verified: `test_01` gives `Observed L=51.0012` (PV)
vs `51.0013` (VV) against `Theoretical L=51.0204`; `test_02`'s rolling orbit is
`1.175e-2` (PV) vs `1.795e-3` (VV) away from the analytic `xt/yt` (6.5x better
for VV); the CPU energy-track run is
bit-identical to plain VV (no atomics on the host); `test_saveload` passes in
all four vv/eb combinations. Cost at 4 threads on a 216-sphere packing: PV
6.96 ms/step vs VV 4.32 ms/step, because the legacy CPU PV force pass
zeroes/reduces per-thread `MTD[i].F/T` arrays while VV accumulates in place; at
`Nproc=1` they are equal (14.76 vs 14.74 ms/step). The coupled CPU loop still
refuses VV.

**PV vs VV accuracy, measured (`test_cu_04`, dt sweep to dt/8).** At the stock
`dt=2e-5` PV's kinematic energy is +1.03 % off its `dt->0` value after the
compression phase and +25.5 % off after the driven phase; VV is dt-flat there
(<=0.15 %, the chaotic-packing noise floor) and its compression-phase limit
matches PV's to 0.004 %. So VV is the better production choice at normal step
sizes, and PV needs about `dt/4` for 0.1 %. The 6x difference in the end-of-run
energy between the schemes is a chaotic rattle attractor, not an integrator
error: height (2e-6 relative), coordination number (~1 %), mean and max normal
overlap (~1 %, 0.1 %) agree, `test_cu_01` gives `Observed L=51.0019` vs theory
`51.0204` for both, and `test_cu_03`'s steady-state energy agrees to 0.15 %.
Cost is within noise (`test_cu_03` +0.6 %, coupled 64^3 -0.7 %, `test_cu_04`
-0.9 %) except the latency-bound single-body `test_cu_01` (+3.7 %, two extra
kernel launches at 73 us/step).

**The PV path was verified bit-identical to the pre-change library in main.**
`tdem/test_cu_04 cu04 1 4.0 0.0` (three runs each; `cu04.inp` is a local
reduced copy of the checked-in `tpercu.inp` -- the driver reads `<key>.inp`, so
keep that fixture beside the run) gives `0 / 0.00914244 /
0.00244635` for both the pre-port build and the ported default; the coupled
Stokes case writes a byte-identical `cur_stokes_force.res`; `tdem/test_01`
prints `Observed L = 51.0012` in both. The dispatch itself costs
**5.3 ns/iteration** (1.78 direct vs 7.10 through member pointers, `-O3`,
2e8 iterations), i.e. 0.007 % of the 73 us `test_cu_01` step, so the three
member-pointer stages stay. Energy tracks are a recording switch: the
`*_EnergyBalance` kernels return the baseline forces, and `CalcEnergy` now
actually sums `Evis`/`Efric` (they previously always read zero). `test_cu_04`
with the tracks matches the default strain to ~1.2 % and VV to ~1.6 %, the
expected chaotic-packing spread, and the coupled terminal velocity agrees to
5e-5. `test_saveload` passes bit-for-bit in all four vv/eb combinations; sum
the domain totals in sorted-key order because `PairtoCInt` is an
`unordered_map`.

**`Domain::Solve` has no integrator test inside the time loop** (done in
`/home/sergio/mechsys-pr2`, uncommitted). Two member function pointers,
`pPreForce` and `pPostForce`, are chosen once at the top of `Solve` from
`UseVelocityVerlet`: VV gets `PreForceVV` (half kick, drift, post-drift list
check, forces with the half-step velocities) plus `PostForceVV` (velocity
finalisation); PV gets `pPreForce = NULL` and `PostForcePV` (one force pass and
the move). Both `#ifdef USE_CUDA` and CPU bodies exist, Allman braces, and the
contact-list statistics (`NumUp`, `IterUp`, `IterUpBase`) had to become members
because the phases touch them. Behaviour-preserving by construction and by
measurement: `test_saveload` PV/VV bit-identical, `test_cu_04` same energies,
coupled Stokes `.res` byte-identical to the pristine tree.
`LBMDEM::Domain::Solve` now drives the same three DEM stages, so the coupled
solver runs either integrator: `DEMDOM.SelectStages()` once before its loop,
then `pPreForce`/`pCalcForce`/`pPostForce` per step, with
`InitVerletAccelerations()` plus `UpLoadDevice(Nproc,true,true)` as the
velocity-Verlet bootstrap and only `CheckpointReady` refused. Two traps found
there: forgetting `SelectStages()` calls a NULL member pointer (segfault at the
first step), and moving the fluid imprint from the top of the step to
"between forces and move" makes the drag diverge (|Flbm| up to 5.7e4).

**PR scaffolding removed** (all uncommitted on `pr2-optional-contact-kernels`):
every PR-added file without a namesake in `/home/sergio/mechsys` is gone --
`lib/dem/checkpoint.h`, `lib/dem/sphere_cube_robust.h`,
`lib/lbmdem/experimental/*`, `scripts/verify_dem_contribution.py`,
`tests/dem_contact/*` (22 files), the four example drivers
(`test_2spheres{,_gpu}.cpp`, `PID_v2.cpp`, `TX_PID_periodic.{cpp,cu}`,
`shearbox_LE_const_vol.cu`) with their `.inp`, and
`tlbmdem/tlbmdem_cu_benchmark.cu`. `test_2spheres*`/`PID_v2`/`TX_PID_periodic*`/
`shearbox_LE_const_vol` were then dropped from `tdem/CMakeLists.txt` and the
`tests/dem_contact` subdirectory from the root one. No `docs/` directory is
added: MechSys keeps its documentation in the root markdown files
(`MECHSYS_HANDBOOK.md` for the feature design and the test switches,
`MEMORY.md` for the distilled facts, the `*_CHANGES.md` reports for
optimisation work), so the ported feature description lives in the handbook.
Net: 46 files, +1477/-6795. Rebuild everything (`make -j8` for the tdem and
tlbmdem targets) after such a sweep; in an in-source build the deleted targets
leave stale untracked executables behind (`tdem/PID_v2`,
`tlbmdem/tlbmdem_cu_benchmark`, ...) that must be removed by hand.

This file is intentionally shorter than the handbook. It records the things
that are easy to get wrong and the shortest path to a trustworthy result.

## 1. Orientation

- The real headers are in `lib/`; `mechsys` is a symlink to `lib/`. Include
  headers as `<mechsys/dem/...>`.
- On Windows, use WSL and the Linux build; MechSys runs there. No native
  Windows port is maintained.
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
- MPM XDMF: grid is `MPM_points` (`Polyvertex`), attributes `Tag`,
  `Velocity`, `Pressure`, `Sigma`, `Strain`; use
  `DefineScalarExpression("speed","magnitude(Velocity)")` and a Pseudocolor
  plot with `pc.pointSizePixels`.
- `SaveWindowAttributes.format` must be `sw.PNG` (not the string `"png"`);
  3.3.2 has no `resolution` member.  A CLI script's `outputDirectory` is not
  honored, so `os.chdir` or move the PNGs afterwards.
- Pseudocolor fixed limits: `limitsMode` only accepts 0/1
  (`OriginalData`/`ActualData`); set `pc.minFlag=1`, `pc.maxFlag=1`,
  `pc.min/max` to force a scale.  For water pressure use
  `AddPlot("Pseudocolor","Pressure")` and a `turbo`/diverging table.
- `SaveWindow` appends a frame number to the requested filename.
- A nonzero CLI exit after a successful save can be a warning; inspect the
  log and PNGs.
- A positive-`Tag` threshold includes any storage particles an application
  keeps with positive tags.  Crop by bounds if the render should show only
  the active region.

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

`Domain::CalcEnergy` loops over **all** `Particle` objects, including free
particles that have left the region of interest.  If only an active subset
matters, filter it explicitly and compute energy from that subset.

A single-smallest-particle energy threshold is also too strict for a large
assembly: many particles at a few mm/s can exceed one particle's energy.  Use
`0.5*M_active*v_settle^2` for a post-check plus a `vmax` bound, but still pass
the strict smallest-particle `minEkin` to `Domain::Solve`; otherwise the solver
can stop while one fast particle is still moving and the post-check will reject
the step.

`Domain::BoundingBox` likewise scans free particles across the whole domain;
inspect and clamp its result before using it as a top/clearance estimate.

### 2b. Periodic boundaries in DEM

`DEM::Domain` supports periodic boxes through `Xmin/Xmax`, `Ymin/Ymax`,
`Zmin/Zmax`.  `Domain::Solve` computes `Per` from those bounds;
`Domain::ResetDisplacements` wraps free particles, and the linked-cell
neighbour search uses the same periodic flags.  `CInteracton` applies `Per`
only when both particles are free, so fixed walls are never periodic.  To make
one direction periodic, remove the fixed planes in that direction and set the
corresponding bounds before `Solve`.

For finite-size particles, choose a separate periodic length several particle
lengths across; do not use a narrow gap as the period.  Fixed base/feed planes
must span the periodic length plus `4*maxd`, and the lost-particle test must
use the same length; otherwise a particle crossing the boundary can hit a
plane edge and register a spurious overlap.  Keep boundary planes non-free
(`FixVeloc` without `FixFree`): `FixFree` changes `IsFree` semantics in
bounding boxes, wrapping, and friction selection.  See
`tdem/test_periodic.cpp` for the reference setup.

### 2c. Save/Load state

`Domain::Save`/`Load` persist both Verlet references, `xb` and `wb`, together
with particle positions, velocities, flags, shapes, interactions and friction
history.

`FreePar`, `NoFreePar` and `Particle::Index` are derived state.  `Solve()`
clears and rebuilds them after `Initialize()`, so `Load()`-then-`Solve()` needs
no rebuild in `Load`.  Code that inspects or extends a loaded domain **before**
the first `Solve()` must rebuild that partition itself.

`Load` sets `Initialized=true` and deliberately preserves the saved `xb` and `wb`, so it
must not reinitialize loaded particles.  Particles added **after** a load are
therefore not initialized by the next `Solve()`; call `InitializeVelocity(dt)`
after adding or repositioning one.

Regression: `tdem/test_saveload` compares all constraint flags, `FixFree`, and
(after the continuation `Solve`) `FreePar`/`NoFreePar` and
`Particle::Index`, in addition to positions, velocities, both `xb` and `wb`,
and orientations.  The 300-step round trip is bit-identical.

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
- `AdaptiveNthread` was implemented for DEM, FLBM and LBMDEM and then rolled
  back on request; no field or sweep remains. The startup controller did not
  produce a reliable win: on the reduced DEM test it selected 32, but the LBM
  cases were flat and the candidate windows were sequential and noisy. Leave
  `Nthread` fixed and benchmark the actual workload.
- `tclbm01` (2400x300 D2Q9, 40k steps, no file output) is flat for
  `Nthread=32..256`: 45.1-45.7 s solve time. Three runs each at 32 and 256 gave
  means 45.11 vs 45.72 s, a 1.3% difference with Welch p=0.11 and a 95% CI of
  [-1.49,+0.26] s, so the apparent 32-thread gain is noise. 320-512 were
  46.4-46.7 s, about 2% slower. Do not quote a tuning gain for this
  memory-bound case.
- An experimental `AdaptiveAlpha` controller was implemented and benchmarked,
  then rolled back on request. Do not assume `Domain::AdaptiveAlpha` exists.
  Any future attempt must change `Alpha` only at a rebuild, rebuild immediately,
  and respect the linked-cell `Beta` cap.

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

Channel-flow obstacle cases add one more: with solid rows 0 and `ny-1` the
walls are at `y=0.5` and `y=ny-1.5`, so the channel width is `L=ny-2` and the
symmetry mirror is `j <-> ny-1-j`.  An off-by-one profile (`L=ny-1`) makes
an otherwise symmetric problem asymmetric and is easy to miss; a no-obstacle
control run (`max|uy| ~ 1e-4 u_max`, symmetric, centreline `uy ~ 1e-16`) is the
cheap test.  `/home/sergio/tflbm_triangle/tclbm_obstacle.cu` holds a validated
obstacle case (`--shape=rect|tri`, Re=40, 20% blockage: steady symmetric wakes,
`Lr/D=1.88` for a 60x40 box and `2.88` for the triangle, both grid- and
Mach-converged for the triangle; still steady at Re=100, i.e. no shedding).
A blunt rectangle gives a *shorter* wake than the apex-upstream triangle.

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

### CUDA MPM (GPU path)

- The GPU path was optimized after the CPU work; see
  [`CUDA_MPM_OPTIMIZATION_CHANGES.md`](CUDA_MPM_OPTIMIZATION_CHANGES.md).
- Hot node fields (`Mass`, `Mn`, `Fn`, `Vn`, `valid`) are separate SoA device
  arrays. `NodeCU` remains only for the cold `Vnf`/fixed flags; do not move the
  hot fields back into it without re-measuring.
- `cudaParticleToNode` builds the touched-node list itself (`atomicExch` +
  `atomicAdd`); `cudaResetNode` resets the previous list first; the device
  counter is ping-ponged. There is no per-step `thrust::fill`, `copy_if`,
  `distance` or count `cudaMemcpy`.
- Reset/solve are launched over `Nnodes` blocks and exit early on the device
  count. A final `cudaDeviceSynchronize()` at the end of `Solve` is required;
  without it `Solve` reports only the queued-work time.
- `ParticleCU` holds 64 nodes per particle (125 before). Cubic B-spline needs
  64; quadratic/GIMP at the usual particle size uses 27. A GIMP particle
  larger than roughly 1.5 cells needs 125 and must be handled explicitly.
  CPDI adds `Corners[4]`, `Gradvec[4]` and `Vt`; the structure is 3712 B
  (3584 B for the grid path).
- `MPM::Domain::Solve` has a CUDA CPDI branch for `BodyMesh!=NULL`:
  `cudaResetNode -> cudaCornerToNode -> cudaParticleToNodeCPDI ->
  cudaSolveNode -> cudaMoveCorners -> cudaNodeToParticleCPDI`.  It uploads
  `CornerCU` plus the CPDI fields of `ParticleCU`; `DnLoadDevice` downloads
  particles and corners and calls `UpdateMesh()` at output.  This is separate
  from LBMMPM, whose `MPMDOM.OneStepCPI` remains CPU.
- Reduced four-frame pre-CPDI `Nthread=256` results:
  `tmpm_cu_01/03/04/05` improved by 4/8/9/10%; `tmpm_cu_02` regressed about
  1% (short steps, frequent re-upload). The GPU uses atomics, so do not claim
  bitwise equality; observed final-frame differences were <=2e-5 absolute in
  `Pressure`/`Sigma` (about 1e-8 relative), and `tmpm_cu_02` happened to be
  bitwise identical.  `tmpm_cu_01` is now the CPDI beam and has its own
  reduced-window comparison above.
- Full `tmpm_cu_04` periodic square duct (`Tf = 200 s`, 216,000 points, 50
  frames, `Nthread=256`) took 22 min 57 s solve / 1378 s wall.  `Vx_max` is
  `7.68731e-4 m/s` against analytical `7.36712e-4 m/s` (+4.35%).  With
  `EOS_LINEAR` and uniform initial density, `Pressure` is exactly 0 in every
  frame; do not try to make a pressure colormap for this case.  Visualize
  `Vx = Velocity[0]` on a cross-section instead.
- Full `tmpm_cu_02` granular column (`Tf = 1.5`, 50 frames, `Nthread=256`)
  took 30 min 37 s Solve / 1837 s wall. The deposit shape is stable from about
  0.6 s (`x,z in [-0.126,0.314]`, `y_max=0.075 m`), but at 1.47 s the basal
  layer still moves: `vmax=1.009 m/s`, `vrms=0.104 m/s`, median `0.017 m/s`.
  It is not a static equilibrium; fixed ground plus `Gn=0.001` leaves a
  persistent sliding layer. Do not claim equilibrium from the 1.5 s state.
- Tried and rejected: packed `NodeHot`, precomputed `Fn/Mass` array, scalar
  shape cache, split gather/update kernel, register caps, block-size changes.

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

After changing `lib/mpm/*`:

```bash
cd tmpm
make -j8 tmpm01 tmpm02 tmpm03 tmpm04 tmpm05
make -j8 tmpm_cu_01 tmpm_cu_02 tmpm_cu_03 tmpm_cu_04 tmpm_cu_05
cd ../tlbmmpm && make -j8 tlbmmpm01 tlbmmpm02   # also compiles the MPM headers
```

The MPM cases are large: `tmpm02`/`tmpm05` are about 26,000 particles and
500,000-580,000 steps, `tmpm03` is about 135,000 particles and 259,000 steps,
and `tmpm04` is 216,000 particles and 42,400 steps. Smoke-test a copy with a
reduced grid and `Tf`; keep `tmpm04` at `tf = 200 s`, because its steady state
is set by viscous diffusion and reducing time hides errors. Reduce only
`particles_per_diameter` there: at 20 ppd (8,000 particles, 14,200 steps,
22.9 s CPU wall) `Vx_max` converged to `8.57e-4` vs analytic `7.37e-4`
(`Rerror = 16.3%`, five cells across the radius). The checked-in 60-ppd case
(216,000 particles, 96 threads) completed the full `tf = 200 s` in 1932 s wall
before the CPU optimization and 1676 s after it, and converged to
`7.68731e-4` vs `7.36712e-4` (`Rerror = 4.35%`).

MPM CPU transfer facts (see `MPM_OPTIMIZATION_CHANGES.md`):

- `Node` is `alignas(64)` with `Mass/Mn/Fn/valid/lck` in the first cache line;
  preserve that layout when adding node fields or the scatter loses its main
  gain.
- `VNodes` is built by first touch during `ParticleToNode`; `Node::valid` is
  also the "already in VNodes" flag.  A later full-grid scan would be slower.
- The shape cache (`FillShapeCache3D`) must keep the same multiplication order
  as `GIMP3D`/`BSpline*3D`; it is bitwise-identical on purpose.
- `Node::Reset` intentionally does not clear `Vn`/`Sn`.  Do not read them for
  invalid nodes.
- Exactness: build the same reduced driver against a copy of pre-change
  `lib/mpm` headers and compare HDF5 datasets at `Nproc=1`; expected result is
  bitwise equality.  At 96 threads expect <=1e-13 relative in
  `PVelocity`/`Sigma` on long runs.
- Measured 96-thread speedups: pre-CPDI `tmpm01` full 92.42 -> 33.77 s
  (2.74x), 20-ppd `tmpm03` 2.03 -> 0.89 s (2.29x), 50-ppd `tmpm03` 0.79 ->
  0.40 s (1.98x), full `tmpm04` 1932 -> 1676 s (1.15x).  `tmpm01` is now a
  CPDI tetrahedral beam, so those pre-CPDI timings are historical.  The
  checked-in refined CPDI case is 44,032 particles, 238,328 nodes and 9,151
  corners.
- Checked-in `tmpm01`/`tmpm_cu_01` now support both paths:
  `./tmpm01 [Nproc] [cpdi|mpm]` and the same for `tmpm_cu_01`; default is
  non-CPDI.  CPDI uses `AddRectangularBeamMesh` with `ndiv=125` and
  `ResizeDomainMesh(...,1.0)` (beta=1), and loads only the boundary corners in
  the tip slice (selected by `x`; the last cross-section is an external face,
  so no `BFac` mask is needed); `ForcePar` is empty.  Non-CPDI
  keeps `ndiv=50`, uses `AddRectangularBeam` + `ResizeDomain(...,Dx)`, and
  loads only `ForcePar`; `ForceCor` is empty.  Both
  report with `2*Fb`; current `Fb=5e-2`.  The CPU/GPU parity smoke at
  `Fb=0.1` (10 frames) matched except one `PVelocity` value (`4.4e-16`) at
  frame 7; solve times 33.88 s CPU vs 5.44 s GPU (`Nproc=1`).
- Halving the tip force (`Fb 0.10 -> 0.05`, both sources) halves the
  deflection but leaves the analytical comparison unchanged: full `Tf=5` GPU
  run tip `1.14806452e-3 m` vs analytical `9.75102376e-4 m`, ratio `1.1774`
  (the `Fb=0.1` ratio was `1.1780`).  The `y/Fb` profiles collapse.  Do not
  use force scaling as an accuracy fix.  Baseline and half-force results are
  `tmpm/tmpm01_Comparison_F01.res` and `tmpm/tmpm01_Comparison_F005.res`.
- CPDI body-mesh refinement with the current corner-only load (`tmpm_cu_01`,
  GPU, `Nproc=1`, `cpdi`, full `Tf=5`, `Fb=0.05`, `beta=1`); refining the body
  mesh also refines the background grid (`Dx=dmax`):
  | `ndiv` | particles | grid nodes | `Dx` (m) | tip (m) | tip/analytic | MAE (m) |
  |---:|---:|---:|---:|---:|---:|---:|
  | 50 | 3,315 | 17,576 | 7.551e-3 | 8.0449e-4 | 0.8250 | 6.842e-5 |
  | 75 | 10,045 | 50,653 | 5.266e-3 | 8.7912e-4 | 0.8878 | 4.554e-5 |
  | 100 | 23,180 | 125,000 | 3.962e-3 | 9.4365e-4 | 0.9486 | 2.065e-5 |
  | 125 | 44,032 | 238,328 | 3.204e-3 | 9.7016e-4 | 0.9694 | 1.217e-5 |
  Refinement moves the tip toward the analytical value and cuts MAE by about
  5.6x from `ndiv=50` to 125.  The checked-in CPDI case is `ndiv=125`;
  non-CPDI remains `ndiv=50`.  2026-11-06 CUDA hot-path optimization: full `ndiv=125` `Solve` fell from
  17m23.33s to 4m02.85s (4.30x); `ndiv=50` fell from 47.53s to 26.74s
  (1.78x).  The final deflection is unchanged (`9.70156508e-4 m`, ratio
  `0.9693702`, MAE `1.2168888e-5 m`; optimized `.res` differs by at most
  `1e-16 m`, `ndiv=50` bitwise).  Changes: `cudaComputeAn` caches
  `An=Fn/Mass`; `ParticleCU::GradvecN[4]` caches `Gradvec/(6*Vt)`; particle
  `v/L/x` and corner `x` are local; CPDI particle scatter no longer rebuilds
  the touched-node list (corner scatter already covers every node); corner
  scatter skips zero-force and zero-shape atomics.  HDF5 fields shift by
  round-off (`~2e-8` relative).  The non-CPDI path has its own pass: separable
  1-D GIMP/B-spline factors (9 evaluations instead of 81 for 3x3x3),
  `cudaComputeAn` for the gather, and register-resident `v/x/L`; full
  `Solve` fell from 1m54.82s to 21.03s (5.46x).  Its final `.res` matches the previous
  GPU baseline to `1e-13 m`; HDF5 shifts by about
  `7e-8` relative in `PVelocity`.  See
  `CUDA_MPM_OPTIMIZATION_CHANGES.md`.  An earlier `Tf=2` test with the old combined
  particle+corner load moved the ratio the other way; that was the load path,
  not evidence that CPDI refinement is harmful.
- `beta` controls the background spacing (`Dx=beta*dmax`): below 1 is finer,
  above 1 is coarser.  At the current `ndiv=125`, `Fb=0.05`, full `Tf=5`:
  | `beta` | grid nodes | `Dx` (m) | tip/analytic | MAE (m) | GPU `Solve` |
  |---:|---:|---:|---:|---:|---:|
  | 1.00 | 238,328 | 3.204e-3 | 0.9694 | 1.217e-5 | 17m23s |
  | 1.25 | 117,649 | 4.005e-3 | 1.1652 | 8.878e-5 | 17m15s |
  | 1.50 | 68,921 | 4.806e-3 | 1.1392 | 7.653e-5 | 17m17s |
  `beta=1.0` is clearly best here; coarsening overshoots the analytical tip.
  The earlier positive impression of `beta=1.5` came from the old combined
  particle+corner load and does not transfer to the boundary-corner path.
  Those beta timings are pre-optimization; at `ndiv=125` the optimized
  `Solve` is about 4m03s, so use beta for accuracy, not as a time knob.
- Final frame at `Fb=0.05`, full `Tf=5`:
  | path / load | tip (m) | tip/analytic | MAE (m) |
  |---|---:|---:|---:|
  | non-CPDI, particles (5,000 particles, 125,000 nodes) | 9.34560567e-4 | 0.9623 | 1.5385e-5 |
  | CPDI, boundary corners, `ndiv=50` (44 corners, 3,315 particles, 17,576 nodes) | 8.04488027e-4 | 0.8250 | 6.8416e-5 |
  | CPDI, boundary corners, `ndiv=75` (99 corners, 10,045 particles, 50,653 nodes) | 8.79124663e-4 | 0.8878 | 4.5539e-5 |
  | CPDI, boundary corners, `ndiv=100` (144 corners, 23,180 particles, 125,000 nodes) | 9.43645930e-4 | 0.9486 | 2.0650e-5 |
  | CPDI, boundary corners, `ndiv=125` (229 corners, 44,032 particles, 238,328 nodes) | 9.70156508e-4 | 0.9694 | 1.2169e-5 |
  | CPDI, old particles + all tip corners, `ndiv=50` | 1.14806452e-3 | 1.1774 | 6.1911e-5 |
  The checked-in CPDI boundary-corner path uses `ndiv=125` and has no particle
  load.  The optimized non-CPDI final GPU `.res` matches the
  previous GPU baseline to `1e-13 m`; CPDI matches to `1e-16`.  CPU/GPU HDF5 payloads are
  not promised bitwise after the separable-shape/`An` caches; the observed
  `PVelocity` shift is at the `1e-8` relative level.  Do not use force-path
  changes alone as an analytical accuracy fix.
- Do not remove `pValidNumberPrev` as a "redundant" counter.  A
  single-`pValidNumber` reset is deterministic but changes the non-CPDI
  `tmpm_cu_01` trajectory (first difference at row 80, final `.res` max
  `8.74e-3 m`); CPDI `ndiv=50` stays exact.  Sync memset, extra device sync,
  and a dummy second allocation all reproduce the difference, and the change
  gave no speedup.  Keep the two-counter ping-pong; validate any redesign on
  the full non-CPDI case.
- `tmpm04` is dense (5.5 particles per node): its scatter is cache-line
  ownership bound.  On a short run 32 threads beat 96 (60.6 vs 71.3 s).
  Eight per-group partial grids gave 60.8 s but regressed sparse cases, so
  they are not in the library.

LBMMPM CUDA facts (`lib/lbmmpm/lbmmpm.cuh`):

- `ForceCorrection` has been removed from the code: CPU `ImprintLattice` and
  CUDA `cudaImprintFacets` now always use the full `-2 rho W Cs C` facet term.
  The former `lbmmpm_aux`/`bMPMAux` upload and its 0.6 default are gone.
  Entries below that mention `FC=0.6` are historical measurements made before
  this removal; they are not current defaults.
- Liu, Ye, Zhang, and Zheng, Int. J. Numer. Methods Eng. 121 (2020)
  5546-5567, doi:10.1002/nme.6501, is the direct source for this term.
  Section 4.2.2, Eq. (37), p. 5557 gives the deformable-solid force as the
  Galilean momentum exchange minus the full reference-pressure term
  `2 rho_f0 omega_h e_h`; Eq. (36) omits that term and is only for
  fluid-rigid interaction.  The paper uses `rho_f0`, the initial fluid
  density, not a fitted 0.6 factor, and its equation uses the actual
  solid-side population, not a copied fluid population.  Thus the code's
  full `-2 rho W Cs C` is paper-consistent up to the code's local `rho`
  (paper says `rho_f0`) and the extra legacy solid-neighbour copy, which is
  not part of Eq. (37).
- Historical full first test (`tlbmmpm_cu_01`, elastic slab/beam, `Re=9`,
  200x200x200, `ndiv=40`, `Tf=1.5`, then-default `ForceCorrection=0.6`,
  `.res`-only, `Nproc=8`):
  solve/wall 643 s.  The 99 reported frames are `t=0...1.47`; after the
  startup transient (`t>=0.75`), mean `Re=8.994` vs target 9, mean tip
  deflection `dx=7.05 mm (dx/Lz=0.0705)`, `dz=-3.42 mm (dz/Lz=-0.0342)`.
  The x deflection still oscillates, about `5.6--8.5 mm` (std 0.94 mm), and
  `dz` only about `-3.53...-3.32 mm`; the last frame is
  `dx=6.34 mm`, `dz=-3.37 mm`.  Peak startup `dx=25.0 mm` at `t=0.21 s`.
  Results: `/tmp/t01_re9/run/tlbmmpm_cu_01time.res`,
  `t01_re9_deflection.png`, `t01_re9_summary.csv`.
- Legacy LBM-MPM history is older than the cube test:
  `rev 2169` (`d629b25b5f52`, 2020-09-17, `MPM: USP integrator implemented`)
  first adds `lib/lbmmpm/Domain.h`, `tlbmmpm/tlbmmpm01.cpp`, and its
  CMakeLists; the commit message does not mention LBM-MPM.  `rev 2172`
  (`4ab7ceb0e0fe`, 2020-09-21, `LBMMPM: started`) is the first dedicated
  LBM-MPM development commit.  `rev 2189` (`440e67b768d4`, 2021-08-04,
  `LBMMPM: First upload`) is the first upload that adds the cube
  `tlbmmpm02.cpp` and moving-cube `tlbmmpm03.cpp`; `tlbmmpm02.cpp` is
  unchanged through the 2024 tip.  So the analyzed legacy cube is the first
  cube test, but not the first LBM-MPM commit.
- Legacy Mercurial first upload (`rev 2189`, `440e67b768d4`, 2021-08-04,
  `LBMMPM: First upload`) has the same `Re=9` t01 setup; `tlbmmpm01.cpp` only
  allocates `dat.Vel` as `new double[ny]` instead of `[nz]` (both 200 here).
  A clean rev-2189 checkout does not compile: `lib/mpm/Domain.h` includes
  `mechsys/mpm/Corner.h`, which was not added until `rev 2200`
  (`723f35c5accc`, 2021-11-24, `LBMMPM: bugs addressed`).  The old coupling
  also uses the full `-2 rho W Cs C` facet term (equivalent to current
  `ForceCorrection=1`), divides refill by `naem` and aborts on negatives, and
  keeps the solid-neighbour copy `F[nxs][nys][nzs][k] = F[nx][ny][nz][k]`
  active; current code uses `fabs(F)/naem`, clamps negative `newF`, and has
  that copy commented out.  The old test is CPU-only, while current
  `tlbmmpm_cu_01` is CUDA.  Rebuilding rev 2189 with a rev-2200 `Corner.h`
  stalls in old mesh/tetgen initialization, so the exact old trajectory is not
  reproducible from the clean repository.  Running current t01 with
  `ForceCorrection=1.0` isolates the force-term change: `t>=0.75` mean
  `dx=4.98 mm`, `dz=-0.12 mm` versus `7.05 mm`, `-3.42 mm` at the default
  0.6 (same `Nproc=8`, wall 660 s).  Results:
  `/tmp/t01_fc1/run/tlbmmpm_cu_01time.res`,
  `/tmp/t01_legacy_compare/fc06_vs_fc1.png`.
- The legacy CPU cube is `tlbmmpm/tlbmmpm02.cpp`; that file is unchanged from
  rev 2189 through the legacy tip.  Running it with 64 OpenMP threads (tip
  legacy headers plus the current gmsh `mesh/unstructured.h` overlay, because
  the clean rev-2189 mesh/tetgen build is broken) gives:
  `b=0.25`: final `Re=38.666`, `CD=4.87784` vs Eq. (38) `1.25021`
  (**+290.2%**), last-10 `CD=4.91114` (+290.4%), wall 32m45s (CPU time
  34h52m, 324 MB RSS); `b=0.84`: final `Re=100.838`, `CD=1.96453` vs
  Eq. (38) `0.580694` (**+238.3%**), last-10 `CD=1.97145` (+238.7%), wall
  32m19s (CPU time 34h28m).  This is the same failure mode as the current
  CUDA full-force + solid-copy runs (`CD=5.174`/+315.7% and
  `CD=1.883`/+229.1%), so the legacy cube's poor accuracy is not a GPU
  artifact: its active solid-neighbour copy and lack of the later empirical
  0.6 facet-force scaling overpredict drag by roughly 3.4--4.0x.  Results:
  `/tmp/legacy_cube/run_b025/tlbmmpm02.res`,
  `/tmp/legacy_cube/run_b084/tlbmmpm02.res`,
  `/tmp/legacy_cube/legacy_vs_current_CD.png`,
  `/tmp/legacy_cube/legacy_vs_current_summary.csv`.
- Reactivating the legacy solid-neighbour copy in the CUDA facet kernel
  (`F[k*Nc+ns]=Fk; Fns=Fk` immediately before the facet force) with
  full facet force (`ForceCorrection=1.0` at the time; now the only behavior)
  changes `t>=0.75` t01 deflection to mean
  `dx=35.58 mm (dx/Lz=0.3558)`, `dz=-6.62 mm`, with range
  `dx=33.88--38.97 mm` and a startup peak `55.52 mm` at `t=0.405 s`; the run
  does not settle by `t=1.47 s` (last `34.28 mm`, `-6.02 mm`, wall 653 s).
  For comparison: `FC=1` without the copy is `4.98 mm/-0.12 mm`, and default
  `FC=0.6` without the copy is `7.05 mm/-3.42 mm`.  The CUDA write is a data
  race because several facet threads can target the same solid `ns`; unlike
  the CPU loop, GPU thread order is not deterministic, so treat the trajectory
  as order-dependent.  Results:
  `/tmp/t01_solidcopy_fc1/run/tlbmmpm_cu_01time.res`,
  `/tmp/t01_solidcopy_fc1/solidcopy_fc1_deflection.png`,
  `/tmp/t01_solidcopy_fc1/solidcopy_fc1_profile.png`,
  `/tmp/t01_solidcopy_fc1/solidcopy_fc1_summary.csv`.
- The t01 bending-line output is `tlbmmpm_cu_01position.res` (not a
  `middle*.res`): it has two normalized columns, `X-position` and
  `Z-position`, and is written only at `idx_out==99` from `middlePar`
  particles (initial x near 0.15).  The time history is
  `tlbmmpm_cu_01time.res`.  The copy is phase-correct for post-streaming
  data: `ImprintLattice` runs at the top of the step after the previous
  `StreamMPM`, so `F[nx]` is post-streaming.  However, the current CUDA copy
  is inside the facet loop *after* `newF`; it therefore does not enter the
  current bounce-back, and later facets/steps can see an order-dependent
  value.  A consistent post-streaming ghost needs one pass after
  `StreamMPM` into a separate `Fghost`, used by both Eq. (30) and Eq. (37).
- The t01 setup is the 10 mm flap of the MDPI experiments: E=1.23 MPa,
  nu=0.3, rho_s=1030, rho_f=1220, mu=1 Pa s, Re=3--12; at Re=9, 0 deg the
  experimental tip deflection is about 39 mm (Appl. Sci. 10(18):6177, 2020,
  doi:10.3390/app10186177).  The copy gives 35.58 mm; no-copy gives about
  5 mm, so the copy is acting as a force amplifier.  A single global
  `ForceCorrection` cannot satisfy both this flap and the cube (the copy
  gives +290--316% cube error).  The likely fix is to replace the
  same-direction copy by a consistent ghost/reflected solid-side population
  used in both Eq. (30) and Eq. (37), or by a one-sided momentum flux/stress
  integration; test `Fns=F[nx][ko]`, `Fns=newF`, and non-equilibrium
  extrapolation before any empirical fit.  Also check the tank: t01 uses a
  0.3 m LBM cube, whereas the MDPI tank is 150w x 20w x 12.5w = 3.0 x 0.4
  x 0.25 m.
- The MDPI simulation is fully **3D**, not 2D.  Its abstract says
  "Two-way coupled three-dimensional transient simulations"; Section 2.2.1
  gives the domain as (150w, 20w, 12.5w) = 3.0 x 0.4 x 0.25 m; Section 2.2.2
  says "For a 3D dynamic mesh domain in Fluent"; and the drag discussion
  refers to "the current 3D case".  The published bending lines are 2D
  projections/extracted centerlines through the 3D flap, not a 2D model.
  The 20 mm `w` is the flap width and the Re length, not evidence of 2D.
  t01 is also 3D (200^3); the VisIt `y=0.15` slice is only a visualization.
- The MDPI paper defines `Re = v w / nu`, where `v` is the constant
  uniform inlet velocity, `w` is the flap width (20 mm), and `nu=mu/rho`.
  It does not use the maximum of a parabolic profile, a hydraulic diameter,
  or the flap length.  At Re=9, `v=9*(1/1220)/0.02=0.36885 m/s`; at Re=12,
  about 0.492 m/s (the paper quotes 0.5).  Current t01 uses this uniform
  inlet at Re=9 (`v=Re*nu/Ly`, `dat.umax=v`), and `Report` prints
  `Re=umax*Ly/nu=9`.  Historical parabolic t01 used
  `umax=1.5*Re*nu/Ly=0.5535 m/s` (mean 0.3689, max-based Re 13.5) and
  varied along z (`dVel[j]`); its no-copy deflection was 4.99 mm versus
  8.78 mm for the current uniform case.
- Historical max-based Re test: `Ly=0.02`, `vmax=0.5535 m/s`, and `Report`
  was temporarily changed to `Re=vmax*Ly/nu=13.5`.  The flow was unchanged
  from the mean-Re=9 case.  Full no-copy run (1.5 s, `Nproc=1`, wall
  719.66 s): `t>=0.75` mean `dx=4.9886 mm` (`dx/Lz=0.049886`), range
  4.234--5.939 mm, `dz=-0.1222 mm`; last-10 mean `4.8307 mm`.  Results:
  `/home/sergio/DATA/t01_remax_w/tlbmmpm_cu_01time.res`,
  `/home/sergio/DATA/t01_remax_w/tlbmmpm_cu_01position.res`.  This diagnostic
  was reverted when the uniform inlet was adopted.
- Uniform-inlet t01 at the MDPI Re=9 (`v=Re*nu/Ly=0.36885 m/s`, all z
  values equal, no copy, 1.5 s, `Nproc=1`, wall 721.17 s): `t>=0.75` mean
  `dx=8.7811 mm` (`dx/Lz=0.087811`, std 0.513 mm, range 8.010--9.704 mm),
  `dz=-0.3666 mm`; last-10 mean `dx=8.6930 mm`; final `dx=8.2306 mm`;
  startup peak 15.026 mm at t=0.195 s.  This is 76% more deflection than the
  parabolic mean-Re=9 run (4.9886 mm), because the parabola was zero near the
  clamp and rose along z, whereas the uniform inlet loads the whole flap.
  Current `tlbmmpm_cu_01.cu` uses this uniform inlet and `Report` writes
  `Re=umax*Ly/nu=9`.  Results:
  `/home/sergio/DATA/t01_uniform_re9/tlbmmpm_cu_01time.res`,
  `/home/sergio/DATA/t01_uniform_re9/tlbmmpm_cu_01position.res`.
- Refined t01 slab mesh: `ndiv=64` (nominal `Lz/64=1.5625 mm`; LBM cell
  `dx=1.5 mm`).  The generated mesh has 1026 vertices, 3586 tets and 1624
  surface triangles; surface edge lengths are min 1.9372, mean 3.0673,
  max 4.4621 mm, so `min/dx=1.291` and every triangle is still larger than
  the LBM cell.  Full uniform-Re=9 no-copy run (1.5 s, `Nproc=1`, wall
  1170.17 s; `dt=4.52154e-5`, C=33.17): `t>=0.75` mean `dx=9.3139 mm`
  (`dx/Lz=0.093139`, std 0.1678 mm, range 9.047--9.606 mm),
  `dz=-0.4182 mm`; last-10 mean `9.3079 mm`; final `9.1766 mm`; startup
  peak 14.312 mm at t=0.180 s.  Compare the old `ndiv=40` result
  (`dx=8.7811 mm`).  Results:
  `/home/sergio/DATA/t01_ndiv64_uniform/tlbmmpm_cu_01time.res`,
  `/home/sergio/DATA/t01_ndiv64_uniform/tlbmmpm_cu_01position.res`.
  VisIt animation of this refined run: 100 states, `speed` on `turbo`
  (0--0.65 m/s), `Slice` at `y=0.15 m`, deformed `MPM_mesh` overlay; 900x900,
  80 ms/frame, wall 135.12 s.  GIF
  `/home/sergio/DATA/t01_ndiv64_uniform/t01_ndiv64_fluid.gif` (6.4 MB),
  poster `t01_ndiv64_fluid_poster.png`, scripts `render_visit_fluid.py`
  and `make_fluid_gif.py`.
  Mesh size sweep: `ndiv=60` min 2.249 mm (1.50 dx), `64` 1.937 (1.29 dx),
  `72` 1.936, `80` 1.646 (1.10 dx), `85` 1.521 (1.01 dx); do not refine
  beyond about `ndiv=80` without rechecking the generated facets.
- Domain-size sweep for the refined (`ndiv=64`) uniform-Re=9 no-copy case.
  `tlbmmpm_cu_01` now takes an optional `argv[2]` domain size `Ldom`; the
  grid uses `n=round(Ldom/dx)` and the body/report centre is explicit `xc`.
  Actual sizes: 0.2 m -> `n=133`, 0.3 m -> `n=200`, 0.5 m -> `n=334`
  (`0.501 m`).  Results for `t>=0.75 s`, `Nproc=1`:
  `L=0.2`: mean `dx=8.8393 mm`, std 0.0477 mm, range 8.763--8.922 mm,
  `dz=-0.4206 mm`, last-10 8.8324 mm; wall 488.53 s, 8.9 GB output.
  `L=0.3`: mean `dx=9.3139 mm`, std 0.1678 mm, range 9.047--9.606 mm,
  `dz=-0.4182 mm`, last-10 9.3079 mm; wall 1170.17 s.
  `L=0.5`: mean `dx=7.8670 mm`, std 0.7678 mm, range 6.736--9.013 mm,
  `dz=-0.3413 mm`, last-10 8.0631 mm; wall 4740.71 s, 138 GB output.
  Increasing to 0.5 m therefore does **not** improve the deflection (15.5%
  lower than 0.3 m) and makes it less steady; decreasing to 0.2 m is 5.1%
  lower than 0.3 m but is much steadier and 2.4x faster.  0.3 m remains the
  largest-deflection and best compromise.  The checked-in default remains
  `Ldom=0.5`; run `./tlbmmpm_cu_01 1 0.3` (or 0.2) to select another size.
  Results:
  `/home/sergio/DATA/t01_L020_uniform/tlbmmpm_cu_01time.res`,
  `/home/sergio/DATA/t01_L050_uniform/tlbmmpm_cu_01time.res`.
  Velocity sweep at `Ldom=0.3` (same mesh, `ndiv=64`, no copy): mean
  `dx` for `t>=0.75 s` is 17.93 mm at v=0.70, 19.23 at 0.80, 30.54 at
  1.20, 36.78 at 1.60, 37.49 at 1.70, 40.58 at 1.815, 41.48 at 1.83,
  45.67 at 2.00.  The transition near 1.81--1.83 is sensitive and
  oscillatory; do not claim a single precise critical velocity.  Current
  default is `Ldom=0.3`, `v=1.815 m/s`, `Re=44.27`; optional
  `argv[3]` overrides v.  Full run (1.5 s, `Nproc=1`, wall 1179.24 s):
  `t>=0.75` mean `dx=40.5825 mm`, std 1.2567, range 38.57--44.56,
  `dz=-8.8224 mm`, last-10 40.2295, final 39.7146, startup peak
  56.5426 mm.  Results:
  `/home/sergio/DATA/t01_L030_v1815_uniform/tlbmmpm_cu_01time.res`.
  Final central centreline (`x~0.15`, `y~0.15`, 60 selected material
  points) has tip `dx=39.61 mm`; plot
  `/home/sergio/DATA/t01_L030_v1815_uniform/t01_v1815_middle_line.png`.
  VisIt animation of this run: 100 states, `turbo` speed range 0--2.0
  m/s on the `y=0.15 m` slice, deformed `MPM_mesh`, 900x900, 80 ms/frame,
  wall 129.33 s.  GIF
  `/home/sergio/DATA/t01_L030_v1815_uniform/t01_v1815_fluid.gif` (6.2 MB),
  poster `t01_v1815_fluid_poster.png`, scripts `render_visit_fluid.py`
  and `make_fluid_gif.py`.  Use the centreline selection from the
  original/final `Position` datasets (fixed material-point ids), not the
  full `middlePar` slice, for a clean line plot.
  The final `tlbmmpm_cu_01position.res` has 731 rows (698 with z>=0), with
  columns `(x-xc)/Lz` and `z/Lz`; its binned free-length tip is
  `X=0.38686` (38.7 mm).  A direct plot of those columns is
  `/home/sergio/DATA/t01_L030_v1815_uniform/t01_v1815_position_file.png`
  (script `plot_position_res.py`); it shows all raw rows as points plus a
  5e-4 height-binned line.
  The same rows plotted as `Delta z` versus `Delta x` (using the matching
  frame-0 particle ids for the initial coordinates) are
  `/home/sergio/DATA/t01_L030_v1815_uniform/t01_v1815_position_delta.png`
  (script `plot_position_delta.py`): 698 free rows, binned tip
  `Delta x=38.84 mm`, `Delta z=-8.23 mm`.
  The requested `z` versus `Delta x` version is
  `/home/sergio/DATA/t01_L030_v1815_uniform/t01_v1815_position_z_dx.png`
  (script `plot_position_z_dx.py`); binned tip `Delta x=38.84 mm`,
  `z=91.04 mm`.

  Half-thickness run (`Lx=0.005 m`, same `v=1.815`, `Re=44.27`,
  `Ldom=0.3 m`, `ndiv=64`, `Nproc=1`, no video) completed in 1081.20 s:
  `/home/sergio/DATA/t01_L030_v1815_Lx005_padLz/`.  It has 1707 material
  points and 643 corners.  For `t>=0.75 s` the mean tip is
  `dx=83.335 mm` (std 0.610), final `83.110 mm`, `dz=-51.81 mm`; mean
  `CD=0.2993`.  The half slab is much more flexible and the reference
  transient peak is about 56.5 mm, so the MPM box must not shrink with
  `Lx`: the original `10*Lx` padding gave only 5 cm and crashed with
  `nhigh out of bounds`; use `PadX=Lz` (10 cm).  `tlbmmpm_cu_01` now
  accepts `argv[5]=Lx`, `argv[6]=ndiv`, `argv[7]=beta` (defaults
  0.01/64/1.5).  Attempts to resolve the 5 mm thickness (`ndiv=80,
  beta=0.75`, and `beta=0.5`) were numerically unstable immediately
  (huge corner forces), so the completed half run uses the same
  under-resolved thickness mesh as the reference and is qualitative for
  the deflection magnitude.
  Comparison plot:
  `/home/sergio/DATA/t01_L030_v1815_uniform/t01_v1815_halfLx_comparison.png`.
  The two middle planes are overlaid in
  `/home/sergio/DATA/t01_L030_v1815_uniform/t01_v1815_middle_planes_10_5mm.png`
  (script `plot_both_middle_planes.py`), plotting `z` against `x-xc` for
  both thicknesses.  Because the 5 mm case was run without HDF5, the plot
  uses `x-xc` consistently for both; the middle-plane selection is within
  1 mm of `xc`.  Binned tips are 10 mm: `dx=38.36 mm, z=91.15 mm`; 5 mm:
  `dx=82.54 mm, z=49.05 mm` (free rows 698 and 549).
  A 100-frame VisIt animation of the half-thickness run is
  `/home/sergio/DATA/t01_L030_v1815_Lx005_padLz_video/t01_Lx005_fluid.gif`
  (5.2 MB; poster `t01_Lx005_fluid_poster.png`; scripts
  `render_visit_fluid.py`, `make_fluid_gif.py`).  The video rerun (same
  parameters, `RenderVideo=1`) took 1137.84 s, wrote 30 GB, and its
  `t>=0.75` mean is 83.30 mm (final 82.94 mm, `CD=0.3001`), consistent
  with the no-video run.  VisIt rendering took 137.70 s; camera focus
  `(0.17,0.15,0.045)`, `parallelScale=0.13`, speed range 0--2.0 m/s.
  The raw `*.h5`/`*.xmf` frames for both the 10 mm reference and the
  half-thickness video run were deleted after producing the GIFs/plots.
  To re-render or post-process fields, rerun with `RenderVideo=1`; the
  retained products are the GIFs, posters, `position.res`, `time.res`, and
  scripts.
  A 100-frame 800x800 bending animation was rendered from the older
  `ndiv=40` MPM HDF5 mesh with VTK 9.5 offscreen + Pillow GIF (12.5 fps,
  dx colour range 0--16 mm, project `turbo` fluid colormap):
  mesh with VTK 9.5 offscreen + Pillow GIF (12.5 fps, dx colour range
  0--16 mm, project `turbo` fluid colormap):
  `/home/sergio/DATA/t01_uniform_re9/t01_uniform_bending.gif`;
  poster `/home/sergio/DATA/t01_uniform_re9/t01_uniform_bending_poster.png`;
  script `/home/sergio/DATA/t01_uniform_re9/render_animation.py`.
- VisIt fluid animation for the uniform t01 slab: 100 LBM/MPM states, a
  `Slice` at `y=0.15 m` (normal `+y`), LBM `speed=magnitude(Velocity_0)` on
  `turbo` with fixed 0--0.65 m/s, and the deformed `MPM_mesh` overlaid.  A
  `Slice` is much cheaper than `Clip` here (~3 s/frame, 128 s wall for 100
  900x900 frames).  The `.visit` playlists are
  `/home/sergio/DATA/t01_uniform_re9/lbm.visit` and `mpm.visit`; render script
  `render_visit_fluid.py`; GIF
  `/home/sergio/DATA/t01_uniform_re9/t01_uniform_fluid.gif` (6.7 MB) and
  poster `t01_uniform_fluid_poster.png`.  `quit()` in the VisIt script is
  needed to avoid a lingering engine after `SaveWindow`.
- Cube drag with the same CUDA solid-neighbour copy and the then-default
  `ForceCorrection=0.6` (full 100 frames, `Nproc=8`; this was before the
  member was removed, so it is historical, not the current full-force cube):
  `b=0.25` reaches `Re=38.884`, `CD=5.4941` vs Eq. (38) `1.2446` (+341.4%);
  no-copy prior is `CD=1.2803` (+1.0%).  `b=0.84` reaches `Re=102.684`,
  `CD=2.0688` vs Eq. (38) `0.5723` (+261.5%); no-copy prior is
  `CD=0.5858` (-0.5%).  The copy roughly triples-to-quadruples cube drag and
  destroys the cube validation.  It was kept active by user request, then
  removed again on 2026-11-08.  Fresh no-copy full-force stock run
  (`tlbmmpm_cu_02.cu` unchanged, `b=0.25`, 20k steps, `Nproc=1`, wall
  360.77 s): final `Re=37.9866`, `CD=0.912101` vs Eq. (38) `1.268068`
  (**-28.07%**), last-10 mean `CD=0.919013`/`Re=37.7113` (-27.95%); Mach is
  only about 0.076.  Removing the copy therefore changes the same stock
  point from `CD=5.174` (+315.7%) to `CD=0.912` (-28.1%), not to the
  historical FC=0.6 value 1.280.  Results:
  `/home/sergio/DATA/cube_nocopy_b025/tlbmmpm_cu_02.res`,
  `/home/sergio/DATA/cube_nocopy_b025/cube_nocopy_vs_copy_CD.png`.
  Historical copy results: `/tmp/cube_solidcopy/b025/run/tlbmmpm_cu_02.res`,
  `/tmp/cube_solidcopy/b084/run/fc06_b2.1e-05.res`,
  `/tmp/cube_solidcopy/cube_solidcopy_CD.png`,
  `/tmp/cube_solidcopy/cube_solidcopy_summary.csv`.
- Current code after removing `ForceCorrection` (full facet force) with the
  solid-neighbour copy still active, full 100 frames, `Nproc=8`:
  `b=0.25` reaches `Re=38.8846`, `CD=5.17416` vs Eq. (38) `1.24459`
  (+315.7%; last-10 +317.3%, wall 354.4 s); `b=0.84` reaches
  `Re=102.686`, `CD=1.88344` vs Eq. (38) `0.572319` (+229.1%; last-10
  +233.1%, wall 355.5 s).  These are slightly lower than the historical
  FC=0.6 + copy values but still 4.0x/3.2x the no-copy drag and far outside
  Eq. (38), so the solid copy remains the dominant error.  Results:
  `/tmp/cube_fc_removed/b025/tlbmmpm_cu_02.res`,
  `/tmp/cube_fc_removed/b084/full_b2.1e-05.res`,
  `/tmp/cube_fc_removed/cube_full_force_CD.png`,
  `/tmp/cube_fc_removed/cube_full_force_summary.csv`.
- Device `F` is direction-major; address it with `FLBM::FIDX(Ncells,Nneigh,
  il,k,ic)` from `lib/flbm/lbm.cuh`.  The boundary kernels in
  `tlbmmpm_cu_01.cu` originally had a local `t01_idx` copy of this formula; it
  was removed.  Pass the already-computed cell index `n`, not x/y/z separately.
- Detailed CUDA LBMMPM architecture (host/device data map, upload/download
  schedule, launch order, kernel responsibilities, extension checklist):
  `CUDA_LBMMPM_STRUCTURE.txt`.  Read it before adding a coupling kernel.
- Do not confuse `MPMDOM.Corners` with `MPM::Particle::Corners[4]`.
  `MPMDOM.Corners` is the body-mesh vertex array and is uploaded flattened as
  `pCornerPos/pCornerVel/pCornerH`; `pCellVerts`/`pFacVerts` are indices into
  that array.  In LBMMPM, `Particle::Corners[4]` is used only by the CPU
  `OneStepCPI` and is not uploaded.  The standalone MPM CUDA CPDI path does
  upload it (as `ParticleCU::Corners`).
- `Mat3::SolAlt` lives with the CUDA `Mat3` class in
  `lib/linalg/matvec.h`.  The LBMMPM kernels call that method; they must not
  reintroduce a local inverse/solver copy.

- D3Q15 only; `MPMDOM.OneStepCPI` stays on the CPU.  Host pair lists cover
  `(body cell, lattice cell)` and `(facet, lattice cell)`; `cudaMarkFacets`
  must run before `cudaImprintFacets`, because the force pass reads
  `Gamma[ns]` and a single mixed kernel races the marking writes.
- Build the CUDA examples with `--fmad=false` (the CMake targets do).  On facet
  edges the CPU `SolAlt` returns exact zeros; FMA contraction in the device
  inverse changes the marked surface cells (145 reduced / 176 full t02 cells).
  With the flag, `Gamma` was bitwise identical at the compared frames.
- 2026-11-06 kernel pass on the RTX A6000 (`Nthread=256`,
  `-maxrregcount=128`, `--fmad=false`): full t02 window (241x61x61, 1001
  steps) 4.993 -> 4.062 s (1.23x); full t01 window (200x200x200, 139 steps,
  `Tf=0.01`) 10.664 -> 8.938 s (1.19x).  The pass fuses Rho/Vel into
  `cudaCollideMPM` (no `cudaVelDen` launch), precomputes `1/Cs`/`1/rho`,
  defers limiter `Ftemp` writes, and uses branch wrap in `cudaStreamMPM`.
  A single-scan `alpha=1` limiter is faster but not equivalent: negative-noneq
  directions can cross the tolerance only after alpha is reduced.  Later-frame
  GPU differences are dominated by baseline run-to-run atomics (`~7.2e-4`
  Density_0, `~2.2e-3` Velocity_0); `frame 0000` remains byte-identical.
- CMAME 2021, Eq. (38): `Cd = 128*Re^-0.8/(1+4.5*beta_s^0.35)`;
  for a cube `beta_s=1`, so the built-in `/5.5` is exact, not a different
  fit.  Reference: P. Zhang, S. Sun, Y. Chen, W. Cui, S.A. Galindo-Torres,
  "Coupled material point Lattice Boltzmann method for modeling
  fluid-structure interactions with large deformations," Computer Methods
  in Applied Mechanics and Engineering 385 (2021) 114040,
  doi:10.1016/j.cma.2021.114040 (hereafter CMAME 2021).  The paper's cube
  validation uses the same 241x61x61 periodic box but `b=2.5e-5`; stock
  `tlbmmpm_cu_02` uses `b=0.25` and reaches `Re~38`, so its high-Re CD
  deviation is not a formula-denominator bug.
- GPU body-force sweep (full 100,000 steps each, same grid, `.res` only):
  `b=2.5e-5` -> `Re=0.00545`, `CD=4657` vs Eq. (38) 1505;
  `b=2.5e-3` -> `Re=0.545`, `CD=46.58` vs 37.82;
  `b=2.5e-2` -> `Re=5.217`, `CD=5.224` vs 6.208;
  `b=0.1` -> `Re=17.8`, `CD=1.796` vs 2.325;
  `b=0.25` -> `Re=37.99`, `CD=0.913` vs 1.268.  At matched Re the CD is
  force-dependent, e.g. Re=2: 13.22 / 15.00 / 19.19 for b=0.025/0.1/0.25;
  Re=5: 5.42 / 5.65 / 6.11.  These constant-force runs are not quasi-steady,
  so instantaneous CD-vs-Re points are not a clean Eq. (38) validation;
  lowering bforce alone does not produce the low-Re agreement.  A fixed-Re
  controller or time-averaged terminal state is needed.
- CMAME 2021 cube geometry: stock `tlbmmpm_cu_02` matches its bounding
  box (`1.45 m` cube in `24.1x6.1x6.1 m`, `241x61x61`), but the stock
  `AddRectangularBeamMesh(...,ndiv=1)` has 24 tetrahedra/14 corners,
  whereas CMAME 2021 states six tetrahedral material points; the body
  force also differs (`0.25` vs `2.5e-5`).  Thus the stock CD is not a
  paper-equivalent cube result.
- Velocity-BC test: one-off GPU D3Q15 case with a Zou-He velocity
  inlet at `x=0` and a zero-gradient outlet at `x=Nx-1` (y,z periodic).
  The six-tet body gives `CD=0.3201` at `Re=100` (`-45.2%` vs Eq. 38);
  the 24-tet body gives `0.3933` (`-32.7%`).
  A mini sweep with the six-tet body gives errors `-26.2, -29.8, -31.8,
  -36.1, -45.2%` at `Re=5,10,20,50,100`, i.e. no match and no full sweep.
  A density outlet was unstable (late force/CD went to zero), so the
  zero-gradient outlet was used.
- Original `tlbmmpm_cu_02` body-force sweep (stock 24-tet cube,
  `nu=0.0020833`, grid 241x61x61, LES Sc=0.17): original `b=0.25` gives
  final Re=37.99; `b=0.80` gives Re=95.54; `b=0.84` gives Re=98.91
  (1.1% below 100).  The `b=0.84` CD is 0.3551 vs Eq. (38) 0.5897
  (-39.8%).  This is the closest match to Re=100 without changing the
  original case values.  Deliverables under `/tmp/lbm_cube_stock/final`:
  `cu02_cube_clip.gif`, `visit_cu02_cube_clip.py`, `lbm.visit`,
  `mpm.visit`, `cu02_bforce_sweep.png`, `cu02_re_vs_b.png`,
  `cu02_CD_vs_formula.png`, `cu02_CD_formula_final.csv`.  Final CD
  errors vs Eq. (38): b=0.25 -28.0%, b=0.80 -39.4%, b=0.84 -39.8%;
  matching Re makes the relative CD error worse.  Historically the facet term
  `-2 rho W Cs C` was multiplied by `LBMMPM::Domain::ForceCorrection`
  (default 0.6, uploaded in `LBMMPM::lbmmpm_aux`); that member and device
  struct were removed, so CPU and CUDA now use the full term.  The CPU
  `ImprintLattice` and CUDA `cudaImprintFacets` both clamp a negative `newF`
  to zero.  MechSys host/CPU code is double; CUDA precision follows the
  user-selectable `USE_GPU_DOUBLE`, which `mechsyscc` and CMake enable by
  default.  Historical `FC=0.6` runs at b=0.84 gave Re=99.160,
  CD=0.585829 vs Eq.38 0.588541 (-0.46%), and at b=0.25 gave CD=1.2803
  vs 1.2681 (+1.0%); these are not the current default.  The VisIt
  animation clips the LBM at y=3.05 and colors `|Velocity_0|`; the cube
  is a filled magenta MPM mesh.  Use `magnitude(Velocity_0)`, not the
  vector variable directly; pseudocolor of the vector produced a nearly
  blank flow field.
- CMAME 2021 low-viscosity values (`nu=2.0803e-6`, source `bforce=2.5e-5`,
  i.e. `fluidforce=6.25e-10`, six-tet, LES Sc=0.17): full 100,000 steps
  ends at `Re=10.553`, `CD=-9170.0` vs Eq. (38) `3.533`.  The CD is
  sign-changing and reaches `-2.4e6`/`+2.2e5` during the run; this is
  numerical instability, not a drag result.  The 24-tet 20k check is also
  negative (`Re=10.63`, `CD=-1823.6`).  Sc=2 gives `CD=161.5` at 20k,
  still not Eq. (38).  Directly setting physical `fluidforce=2.5e-5`
  gives `Re=4251` and `CD=-1.53` after 1,000 steps.  That viscosity and
  body force therefore cannot be validated with the present low-Tau
  SRT/MPM coupling; no sweep.
- LES is **implemented and active** in the GPU LBMMPM path:
  `cudaCollideMPM` (`lib/lbmmpm/lbmmpm.cuh:322`) computes the non-equilibrium
  `Q` and uses `tau=0.5*(tau+sqrt(tau*tau+6*Q*Sc/rho))`; the caller passes
  `LBMDOM.Sc` (`lib/lbmmpm/Domain.h:781`) and `Sc=0.17` is uploaded in
  `FLBM::Domain::UpLoadDevice` (`lib/flbm/Domain.h:2227`).  A runtime check
  with `Sc=0` versus `Sc=0.17` at Re=100 changes the 1,000-step CD from
  `0.2682` to `0.3117`, so it is not dead code.  The GPU formula mirrors the
  CPU `CollideMPM`, but uses the first power of `Sc`; `cudaCollideSW` uses
  `Sc*Sc`.  All velocity-BC results above already had LES on (Sc=0.17).
  The low-viscosity run is still wild with LES on (and with Sc=2), so LES
  is not the missing piece and does not by itself stabilize Tau=0.500125.
- Major separate discrepancy: the current cube uses
  `nu=dx^2/(dt*24)=0.0020833 m^2/s`, **~1000x** the CMAME 2021 low-viscosity
  value `nu_f=2.0803e-6 m^2/s`.  At `Re=100` that means inlet
  `u=0.1158` (Mach 0.116) and lattice `Tau=0.625`, versus `u=1.16e-4`
  and `Tau=0.500125`.  Using the low viscosity made the current
  SRT/coupling unstable (CD of order +/-1e3); an intermediate `nu=1e-4`
  was also unstable (negative CD).  The velocity BC therefore does not
  make the current case equivalent to that low-viscosity setup.
- Six-tet cube run (manual 8-corner, 6-tet mesh; BFac=12, total particle
  volume `3.048625 = 1.45^3`): with source `bforce=0.843` the full
  100,000-step GPU run reaches **Re=100.03** at `t=19800`; `CD=0.30258`
  versus Eq. (38) `0.58444`, **-48.2%**.  At matched Re the 6-tet body is
  *worse* than the stock 24-tet body (Re=38: -35.0% vs -28.0%; Re=20:
  -28.1% vs -22.1%).  Changing only the mesh to the six tets stated in
  CMAME 2021 therefore does not recover the match; the remaining
  difference is in the force/coupling normalization or another model detail.
- Full settled cube (`tlbmmpm_cu_02`, 100,000 steps (`Tf=20,000`,
  `dt=0.2`), `.res`-only output, no HDF5): `Solve=5m54.3s`; 99 frames, `Re` 1.06 -> 37.99, `CD` 41.11 -> 0.913,
  built-in `128*Re^-0.8/5.5` 22.25 -> 1.268.  Final relative error -28.0%;
  last-20% mean `CD=0.928` vs formula `1.286` (-27.8%).  The formula is the
  case's empirical correlation, not exact cube theory.  `tlbmmpm_cu_01` is
  the elastic beam; `tlbmmpm_cu_02` is the cube.
- Cube domain size: length and width behave oppositely.  At matched `Re=30`,
  the formula error is -25.5% for 241x61x61, -24.1% for 361x61x61, -23.2%
  for 481x61x61, but -39.2% for 241x81x81 and -41.5% for 241x91x91.  At
  `Re=35`: -27.0%, -25.3%, -24.1%, -40.7%, -43.1%.  A longer streamwise box
  gives a modest real improvement (up to ~3.3 percentage points); a wider
  transverse box lowers blockage/`CD` and makes the match much worse.  The
  361/481/81/91 runs stopped at `t=9800` (50,000 steps) with `CD` still
  drifting; only the 61 run is full.  Test a smaller transverse box or a
  different correlation/normalization next.
- Legacy `ImprintLattice` check: `/home/sergio/mechsysv1` (2020-source copy)
  has one extra active assignment before the facet force:
  `F[nxs][nys][nzs][k] = F[nx][ny][nz][k]`.  It was commented out in
  `e1ad057` (2024-07-23) and is absent from `cudaImprintFacets`.  Restoring
  it on the current CPU t02 cube (5,000 steps to t=1000) gives
  `CD=31.52` at `Re=5.02` vs current `6.83` and formula `6.40`; a
  force-only `Fk` variant gives 28.59.  The line corrupts the
  solid-neighbour population and worsens the match; do not re-enable it.
  If the old code matched better, its cause is elsewhere (old MPM/LBM force
  normalization or reporting).
- Older measurements: reduced t02 (61x21x21, 500 steps) 0.937 s CPU solve vs
  0.617 s GPU; full t02 window (241x61x61, 1000 steps) 39.62 s vs 5.80 s
  (6.8x); reduced t01 beam (40x40x40, 277 steps) 1.285 s vs 0.806 s;
  full-grid t01 startup window (200x200x200, 28 steps) 14.69 s vs 6.35 s
  (2.31x).
- `Re`/velocity agree closely, but GPU `CD` is lower by 8.5% at reduced t02
  and 24% at the end of the full t02 window; corner `h` is accumulated with
  `atomicAdd` and the CPU/GPU update order differs.  Do not quote GPU `CD` as
  validated.  The t01 deflection differs by 4.1% at `t=0.015`.
- The GPU `Solve` must imprint once after `UpLoadDevice` and before the first
  output; otherwise `frame 0000` has an empty initial body `Gamma`.  With that
  call, all common HDF5 dataset payloads at `frame 0000` are byte-identical.
  Whole HDF5 files still differ by metadata/timestamps (12 bytes LBM, 22 bytes
  MPM); compare datasets, not file MD5.  After one step the payloads diverge,
  so the CUDA trajectory is not bitwise.
- A byte comparison of the current benchmark outputs: reduced t02 `frame 0001`
  has exact `Gamma` but 8706/28862 differing `Density_0`/`Velocity_0` bytes;
  full t02 `frame 0001` has 2.95M/13.4M; t01 reduced `frame 0001` already has
  MPM `Sigma/Strain/PVelocity/Position/Vertex` differences.

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
