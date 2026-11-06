# List of changes:

- Introduced velocity Verlet as integration scheme
- Added contact outputs (Normal and Tangential Overlap, Rolling, and Dashpot components)
- Added four mode of Coulomb criterion. The default mode (0) is the Legacy one (i.e. elastic threshold and persistent dashpot), mode 1 applies the threshold on the total force (elastic + dashpot) and updates the tangential overlap to reduce the elastic force so that the sum is equal to the cap. Mode 2 uses an elastic threshold and deactivates the tangential dashpot on sliding (as done by PFC on dashpot mode 2/3). Mode 3 uses an hybrid threshold (elastic in tangential, total in normal direction) and then returns the total tangential force to the threshold (legacy Yade model). 

I assume that any 'serious' simulation will use either mode 0 or 1, so only these are implemented in the GPU path to reduce cluttering (and improve efficiency). The latter two are only kept in the CPU path for illustrative purposes.

I implemented the translational part of velocity Verlet and the CPU version of the rotational integration. I made a mess with the rotational one on the GPU so I used AI (GLM-5.3 iirc) to fix and test it. Likewise, the Lees-Edwards implementation was first done manually (copying the one in MercuryDPM) and then fixed/improved with AI (my own implementation triggered Verlet list rebuilding at every step).

Because the DEM-LBM coupling was carried out with Position Verlet in mind, the scheme does not work with velocity Verlet. A mix of GLM, Deepseek and Codex was used to implement an implicit coupling with velocity verlet, which has only been tested on a few benchmarks and so it is still tentative/experimental.


# AI Review:
The changes are reviewed by a deepseek agent in the following:

- **Rotation.** The registered regressions run with zero angular velocity, so
  `w` stays zero there and the vertex-rotation branch of `OrientationUpdate`,
  which requires a non-zero `w_half`, never fires; they assert only that a
  pre-force rebuild preserves `Wdot`, and `vertex_rotation()` separately feeds a
  hand-rotated cube to `MaxD` requiring `sqrt(2)*sin(theta/2)` with no centre
  motion. The **analytic check of the rotational path is `tdem/test_cu_02`**:
  a sphere on a spinning turntable, both bodies given spin so the contact starts
  without sliding, follows the closed-form orbit
  `x = x0 + vx/wr*sin(wr*t)`, `y = vx/wr - (vx/wr)*cos(wr*t)` with
  `wr = (2/7)*omega`. That solution exists only because the sphere's spin
  couples into the contact through friction, which makes it a sharp test of the
  rotational integration. The bounded verifier runs to approximately 0.1 s
  and checks the reported endpoint against a `1e-3` m tolerance for both
  integrators through `tests/dem_contact/check_turntable.py`. This is an
  analytic endpoint check of rotation acting through contact friction.
  What is still **not** covered: a free spinning body with no contact, so there is
  no direct check of angular-momentum conservation in the absence of torque, and
  the quaternion half kick is exercised only through its effect on the contact.
- **Coupled VV–LBM.** The path is opt-in, default-off, described as **first
  order**, and its guards reject configurations it cannot represent. It is not
  claimed to be accurate or fast, and the legacy coupled solver is unchanged.
  `dem_coupling_smoke` runs three tiny `Selection::ImplicitVV` fixtures and
  checks the combined mass and momentum against fixed budgets, so the exchange
  arithmetic does have a numeric guard rail. It is an endpoint check on three
  fixtures, not a physical validation, and the rejection paths themselves are
  only exercised by `tlbmdem/tlbmdem_cu_benchmark --guards`, which is not part
  of the registered regression suite.

## Sphere contact controls

These controls apply to linear and Hertz sphere contacts. All four Coulomb
modes are available on CPU; CUDA supports modes 0 and 1. Tensile cutoff and
first-contact correction are available on both backends. These settings do not
select new polyhedral contact laws.

| Domain member | Default | Meaning |
| --- | --- | --- |
| `SphereCoulombMode` | `0` | Mode 0, the legacy law: cap on the elastic tangential force using the elastic normal force, tangential dashpot retained permanently. |
| `SphereCoulombMode = 1` | opt-in | Mode 1, LAMMPS: cap on the total tangential force using the total normal force. |
| `SphereCoulombMode = 2` | opt-in | Mode 2, PFC with dashpot mode > 2: elastic threshold, tangential dashpot disabled once sliding. |
| `SphereCoulombMode = 3` | opt-in | Mode 3, Yade's legacy law: `mu` times the total normal force as the threshold, compared with the *elastic* tangential force, then the total is projected onto that limit. |
| `SphereTensileCutoff` | `false` | Clamps a tensile elastic-plus-dashpot normal result to zero. |
| `SphereFirstContactCorrection` | `false` | Scales the first tangential increment by the fraction of the step spent in contact, interpolating between previous and current signed overlap. |

The four modes are three independent choices crossed, and the numbering follows
the choice of trigger:

| Mode | Normal force in the `mu` threshold | Tangential force compared with it | What is then capped | Dashpot while sliding | Reference |
| --- | --- | --- | --- | --- | --- |
| 0 | elastic | elastic | elastic | kept | legacy |
| 1 | total | total | total | kept inside the cap | LAMMPS |
| 2 | elastic | elastic | elastic | dropped | PFC, dashpot mode > 2 |
| 3 | total | elastic | total | kept inside the cap | Yade (legacy) |

Modes 1 and 3 use the same total-force projection when their respective
triggers fire. With `L = mu*|Fn_total|` and trial total `T`, the projected force
is `T` if `|T| <= L`, otherwise `L*T/|T|`. The zero vector stays zero. The
elastic history is recovered as `Ft_elastic = Ft_projected - Ft_dashpot`.
This operation never increases the total-force magnitude.

Mode 1 triggers on the total tangential force; mode 3 triggers on its elastic
part. Neither comparison is universally stricter because the spring and
dashpot vectors can oppose each other. If mode 3 triggers while the total is
already inside the limit, the total remains unchanged. If its elastic trigger
does not fire, the total can exceed the limit. Modes 0 and 2 can also exceed a
total-force threshold while their dashpots are present. Mode 2 removes the
dashpot only after its elastic sliding trigger fires.

Measured on a single contact with `Kn=100`, `Gn=2`, `Kt=20`, `Gt=3`,
`Mu=0.5` and a normal approach giving `|Fn_elastic|=10` and
`|Fn_total|=12` (so the two limits are 5 and 6):

```
mode 0: Ft_elastic=+5.000 Ft_dashpot=+12.000 total=+17.000   (limit 5, caps elastic)
mode 1: Ft_elastic=-6.000 Ft_dashpot=+12.000 total= +6.000   (limit 6, caps total)
mode 2: Ft_elastic=+5.000 Ft_dashpot= +0.000 total= +5.000   (limit 5, caps elastic, no dashpot)
mode 3: Ft_elastic=-6.000 Ft_dashpot=+12.000 total= +6.000   (limit 6, caps total)
```

Modes 1 and 3 give the identical total and the identical decomposition. The
negative elastic component in those rows is not a sign error: the elastic
history is defined as whatever remains after the dashpot is subtracted from the
cap.

All four modes exist on the CPU path (`CInteractonSphere::CalcForce`, both the
linear and the Hertz law) and in the fast `CalcForceSphere` path. **Only modes
0 and 1 exist on CUDA**: `CalcForceVV` and `CalcForceVV_Hertz` clamp anything
above 1 to 0 and the host warns once when it does so, because uploading mode 2
or 3 while running a different law would be silently wrong. Any value above 3
falls back to mode 0 everywhere.

The cutoff's reported normal dashpot component includes the clamp adjustment:
elastic plus dashpot reconstructs the normal force actually applied. In modes
1 and 3 the elastic tangential component can oppose the dashpot; their sum is
the projected total force when the cap is triggered. The optional onset
correction estimates an activation fraction; it is not an exact collision-time
solver. Energy diagnostics follow the
selected force decomposition. The two-sphere example accepts these settings
in `tdem/test_2spheres.inp`; its CPU and CUDA targets share one source.

## HDF5 and XDMF contact output

`WriteBF` exports a particle-pair row when there is an active overlapping
feature **or any nonzero contact-force component**, including tangential,
dashpot, or custom thermostat force. An overlapping contact remains visible
when its resultant force is zero. A separated candidate without a force is
omitted. Existing dataset names and float32 visualization precision remain
unchanged; solver precision is independent of visualization precision.

| Dataset | Meaning |
| --- | --- |
| `Normal`, `Tangential` | Elastic normal and tangential contact force, summed over active features. |
| `DashpotForceNormal`, `DashpotForceTangential` | Normal and tangential dashpot contributions. |
| `Rolling` | Sphere rolling-resistance force history; zero for nonspherical rows when a mixed output contains this dataset. |
| `NormalOverlap` | Maximum positive active feature overlap for the pair; sphere overlap is radius sum minus separation. Zero for force-only rows. |
| `TangentialOverlap` | Sphere tangential displacement history; zero for nonspherical pairs, which have multiple feature histories. |
| `ID1`, `ID2`, `Position`, `Branch` | Existing pair identification and display geometry. |

All new fields have matching XDMF dataset references. Polyhedral
`NormalOverlap` is a maximum, not a center-distance estimate or an average.
Custom CUDA force kernels that replace the built-in kernels must maintain
the common contact count and overlap fields to supply geometric-only rows.
`RecordContactOverlap` provides that bookkeeping without changing forces.

## Coupled VV: existing experimental API

The separate `LBMDEM::Experimental::Benchmark::Session` API is opt-in through
`MECHSYS_EXPERIMENTAL_BENCHMARK_COUPLING`; the example target is built only
with `A_WITH_EXPERIMENTAL_VV_COUPLING=ON`. Include
`mechsys/lbmdem/experimental/benchmark_coupling.h`, select
`Selection::ImplicitVV`, set `DEMDOM.UseVelocityVerlet=true`, and supply a
finite `max_intervals` before initializing and advancing the session.
`tlbmdem/tlbmdem_cu_benchmark.cu` is the complete usage example.

This implementation combines mechanical VV with implicit rigid-body
momentum exchange, retaining physical particle mass and inertia. Opposite
particle impulses are computed from the actual fluid population increments.
The supported configuration is FP64 CUDA, smooth single-fluid D3Q15,
canonical spheres, free or fully prescribed bodies, and the original sphere
contact-law subset. Its guards reject unsupported combinations. It includes
the previously developed sphere/grid occupancy calculation and optional
compensated density summation. Coupled restart, partial constraints,
prescribed spin, multiple fluids, and general polyhedral coupling are not
supported by this API. Ordinary `LBMDEM::Domain::Solve` retains its upstream
legacy algorithm. The experimental session downloads particle/fluid snapshots;
call `DEMDOM.DnLoadDevice(1,true)` before `DEMDOM.WriteBF(...)` when current
contact diagnostics are required from that session.

Incoming fluid populations and entry geometry are still used. Mechanical VV
does **not** make the complete coupling second order: the composition is
conservatively described as first order. Prior restricted validation does
not establish general continuum accuracy, universal stability, or a speedup.
The publication integration does not reopen or extend that physical study.

## Additional integration support

Quaternion helpers, Hertz force/displacement history conversion, preserved
prescribed-velocity position history, and pre-force rebuild state handling
are required by these paths. The explicit, versioned standalone VV
`SaveCheckpoint`/`LoadCheckpoint` API and persistent particle IDs accompany
the imported validated VV implementation. Legacy HDF5 `Save`/`Load` retains
upstream's history serialization; coupled checkpoint continuation is rejected.
No performance experiment, mixed-precision path, unfinished application
driver, or private evidence package is part of this contribution.

## Example applications

Three complete drivers exercise the optional integration and contact paths.
They are ordinary example programs, not regression tests, and they are built by
default with the other `tdem` targets.

| Target | Source / input | Purpose |
| --- | --- | --- |
| `PID_v2` | `tdem/PID_v2.cpp`, `tdem/test_pid_v2.inp` | Wall-particle servo that compresses a box to a target pressure and then shears it at constant axial wall velocity. Uses velocity Verlet and Coulomb mode 1. |
| `TX_PID_periodic` / `TX_PID_periodic_gpu` | `tdem/TX_PID_periodic.{cpp,cu}`, `tdem/test_tx_pid_periodic.inp` | Same PID workflow driven through the periodic cell instead of walls: the stress is measured from the contact virial and the servo resizes `Xmin`/`Xmax`/`Per`. The GPU variant reduces the virial on the device. |
| `shearbox_LE_const_vol` | `tdem/shearbox_LE_const_vol.cu`, `tdem/test_shearbox_LE.inp` | Constant-volume three-phase periodic shearbox with Lees-Edwards kinematics. It owns its step loop and supplies the Lees-Edwards branch vector, force kernels and reset through `Domain`'s existing function pointers, so no Lees-Edwards support is required inside the library. |

All PID drivers select velocity Verlet explicitly. The shearbox uses legacy
position Verlet during compression and velocity Verlet during shear; its local
shear loop preserves contact histories and checks candidate validity after the
drift and before forces. Nothing here changes a library default. `TX_PID_periodic` and
`TX_PID_periodic_gpu` share the same physics but are separate sources, so the
CPU and GPU builds can be compared directly. The shearbox example in
`tdem/test_shearbox_LE.inp` is a small, quick configuration (1188 spheres); the
solvers read every physical parameter from the input file, so larger
production boxes only need a different `.inp`.

## Selecting the integrator in the CUDA examples


The four CUDA example drivers take the integrator as an optional trailing
argument, so the same binary runs both schemes and prints the selection it
actually used.  Position Verlet remains the default, matching the library.

| Driver | Invocation |
| --- | --- |
| `tdem/test_cu_01` | `test_cu_01 [Nproc] [final time] [integrator]` |
| `tdem/test_cu_02` | `test_cu_02 [Nproc] [final time] [integrator]` |
| `tdem/test_cu_03` | `test_cu_03 [Nproc] [Alpha] [Nthread] [final time] [integrator]` |
| `tdem/test_cu_04` | `test_cu_04 <input> [Nproc] [dthoot] [qin] [integrator]` |

`integrator` is `0` for position Verlet and `1` for velocity Verlet, and every
run prints `DEM_BENCHMARK_INTEGRATOR velocity_verlet=<0|1>` before solving.  The
existing arguments keep their positions and meanings, so previous command lines
behave exactly as before.  Tests 1 and 2 also accept a final time so a
comparison can be shortened without editing the source.

## Bounded regression checks

Nothing in this contribution needs a build option. A plain
`cmake -S . -B build` produces the same library, the same examples and all four
CUDA demos; `A_WITH_TDEM` and `A_MAKE_USE_GPU_DOUBLE` are existing MechSys
options that already default to `ON`, so passing them changes nothing. The
regression suite below is separate tooling and is the only opt-in part, so that a
normal MechSys build does not pay for thirteen extra test binaries:

```sh
cmake -S . -B build \
  -DA_WITH_DEM_CONTACT_REGRESSIONS=ON \
  -DA_WITH_EXPERIMENTAL_VV_COUPLING=ON
```

`A_WITH_DEM_CONTACT_REGRESSIONS` adds `tests/dem_contact`;
`A_WITH_EXPERIMENTAL_VV_COUPLING` additionally builds the coupling smoke test
and `tlbmdem/tlbmdem_cu_benchmark`. Both default to `OFF` and neither changes
the library. Build the `dem_test_*` targets, then run:

```sh
ctest --test-dir build -R '^dem_' --output-on-failure
python3 tests/dem_contact/check_exports.py build/tests/dem_contact/contact_cpu
python3 tests/dem_contact/check_exports.py build/tests/dem_contact/contact_gpu
python3 tests/dem_contact/check_histories.py build/tests/dem_contact fresh-history-output
# Analytic rotational check: run tdem/test_cu_02 first, then
python3 tests/dem_contact/check_turntable.py test_cu_02_trayectory.res
```

The analytic contact test covers both laws and all four Coulomb modes on CPU
(modes 0 and 1 additionally on CUDA),
tensile cutoff, activation correction, force decomposition, action/reaction,
and CPU/CUDA agreement. The VV test checks fast contact formation, repeated
rebuild reference resets, preserved GPU state/history, and vertex rotation.
The coupling smoke checks three tiny fixtures against fixed mass/momentum
budgets. Two-step CPU/CUDA checks also run ordinary legacy coupled Solve,
including its entry without standalone VV acceleration initialization. The history probe supports short CPU/CUDA comparisons, feature
oracles, and full-step checkpoint continuation. The turntable check asserts the
closed-form orbit of the spinning sphere in `tdem/test_cu_02`, which is the
analytic test of the rotational integration path. Export checks require Python
with `h5py` and `numpy`. These are bounded regressions, not a full workload
or scientific-accuracy certification.

## Custom CPU contact compatibility and example verification

The original four-argument virtual `Interacton::CalcForce(dt, Per, iter, law)`
remains the extension point for custom CPU interactions. The overload accepting
sphere options sets a scoped per-interaction context and dispatches through that
original virtual method. Existing four-argument overrides keep receiving solver
calls; options are restored after a return or exception. Builtin sphere contacts
consume the options, while other interactions can retain their existing laws.
This preserves source overrides; CUDA structures have changed, so rebuild all
code sharing these headers and device buffers together.

The optional regression suite also exercises the real periodic PID callbacks
and the shearbox's custom force kernels and physical-step helper. These checks
cover total contact stress, zero-duration initialization, history-preserving
rebuilds, fast post-drift contact formation, shear-image motion and periodic
wrapping. Both custom shearbox laws export overlap/count bookkeeping and force
components that reconstruct the applied force, including the linear tensile
clamp. The shearbox retains its existing custom friction laws, including their
retained tangential/rolling force histories during a temporary gap. Rebuilds
migrate those forces by pair identity; the Hertz case does not convert through
zero stiffness when separated.

After configuring a build with the two regression options above, run:

```sh
python3 scripts/verify_dem_contribution.py --build build --output fresh-verification
```

The verifier builds the required targets, rejects a build configured for another
checkout, checks each return code and integrator banner, and records source,
input and binary hashes with individual logs and a final `SUMMARY.json` verdict.

It runs all bounded regressions, exports, history/restart comparisons, legacy
save/load, all four CUDA examples in both integration modes, and all four new
application targets. Default application inputs use documented smaller windows
or boxes, recorded alongside their hashes; `--full-examples` uses the checked-in
application configurations. Output directories must be new and existing evidence
is never deleted. Configure the normal library dependencies beforehand; the
verifier does not borrow paths from unrelated projects or make timing claims.

The contact regression evaluates 88 CPU linear/Hertz analytic cases and 32
CPU fast-sphere cases. The CUDA executable repeats these host checks and adds
44 device contact comparisons for modes 0 and 1. The log labels host and device
counts separately. Cases include different mode-1/mode-3 trigger outcomes,
opposing spring and dashpot forces, exact and near cancellation, reversal, and
two-dimensional tangential directions.
