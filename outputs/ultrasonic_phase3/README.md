# Ultrasonic control: executed Phase-3 study

Start with [phase3_report.md](phase3_report.md). The report separates stationary
stiffness, nonlinear periodic stability, finite independent motion and computer
deadlines. Nothing in this package establishes a physical particle-count record.

## Run

Python 3.13 was used with the versions in `requirements.txt`. PyTorch uses CUDA
when available and falls back to CPU; the compiled nonlinear integrator is CPU.
All seeds and resource limits are recorded in each `case.json`. Optimization is
local and its unsuccessful attempts are not impossibility certificates.

```powershell
python -m pip install -r requirements.txt
python -B tests.py
python -B tests_dynamics.py
python -B experiments.py --group comparisons
python -B experiments.py --group temporal
python -B experiments.py --group refinements
python -B experiments.py --group scaling
python -B experiments.py --group power_density
python -B convergence.py --nonlinear
python -B convergence_high.py
python -B power_match.py
python -B sdp_recovery.py
python -B validate_cases.py compare_2_static compare_2_covariance --samples 64
python -B build_report.py
python -B audit_package.py
```

Run these commands from this directory. A case with existing `case.json` is
skipped by default; `--rerun` deliberately replaces it. Use distinct validation
tags to retain different motion resolutions or controller versions. Do not
replace the archived raw data when testing a new design.

`run_all.py` provides automation for the complete pipeline. Full covariance
trajectories are substantially slower than stationary screens. Memory limits in
`Config` bound specific caches, not the entire Python process. High-count screens
are explicitly distinguished from completed independent-motion tests.

If pip temporary-directory permissions fail on this recorded Windows runtime:

```powershell
python -B install_task_dependencies.py C:/path/to/task-deps
$env:PYTHONPATH = 'C:/path/to/task-deps'
```

This alternative verifies every wheel against `dependency_wheels.json`. Install
PyTorch and NumPy separately if unavailable; they were already installed on the
execution host. On other operating systems use normal pip-compatible wheels.

## Contents and units

- `physics.py`: finite disks; pressure derivatives; filtered spectra and beats;
  Born/Foldy Rayleigh scattering; full coupled Jacobians/Floquet matrices.
- `optimization.py`: IB, force, covariance, support, two-tone, schedule/gauge,
  order, PSD and nonlinear-flow optimization.
- `dynamics.py`, `validation.py`: converged local pressure surrogates, source
  envelope delays, nonlinear drag and contacts, noisy delayed camera feedback,
  every signed move and simultaneous commands, resource/deadline checks.
- `experiments.py`: automated attempts with raw coefficients and histories.
- `convergence*.py`: spatial/spectral/time and nonlinear shooting checks.
- `reproduce.py`, `independent_legacy.py`: legacy reproduction and separate
  NumPy/SciPy physics recalculation. Verified inputs are in `inputs/legacy`.
- `benchmark_controller.py`: emitted-command timing on saved camera inputs;
  every control coordinate is checked against the saved nonlinear trajectory.
  This separates controller cost from offline verification/simulation oracles.
- `results/`: every executed attempt, arrays, trajectories, convergence and
  source hashes; `all_experiments.csv` is the compact numerical inventory.
- `figures/`: plots generated only from those saved calculations.
- `provenance/`: execution logs, dependency/runtime details and source snapshots.

SI units are used internally: metres, seconds, Pa peak phasors, m/s velocity,
newtons, joules, watts and N/m. Plot/report conversions are explicit. Power is an
isolated-piston **acoustic model**, not measured electrical input or true coupled
array radiation power. Source drive `1` is a normalized peak surface-velocity
command, not a voltage specification.

## How to interpret outcomes

`primary_screen_rejection` tests uncoupled target forces/stiffness only.
`stationary_screen_rejection` includes a coupled equilibrium and periodic test.
`stationary_candidate_not_motion_validated` is never an independent count.
`strict_pass` requires all signed/simultaneous moves, ≤0.10-mm target/held error,
≥17.441-µN/m weakest stiffness, ≥60 simulated seconds, no domain/contact failure,
resource compliance and ≤50-ms measured emitted-drive controller updates.
Separate slow-threshold SDP runs do not meet the main benchmark automatically.

A surrogate domain exit means numerical fidelity is no longer established. It
must not be relabeled as a physical particle escape. A finite polynomial fit is
checked against fresh finite-disk derivatives before trajectory acceptance.

The default 1-mm ramp lasts two seconds; the 0.25-mm ramp lasts 0.4 seconds. These
are slow spatial motions with 20-Hz command frames. Reported controller replay
passes have identical saved nonlinear motion and identical emitted commands;
their timing is measured separately and is not an OS hard real-time guarantee.

Streaming, boundaries, thermoviscous/finite-size scattering, particle absorption,
real source impulse variation and mutual radiation impedance remain excluded.
Measure these before treating any modeled count, stiffness or watt value as a
physical display specification.
