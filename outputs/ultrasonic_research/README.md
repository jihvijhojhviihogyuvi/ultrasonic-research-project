# Ultrasonic particle-control research package

Start with **research_report.md**. This package tests hypotheses; it does not establish a hardware record or a maximum particle capacity.

## Run

Python 3.11 or newer is recommended. Install the two requirements in a virtual environment:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe selftest.py
.\.venv\Scripts\python.exe simulator.py --n 128 --counts 2,4 --iterations 300 --probe-iterations 150 --out my_results
```

The supplied measured run used Python 3.13 and the installed PyTorch/NumPy versions recorded in each results folder's metadata.json. The simulator limits its own CPU threading to one thread. Importing PyTorch can take several seconds.

The full comparison evaluates eight configurations:

| CLI name | Interpretation |
|---|---|
| ib | Phase-only IB-style pressure holography with top-array phase reversal; uncorrected gravity offset |
| static_force | Single coherent field, jointly optimized complex drives for force balance and stiffness |
| multifrequency | Joint optimization at 39.6 and 40.4 kHz; both tones act on every particle; common crest constraint |
| global_local | Fixed common axial reference plus sequential single-target traps; a specific unoptimized hybrid baseline |
| time_single | Sequential single-target IB traps; includes envelope settling |
| covariance | Up to four jointly optimized temporal field modes; includes settling and cross terms during transitions |
| support_residual | One coherent support field plus bounded optimized residual; includes coherent interference |
| aperiodic | Static force optimization on a fixed slightly jittered aperture; a geometry candidate, not a layout search |

Global-local, covariance, and aperiodic code are deliberately concrete implementations, not assertions about every possible architecture with those names. Static force and covariance have the same amplitude/phase actuator capability. Comparing covariance only to phase-only IB would confound a control-algorithm improvement with a scheduling improvement.

## Useful commands

```powershell
# Canonical 512-source screening: no independence certificate without finite probes.
python simulator.py --n 512 --counts 8,16 --iterations 400 --no-probes --out screen_512

# Full local finite-motion tests; considerably more work at large M.
python simulator.py --n 512 --counts 8 --iterations 600 --probe-iterations 300 --out probes_512

# Ideal instantaneous switching versus measured-envelope hypothesis.
python simulator.py --n 128 --counts 4 --architectures static_force,covariance --tau-ms 0 --out ideal
python simulator.py --n 128 --counts 4 --architectures static_force,covariance --tau-ms 0.4 --out slow

# Separate low-bandwidth two-particle prototype; not the 15-Hz main benchmark.
python simulator.py --n 32 --counts 2 --architectures static_force,covariance --mechanical-hz 5 --iterations 500 --probe-iterations 200 --out slow_prototype

# Leading particle scattering, including a coupled force Jacobian.
python scattering_check.py --selftest
python scattering_check.py --results results_512_screen --out my_scattering.json
python convergence_check.py

# Spatial finite-difference convergence, keeping every other parameter unchanged.
python simulator.py --n 128 --counts 2 --architectures static_force --fd-um 25 --out stencil_25

# Independent motion of particle 0, while the other particle is held.
python tracking.py --architecture static_force --particles 2 --motion-hz 2 --seconds 0.5 --out motion.json
python tracking.py --architecture covariance --particles 2 --motion-hz 2 --seconds 0.5 --out motion_covariance.json
```

Use independent output directories for each run. To reproduce an archived run with the exact implementation, execute that directory's **simulator_snapshot.py** with the recorded arguments; the SHA-256 is in metadata.json. The tracking script uses the current simulator module. All seeds are explicit. Target coordinates and complex command coefficients are saved in fields.json.

## Interpretation of outputs

`results.json` and `results.csv` contain every attempted case, including failures. `stable_local_traps` counts sites in the complete attempted configuration meeting both stiffness and accuracy at converged equilibrium. It does not prove a stable subset exists after removing other sites. `requested_particles` is a target count, not an achieved count.

`tier1_local_pass` requires a converged accurate equilibrium, sufficient stiffness, full local force-Jacobian row rank, zero failed positive-axis finite probes, and stable local periodic dynamics with bounded sampled ripple. Each particle is moved by 0.25 mm along each axis in turn. The other particles must drift less than 10% of that displacement. This is a **local** independence screen. Negative directions, arbitrary trajectories, simultaneous random commands, lifetime, collisions, measurement noise, and 1-mm displacement certification remain additional tests.

The multifrequency configuration intentionally cannot receive `tier1_local_pass`: beat-driven dynamics have not been implemented. Its static mean-force results are useful, but incomplete. `floquet_radius: null` or `ripple_mm: null` denotes an unperformed or unsuccessful test, never zero ripple.

`force_rank` is a diagnostic of the available continuous complex drive parameterization, not a proof that the chosen IB algorithm can exploit it or that a full-rank system has sufficient power. The finite-motion tests are the stronger check. The independent-particle dynamics omit off-diagonal interaction blocks.

`drive_effort` is the nominal sum of squared emitter commands, averaged over temporal slots or summed over simultaneous tones. `acoustic_power_estimate_W` maps that effort to nominal isolated-piston acoustic watts at 40 kHz. It is a common conservative reference budget, **not measured electrical watts**. Real emitter impedance, mutual radiation impedance, and tone-dependent conversion losses require calibration. The pressure gain, compressibility ratio, Q, and settling constant are explicitly assumed in Config.

`optimization_s` is wall time for the offline gradient solve on this host. A null or zero does not imply zero total control cost: IB, field construction, equilibrium search, I/O, and diagnostics have separate costs. The tracking script reports whether its offline optimizer actually meets the requested frame interval. It does not claim real-time execution when the simulated trajectory was computed offline.

## Model validity and missing physics

- Linear homogeneous air and ideal baffled-piston disk quadrature, 24 source patches per emitter.
- Peak phasors and the correct negative velocity term in Gor'kov's potential.
- Small spherical nonabsorbing particles. The supplied radius gives ka about 0.183.
- Gravity and Stokes drag. Larger excursion speeds need finite-Re drag and streaming.
- Central spatial finite differences (50 micrometers by default), analytic quadrature pressure gradients, differentiable complex drives.
- Periodic first-order emitter envelopes and local Floquet checks. Sampled ripple may miss a maximum between samples; increase envelope_samples and time resolution before accepting an experiment.
- No particle scattering, acoustic streaming, nonlinear propagation, enclosure reflections, transducer-to-transducer mechanical coupling, thermal drift, sensing noise, or human intrusion.

The optional scattering_check.py adds leading Rayleigh/Born particle interactions and a coupled mean-force stability check. It has a separate r^-4 attraction control. It does not implement full multiple scattering, finite-size stress integration, or coupled temporal Floquet dynamics. The temporal envelope model also omits spatially retarded envelopes and frequency-resolved modulation sidebands; these are necessary before relying on short-slot predictions.

These omissions prevent a reliable maximum particle count or a claimed high-density breakthrough. The report specifies the interaction, boundary, and robust-control extensions needed before hardware investment.

## Files

- research_report.md: full eleven-section research deliverable.
- simulator.py: eight-way differentiable screening.
- tracking.py: local affine transient tracking test.
- selftest.py: analytic and numerical physics controls, including known failure cases.
- scattering_check.py: leading scattering sensitivity and coupled mean-field dynamics.
- convergence_check.py: fixed-field spatial stencil and source quadrature refinement.
- results_128/: executed local-screen results and source snapshot.
- results_512_screen/: larger executed mean-field screen; no finite-motion independence claim.
- tracking_static.json: executed local trajectory test.
- results_32_slow/: separately executed 5-Hz prototype-class study.
- results_cycle_1ms/: executed sequence-period sensitivity study.
- scattering_512.json and convergence_results.json: executed sensitivity checks.
- requirements.txt: minimal dependencies.

Read the report's quantitative comparison for the actual run outcomes and limitations.
