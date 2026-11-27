# Phase 4: independent motion of ultrasonically levitated particles at 99% modeled power

26 November 2026 · Numerical simulation · 512 emitters · 16 particles · No physical experiment

**Decision: the tested controllers did not demonstrate independent control of all 16 particles under the complete benchmark.** Phase retuning recovered stationary stiffness without increasing steady modeled power. Short time-resolved simulations demonstrated selected motion with other beads held, but failed the project’s conservative frame-power gate. All five controller designs also missed both coordinated all-16 moves in a local screening model. This identifies failures in the tested designs and qualification procedure, not a physical upper bound on the array.

## 1. Reproduction without the missing archive

The saved [design](https://github.com/jihvijhojhviihogyuvi/ultrasonic-research-project/blob/90c4651539173cc3fc86adf9aa8c54d7b57413ce/outputs/ultrasonic_phase3/results/experiments/compare_16_covariance/design.npz), [case configuration](https://github.com/jihvijhojhviihogyuvi/ultrasonic-research-project/blob/90c4651539173cc3fc86adf9aa8c54d7b57413ce/outputs/ultrasonic_phase3/results/experiments/compare_16_covariance/case.json), [equations](https://github.com/jihvijhojhviihogyuvi/ultrasonic-research-project/blob/90c4651539173cc3fc86adf9aa8c54d7b57413ce/outputs/ultrasonic_phase3/physics.py) and [published three-power comparison](https://github.com/jihvijhojhviihogyuvi/ultrasonic-research-project/blob/90c4651539173cc3fc86adf9aa8c54d7b57413ce/outputs/ultrasonic_phase3/16_particle_power_margin.md) were sufficient. Inputs were pinned to `90c4651539173cc3fc86adf9aa8c54d7b57413ce` and their Git blob and SHA-256 digests checked. Existing repository files were not modified. The former ZIP was not used.

| Parked configuration | Steady power (W) | Weakest stiffness (µN/m) | Worst equilibrium error (mm) | Largest linear multiplier |
|---|---:|---:|---:|---:|
| Original, 100% | 1.492952 | 17.586059 | 0.041533 | 0.897091 |
| Original, 99.5% | 1.485487 | 17.496744 | 0.041468 | 0.897099 |
| Original, 99% | 1.478022 | 17.407414 | 0.041413 | 0.897090 |
| Phase retune, 99% | 1.478022 | 18.960238 | 0.041351 | 0.897089 |
| Retune attenuated to satisfy parked bound | 1.356618 | 17.365207 | 0.042350 | 0.897089 |

Requirement: 17.441031 µN/m. Each of the three original stiffness values reproduced to within 0.000001 µN/m of its published rounded value. The independent ±8-harmonic fixed-drive regression also passed; maximum discrepancy on the archived complex emitter coefficients was below 4 × 10⁻¹⁶. At 99% power the original design failed stiffness, not position accuracy or linear periodic stability.

## 2. What the input Jacobian and phase retuning establish

At fixed 99%-power equilibrium we evaluated J = ∂⟨F⟩/∂θ for the 48 real force coordinates. Each of the 512 phase variables rotates all four slots of one emitter together, preserving the magnitudes of that channel’s steady spectrum and periodic envelope. **Rank: 48/48; condition number: 21.85.** Largest and smallest singular values: 68.17 and 3.12 nN/rad, with a relative cutoff of 10⁻⁷. A centered random-direction force derivative agreed within 2e-09 relative error. A uniform phase rotation is a null direction, as expected.

The off-diagonal particle blocks account for 0.0241% of the mean stiffness matrix’s Frobenius norm at the original 99% equilibrium. This refers to *local particle–particle force derivatives*, not all shared-field interference. Full row rank rules out shortage of local force directions in these controls at this point. It does not imply a large capture range, sufficient stiffness along a finite trajectory or controllability of the full nonlinear system.

The retune optimized relative phases, not overall amplitude. We differentiated the weakest coupled stiffness directions, projected their gradient into the approximate null space of mean force, tested bounded phase steps and re-solved gravity balance. After three iterations an independent field rebuild gave **18.960238 µN/m**, a **8.71%** stationary margin at unchanged steady power. This is one saved geometry and one target cloud, not a robustness claim across manufacturing or calibration variation.

## 3. Nonlinear motion: measure the held particles, not just the commanded one

All 16 particles began at their retuned mean equilibria with zero initial velocities. The command ramp lasted 0.4 s, the displaced hold 0.2 s, the return 0.4 s and the final hold 0.2 s. Position errors below are Euclidean maxima over modeled substeps relative to the commanded paths from these starting positions; held drift is relative to each held particle’s own start. Final error is the largest 3D error among all particles at the last sample. The motions were preplanned; no camera feedback was used in these nonlinear pilots.

| Maneuver | Worst moving error (mm) | Worst held drift (mm) | Final max error (mm) | Minimum cycle-mean stiffness (µN/m) | Max frame-power bound (W) |
|---|---:|---:|---:|---:|---:|
| Particle 0: −x; 15 held | 0.06744 | 0.08326 | 0.02839 | 18.61522 | 1.610294 |
| Particle 0: −x and particle 7: +y; 14 held | 0.08326 | 0.06297 | 0.02839 | 18.29872 | 1.610292 |

At the end of the outbound hold, particle 0 reached **-0.251238 mm** on its commanded −0.25 mm x move; its final x displacement was -0.001361 mm. In the paired trial particle 0 reached -0.251282 mm in x and particle 7 reached 0.254401 mm in y (commands −0.25 and +0.25 mm). No contacts were recorded in either trial. Largest frozen-field multiplier magnitudes were 0.897091 and 0.897092, respectively; these are not stability certificates for arbitrary switching.

The original benchmark thresholds cycle-averaged stiffness. Instantaneous stiffness is a different metric: its minima were **-18.404** and **-19.006 µN/m**, respectively. If positive stiffness above threshold is required at *every* instant, neither pilot passes. Stiffness was evaluated at the commanded sites; a continuous stiffness guarantee at every actual bead coordinate was not established.

**Per-source checks.** Single-pilot maximum slot amplitude was 0.965168 (normalized limit 1); maximum sampled emitted peak was 0.965166. The independent upper bound on the steady periodic command peak was 0.965640. A sampled peak is not a continuous-time certificate for the transition. Moving and held errors against the original specified target coordinates can also be recovered with `verify_results.py`.

## 4. Five control policies at equal modeled steady power

These comparisons are *trajectory screens in a periodically varying linearization*, not full nonlinear qualifications. The state has 48 position and 48 velocity coordinates. For each policy we ran 96 reset single-particle trials (16 × 3 axes × 2 signs) and two coordinated all-16 signed trials, both at 0.25 mm. Each trial lasted 1.2 s. A screen pass requires both worst moving error and worst held drift ≤0.10 mm; the all-16 commands have no held beads. All designs use the same retuned state, 512 sources, 0.8-rad local phase trust region and 99% steady ceiling.

| Policy | Single passes / 96 | Beads passing all six directions (screen only) | All-16 passes / 2 | Worst moving error (mm) | Worst held drift (mm) |
|---|---:|---:|---:|---:|---:|
| Direct force | 89 | 12 | 0 | 0.1068 | 0.0831 |
| Support plus correction | 89 | 12 | 0 | 0.1784 | 0.0834 |
| Temporal multiplex | 64 | 0 | 0 | 0.3924 | 0.0838 |
| One-frame predictive | 84 | 9 | 0 | 0.1319 | 0.0957 |
| Joint two-frame | 81 | 7 | 0 | 0.1423 | 0.1001 |
| No movement command (control) | 0 | 0 | 0 | 0.2852 | 0.0824 |

**No architecture passed either all-16 move.** The two simultaneous direct-force errors were 0.1068 and 0.1049 mm; both commands reached the tested phase trust-region limit. Direct force and support-plus-correction passed 89 signed moves each, but that is not 89 particles or a qualified count. The parked control passed no commanded motion. The original controller was compared at rest; its full nonlinear closed-loop 16-particle benchmark was not run.

**Implemented mechanisms.** Direct force uses Δθ = J⁺KΔq, with zero desired displacement at held sites. The cached dominant online operation is a 512 × 48 synthesis; the selected nonlinear pilots also require an offline force trim. Support-plus-correction leaves half the emitters providing the support waveform and uses 256 for correction. Its restricted Jacobian still has rank 48, but more phase change is required; all 512 sources still emit. Temporal multiplexing alternates twice the correction with an uncorrected cycle: local clipping, switching transients and micromotion can defeat the time-average argument. It requires switching at the acoustic-cycle level, not simply changing a 20 Hz camera command.

One-frame predictive control minimizes regularized error on a one-frame propagation of the coupled state. Joint control instead fits the next two frames for all coordinates with a penalty on input changes and threefold residual weights on held coordinates. The horizon has 96 control coefficients for joint fitting; gains are precomputed. Predictive and joint screens include 20 ms modeled observation delay and 10 µm noise (seed 311); open-loop designs do not use the observations. All designs are compatible with the modeled phase channels. Local-model error, saturation and delay limit these gains, and none has demonstrated hardware performance.

## 5. Exact rejection criterion, numerical checks and limits

**The power definitions matter.** The steady modeled spectral sum is Σₗₙ wₗₙ|aₗₙ|² = 1.478022 W, meeting the 99% ceiling. The [original trajectory validator](https://github.com/jihvijhojhviihogyuvi/ultrasonic-research-project/blob/90c4651539173cc3fc86adf9aa8c54d7b57413ce/outputs/ultrasonic_phase3/validation.py) uses a conservative frame-energy bound based on each emitter’s maximum radiation weight over sidebands. Even at rest this gives Σₙ maxₗ(wₗₙ) Σₗ|aₗₙ|² = **1.610291 W**, 8.95% above **1.478022 W**. This explains a rejection without claiming that physical power exceeded the ceiling. We deliberately completed the simulated moves after this diagnostic failure to learn about motion; the original validator would stop and would not credit a pass.

**Derating alone did not work.** Uniform amplitude scaling by 0.958050 meets this parked bound. Re-solving gravity balance gave 17.365207 µN/m at 1.356618 W steady: below 17.441031. Thus the phase-only retune solved one constraint, but not stiffness and conservative power together. No physical capacity theorem follows.

**Fidelity and timing.** Baseline runs retained first-Born coupling, gravity, finite disks, source spectral filtering and all 48 force coordinates. The full-band checks used ±40 harmonics, 216 disk points per emitter and 256 phases per cycle. For motion, finite-disk jets were rebuilt at commanded sites; a local cubic pressure approximation, a first-order spatial correction to interaction forces, nonlinear drag, filtered source transitions and eight propagation-delay bins supported time integration at 128 steps per 5 ms cycle. In 21 single-pilot spot checks the maximum force discrepancy against freshly evaluated disk/Born forces was 0.0681% of particle weight. Checks are sampled after transitions, not a uniform error bound. Doubling both phase resolution and delay bins for the highest-error first frame changed matched-time positions by at most 0.00002331 mm; this was not a full-trajectory refinement.

The single pilot issued 24 preplanned commands on a 20 Hz simulated schedule. Cached slot synthesis took at most 0.078 ms, but offline force correction alone took up to 124.5 ms. Computing steady/transition fields for validation is also offline; per-frame times are in the raw records. The single pilot used 147.5 s of CPU-environment wall time and peak process RSS of 407 MiB. These are not measured end-to-end camera-to-driver update times. Blocked integration and a compact derivative cache avoided the original degree-8 transfer array’s full-band 1,670.625 MiB estimate, above its 665.6 MiB guard. No GPU performance was measured.

**Scope.** One saved target cloud in idealized air, fixed aperture, bead properties and frequency. No streaming, boundary effects, source-calibration variation, mutual radiation-impedance measurement or physical experiments. The tangent comparison neglects finite-command transition nonlinearities. Each policy had one 60 s parked tangent continuation, not 60 s of nonlinear closed-loop retention. Other transducer/particle counts were not attempted because the all-16 prerequisite failed. Scaling to 25–100 cannot be inferred from this study.

## 6. Acceptance decision and reproduction

| Engineering target | Status and scope |
|---|---|
| 16 particles, 512 emitters | **Pass** for the simulated geometry. |
| Independent control of all 16 | **Not demonstrated.** Signed nonlinear trajectories for every particle and a continuous full protocol were not tested. |
| Moving error and held drift ≤0.10 mm | See the full-model pilot maxima in §3; **not established for all particles**. |
| Simultaneous independent moves | Selected pair tested in §3; **fail** for both all-16 patterns in every local controller screen. No full nonlinear all-16 qualifying run. |
| Weakest stiffness ≥17.441031 µN/m | **Pass** for the retuned parked state and minimum frozen cycle means at the pilot target sites. **Fail** if required at every instant. |
| Modeled power ≤99% of nominal | **Pass** for frozen spectral watts; **fail** under the existing conservative trajectory bound. |
| 60 s retention | **Not tested** by a nonlinear closed-loop run. |
| ≥20 Hz end-to-end updates | **Not tested.** Timed playback is not a closed-loop deployment. |
| Hardware validation | **Not tested.** |

**Next discriminating experiment.** Jointly optimize slot amplitudes, timing and phases against the original conservative frame bound, source-peak limits and minimum stiffness at parked and displaced waypoints. Compare with a single-carrier support design using the same ceilings. Reject a candidate at the stationary gate unless it passes both constraints; do not loosen the bound without separately validating a better power model. For a surviving candidate, run every signed move, mixed simultaneous moves with held beads, then a continuous 60 s nonlinear retention and latency test. The forecast is a testable proposal, not a simulated achievement.

**Reproduce.** Use the companion source and data package; the old ZIP is not required. CPU environment: Python 3.13.18, NumPy 2.3.5, SciPy 1.17.1, Numba 0.66.0, and g++ with C++17. Original inputs are in `inputs/`, original sources in `baseline/`, and Phase 4 settings in `config_phase4.json`. Use a new output path. All raw trajectory files are SI units; tables above convert to mm and µN/m.

```sh
python -m pip install -r requirements.txt
python verify_results.py
python run_phase4.py --out ../phase4-replay --stage all
```

`results/` records stationary fields, stiffness matrices, force Jacobian and singular values, phase-search history, per-frame commands, nonlinear substep positions and velocities, and all 98 trials per controller. `report/` contains the numerical tables and figures. `SHA256SUMS` lists package file hashes. The earlier file-write problem was addressed by writing nonlinear trajectories per frame; completed reruns, not incomplete records, are used for this report. The [Phase 3 report](https://github.com/jihvijhojhviihogyuvi/ultrasonic-research-project/blob/90c4651539173cc3fc86adf9aa8c54d7b57413ce/outputs/ultrasonic_phase3/phase3_report.md) remains the reference for the previously qualified two-particle simulations.
