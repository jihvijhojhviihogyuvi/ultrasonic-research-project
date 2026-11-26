# How much power margin does the 16-particle candidate have?

At the original modeled acoustic-power ceiling, the archived four-field design meets the stationary stiffness requirement. Reducing the modeled acoustic power by 0.5% still passes; reducing it by 1.0% drops the weakest stiffness below the requirement. The mean position error and linear periodic stability still pass at all three powers. The failure is a stiffness shortfall, not a simulated loss of particles.

## Results on the finest tested grid

All three runs use the same saved 512-source, 16-particle configuration, controller projection, target cloud and numerical resolution. Each solves the coupled gravity-balanced equilibrium afresh. The required weakest stiffness is **17.441031 µN/m**, and the allowed worst equilibrium error is **0.10 mm**.

| Nominal acoustic power | Modeled power (W) | Weakest stiffness (µN/m) | Stiffness margin | Worst equilibrium error (mm) | Largest linear multiplier magnitude | Stationary result |
|---:|---:|---:|---:|---:|---:|---|
| 100.0% | 1.492952 | 17.586059 | +0.8315% | 0.041533 | 0.897091 | Pass |
| 99.5% | 1.485487 | 17.496744 | +0.3194% | 0.041468 | 0.897099 | Pass |
| 99.0% | 1.478022 | 17.407414 | −0.1927% | 0.041413 | 0.897090 | Below stiffness threshold |

The two reductions are in **power**, not amplitude. The test holds relative drives fixed; it does not reoptimize them at lower power. These points bracket a pass and a failure for one design and one seeded target cloud, not a universal tolerance. No electrical power, robustness to arbitrary channel errors, nonlinear controlled motion or hardware performance was measured. The largest independently moved count qualified in the [original Phase 3 simulations](phase3_report.md) remains **two**.

## Cross-checks

- **Two fixed-drive controls.** The ±8- and ±16-harmonic rows in the [archived convergence results](results/convergence/sixteen_particle.json) were reproduced for stiffness, linear multiplier magnitude, power and gravity residual. The unchanged ±16 drive used 1.494590 W, above the 1.492952 W ceiling; the refined model was not in force balance at the old coordinates.
- **Two ways to enforce drive limits.** At ±24 harmonics and 96 disk integration points per emitter, scaling the archived projected slots together produced 17.679338 µN/m. Reapplying the [controller’s per-source projection](validation.py) to its saved raw variables produced 17.564687 µN/m at the same modeled power. The latter method was used for the power-margin tests; it is not a closed-loop run.
- **Spatial sensitivity.** At full power and ±40 harmonics, the 24-point disk calculation had only a +0.00064% stiffness margin; the 216-point calculation gave +0.832%. The coarser row uses 128 time samples per cycle and a 25 µm finite-difference step; the finest uses 256 and 12.5 µm. These are numerical sensitivities, not experimental error bars. Only like-for-like finest runs were used to compare power.

## Method

An independent NumPy/SciPy implementation rebuilt the frequency-dependent finite-disk field from the [archived design](results/experiments/compare_16_covariance/design.npz) and [case settings](results/experiments/compare_16_covariance/case.json). It reconstructed all archived ±8-harmonic emitted coefficients with maximum absolute discrepancy below 5 × 10⁻¹². The calculation includes first-Born particle coupling, gravity and all 48 position coordinates. The per-source peak bounds and quadratic modeled acoustic-power cap are applied before equilibrium solving. The ±40 harmonics span the configured finite bandwidth of the model, not all real acoustic frequencies.

The Newton solve rebuilds pressure derivatives at updated coordinates. The equilibrium is accepted only when the maximum force residual is below one millionth of a particle’s weight. The weakest stiffness is the smallest eigenvalue of the symmetric part of the cycle-averaged coupled stiffness matrix. The independently sampled phase average and frequency-domain mean must agree. Linear Floquet multipliers describe small disturbances about the mean equilibrium; the affine displacements are not closed-loop trajectories. See the [original model](physics.py) for the underlying assumptions.

## Next experiment

Retune the controller to pass the same stationary test at 99% of nominal modeled power. Then command particle 0 to move −0.25 mm in x and return while the other 15 targets stay fixed. Record moving and held errors, weakest stiffness, channel peaks, power and frame deadlines. One passing move would be a pilot, not completion of the signed, simultaneous and retention protocol.

The [existing degree-8 motion surrogate](dynamics.py) estimates a 1,670.625 MiB dense transfer array at ±40 harmonics, above its configured 665.6 MiB guard. These are calculated array sizes, not measured process memory. Stream or block the transfer and compare its results with direct-force calculations at the same disk resolution before relying on a trajectory.

## Reproduction

[Download the separate verification package](https://backend.composio.dev/api/v3/sl/4XOY5DfZL8) for the standalone implementation, inputs, raw arrays and eight run records. The inputs were checked against [commit `58d087b`](https://github.com/jihvijhojhviihogyuvi/ultrasonic-research-project/commit/58d087ba42cf12a4d15311de3202fa929c0a45bf). The package was tested with Python 3.13.18, NumPy 2.3.5 and SciPy 1.17.1. No network or GPU is needed once dependencies are installed. Run from the extracted package directory:

```sh
python verify.py
python run.py --out reproduced --only fixed-8 fixed-16 uniform-24-r2 native-24-r2 native-40-r1 native-40-r3 native-40-r3-P995 native-40-r3-P990
python verify.py --results reproduced
```

The runner refuses to overwrite an existing run. `summary.csv` covers the eight included runs. Each `results/<run>/` contains thresholds and a compressed array file with positions, forces, mean coupled stiffness, multipliers and affine displacements in SI units. `provenance.json` records source links and hashes; `SHA256SUMS` covers package files except the manifest itself.
