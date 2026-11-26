# Phase 3 computational validation contract

The initial supplied report is in the conversation; the complete follow-up report and original implementation/results are in outputs/ultrasonic_research. They are inputs, not accepted truth. Preserve them.

All new deliverables live in outputs/ultrasonic_phase3. Intermediate dependency packages, scripts, and logs live in work/phase3. No physical capacity or hardware accuracy claim follows from model-only results.

## Ordered work

1. Reproduce the prior two-, eight-, and sixteen-particle runs from inspected archived source snapshots in fresh runs; compare all metrics and save coefficients. Independently recalculate selected forces/Hessians and Floquet matrices using a separate numerical implementation. Investigate convergence and limitations.
2. Implement frequency-resolved finite-aperture pressure/velocity/derivatives, filtered temporal command spectra, propagation delays, beats, source/material configuration, consistent acoustic power/peak constraints, and CPU/GPU paths.
3. Implement particle interactions/scattering, full coupled stiffness and dynamics, nonlinear time integration, collision/escape handling, camera-like feedback, and signed/simultaneous motion tests with retention >=60 s for candidates passing prerequisites.
4. Implement mean-force and direct periodic-stability optimization; vary field counts 2,3,4,8,16, timing, unequal duties, ordering, transitions, common support, gauge/factorization, and control-Jacobian objectives. Implement an SDP bound/initializer on a tractable source count and incremental predictive feedback.
5. Execute matched comparisons, numerical/physics tests, timing/scheduling ablations, systematic source/particle/power/density scaling, retaining failed and infeasible attempts. Large runs may use explicit fail-fast physical prerequisites; do not call them completed independence tests or global capacity bounds.
6. Generate graphs only from saved runs. Audit resource comparability, units, convergence, provenance, and pass labels. Produce nine-part report, complete code/dependencies/raw tables/graphs, and a small falsifiable prototype with clearly labeled approximate budget.

## Acceptance inventory

- [ ] G1: both prior reports and all legacy code assumptions are reviewed; originals remain preserved.
- [ ] G2: Experiments A/B/C have fresh numerical runs, comparison tables, and independent selected physics/Floquet recalculations.
- [ ] G3: the new acoustic engine covers configurable geometry/count/material/frequency, finite aperture, complex pressure and velocity, standing/traveling/multiple fields, and transitions.
- [ ] G4: filtered frequency-dependent propagation, retarded envelopes, beats, peak drives, and consistent power are implemented and tested.
- [ ] G5: the particle engine includes gravity, drag, scattering/interactions, coupled stiffness/dynamics, collisions, escapes, and trajectories.
- [ ] G6: all required controllers, camera-like feedback, deterministic updates, and bounded incremental control are implemented.
- [ ] G7: periodic dynamics and Floquet stability enter optimization rather than only mean stiffness.
- [ ] G8: S=2,3,4,8,16, unequal timing/duties, ordering, smoothing, and common support are executed and compared.
- [ ] G9: differentiable optimization, covariance factorization, a tractable SDP experiment, control-Jacobian conditioning, and geometry/timing refinements have mathematical and executed tests.
- [ ] G10: positive/negative 0.25- and 1-mm moves, simultaneous motion, held/target tolerance, stiffness, >=60-s retention, and >=20-Hz updates are explicitly checked for candidates; fail-fast rejections remain distinct.
- [ ] G11: source counts 32,64,128,256,512 and particle counts 2,4,8,16,25,32,50 plus 100 where feasible are attempted under a declared fair protocol, with power/density/field-count/cost curves from actual runs.
- [ ] G12: complete runnable Python, dependencies, seeds, CPU fallback/GPU support, memory limits, automation, raw results, graphs, meaningful tests, convergence checks, and error handling are delivered.
- [ ] G13: scientific audit distinguishes analytic controls, model results, published benchmarks, predictions, and unverified hypotheses; limitations remain explicit.
- [ ] G14: the recommended architecture follows the evidence and its smallest physical test has components, approximate costs, measurements, and success/failure criteria.
- [ ] G15: the final nine-part report reconciles reproduction, framework, every executed result, temporal analysis, discoveries, scaling, architecture, prototype, and remaining research.

Scientific outcomes are manually reviewed against source, tests, and saved runs; no fixed expected performance number is its own acceptance oracle. A reproducible failure is useful evidence, not architecture success.
