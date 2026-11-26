# Gates: ultrasonic architecture research

OWNS: outputs/**, work/**

Scope: advance the supplied report with a source-checked critique, distinct architectures, governing equations, optimization, runnable comparative screening, and falsifiable hardware experiments.

- [x] G1: critique distinguishes errors, empirical limits, and unverified citations in the supplied report
  EVIDENCE: research_report.md section 1 distinguishes sign/dimension errors, qualified HAT results, the Wu spacing intervention, and unverifiable citations.
- [x] G2: original publications support the central experimental claims, with air and liquid results separated
  EVIDENCE: primary HAT full text, Wu full text, Puranen original manuscript, IEEE time-sharing abstract plus water-specific conference record, and original conference abstracts were inspected; citations are in sections 1-2.
- [x] G3: at least five technically distinct new architectures have mechanisms, constraints, and falsification tests
  EVIDENCE: section 3 specifies eight designs, their force mechanisms, limitations, and distinct falsification tests; no priority claim.
- [x] G4: the physics model handles gravity, positive-definite stiffness, shared-field interference, degrees of freedom, and interactions correctly
  EVIDENCE: section 4, selftest.py success, and scattering_check.py analytic control; retarded envelopes and finite-size/full scattering limitations are explicit.
- [x] G5: the optimization explicitly tests independent control, power, stiffness, separation, accuracy, and update rate
  EVIDENCE: section 5 specifies operational independence, lifted SDP, full dynamics, real-power and peak constraints, robust paths, and computational latency; runnable code states its subset.
- [x] G6: runnable screening covers holography, multiple frequencies, global-local support, and at least three new approaches with fair assumptions
  EVIDENCE: 37 saved cases include all eight implementations at common counts/targets and budgets internally; N1 static force, N2 covariance, and N3 support residual are executed.
- [x] G7: the implementation is scientifically reviewed against analytic controls, derivative checks, equilibrium checks, and reproducible run records
  EVIDENCE: physics controls passed, scattering asymptotic ratio 1.008187, limited convergence test, and package provenance audit passed; verification_summary.md records scope.
- [x] G8: numerical comparisons identify measured simulations, analytic calculations, and untested targets without claiming validated particle capacity
  EVIDENCE: section 7 derives both stiffness gains from JSON; no 512-source independence certificate; covariance instability and incomplete beat dynamics are reported.
- [x] G9: the recommended architecture has a falsifiable advantage and a software-first roadmap with success and failure criteria
  EVIDENCE: sections 8-9 recommend the simpler force-based controller, make covariance conditional, and define slow two-particle and later power-matched capacity experiments.
- [x] G10: the final deliverables cover all eleven required sections, working source links, unresolved questions, and specific next experiments
  EVIDENCE: package audit confirms eleven ordered sections; sections 10-11 list unresolved questions and specific experiments; README documents all runnable scripts and archived results.

These are manual scientific acceptance gates: command success alone cannot establish the corresponding physical or research conclusions. Script checks and run logs are supporting evidence, not automated certification of the scientific model.
