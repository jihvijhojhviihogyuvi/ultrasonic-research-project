# Verification record

The research deliverable covers all eleven requested sections and was reviewed against ten scientific acceptance gates: **10 met, 0 unmet, 0 abandoned**. These are deliverable checks, not certification of a physical display or an optimal architecture.

The saved executions comprise **37 architecture cases**. Five receive the limited `tier1_local_pass` label: three two-particle cases at 128 sources, and two two-particle cases at 32 sources under the separately specified 5-Hz prototype requirement. Seventeen cases did not run finite independence probes. No case is described as a validated high-density hardware result.

Checks completed:

- Analytic pressure-gradient and standing-wave stiffness controls, including gravity sign.
- A known saddle confirming that positive Laplacian alone does not establish stability.
- Coherent interference versus incoherent slot addition.
- Drive-gradient finite differences, away from a projection kink.
- Periodic envelope filtering, including a case where fast opposite-phase commands are strongly suppressed.
- Peak and total drive-effort constraints for all eight implementations.
- Source propagation batching consistency.
- Leading Rayleigh scattering checked against near-node r^-4 attraction.
- Fixed-field stencil and source disk-quadrature refinement.
- Artifact provenance: source-snapshot SHA-256, identical target clouds within each comparison, saved command power/peak limits, and honest pass labels.
- Reported 43.2% and 87.1% stiffness gains recalculated directly from saved results.
- The local trajectory test reports its actual offline optimizer latency and marks the requested 50-Hz real-time rate as unmet by that implementation.

Success output from the current scripts:

```
PHYSICS_CONTROLS_PASSED: derivatives, standing wave, gravity, saddle, interference, envelope, power, rank
SCATTERING_CONTROL_PASSED: near-node attraction / asymptotic = 1.008187
PACKAGE_AUDIT_PASSED: 37 cases; 5 local passes; 17 without finite independence probes; 11 sections; gain percentages [43.2100948736099, 87.11686261322737]
```

Material limitations are explicit in the report and README: assumed source/particle properties, nominal power conversion, limited finite-motion directions, leading rather than full scattering, missing streaming/nonlinear propagation/enclosure physics/sensing latency, and omitted spatially retarded temporal envelopes and spectral beat dynamics. The model does not establish maximum particle capacity, an optimal configuration, or a physical accuracy specification.
