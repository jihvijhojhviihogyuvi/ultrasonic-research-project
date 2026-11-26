# Source and execution provenance

Archived Phase-2 snapshots are preserved byte for byte under `inputs/legacy`;
their original SHA-256 values are checked by `reproduce.py`.

The Phase-3 pre-review core files were reconstructed by reversing two inspected
patches, then verified against the hashes recorded in the original executions:

- `physics.py`: `92a0e956e4bff2931aa34265435e683b5abf64d761a97c352ab85102cf868785`
- `optimization.py`: `872b178b0de493e592e2155758fea668dc08551030f51ae7e5afdcb91a644c9a`

These exact versions are in `source_versions/reconstructed_pre_review_*`.
The final physics adds configuration validation without changing numerical
equations. The final nonlinear optimizer replaces velocity norm at zero with a
negligibly regularized norm to make its second derivative finite. The original
failed nonlinear pilot and corrected pilot are both retained.

The trajectory harness evolved during testing. Its early power bound used a
global peak for every source and was overly conservative; the corrected bound
weights each actual source amplitude. Early timing included offline stiffness
and acoustic-polynomial oracles. The emitted-drive controller is timed
separately in the final harness. Old timing failures remain in the raw records.
Controller replay verifies every saved control coordinate before combining its
timing with an already completed nonlinear trajectory; it does not fabricate
another trajectory. The first faster-ramp pilot failed the accuracy contract;
the final two-second 1-mm ramps are explicitly identified in the report.

Logs preserve the rejected degree-6 surrogate, failed early harness runs and
successful/unsuccessful final cases. Every new architecture/scaling attempt has
its full configuration, seed, fields, optimizer history and source hashes.
Numerical model changes are not silently represented as the same experiment.
