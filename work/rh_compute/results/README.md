# Result artifact layout

This directory contains reproducible finite certificates, manifests, and
diagnostic caches used by the RH research programme.

High-cardinality historical families are grouped so GitHub can display every
directory without truncating its file list:

- `arb_edrei/log_sign/`: finite Edrei logarithmic-sign diagnostics;
- `arb_edrei/power_hankel/lam*/`: Edrei power-Hankel diagnostics grouped by
  lambda label;
- `arb_toeplitz/`: Arb Toeplitz finite-minor certificates and logs.

The remaining result families stay at this level to preserve the established
paths used by current certificate builders and checkers.

## Local-only large checkpoints

Raw resumable checkpoints at or above 100,000,000 bytes are kept out of new
Git objects. The repository carries their compact run contracts, pilots,
progress reports, and completed certificates where available. When an
incomplete interval claim still depends on the full prefix, a lossless gzip
archive and an archive-backed checker are published instead.

The current local-only raw checkpoints are:

- `jensen_window_pf_compound_order11_compact_adaptive_h23_segments.jsonl`;
- `jensen_window_pf_compound_order12_sparse_h23_lower_bridge_segments.jsonl`.

The second checkpoint is published losslessly as
`jensen_window_pf_compound_order12_sparse_h23_lower_bridge_segments.jsonl.gz`.
The first has a completed compact certificate and remains local rather than
placing a multi-gigabyte cache in repository history.

Earlier large caches already present in repository history remain at their
published revisions when the live resumable file has subsequently grown past
that threshold.
