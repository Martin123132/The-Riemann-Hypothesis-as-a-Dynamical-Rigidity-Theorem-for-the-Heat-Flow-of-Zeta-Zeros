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
