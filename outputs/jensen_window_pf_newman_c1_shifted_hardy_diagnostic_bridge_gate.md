# Newman C1 Shifted-Hardy Diagnostic Bridge Gate

Date: 2026-08-05

Status: pinned external evaluator, low-height calibration, and physical kernel interpolation validated; physical carrier values and proof error bounds remain open; not a proof of RH.

## Structural decision

The exact ridge identity 2 M_N-F_N=M_N-tau_band+E_band leaves one raw carrier. The diagonal-plus-adjacent collar is locally exact but does not by itself give an O(polylog N) physical evaluator.

The alternate diagnostic bridge uses

```text
If K_j(lambda)=2 exp[-i theta_RS(T+delta_j)]exp(i delta_j lambda) and g(lambda) is fitted as sum_j c_jK_j(lambda) with real c_j, then sum_j c_j H_main(T+delta_j)=Re sum_(n<=N)n^(-1/2+iT)g(log n). The half-endpoint carrier convention is restored by explicit endpoint subtraction.
```

The common phase is essential. Without the physical normalizer, the stationary endpoint loses one quadrature. The retained mismatch is

```text
At the first physical root omega_c/exp(-i theta_RS(T_0))=-1-4.83903158037279e-17i, so the retained perpendicular endpoint fraction is 4.83903158037e-17.
```

## External calibration

- Upstream: `https://github.com/dml2391/Hardy-function-fastcodes`.
- Commit: `2e16dac3206b707052c3ac4cacdf3d1a2325e636`.
- License: `GPL-3.0`.
- Single T=10^10 absolute error: `0.0018623960398`.
- Accepted 0.01-grid maximum main error at T=10^10: `0.00724738923861`.
- Accepted 0.01-grid maximum main error at T=10^12: `0.00288990976752`.
- Rejected 0.02-grid maximum main error: `1.43906249489`.
- Rejected 0.04-grid maximum main error: `4.0876351952`.

The wider grids are falsification controls. They are not admissible physical samplers.

## Physical kernel

No terminal band peel is required. Across `3` roots and `24` base/tangent rows, including both interval endpoints, the maximum held-out relative error is `3.80043545001e-11`. The maximum coefficient l1 norm is `2311004.4995`.

The largest physical coefficient l1 norm is 2311004.4995. A worst-case independent sample-error estimate is therefore amplified by that factor; only empirical low-height correlation currently reduces it.

## Runtime boundary

The upstream multi evaluator has no restart facility. Measured one-rank runtime rises from about 9 seconds at T=10^10 to 153 seconds at T=10^12, so an uncheckpointed physical T about 3.26e22 run cannot be assumed to finish inside one four-hour cycle even with four night workers.

## Proof boundary

The external algorithm documents only an asymptotic relative-error order, with no validated constant or interval enclosure here. Kernel interpolation, raw ten-decimal outputs, and correlated low-height errors prove no physical observation, determinant sign, flow inequality, Lambda<=0, PF-infinity, RH, or prize-level conclusion.

Add append-only block checkpoints and deterministic resume to the accepted shift_step=0.01 multi evaluator, preserve each of the 15 partial sums and the RS-tail state, and validate stop/resume equivalence at T=10^10 and T=10^12. Then formulate a joint correlated-error diagnostic; do not launch the physical job or call it proof evidence before both gates close.

built Newman C1 shifted-Hardy diagnostic bridge gate: 15 rows, 0 issues, 4 portability patches, 6 raw fixtures, 5 direct-main calibrations, 3 low-height reconstructions, 24 physical kernel fits, shift_step 0.01 retained, 0.02 and 0.04 rejected, 0 physical values and 1 open resumable-evaluator row
