# Jensen-Window PF Suzuki K2 Residual Scout

Date: 2026-07-23

Status: finite double-precision reconnaissance. This is not a proof
of an L2 tail, fixed-shift innerness, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_suzuki_k2_residual_scout.json
python work/rh_compute/scripts/jensen_window_pf_suzuki_k2_residual_scout.py
python work/rh_compute/scripts/jensen_window_pf_suzuki_k2_residual_scout.py --resume
python work/rh_compute/scripts/check_jensen_window_pf_suzuki_k2_residual_scout.py
```

## Target

For `t=log x`, the exact level-two residual is

```text
r_(omega,2)(t)
 =-q_(omega,1)
  +integral_0^t[H_(omega,1)(exp(u))-1]du,
q_(omega,1)=-2*xi'(1/2+omega)/xi(1/2+omega).
```

A cofinal family of genuine `L2(0,infinity)` estimates would be
RH-equivalent. This scout integrates only through `x=5000`.

## Finite Results

| omega | q1 | r2(10) | r2(100) | r2(1000) | r2(5000) | finite energy | coarse/fine rel. |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1/4 | -0.02310267 | 0.0024403532 | -0.00078215859 | -0.00067013994 | 9.6809901e-05 | 4.0777474e-05 | 0.00267 |

The fine grid is uniform in `t` with `2401` points; its even-index
subgrid provides the independent coarse quadrature comparison.

## Interpretation

These values test normalization and reveal finite-range shape.
They cannot establish convergence, square integrability,
positive-time Hardy support, uniformity as `omega->0`, or any
unbounded arithmetic estimate. Boundary spectral energy is
automatic and is not what this time-domain scout certifies.

## Runtime

Mode: `daytime_one_worker`; one below-normal
worker, numerical thread pools limited to one.
Completed shifts: `1/4`; parked:
`true`.

## Proof Boundary

Finite double-precision time-domain quadrature only. It does not prove convergence as x tends to infinity, an L2 residual, Hardy causality, fixed-shift innerness, RH, or Lambda<=0.
