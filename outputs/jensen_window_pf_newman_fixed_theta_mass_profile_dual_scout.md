# Fixed Theta Mass-Profile Dual Scout

Date: 2026-08-03

Status: finite primal-dual diagnostic. This is not a proof of contact exclusion or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_fixed_theta_mass_profile_dual_scout.json
python work/rh_compute/scripts/jensen_window_pf_newman_fixed_theta_mass_profile_dual_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_fixed_theta_mass_profile_dual_scout.py
```

## Fixed-Coefficient Question

The true unit-coefficient theta sum induces the normalized mass profile

```text
p_(n,t)=M_(n,t)/sum_k M_(k,t).
```

At each sampled `(t,x)` this scout solves the finite-dimensional
total-variation projection from `p` to the positive contact slice and the
closed-form affine Fisher projection. The latter returns a two-coordinate
dual observable that can be tested for an arithmetic pattern in `n`.
The smallest double-precision objective is below the conservative
coarse/fine quadrature resolution gate, so the table is a locator rather
than a certified distance report; the independent checker resolves the
selected face at 80 decimal places.

## Depth Stability

| depth | locator TV objective | t | x | normalized omitted-mass bound | TV/tail | Fisher d^2 at TV point |
|---:|---:|---:|---:|---:|---:|---:|
| 4 | 1.289236e-14 | 0.2 | 133.675 | 1.010e-30 | 1.277e+16 | 9.474159e-19 |
| 8 | 1.233822e-14 | 0.2 | 133.675 | 1.270e-106 | 9.715e+91 | 9.474068e-19 |
| 12 | 1.232082e-14 | 0.2 | 133.675 | 2.274e-226 | 5.419e+211 | 9.474119e-19 |
| 12 + infinity endpoint | 1.231660e-14 | 0.2 | 133.675 | 2.274e-226 | 5.417e+211 | n/a |

## Dual Audit

At the deepest nearest point, the TV primal-dual gap is `1.578e-30` and the scaled contact residual is `1.110e-16`.

The Fisher feature dual in original `(h,h')` coordinates is

```text
[0.0006624225143609716, 0.04717828781428395]
```

Its affine projection has minimum weight `1.229e-193`; positivity-active is `False`.

## Decision Gate

A depth-stable numerical gap is not itself a theorem. A useful next step
requires the dual direction to simplify into an identity or inequality
that explicitly uses the unit theta coefficients and controls every
omitted component. If the direction wanders with `(t,x)`, this local
projection route should be deprioritized in favor of the global flow and
degree-uniform Jensen branches.

## Proof Boundary

The LPs and Fisher projections are finite high-accuracy diagnostics. The selected LP objective lies below the conservative quadrature resolution gate and is a locator, not a certified distance. The stored tail envelopes are analytic, but this scout does not interval-certify the quadrature, optimize over continuous (t,x), or turn a sampled dual into an infinite arithmetic separator.
