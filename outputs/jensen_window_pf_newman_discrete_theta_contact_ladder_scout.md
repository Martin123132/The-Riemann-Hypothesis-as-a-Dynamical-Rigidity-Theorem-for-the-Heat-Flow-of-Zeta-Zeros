# Discrete Theta Contact Ladder Scout

Date: 2026-08-03

Status: finite adversarial diagnostic. This is not a proof of `Lambda<=0` or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_discrete_theta_contact_ladder_scout.json
python work/rh_compute/scripts/jensen_window_pf_newman_discrete_theta_contact_ladder_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_discrete_theta_contact_ladder_scout.py
```

## Division-Free Test

At multiplicity two, normalized positive component weights satisfy

```text
sum p_a h_a=0,  sum p_a h_a'=0,
c sum p_a h_a'''-6(ell+1/2)sum p_a h_a''=0.
```

Together with `sum p_a=1`, this is a four-equation linear witness problem.
The adjacent-jet ratio is never divided during the search.

## Ladder Results

| roster | targets realized | unbounded guard | sampled definite field span | max jet drift |
|:---|---:|:---:|:---|---:|
| continuum_log_mesh | 7/7 | True | [-28.2, 12.9] | 1.257e-12 |
| geometric_sqrt2 | 7/7 | True | [-196, 12.9] | 1.224e-12 |
| integer_1_8 | 5/7 | False | [-957, 482] | 1.225e-12 |
| integer_1_12 | 5/7 | False | [-1.06e+03, 482] | 1.225e-12 |

The target set is `-20,-5,-1,0,1,5,20`. Each realized value has an explicit
positive four-component witness in the JSON result.

## Integer Gap Refinement

- `ell=-1`: target witness `True`; nearest definite interval distance `0.000000e+00`; pole candidate `False` over `729` refined points.
- `ell=20`: target witness `True`; nearest definite interval distance `0.000000e+00`; pole candidate `True` over `1449` refined points.

These are local finite refinements. Failure to find a witness is not
a continuum exclusion theorem.

## Fixed Theta Rung

The smallest sampled normalized fixed-weight contact residual was `8.030338e-07` at `t=0.0, x=77.0`.

The nearest stored sparse free-weight contact has total-variation distance at most `1.184985e-10` from the retained fixed mixture. This is an upper bound from explicit contact vertices,
not the exact distance to the full contact polytope.

The result stores analytic envelopes for both the omitted `u>2` integral
and the fixed-weight integer tail `n>=9`. Quadrature agreement and finite
grid coverage remain diagnostic rather than interval-certified.

## Decision Gate

After local refinement, the free integer roster realizes all seven tested
fields. Thus the sampled evidence rejects the integer shift locations alone
as the missing rigidity. This does not prove realization of every real field.
The next local theorem must use the exact unit theta coefficients, or the
equivalent time-dependent normalized mass profile, jointly with the contact
equations and the infinite roster. A global zero-flow bypass remains separate.

## Proof Boundary

Positive point witnesses are high-accuracy quadrature diagnostics, not interval certificates. No fixed-weight contact exclusion, Xi field bound, Lambda<=0, RH, or prize-level conclusion is proved.
