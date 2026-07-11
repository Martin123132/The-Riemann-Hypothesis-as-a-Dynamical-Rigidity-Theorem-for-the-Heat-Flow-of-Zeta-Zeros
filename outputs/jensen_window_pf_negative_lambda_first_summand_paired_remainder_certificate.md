# Jensen-Window PF Negative-Lambda First-Summand Paired-Remainder Certificate

Date: 2026-07-10

Status: interval compact paired-remainder and negative-skewness theorem with open asymptotic ray.
This is not a proof of the global cumulant wall, cone entry, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate`.

Machine-readable result:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.json
```

Generator and checker:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.py
```

## Exact Paired Form

At the log-variable mode put `y=sqrt(a_t)*(x-x_t)` and pair `y` with
`-y`. The four raw moments are integrated together, then normalized and
centralized before evaluating

```text
H_t=(t^2/a_t^(3/2))*(kappa_3(Y_t)+alpha_t+C3_t+C5_t).
```

The target is `H_t>=-79/1000`.

## Taylor And Tails

Arb certifies `42800` positive `V^(9)` intervals and
`40736` normalized eighth-derivative
envelope intervals, proving

```text
sup_|y|<=6 V^(8)(x_t+y/sqrt(a_t))/a_t^4<=1/50000.
```

On `|y|<=6`, a paired midpoint rule uses interval second derivatives.
Beyond each endpoint, the certified outward slope supplies an exact
exponential moment tail. The left-tail use is legitimate because `V'` is
increasing for `u>=1/100` and negative below that threshold.

## Compact Theorem

The adaptive cover proves

```text
H_t>=-79/1000 for every real mode parameter 0.9264<=u_t<=5.
kappa_3,t(2*log(U))<0 for every real mode parameter 0.9264<=u_t<=5.
```

- Base blocks: `4074`.
- Accepted adaptive blocks: `4074`.
- Maximum refinement depth: `0`.
- Minimum outward-rounded margin: `6.39581237362933276527177747885258985212094684222084936444645E-5`.
- Negative third-cumulant blocks: `4074`.
- Maximum outward-rounded third-cumulant upper bound: `-9.13150910388589240652708356293943489101299308784061615122009E-6`.
- Upper endpoint `t=V'(x(5))`: `[15241916613.768536123445468287824895052864753573752 +/- 3.15E-40]`.

## Remaining Ray

The only unproved row in this artifact is

```text
H_t>=-79/1000 for u_t>=5.
```

Here `q>=pi*exp(20)>10^9` and every normalized derivative is very small.
The next step is an analytic paired perturbation bound on this ray; finite
compact quadrature is not promoted across the unbounded parameter range.

```text
outputs/jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate.md
outputs/jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.md
outputs/jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

Paired interval quadrature, a global eighth-derivative Taylor envelope, and explicit two-sided tail budgets prove the seventh-order normalized remainder floor -79/1000 and strict negative third-cumulant sign for every real mode 0.9264<=u<=5. The asymptotic mode ray u>=5 remains open.
