# Jensen-Window PF Negative-Lambda First-Summand Paired-Remainder Ray Certificate

Date: 2026-07-10

Status: analytic asymptotic-ray theorem and global paired-remainder closure.
This is not a proof of the remaining all-order bridge, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate`.

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate.json
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate.py
```

## Ray Geometry

For `u>=5`, put `q=pi*exp(4*u)`, `a=V''(x_t)`, and
`Y=sqrt(8*log(q))`. Exact component bounds give

```text
|alpha|<=2/sqrt(q),
sup_|y|<=Y |V^(4)(x_t+y/sqrt(a))|/a^2<=3/q,
P=t^2/a^(3/2)<=6*sqrt(q)/(5*u).
```

At the worst endpoint, the central perturbation is below `2.25043943534960228992130917200459785961038038747658954735372E-2`.
The adaptive window has outward slopes at least `9*Y/10`, and both
actual and Gaussian-linear model tails contribute at most `q^(-2)`
to each normalized raw-moment error.

## Moment And Cumulant Bounds

The first-order Gaussian model has moments

```text
p0=1, p1=-alpha/2, p2=1, p3=-5*alpha/2.
```

The exact comparison proves

```text
|m1+alpha/2|<=13/q,
|m2-1|<=36/q,
|m3+5*alpha/2|<=80/q,
|kappa_3(Y)+alpha|<=120/q.
```

After the prefactor, this is below `1/1000`. Adding the already
certified cubic and fifth-order caps gives

```text
H_t>=-299/25000>-3/250>-79/1000, u>=5.
```

## Global Closure

The compact paired theorem covers `0.9264<=u<=5`; the mode at `t=318`
lies to the right of `0.9264`. Therefore

```text
H_t>=-79/1000 for every real t>=318,
kappa_3,t(2*log(U))>=-37/(50*t^2) for every real t>=318.
```

This closes the analytic hypothesis in the exact first-summand
cumulant/Gamma bridge. It does not by itself prove the remaining
cone-propagation or all-order Newman-direction bridge.

```text
outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.md
outputs/jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate.md
outputs/jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.md
outputs/jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```
