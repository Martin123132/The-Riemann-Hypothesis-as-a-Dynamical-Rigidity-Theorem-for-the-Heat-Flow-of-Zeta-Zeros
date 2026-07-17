# Jensen-Window PF Multiplier Leading-Atom Bound Certificate

Date: 2026-07-11

Status: conditional interval leading-atom bounds; not a proof of the
multiplier factorization, PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_multiplier_leading_atom_bound_certificate`.

```text
work/rh_compute/results/jensen_window_pf_multiplier_leading_atom_bound_certificate.json
python work/rh_compute/scripts/jensen_window_pf_multiplier_leading_atom_bound_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_multiplier_leading_atom_bound_certificate.py
```

Current result:

```text
validated Jensen-window PF multiplier leading-atom bound certificate: 8 rows, 0 issues, 56 difference orders, beta6 in (4.863538496,4.863538497), alpha_min>4.863538496, N(11/2)<=1, 1 open existence handoff
```

## Conditional Kernel

Under the open unit-atomic multiplier target, define

```text
a_m=(-1)^m*Delta^m y_1,
g_m(alpha)=integral_0^1 r^(alpha-1)*(1-r)^(m+2)/(-log r)dr.
```

Then

```text
a_m=sum_j g_m(alpha_j),
partial_alpha g_m(alpha)=-Beta(alpha,m+3)<0.
```

If `beta_m` is the unique solution of `g_m(beta_m)=a_m`, monotonicity
gives `alpha_min>=beta_m`.

## Interval Bound

At order 6, 250-digit Arb arithmetic proves

```text
g_6(4.863538496)-a_6=[3.6340198025246959890220872597646469641634080905380338232387292267659619923440835845955394274199001689031528386553693e-14 +/- 5.90e-130]>0,
g_6(4.863538497)-a_6=[-1.44400690245213671439729637962228511397239566746416892150369167711981335757517927319511474915782715131436739220893555e-13 +/- 4.27e-130]<0.
```

Therefore

```text
4.863538496<beta_6<4.863538497,
alpha_min>4.863538496.
```

All other orders `0..55` are interval-certified weaker at the lower
endpoint, so order 6 gives the strongest bound in the available triangle.

## Conditional Count

For any cutoff `A`, positivity and monotonicity give

```text
N(A)*g_m(A)<=a_m.
```

At `A=11/2`, order 3 gives

```text
a_3/g_3(11/2)=[1.63178855134584329177740541990111930684216062759575212148207323650049948826802369155744924183606733607412660456855010332782239084 +/- 9.05e-129]<2,
N(11/2)<=1.
```

This does not prove that an atom lies below `11/2`. It says that if one
does, it is the unique atom in that interval. The subsequent certificate
`outputs/jensen_window_pf_multiplier_unit_atomic_obstruction_certificate.md`
composes the cutoff with a consecutive-moment ratio and rules out the measure
entirely; the conditional bounds here are retained as its first step.

```text
outputs/jensen_window_pf_multiplier_hausdorff_uniqueness_bridge.md
outputs/jensen_window_pf_multiplier_complete_monotonicity_frontier_scout.md
outputs/jensen_window_pf_multiplier_unit_atomic_obstruction_certificate.md
outputs/jensen_window_pf_multiplier_counting_measure_target.md
```
