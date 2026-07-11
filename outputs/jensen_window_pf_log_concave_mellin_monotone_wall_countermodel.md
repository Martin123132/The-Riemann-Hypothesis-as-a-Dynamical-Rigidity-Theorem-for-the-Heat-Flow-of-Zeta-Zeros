# Jensen-Window PF Log-Concave Mellin Monotone-Wall Countermodel

Date: 2026-07-10

Status: exact interval countermodel gate. This is not a proof or disproof
of zeta monotone contractions, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_log_concave_mellin_monotone_wall_countermodel`.

Machine-readable result:

```text
work/rh_compute/results/jensen_window_pf_log_concave_mellin_monotone_wall_countermodel.json
```

Generator and checker:

```text
python work/rh_compute/scripts/jensen_window_pf_log_concave_mellin_monotone_wall_countermodel.py
python work/rh_compute/scripts/check_jensen_window_pf_log_concave_mellin_monotone_wall_countermodel.py
```

Current result:

```text
validated Jensen-window PF log-concave Mellin monotone-wall countermodel: 6 rows, 0 issues, 2 upper-wall contractions, 1 monotone-wall violation
```

## Exact Witness

Take the log-concave density

```text
f(y)=exp(-5*y)*1_[0,1](y).
```

Its Gamma-normalized Mellin transform is

```text
H(p)=integral_0^1 y^(p-1)*exp(-5*y)dy/Gamma(p)
    =gamma_lower(p,5)/(5^p*Gamma(p)).
x_k=H(k+3/2)*H(k-1/2)/H(k+1/2)^2.
```

Arb certifies:

- `x_1` lies between `9.585809898433196985470820434295664830224263410373725967025382207166572E-1` and `9.585809898433196985470820434295664830224263410373725967025382207166574E-1`.
- `x_2` lies between `9.312101606377941403101410236154778849453880349801908190968314813658825E-1` and `9.312101606377941403101410236154778849453880349801908190968314813658827E-1`.
- `x_2-x_1 < -2.737082920552555823694101981408859807703830605718177760570673935077461E-2`.

Both contractions satisfy the Berwald-Borell upper wall `0<x_k<1`,
but `x_2<x_1`. Thus log-concavity of the measure, even with compact
support and a log-affine density, does not imply `x_(k+1)>=x_k`.

## Proof Boundary

This witness blocks only a generic theorem. It does not touch the strict
zeta-kernel log-concavity certificate and does not rule out a zeta-specific
ratio recurrence, determinant identity, or uniform saddle theorem.

```text
outputs/jensen_window_pf_kernel_mellin_upper_wall_certificate.md
outputs/jensen_window_pf_monotone_contraction_theorem_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

The compactly supported log-concave density exp(-5y) on [0,1] has Gamma-normalized Mellin contractions x_1 and x_2 inside the Berwald-Borell upper wall but x_2<x_1. Therefore the remaining adjacent-k wall is not a generic consequence of log-concavity and needs zeta-specific structure.
