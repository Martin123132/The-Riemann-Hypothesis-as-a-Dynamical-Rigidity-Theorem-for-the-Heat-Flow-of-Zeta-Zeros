# Jensen-Window PF Multiplier Unit-Atomic Obstruction Certificate

Date: 2026-07-11

Status: interval-certified unit-atomic obstruction; not a proof of
PF-infinity, the all-order Jensen bridge, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_multiplier_unit_atomic_obstruction_certificate`.

```text
work/rh_compute/results/jensen_window_pf_multiplier_unit_atomic_obstruction_certificate.json
python work/rh_compute/scripts/jensen_window_pf_multiplier_unit_atomic_obstruction_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_multiplier_unit_atomic_obstruction_certificate.py
```

Current result:

```text
validated Jensen-window PF multiplier unit-atomic obstruction certificate: 8 rows, 0 issues, 1 atom cutoff, 1 ratio cap, 1 unit-atomic route ruled out, 0 open requirements
```

## Target Consequences

If the proposed unit counting multiset existed, then

```text
a_m=(-1)^m*Delta^m y_1=sum_j g_m(alpha_j),
g_m(alpha)=integral_0^1 r^(alpha-1)*(1-r)^(m+2)/(-log r)dr.
```

The kernel is positive and strictly decreasing in `alpha`. Arb proves

```text
g_6(4.863538496)-a_6=[3.6340198025246959890220872597646469641634080905380338232387292267659619923440835845955394274199001689031528386553693e-14 +/- 5.90e-130]>0.
```

Because every target atom has unit weight, this forces `alpha_j>4.863538496`
for every `j`.

## Ratio Contradiction

Put `R_m(alpha)=g_(m+1)(alpha)/g_m(alpha)`. Under the probability
density proportional to the integrand defining `g_m`,

```text
R_m(alpha)=E_alpha[1-r],
R_m'(alpha)=Cov_alpha(1-r,log r)<0.
```

The covariance is strictly negative because `1-r` decreases and `log r`
increases. Therefore the target would require

```text
a_7/a_6<R_6(4.863538496).
```

The actual intervals give

```text
a_7/a_6=[0.60329213501285476539895444721995922153553918551423103414127248752101292774130285953284935881402312478229139424891605140394779 +/- 4.76e-126],
R_6(4.863538496)=[0.6025238849153418649240749307941747635646417680649518410079275535058583119279131726463225567787733712573160877741236982418582741262310526346166283516150332751712563336295821838120906942619746924833555077744248045915567093599748222065622108860579 +/- 4.06e-245],
gap=[0.00076825009751290047487951642578445797089741744927919313334493401515461581338968688652680203524975352497530647479235316208951 +/- 8.00e-126]>0.
```

This is incompatible with the required weighted-average cap. Hence the
normalized zeta coefficients do not admit the proposed convergent
unit-atomic elementary multiplier product.

## Boundary

The contradiction uses unit multiplicity in the cutoff step. It does not rule
out arbitrary positive subunit weights, but such fractional weights were
already insufficient for multiplier preservation. More general multiplier
sequences and the other all-order Jensen/PF routes remain open.

```text
outputs/jensen_window_pf_multiplier_leading_atom_bound_certificate.md
outputs/jensen_window_pf_multiplier_hausdorff_uniqueness_bridge.md
outputs/jensen_window_pf_multiplier_counting_measure_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```
