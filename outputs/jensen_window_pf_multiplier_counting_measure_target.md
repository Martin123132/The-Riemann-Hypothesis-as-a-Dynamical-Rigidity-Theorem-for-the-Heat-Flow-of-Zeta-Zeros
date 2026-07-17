# Jensen-Window PF Multiplier Counting-Measure Target

Date: 2026-07-10

Status: open theorem target. This is not a proof of PF-infinity, Jensen
hyperbolicity for zeta, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_multiplier_counting_measure_target`.

```text
work/rh_compute/results/jensen_window_pf_multiplier_counting_measure_target.json
python work/rh_compute/scripts/jensen_window_pf_multiplier_counting_measure_target.py
python work/rh_compute/scripts/check_jensen_window_pf_multiplier_counting_measure_target.py
```

Current result:

```text
validated Jensen-window PF multiplier counting-measure target: 10 rows, 0 issues, 4 exact rows, 1 finite evidence row, 3 countermodel rows, 1 live route, 0 ready-to-apply rows
```

## 2026-07-11 Unit-Atomic Obstruction

The interval certificate
`outputs/jensen_window_pf_multiplier_unit_atomic_obstruction_certificate.md`
rules out the sufficient theorem stated below for the actual lambda-zero zeta
coefficients. The order-6 magnitude bound forces every unit atom above
`4.863538496`, while the consecutive-moment ratio is strictly larger than the
maximum ratio available above that cutoff. The certified contradiction gap is
greater than `7.68e-4`.

The target is retained below as theorem-search history and as a precise record
of the rejected subclass. This does not rule out every multiplier sequence or
the other all-order Jensen/PF routes.

## 2026-07-11 High-Precision Frontier

The interval certificate
`outputs/jensen_window_pf_multiplier_complete_monotonicity_frontier_scout.md`
raises the necessary complete-monotonicity check from order 8 to the complete
finite triangle through order 55:

```text
(-1)^m*Delta^m[-log(x_k)]>0,
0<=m<=55, 1<=k<=56-m,
lambda in {0,1e-6,1e-4,1e-2,1e-1}.
```

All 7980 Arb intervals are strictly positive. Thus finite complete
monotonicity does not falsify the ansatz on the available frontier. The live
burden is now the stronger unit-atomic recovery or an interpolation-uniqueness
theorem that can connect the integer sequence to the certified continuous
Mellin obstruction.

## 2026-07-11 Hausdorff-Uniqueness Handoff

The exact bridge
`outputs/jensen_window_pf_multiplier_hausdorff_uniqueness_bridge.md` proves
that any integer-only product produces one uniquely determined finite
Hausdorff measure:

```text
y_k=integral_[0,1] r^(k-1)dnu(r),
dnu/dr=q(-log r)*sum_j r^alpha_j.
```

Consequently, if the unit-atomic representation exists, both `nu` and the
multiset `{alpha_j}` are unique. The remaining task is no longer an arbitrary
fit: recover this unique measure, prove the displayed density form with unit
integer multiplicities, or rule it out.

The same bridge gives an exact non-promotion guard. A periodic factor such as
`sin(2*pi*s)` vanishes at every integer and has zero centered second
difference, so integer and curvature agreement still do not identify the
natural Mellin interpolation. A separate growth or canonical-interpolation
theorem is required before importing the continuous power-sum obstruction.

## 2026-07-11 Conditional Leading-Atom Bounds

The interval certificate
`outputs/jensen_window_pf_multiplier_leading_atom_bound_certificate.md`
extracts the first rigorous consequences for any putative unit counting
measure. Writing

```text
a_m=(-1)^m*Delta^m y_1=sum_j g_m(alpha_j),
g_m(alpha)=integral_0^1 r^(alpha-1)*(1-r)^(m+2)/(-log r)dr,
```

the kernel is strictly decreasing in `alpha`. Order 6 gives

```text
4.863538496<beta_6<4.863538497,
alpha_min>=beta_6>4.863538496.
```

All other available orders `0..55` give weaker lower bounds. A separate
order-3 count inequality proves `N(11/2)<=1`. These conclusions remain
conditional on existence of the target measure; in particular, they do not
prove that an atom occurs below `11/2`.

## Sufficient Theorem

Normalize the actual coefficients by

```text
B_k=A_k*A_0^(k-1)/A_1^k, B_0=B_1=1.
```

For `alpha>0`, define

```text
M_k^(alpha)=(alpha/(alpha+1))^(k-1)*(alpha+k)/(alpha+1),
EGF(M^(alpha))=(1+z/(alpha+1))*exp(alpha*z/(alpha+1)).
```

The EGF has one negative zero and is Laguerre-Polya type I, so the
Pólya-Schur theorem makes `M^(alpha)` a multiplier sequence. Finite
pointwise products preserve real-rootedness by composition. Therefore the
following would be sufficient:

```text
B_k=product_j M_k^(alpha_j),
alpha_j>0 with unit integer multiplicity,
sum_j alpha_j^(-2)<infinity.
```

Coefficientwise limits then preserve every fixed Jensen polynomial. This is a
sufficient subclass theorem, not a claimed characterization of every multiplier
sequence.

## Ratio Form

The same target is

```text
-log x_k=sum_j -log(1-1/(k+alpha_j)^2).
```

The elementary kernel has the positive Laplace representation

```text
-log(1-1/(k+alpha)^2)
  =integral_0^infty exp(-(k+alpha)t)*2*(cosh(t)-1)/t dt.
```

Thus complete monotonicity of `y_k=-log x_k` is necessary. The Arb scout
certifies that sign pattern through order 8 on all five cached heat times,
but finite evidence does not construct a counting measure.

## Failure Gates

The atoms cannot simply be assigned arbitrary positive weights. The exact
fractional-power cubic guard has discriminant `<-27/125`, and the equal
positive mixture of the `u=3/5,2/3` boundary families has discriminant `-937/3456`.
A continuous positive fit is therefore not enough.

The separate Arb Mellin-interpolation certificate constructs the continuous
power-sum candidates from log moments of `Phi` and finds a strictly negative
shift-2 size-4 Hankel determinant. Therefore the natural Gamma-normalized
Mellin interpolation cannot be the elementary multiplier product. This does not
rule out equality asserted only for integer `k`; such a promotion would need a
separate interpolation-uniqueness theorem.

## Acceptance

A closing proof must construct the actual discrete multiset from Phi/Newman
data, prove convergence and equality for every k, and avoid endpoint
hyperbolicity or RH as an input. No such construction is currently known.

Primary theorem anchors:

```text
https://annals.math.princeton.edu/wp-content/uploads/annals-v170-n1-p14-p.pdf
https://arxiv.org/abs/math/0606360
```

```text
outputs/jensen_window_pf_rank_two_boundary_family_lemma.md
outputs/jensen_window_pf_defect_complete_monotonicity_scout.md
outputs/jensen_window_pf_mellin_multiplier_power_sum_obstruction.md
outputs/jensen_window_pf_bridge_target.md
```
