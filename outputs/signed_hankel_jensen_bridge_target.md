# Signed-Hankel/Jensen Bridge Target

Date: 2026-07-16

Status: theorem target. This is not a proof of RH or `Lambda <= 0`; it is the
surviving Xi/Phi-specific Jensen question after the proposed all-order
signed-Hankel antecedent was rigorously rejected at order ten.

## Purpose

The rigorous results support a signed-Hankel/sign-consistency pattern through
order nine for the exponential/Jensen coefficients:

```text
A_k(lambda) = mu_{2k}(lambda) k! / (2k)!
```

The pattern does not continue at all shifts: four required order-ten endpoint
determinants are strictly negative. The missing proof step is therefore not
another finite grid or completion of that false hierarchy. It is a weaker
Xi/Phi-specific theorem that gives all-degree/all-shift Jensen hyperbolicity,
Laguerre-Polya membership, or directly excludes positive Newman birth.

## Objects

For each `lambda`, define:

```text
P_{d,n,lambda}(x)
  = sum_{j=0}^d binom(d,j) A_{n+j}(lambda) x^j
```

The Jensen target is:

```text
for every d >= 1 and n >= 0,
P_{d,n,0}(x) has only real nonpositive zeros.
```

Equivalently, `Q_{d,n,0}(y)=P_{d,n,0}(-y)` has `d` positive real zeros.

The total-positivity reformulation is recorded in:

```text
outputs/jensen_window_pf_bridge_target.md
outputs/jensen_window_pf_obligation_algebra.md
outputs/arb_jensen_window_pf_obligation_diagnostic.md
outputs/arb_jensen_window_sturm_hyperbolicity_diagnostic.md
python work/rh_compute/scripts/check_jensen_window_pf_bridge_target.py
python work/rh_compute/scripts/check_jensen_window_pf_obligation_algebra.py
python work/rh_compute/scripts/check_arb_jensen_window_pf_obligation_manifest.py
python work/rh_compute/scripts/check_arb_jensen_window_sturm_manifest.py
```

It asks for every binomially weighted Jensen window:

```text
B^{d,n,0}_j = binom(d,j) A_{n+j}(0)
```

to be a finite PF-infinity sequence. This is an equivalent Jensen-window
target, not a proof that the target holds.

The obligation algebra gate records that degree 2 is the exact signed-Hankel
contact point, while degree 3 and degree 4 introduce additional
Jensen-window Toeplitz obligations and exact low-order countermodel failures.
The Arb diagnostic validates `1470/1470` selected degree-3/4
Jensen-window interval determinants across the five-lambda grid and shifts
`n=0..20`; this remains finite evidence only.
The companion Arb/Sturm diagnostic validates `210/210` degree-3/4 and
`105/105` degree-5 positive-root counts for `Q_{d,n,lambda}(y)=P_{d,n,lambda}(-y)`
on the same lambda and shift grid; this also remains finite evidence only.

For the reshaped-Hankel sign-consistency route, define the shifted row block:

```text
R_{k,n}(j_1,...,j_k)
  = det(A_{n+i+j_l})_{i=0..k-1, l=1..k}
```

where:

```text
0 <= j_1 < ... < j_k.
```

The all-order sign-consistency hypothesis would be:

```text
(-1)^(k(k-1)/2) R_{k,n}(j_1,...,j_k) > 0
```

for every `k >= 1`, `n >= 0`, and every strictly increasing column set, with
the correct weak/zero clauses if structural degeneracies are proved.

The first Arb certificate verifies only a finite unshifted frontier:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
n = 0
k = 2..7
N = 20 columns
689,795/689,795 finite minors positive
```

The shifted Arb certificate extends the finite evidence to a bounded shift
range:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
n = 0..20
k = 2..5
N = 18 columns
1,322,685/1,322,685 finite minors positive
```

The order-6 shifted Arb certificate adds the next-order finite slice:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
n = 0..20
k = 6
N = 16 columns
840,840/840,840 finite minors positive
```

The order-7 shifted Arb certificate adds another bounded next-order slice:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
n = 0..20
k = 7
N = 15 columns
675,675/675,675 finite minors positive
```

The order-8 shifted Arb certificate adds the next bounded slice:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
n = 0..20
k = 8
N = 14 columns
315,315/315,315 finite minors positive
```

The consolidated shifted staircase checker validates:

```text
3,154,515/3,154,515 finite shifted minors positive
```

All certificates remain finite. They are evidence for the theorem target, not
the all-shift/all-order theorem itself.

## Candidate Theorem B-Star (Actual Antecedent Rejected)

Prove a theorem of the following shape:

```text
Assume:
  1. A_k(lambda) are the Xi/Phi heat-flow coefficients, not an arbitrary
     positive sequence.
  2. A_k(lambda) satisfy the all-order shifted reshaped-Hankel
     sign-consistency condition above.
  3. The associated entire functions have the required convergence, parity,
     growth, and heat-flow normalization.

Then:
  for lambda = 0, all Jensen polynomials P_{d,n,0} are hyperbolic;
  equivalently H_0 belongs to the required Laguerre-Polya class;
  hence Lambda <= 0.
```

Stronger acceptable version:

```text
the same conclusion holds uniformly for lambda in [0, epsilon]
```

or:

```text
the sign-consistency structure excludes positive Newman boundary birth
directly.
```

As an abstract implication this remains a legitimate theorem question, but it
cannot be applied to the actual endpoint sequence. Assumption 2 is false:
rigorous order-ten determinants give `Q_(10,n)(-100)<0` for `n=0,1,2,3`.
Any theorem intended to close the Xi/Phi problem must replace that assumption
with a strictly weaker condition actually satisfied by the coefficients.

## Exact Low-Degree Gates

Degree 1 is automatic from positivity.

Degree 2 is exact:

```text
P_{2,n}(x) = A_n + 2 A_{n+1} x + A_{n+2} x^2
Delta_2 = 4(A_{n+1}^2 - A_n A_{n+2})
Delta_2 = -4 det [[A_n, A_{n+1}], [A_{n+1}, A_{n+2}]]
```

Thus the `m = 1` signed-Hankel condition is exactly degree-2 Jensen
hyperbolicity for positive coefficients.

The reciprocal-gamma source of this sign has now been isolated exactly:

```text
outputs/jensen_window_pf_reciprocal_gamma_mixture_sign_gate.md
work/rh_compute/results/jensen_window_pf_reciprocal_gamma_mixture_sign_gate.json
python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_gamma_mixture_sign_gate.py
```

Karlin's 1964 reciprocal-gamma theorem proves that every fixed-scale kernel

```text
gamma_(n+i+j)*t^(n+i+j), gamma_k=k!/(2k)!, t>0,
```

has the required strict all-order signature for arbitrary row and column
sets. The Xi coefficients are positive scale mixtures of these kernels, but
the exact row-wise determinant integral has independent scale variables and a
sign-changing integrand. A positive three-atom measure on strictly positive
scales already reverses the required order-two sign, so positive mixing is not
a closure theorem.

For the actual Xi measure, the established all-shift order-two cone is exactly
the tilted concentration estimate

```text
CV_n^2<=2/(2n+1).
```

Thus the all-order antecedent has been sharpened: prove a higher compound
concentration or positive coupling theorem for the actual Xi mixing measure.
Fixed-scale reciprocal-gamma sign regularity by itself is no longer an open
step and cannot be promoted through mixing.

The first higher compound coordinate is now exact:

```text
outputs/jensen_window_pf_reciprocal_defect_compound_order3_gate.md
work/rh_compute/results/jensen_window_pf_reciprocal_defect_compound_order3_gate.json
python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_defect_compound_order3_gate.py
```

For the contiguous `3` by `3` minor, set `d_k=1-x_k` and
`q_k=d_k^(-1/2)`. Its required negative sign is equivalent to

```text
C_n=q_(n+1)*q_(n+3)-q_(n+2)^2+1>0.
```

This is unit-buffered log-concavity of the reciprocal defects. A strict
rational countermodel satisfies the full ratio cone, decreasing defect,
increasing scaled defect, reciprocal-defect increments below one, and both
neighboring cubic Jensen inequalities, yet has `C_n<0` and the wrong positive
order-three Hankel sign. Therefore the new curvature margin is not a
consequence of any currently proved scalar or cubic cone.

The actual lambda=-100 entry problem for this contiguous family is now closed:

```text
outputs/jensen_window_pf_negative_lambda_m100_compound_order3_entry_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_m100_compound_order3_entry_certificate.json
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_m100_compound_order3_entry_certificate.py
```

A repaired Arb prefix proves `C_n>0` through `n=317`. The exact all-index
scaled-defect theorem, together with the certified anchor `s_319>251/500`,
gives the uniform analytic tail bound

```text
C_n(-100)>57613471/66107054971>0, n>=318.
```

Hence every shifted contiguous `3` by `3` signed-Hankel minor has the required
negative sign at `lambda=-100`. Its forward propagation is now also closed:

```text
outputs/jensen_window_pf_compound_order3_forward_invariance_certificate.md
work/rh_compute/results/jensen_window_pf_compound_order3_forward_invariance_certificate.json
python work/rh_compute/scripts/check_jensen_window_pf_compound_order3_forward_invariance_certificate.py
```

The exact heat derivative factors into the cooperative chain

```text
C_n'/r_(n+2)=alpha_n*C_(n+1)+beta_n*C_n, alpha_n>0.
```

The proved cubic forward-uniform tail makes a weighted infinite maximum
principle legitimate. Therefore

```text
C_n(lambda)>0 for every n>=0 and finite lambda>=-100,
D_(3,n)(0)<0 for every n>=0.
```

Noncontiguous order-three minors now follow from exact secant geometry:

```text
outputs/jensen_window_pf_order3_noncontiguous_secant_transfer_lemma.md
work/rh_compute/results/jensen_window_pf_order3_noncontiguous_secant_transfer_lemma.json
python work/rh_compute/scripts/check_jensen_window_pf_order3_noncontiguous_secant_transfer_lemma.py
```

After positive column scaling, the three-row Hankel columns are planar points
with strictly decreasing abscissas. Contiguous negativity is exactly strict
decrease of their edge slopes; every longer secant is a positive weighted
average of those slopes. Consequently

```text
R_(3,n)(j_1,j_2,j_3)<0
for every n>=0 and 0<=j_1<j_2<j_3 at lambda=0.
```

The arbitrary-column order-two sign follows from the same strict ratio order.
Thus reshaped-Hankel orders two and three are complete. Compound order four is
the first genuinely new signed-Hankel layer. Its contiguous part is now
complete, and the downstream finite-block theorem below also closes its
arbitrary-column part.

Its contiguous coordinate is now exact:

```text
outputs/jensen_window_pf_compound_order4_condensation_gate.md
work/rh_compute/results/jensen_window_pf_compound_order4_condensation_gate.json
python work/rh_compute/scripts/check_jensen_window_pf_compound_order4_condensation_gate.py
```

Desnanot-Jacobi condensation shows that the required positive `4` by `4`
sign is strict log-concavity of `T_n=-H_(3,n)>0`. In defect-gap coordinates,
the all-shift target is

```text
G_(n+1)^2>x_(n+3)^3*G_n*G_(n+2).
```

The repaired lambda=-100 source certifies all 317 available margins through
`n=316`, including the stronger finite cap `P_n<2/(n+3)^2` for the gap
log-curvature penalty. The exact sufficient tail handoff is

```text
P_n<=4/(n+3)^2, n>=317.
```

The scaled-defect anchor makes this `O(1/n^2)` cost strictly smaller than the
available `-3 log x_(n+3)=O(1/n)` buffer. This tail bound is now reduced by

```text
outputs/jensen_window_pf_compound_order4_first_summand_curvature_bridge.md
work/rh_compute/results/jensen_window_pf_compound_order4_first_summand_curvature_bridge.json
python work/rh_compute/scripts/check_jensen_window_pf_compound_order4_first_summand_curvature_bridge.py
```

to the single continuous first-summand ceiling

```text
K_1(t)=(log g_1(t))''<=7/(2t^2), t>=319.
```

The bridge proves `J_1(t)>=1/(7t)` for every real `t>=319`, gives
`P_n^(1)<=18/(5k^2)` conditionally on the displayed ceiling, and proves the
unconditional full-kernel perturbation budget `|P_n-P_n^(1)|<=2/(5k^2)`.
It also localizes the ceiling to same-point `B_1` derivatives and explicit
fourth- through sixth-derivative envelopes, avoiding adjacent-moment interval
cancellation.
A new `107,452`-tile Arb certificate proves the displayed ceiling on the
full real-parameter compact range `319<=t<=V'(2)`, using `1,073` positive
central blocks and a weakest margin above `1.14e-10`:

```text
outputs/jensen_window_pf_compound_order4_localized_curvature_compact_certificate.md
work/rh_compute/results/jensen_window_pf_compound_order4_localized_curvature_compact_certificate.json
python work/rh_compute/scripts/check_jensen_window_pf_compound_order4_localized_curvature_compact_certificate.py
```

At this intermediate stage only the analytic mode ray `u>=2` remained open
for entry. A rationalized theorem contract for that ray is recorded in
`outputs/jensen_window_pf_compound_order4_gaussian_cumulant_ray_target.md`:
the exact epsilon-six cumulant algebra and seven conditional full-collar
checks pass. The companion 1,800,000-block finite theorem and
coefficient-positive asymptotic theorem prove the formal corridors for every
`u>=2`. The first omitted epsilon-eight coefficient layer is also proved
globally, reducing the remaining exact-density central and two-tail errors to
the sufficient scaled budgets `1/100` on `2<=u<=20` and `1/(50u)` on
`u>=20`:

```text
outputs/jensen_window_pf_compound_order4_formal_cumulant_corridor_certificate.md
outputs/jensen_window_pf_compound_order4_formal_cumulant_asymptotic_ray_certificate.md
```

A
rational sequence with `s=(3/10,2/5,41/100,49/100,11/20)` satisfies the
strict ratio, decreasing-defect, increasing-scaled-defect, and cubic cones and
every available arbitrary-column order-three sign while having `H_(4,0)<0`.
Therefore order four is a genuinely new compound inequality, not a formal
consequence of the completed lower layers.

The recent positive order-moment transport theorem has also been checked in
`outputs/jensen_window_pf_order_moment_transport_fit_gate.md`. Its complete
monotonicity hypothesis fails for the transformed first Newman summand at the
origin, and its positive Hankel orientation is not the reciprocal
sign-regular orientation needed here. It does not discharge the curvature target.

The downstream exact-corridor, compact-heat asymptotic, and finite-confinement
theorems now close the complete contiguous order-four layer:

```text
outputs/jensen_window_pf_compound_order4_uniform_heat_forward_invariance_certificate.md
H_(4,n)(lambda)>0 for every integer n>=0 and every lambda in [-100,0]
```

The proof combines all-shift entry at `lambda=-100`, a compact-uniform
eventual positive tail, and backward variation of constants in the
cooperative flow. Reverse the first `N+1` columns of a finite four-row Hankel
block. Every solid minor of the reversed block is a completed signed
contiguous minor, so every initial minor is positive. The rectangular
Gasca-Pena criterion then makes the block strictly totally positive:

```text
outputs/jensen_window_pf_order4_noncontiguous_total_positivity_transfer.md
R_(4,n)(j_1,j_2,j_3,j_4)>0
for every n>=0 and 0<=j_1<j_2<j_3<j_4.
```

The same argument transfers contiguous layers through any fixed order `m` to
all arbitrary columns through `m`. It is the transfer used below after the
now-completed contiguous order-five theorem. PF-infinity, the all-degree
Jensen bridge, RH, and `Lambda<=0` remain open.

Its uniform tail and dynamical propagation are already reduced exactly:

```text
outputs/jensen_window_pf_compound_order5_uniform_tail_flow_reduction.md
[h^0,...,h^10] det K=[0,0,0,0,0,0,0,0,0,0,294912*G_2^10]
Q_n'=a_n*Q_(n+1)+b_n*Q_n,
a_n=(4*n+34)*H_(4,n)/H_(4,n+1)>0.
```

The determinant identity gives one compact-uniform eventual positive tail,
and finite confinement reduced the then-open endpoint problem to

```text
target_compound_order5_m100_entry:
H_(4,n+1)(-100)^2>H_(4,n)(-100)*H_(4,n+2)(-100), every n>=0,
equivalently H_(5,n)(-100)>0 for every n>=0.
```

An exact rational countermodel passes every available signed contiguous layer
through order four while having `H_(5,0)<0`, so this endpoint theorem cannot
be promoted from the completed lower layers.

The exact stable factorization and a new rigorous Arb prefix now close the
endpoint through `n=316`:

```text
outputs/jensen_window_pf_compound_order5_m100_prefix_certificate.md
H_(5,n)(-100)=W_n*J_n(-100)>0 for every 0<=n<=316,
relative_316>0.006269.
```

Thus only `J_n(-100)>0` on the analytic tail `n>=317` remained at that
stage.

Exact logarithmic reduction sharpens the tail once more. With `k=n+4`,

```text
outputs/jensen_window_pf_compound_order5_m100_tail_curvature_reduction.md
C_n=Delta^2 log(F_n)-Delta^2 log(d_(n+3)),
J_n>0 iff C_n<-4log(x_k).
```

The proved defect anchor makes the deliberately loose ceiling

```text
C_n<=100/(n+4)^2, every n>=317,
```

sufficient for the complete tail.

The cancellation-preserving first/full bridge is

```text
outputs/jensen_window_pf_compound_order5_first_summand_curvature_bridge.md
work/rh_compute/results/jensen_window_pf_compound_order5_first_summand_curvature_bridge.json
python work/rh_compute/scripts/check_jensen_window_pf_compound_order5_first_summand_curvature_bridge.py
```

It proves positive floors at both nested logarithmic layers and the
unconditional perturbation theorem

```text
|C_n-C_n^(1)|<=37/k^2, k=n+4>=321.
```

The exact tent identity then leaves the continuous first-summand target

```text
q_1''(t)<=60/t^2 for every real t>=320,
```

whose `63/k^2` discrete budget combines with `37/k^2` to recover the original
`100/k^2` ceiling.

That continuous theorem is now closed on three overlapping ranges:

```text
outputs/jensen_window_pf_compound_order5_nested_curvature_compact_certificate.md
outputs/jensen_window_pf_compound_order5_nested_curvature_finite_ray_certificate.md
outputs/jensen_window_pf_compound_order5_nested_curvature_asymptotic_ray_certificate.md
```

The hashed compact cover proves the ceiling on `320<=t<=V'(2)`. A rigorous
mode-two collar and 1850 exact-cumulant blocks prove it for `2<=u<=20`.
Normalized H-derivative boxes, analytic stable logarithms, and a single
dimensionless interval prove

```text
t^2*q_1''(t)<10<60 for every mode u>=20.
```

The endpoint composition

```text
outputs/jensen_window_pf_compound_order5_m100_entry_certificate.md
H_(5,n)(-100)>0 for every integer n>=0
```

therefore discharges `target_compound_order5_m100_entry`. Finally,

```text
outputs/jensen_window_pf_compound_order5_uniform_heat_forward_invariance_certificate.md
H_(5,n)(lambda)>0 for every n>=0 and -100<=lambda<=0
R_(5,n)(j_1,j_2,j_3,j_4,j_5;lambda)>0
```

for every increasing column quintuple on the same heat interval. Contiguous
and arbitrary-column order five are complete.

The next stable layer is now complete as well. The exact order-six flow
reduction and prefix are recorded in

```text
outputs/jensen_window_pf_compound_order6_uniform_tail_flow_reduction.md
outputs/jensen_window_pf_compound_order6_m100_prefix_certificate.md
```

The reduction gives a uniform eventual signed tail and a cooperative system
for `Q_(6,n)=-H_(6,n)`. The 1024-bit prefix proves entry through `n=316`.
The canonical factorization

```text
H_(5,n)=A_(n+4)^5*exp(p(n+4))
```

reduces the remaining tail to `P_k<=320/k^2` for `k=n+5>=322`. The
inverse-seventh-power first-summand theorem and a degree-75
coefficient-positive transfer reduce this to

```text
p_1''(t)<=200/t^2 for every real t>=321.
```

That continuous theorem is closed on three ranges:

```text
outputs/jensen_window_pf_compound_order6_nested_curvature_compact_certificate.md
outputs/jensen_window_pf_compound_order6_nested_curvature_finite_ray_certificate.md
outputs/jensen_window_pf_compound_order6_nested_curvature_asymptotic_ray_certificate.md
```

The compact H2-H10 cover proves it on `321<=t<=V'(2)`. A mode-two collar and
17,999 exact-corridor blocks prove it for `2<=u<=20`. Coarse exact ninth- and
tenth-cumulant corridors, normalized H boxes, and one dimensionless stable-log
interval prove

```text
t^2*p_1''(t)<22.769<200 for every mode u>=20.
```

The endpoint and forward compositions then prove

```text
outputs/jensen_window_pf_compound_order6_m100_entry_certificate.md
Q_(6,n)(-100)=-H_(6,n)(-100)>0 for every n>=0

outputs/jensen_window_pf_compound_order6_uniform_heat_forward_invariance_certificate.md
Q_(6,n)(lambda)=-H_(6,n)(lambda)>0
epsilon_6*R_(6,n)(j_1,...,j_6;lambda)>0
```

for every shift, every increasing column sextuple, and every
`-100<=lambda<=0`. Contiguous and arbitrary-column order six are complete;
the first unclosed all-shift compound layer is order seven.

The universal fixed-order tail theorem is now explicit:

```text
outputs/jensen_window_pf_graded_kernel_vandermonde_all_order_lemma.md
for every fixed m there exists N_m such that Q_(m,n)(lambda)>0 for n>=N_m
```

uniformly on `-100<=lambda<=0`. The threshold may depend on `m` and is
non-effective, so this does not reverse the quantifiers or fill any finite
prefix.

At order seven the exact specialization and flow reduction are

```text
outputs/jensen_window_pf_compound_order7_uniform_tail_flow_reduction.md
52183852646400*G_2^21*h^21
Q_(7,n)*H_(5,n+2)=Q_(6,n+1)^2-Q_(6,n)*Q_(6,n+2)
```

The completed order-six cone makes the order-seven heat flow cooperative.
Thus an all-shift `lambda=-100` entry theorem would propagate to the full heat
interval and then to arbitrary columns. The reduction also contains an exact
lower-cone countermodel, so orders at most six cannot supply that entry by
formal promotion.

The rigorous endpoint prefix is

```text
outputs/jensen_window_pf_compound_order7_m100_prefix_certificate.md
Q_(7,n)(-100)>0 for every 0<=n<=314
L_314>9/1000
```

and the missing `n>=315` tail has the canonical reduction

```text
outputs/jensen_window_pf_compound_order7_m100_tail_curvature_reduction.md
R_k<=900/k^2 for every k>=321
```

The complete-to-first-summand transfer is now closed by

```text
outputs/jensen_window_pf_negative_lambda_first_summand_power8_rebalanced_dominance_extension.md
0<=delta_k<2/k^8 for every k>=300

outputs/jensen_window_pf_compound_order7_first_summand_curvature_bridge.md
min(T_j,T_j^(1))>=3/(2*j), j>=320
|R_k-R_k^(1)|<262/k^2, k>=321
r_1''(t)<=600/t^2 => R_k<863/k^2<900/k^2
```

The remaining fixed-order-seven work is now only the continuous fourth-nested curvature theorem
`r_1''(t)<=600/t^2` on `t>=320`. All-shift order seven remains open. These
advances are fixed-order and eventual-tail progress, not
PF-infinity, the all-degree Jensen bridge, RH, or `Lambda<=0`.

### 2026-07-16 Endpoint Update

The historical order-seven frontier above is now superseded. Signed
contiguous and arbitrary-column layers are rigorous through order nine on
the full heat interval. Moreover, the arbitrary-order heat step is closed:

```text
outputs/jensen_window_pf_all_order_endpoint_heat_reduction.md
[Q_(m,n)(-100)>0 for every m>=10,n>=0]
iff
[Q_(m,n)(lambda)>0 for every m>=10,n>=0 and -100<=lambda<=0].
```

The former endpoint antecedent has the normalized Jacobi-Trudi coordinate

```text
outputs/jensen_window_pf_endpoint_deep_schur_coordinate.md
h_k=A_k(-100)/A_0(-100),
Q_(m,n)(-100)=A_0(-100)^m s_((n+m-1)^m)(h).
```

Arbitrary columns are in bijection with partitions
`lambda_1>=...>=lambda_m>=m-1`. Therefore the proposed all-order
signed-Hankel input was exactly the deep rectangular hierarchy

```text
s_((N^m))(h)>0 for every m>=10 and N>=m-1.
```

That hierarchy is false for the actual endpoint sequence. The checked artifact

```text
outputs/jensen_window_pf_endpoint_order10_counterexample.md
```

uses independent `4096`-bit Arb Hankel, Jacobi-Trudi, and Toda evaluations to
certify

```text
s_((N^10))(h)<0 for N=9,10,11,12,
s_((13^10))(h)>0.
```

The four negative shapes lie inside the required deep cone. The fixed-order
initial-minor theorem and heat equivalence remain valid conditional results,
but their all-order endpoint antecedent cannot be supplied. Full PF-infinity
was already unavailable because rigorous endpoint balls prove
`s_(1,1,1)(h)<0`; the direct order-ten result now rejects the narrower deep
hierarchy itself.

The endpoint hierarchy also has the exact rectangular Toda coordinate

```text
outputs/jensen_window_pf_deep_schur_toda_boundary_gate.md
tau_(m+1,N) tau_(m-1,N)
 =tau_(m,N)^2-tau_(m,N-1) tau_(m,N+1).
```

With a positive lower row, the next rectangle is positive exactly when the
current row is strictly log-concave in its width. This identity is
subtractive, and the first previously open row fails that inequality at
widths `9,10,11,12`. It explains the failure; it is not a propagation theorem.

The same gate blocks two overstrong bridge routes. The ordinary shifted tail
with negative indices reset to zero does not reproduce shallow translated
Jacobi-Trudi determinants: the smallest defect is
`h_0 h_2/h_1^2>0`. More decisively, the Edrei specialization

```text
H(z)=exp(z/100)/((1-z)(1-2z))
```

is strictly Schur-positive for every partition, yet its shift-zero cubic
Jensen polynomial has discriminant
`-222484532394597/2000000000000<0`. Hence even strict full unweighted
Schur/PF positivity does not imply Jensen hyperbolicity. Candidate Theorem
B-Star cannot be applied because its actual signed-Hankel antecedent fails;
the surviving theorem must use weaker properties specific to the actual
Xi/Phi coefficients and cannot be about Schur-positive sequences alone.

Degree 3 is the first nontrivial bridge obstruction. The exact algebra gate in:

```text
outputs/jensen_hankel_bridge_algebra.md
work/rh_compute/results/jensen_hankel_bridge_algebra.json
```

gives a positive rational sequence that passes finite reshaped-Hankel signs for
`k = 2,3`, `N = 3`, but has negative degree-3 Jensen discriminant. Therefore
low-order finite sign checks cannot be promoted into Jensen hyperbolicity.

## Proof Obligations

1. Replace the rejected all-order shifted signed-Hankel antecedent with a
   weaker condition that is provably satisfied by the actual Xi/Phi
   coefficients. The condition must tolerate the four negative order-ten
   endpoint rectangles while retaining enough structure to control Jensen
   windows.
2. Prove an Xi/Phi-specific theorem that converts that weaker structure into
   hyperbolicity of every Jensen polynomial, not just degree 2. Strict
   all-shape Schur positivity alone has an exact cubic counterexample, and the
   actual deep rectangular cone is false.
3. Equivalently, prove that every binomially weighted Jensen window
   `B^{d,n,0}_j = binom(d,j) A_{n+j}(0)` is PF-infinity, using hypotheses
   weaker than all-shift signed-Hankel positivity, or prove hyperbolicity
   directly without this PF reformulation.
4. Identify the exact theorem ecosystem: signed total positivity,
   sign-regular kernels, compound matrices, variation-diminishing transforms,
   multiplier sequences, or a new Xi-specific determinant identity.
   The reciprocal-gamma audit reduces this to higher compound concentration
   or a positive coupling theorem for the Xi scale mixture; common-scale sign
   regularity alone is already proved and positive-mixture closure is false.
5. Prove the limiting passage from all Jensen polynomials to the required
   Laguerre-Polya or Newman conclusion without assuming RH.
6. Preserve the heat-flow parameter discipline: finite positive-lambda
   certificates can guide the theorem, but they cannot replace a lambda-zero
   proof or a uniform interval theorem.

## Kill Gates

Reject a proposed proof if it:

```text
uses only finitely many A_k(lambda);
uses only the unshifted n=0 reshaped-Hankel block;
treats degree-2 Jensen as representative of all degrees;
assumes Jensen hyperbolicity, Laguerre-Polya membership, RH, or Lambda <= 0;
uses asymptotic fixed-degree Jensen results as all-shift Jensen results;
depends only on local zero repulsion;
confuses ordinary coefficient PF with signed-Hankel sign consistency.
assumes the rejected all-shift signed-Hankel/deep-rectangle hierarchy.
```

## Current Status

This theorem target is open, but the former all-order signed-Hankel route to
it is closed by counterexample. The fixed-order certificates and exact
low-degree gates constrain a replacement theorem; they do not prove the
Xi/Phi Jensen bridge.
