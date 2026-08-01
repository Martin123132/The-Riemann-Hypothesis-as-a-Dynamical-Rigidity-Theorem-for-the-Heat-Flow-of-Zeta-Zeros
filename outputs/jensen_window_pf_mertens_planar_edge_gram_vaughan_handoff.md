# Jensen-Window PF Mertens Planar Edge-Gram/Vaughan Handoff

Date: 2026-07-24

Status: exact edge-Gram and Vaughan-symmetrization reduction with
one open signed edge-pair gate; not a proof of RH, PF-infinity,
or `Lambda <= 0`.

```text
work/rh_compute/results/
jensen_window_pf_mertens_planar_edge_gram_vaughan_handoff.json
python work/rh_compute/scripts/
check_jensen_window_pf_mertens_planar_edge_gram_vaughan_handoff.py
```

## Exact Edge Lift

Let

```text
E_K={(p,q):K<p<q<=4K-1},
c_(p,q)=mu(p)mu(q),
Phi_(x,H)(p,q)=1_(p<=x)1_(q-p<=H).       (PEGV.1)
```

Then the planar prefix from Corollary 11.22Z.5 is

```text
S_K(x,H)=sum_(e in E_K)c_e Phi_(x,H)(e). (PEGV.2)
```

Put `nu_(x,H)=|Delta Wtilde_K(x,H)|` and define the edge
feature

```text
v_e(x,H)=sqrt(nu_(x,H))Phi_(x,H)(e).      (PEGV.3)
```

Its Gram kernel is the explicit tail mass

```text
G_K(e,e')
 =<v_e,v_(e')>
 =sum_(x>=max(p,p'))
   sum_(H>=max(q-p,q'-p'))nu_(x,H).       (PEGV.4)
```

Consequently

```text
E_(alpha,K)
 =||sum_e c_e v_e||_2^2
 =sum_(e,e')c_e c_(e')G_K(e,e').         (PEGV.5)
```

`G_K` is positive semidefinite and entrywise nonnegative. This
positivity lives in edge-pair space; it does not control the sign
of its Mobius-labelled off-diagonal.

## The Entire Diagonal Closes

The edge diagonal is

```text
D_(alpha,K)
 =sum_e mu(p)^2mu(q)^2G_K(e,e)
 =sum_(x,H)nu_(x,H)
   sum_e mu(p)^2mu(q)^2Phi_(x,H)(e).      (PEGV.6)
```

There are `3K-1` available vertices, so every threshold contains
at most

```text
binom(3K-1,2)<(9/2)K^2
```

edges. Corollary 11.22Z.5 therefore gives

```text
D_(alpha,K)
 <(9/2)K^2 V_(alpha,K)
 <(27/2)pi^2 K(1+log(2R)).               (PEGV.7)
```

Thus `D_(alpha,K)=O_epsilon(K^(1+epsilon))` for every
`epsilon>0`. The energy obstruction is purely the signed
off-diagonal edge-pair form

```text
C_(alpha,K)
 :=sum_(e!=e')c_e c_(e')G_K(e,e').       (PEGV.8)
```

Since `E=D+C` and `D` is already at the target scale,

```text
E_(alpha,K)=O_epsilon(K^(1+epsilon))
 iff
|C_(alpha,K)|=O_epsilon(K^(1+epsilon)).  (PEGV.9)
```

The off-diagonal splits exactly into pairs of distinct edges
sharing one vertex and pairs with four distinct endpoints. This
is a weighted collision/disjoint three- and four-point Mobius
problem. Estimating the two strata separately is sufficient but
not equivalent, because they may cancel.

For example, the same-base collision at a fixed threshold is

```text
sum_p mu(p)^2
 [(sum_q mu(q))^2-sum_q mu(q)^2],         (PEGV.10)
```

with the admissible tip interval determined by `(x,H)`. This
exposes a local Mertens-square input; it does not prove one.

## Exact Vaughan Symmetrization

For real sequences define

```text
B_t(f,g)=sum_(p,q)Phi_t(p,q)f(p)g(q),
A_t=(B_t+B_t^*)/2.                       (PEGV.11)
```

Although `B_t` is directed,

```text
S_t=B_t(mu,mu)=<mu,A_t mu>.              (PEGV.12)
```

For `1<=U,V` with `UV<=K` and for the present support `n>K`,
the pointwise finite Vaughan identity is

```text
mu=-u_I+u_II,

u_I(n)=sum_(b<=U,c<=V,bc|n)mu(b)mu(c),
u_II(n)=sum_(b>U,c>V,bc|n)mu(b)mu(c).    (PEGV.13)
```

Equivalently, with

```text
a_d=sum_(bc=d,b<=U,c<=V)mu(b)mu(c),
b_d=sum_(c|d,c>V)mu(c),

u_I(n)=sum_(d|n,d<=UV)a_d,
u_II(n)=sum_(d|n,d>V,n/d>U)b_d mu(n/d).  (PEGV.14)
```

The directed four-term form is

```text
S_t
 =B_t(u_I,u_I)+B_t(u_II,u_II)
  -B_t(u_I,u_II)-B_t(u_II,u_I).          (PEGV.15)
```

Both cross terms remain. Exact self-adjoint symmetrization
combines them, without dropping either one:

```text
S_t
 =<u_I,A_tu_I>+<u_II,A_tu_II>
  -2<u_I,A_tu_II>.                       (PEGV.16)
```

Equation (PEGV.16) holds as a vector identity after multiplying
each threshold by `sqrt(nu_t)`. The Hilbert norm must be taken
after the signed combination.

The exact specialization `U=V=1` is a cancellation guard. Then
`u_I=1` and `u_II=1+mu`. If

```text
N_t=B_t(1,1), L_t=B_t(1,mu),
R_t=B_t(mu,1),
```

the four terms in (PEGV.15) are

```text
N_t, N_t+L_t+R_t+S_t, N_t+L_t, N_t+R_t.
```

Their signed sum is exactly `S_t`. Thus deterministic
edge-count pieces cancel only across the complete Vaughan
combination. Also, `L_t` and `R_t` need not coincide for a
truncated directed band.

`A_t` is not positive semidefinite. Any included edge gives a
principal minor

```text
[[0,1/2],[1/2,0]]
```

with determinant `-1/4`. Positivity-only coercivity is therefore
not available.

## Fourier And Theorem-Fit Audit

At the full base prefix, let

```text
F_K(theta)=sum_(K<n<4K)mu(n)e(ntheta).
```

Then

```text
S_K(4K-2,H)
 =integral_0^1 |F_K(theta)|^2
   sum_(h<=H)cos(2pi h theta)dtheta.      (PEGV.17)
```

This confirms that the target is quadratic in one Mobius
exponential sum before the curvature energy squares it.

- Green and Tao, Lemma 4.1, supplies the exact finite Vaughan
  identity used in (PEGV.13)-(PEGV.16).
- Lewko and Lewko, Lemma 23, controls squared variations of
  linear interval sums after averaging over primitive
  characters. It supplies neither the quadratic edge lift nor
  the missing arithmetic cancellation in (PEGV.9).
- Tao and Teravainen prove quantitative `U^k` uniformity of
  Mobius with iterated-logarithmic decay. Even a formal transfer
  to this weighted four-linear form is power-short.
- Matomaki, Radziwill, and Tao prove `o(H^k X)` averaged Chowla.
  Even granting the required weighted-cutoff transfer, its
  natural `o(K^3)` scale after the `1/K` curvature normalization
  remains two powers above `K^(1+epsilon)`.

The generic geometry-only estimate is

```text
E_(alpha,K)=O(K^3(1+log(2R))),            (PEGV.18)
```

which records the same two-power deficit.

## Open Gate And Proof Boundary

For every member of one fixed cofinal positive-`alpha` sequence,
the sufficient energy route now asks for

```text
|C_(alpha,K)|=O_epsilon(K^(1+epsilon))   (PEGV.19)
```

on dyadic `K`, for every `epsilon>0`. Alternatively, one may
bypass this stronger three/four-point energy condition by
directly proving the signed planar criterion (11.22Z.99).

The edge-Gram identity, diagonal estimate, pointwise Vaughan
regrouping, self-adjoint symmetrization, unit-cutoff cancellation,
and theorem-fit exclusions are exact. No estimate of (PEGV.19),
curvature-energy gain, full Burnol bound, RH, PF-infinity,
`Lambda<=0`, or Clay-prize result is proved.

Primary sources:

- Green and Tao, *Quadratic Uniformity of the Mobius Function*,
  Lemma 4.1: https://doi.org/10.5802/aif.2401
- A. Lewko and M. Lewko, *A Variational
  Barban-Davenport-Halberstam Theorem*, Lemma 23:
  https://arxiv.org/abs/1111.6190
- Tao and Teravainen, *Quantitative Bounds for Gowers Uniformity
  of the Mobius and von Mangoldt Functions*:
  https://doi.org/10.4171/JEMS/1404
- Matomaki, Radziwill, and Tao, *An Averaged Form of Chowla's
  Conjecture*: https://doi.org/10.2140/ant.2015.9.2167

## Claim Ledger

| ID | Role | Status | Statement | Boundary |
|---|---|---|---|---|
| `pegv_01_inherited_energy` | `exact_equivalence` | `available_exact` | Corollary 11.22Z.5 isolates E_(alpha,K)=sum_t nu_t\|S_t\|^2 with nu_t=\|Delta Wtilde_K(t)\|. | The curvature-energy estimate itself remains open. |
| `pegv_02_edge_set` | `exact_definition` | `available_exact` | E_K={(p,q):K<p<q<=4K-1}; c_(p,q)=mu(p)mu(q). | The edge orientation records the positive shift q-p. |
| `pegv_03_threshold_feature` | `exact_definition` | `available_exact` | Phi_(x,H)(p,q)=1_(p<=x)1_(q-p<=H). | The global edge set already enforces q<=4K-1. |
| `pegv_04_prefix_edge_pairing` | `exact_identity` | `available_exact` | S_K(x,H)=sum_(e in E_K)c_e Phi_(x,H)(e). | This is the zero-extended planar prefix from Corollary 11.22Z.5. |
| `pegv_05_curvature_measure` | `exact_definition` | `available_exact` | nu_(x,H)=\|Delta Wtilde_K(x,H)\| and V_(alpha,K)=sum_(x,H)nu_(x,H). | The absolute value belongs to the sufficient energy route, not to the weaker direct signed pairing. |
| `pegv_06_edge_feature_vector` | `exact_definition` | `available_exact` | v_e(x,H)=sqrt(nu_(x,H))Phi_(x,H)(e). | Zero-curvature cells contribute a zero coordinate. |
| `pegv_07_tail_gram` | `exact_identity` | `available_exact` | G_K(e,e')=<v_e,v_e'>=sum_(x>=max(p,p')) sum_(H>=max(q-p,q'-p'))nu_(x,H). | All sums remain inside the finite planar rectangle. |
| `pegv_08_gram_psd` | `exact_property` | `available_exact` | G_K is positive semidefinite and entrywise nonnegative. | Positivity is in edge-pair space and does not determine the sign of the Mobius-labelled off-diagonal. |
| `pegv_09_energy_gram` | `exact_identity` | `available_exact` | E_(alpha,K)=sum_(e,e')c_e c_(e')G_K(e,e') =\|\|sum_e c_e v_e\|\|_2^2. | This is a quartic Mobius identity after expanding edge labels. |
| `pegv_10_diagonal` | `exact_identity` | `available_exact` | D_(alpha,K)=sum_e mu(p)^2mu(q)^2G_K(e,e) =sum_t nu_t sum_e mu(p)^2mu(q)^2Phi_t(e). | This is the diagonal e=e' in edge-pair space. |
| `pegv_11_edge_count` | `exact_bound` | `available_exact` | sum_e Phi_t(e)<=binom(3K-1,2)<(9/2)K^2. | There are 3K-1 possible ordinary-Mobius vertices. |
| `pegv_12_diagonal_bound` | `exact_bound` | `available_exact` | D_(alpha,K)<(27/2)pi^2 K(1+log(2R)). | Combine row 11 with V_(alpha,K)<3pi^2(1+log(2R))/K. |
| `pegv_13_diagonal_scale` | `exact_asymptotic` | `available_exact` | For every epsilon>0, D_(alpha,K)=O_epsilon(K^(1+epsilon)). | The implied constant may depend on epsilon; R<=K. |
| `pegv_14_offdiagonal_edge_form` | `exact_definition` | `available_exact` | C_(alpha,K)=sum_(e!=e')c_e c_(e')G_K(e,e'). | This ordered sum equals twice the unordered off-diagonal. |
| `pegv_15_energy_equivalence` | `exact_equivalence` | `available_exact` | E_(alpha,K)=O_epsilon(K^(1+epsilon)) iff \|C_(alpha,K)\|=O_epsilon(K^(1+epsilon)). | Use E=D+C and the already closed nonnegative diagonal D. |
| `pegv_16_collision_partition` | `exact_partition` | `available_exact` | C=C_collision+C_disjoint according as two distinct edges share one vertex or have four distinct endpoints. | Separate estimates are sufficient but are not asserted to be necessary because the two strata may cancel. |
| `pegv_17_same_base_collision` | `exact_identity` | `available_exact` | At a fixed threshold, same-base edge pairs contribute sum_p mu(p)^2[(sum_q mu(q))^2-sum_q mu(q)^2] over the admissible tips q. | The inner ranges depend on (x,H); this exposes a local Mertens-square obligation rather than closing it. |
| `pegv_18_band_operator` | `exact_definition` | `available_exact` | B_t(f,g)=sum_(p,q)Phi_t(p,q)f(p)g(q), and A_t=(B_t+B_t^*)/2. | B_t is directed while A_t is self-adjoint. |
| `pegv_19_symmetrization` | `exact_identity` | `available_exact` | S_t=B_t(mu,mu)=<mu,A_t mu>. | Every real quadratic form sees only the symmetric part. |
| `pegv_20_pointwise_vaughan` | `exact_identity` | `available_exact` | For 1<=U,V with UV<=K and n>K, define u_I(n)=sum_(b<=U,c<=V,bc\|n)mu(b)mu(c) and u_II(n)=sum_(b>U,c>V,bc\|n)mu(b)mu(c); then mu(n)=-u_I(n)+u_II(n). | Here n>max(U,V), the support condition underlying Green-Tao Lemma 4.1. |
| `pegv_21_grouped_vaughan` | `exact_identity` | `available_exact` | u_I(n)=sum_(d\|n,d<=UV)a_d and u_II(n)=sum_(d\|n,d>V,n/d>U)b_d mu(n/d), with a_d=sum_(bc=d,b<=U,c<=V)mu(b)mu(c) and b_d=sum_(c\|d,c>V)mu(c). | These are exact finite divisor regroupings. |
| `pegv_22_directed_four_term` | `exact_identity` | `available_exact` | S_t=B_t(u_I,u_I)+B_t(u_II,u_II)-B_t(u_I,u_II)-B_t(u_II,u_I). | The two directed cross terms must both be retained. |
| `pegv_23_symmetric_three_term` | `exact_identity` | `available_exact` | S_t=<u_I,A_tu_I>+<u_II,A_tu_II>-2<u_I,A_tu_II>. | This combines, but does not discard, the two directed cross terms. |
| `pegv_24_hilbert_vaughan` | `exact_identity` | `available_exact` | The vector (sqrt(nu_t)S_t)_t equals the corresponding signed three-term A_t combination before its Hilbert norm is squared. | Taking norms term by term would destroy the exact Vaughan cancellation. |
| `pegv_25_unit_cutoff` | `exact_identity` | `available_exact` | For U=V=1, u_I=1 and u_II=1+mu. | This is an exact actual-Mobius specialization, not a synthetic model. |
| `pegv_26_unit_cutoff_cancellation` | `exact_identity` | `available_exact` | Writing N_t=B_t(1,1), L_t=B_t(1,mu), R_t=B_t(mu,1), the four Vaughan terms are N_t, N_t+L_t+R_t+S_t, N_t+L_t, and N_t+R_t, whose signed sum is S_t. | Deterministic edge-count pieces cancel only in the full four-term combination. |
| `pegv_27_cross_term_guard` | `proof_guard` | `guard_active` | L_t and R_t need not agree for a directed truncated band. | Replacing the two cross terms by one directed term is invalid; only self-adjoint symmetrization combines them. |
| `pegv_28_indefinite_band_guard` | `proof_guard` | `guard_active` | Whenever a threshold contains an edge, A_t has a 2 by 2 principal minor [[0,1/2],[1/2,0]] with determinant -1/4. | The symmetrized band is indefinite, so positivity-only coercivity is unavailable. |
| `pegv_29_generic_bound` | `exact_bound` | `available_exact` | \|S_t\|<=binom(3K-1,2) and hence E_(alpha,K)=O(K^3(1+log(2R))) without arithmetic input. | This uses only \|mu\|<=1 and the total curvature variation. |
| `pegv_30_generic_power_guard` | `proof_guard` | `guard_active` | The generic edge-Gram estimate is two powers above K^(1+epsilon). | A geometry-only Bessel, trace, or operator-norm argument cannot be promoted to the arithmetic target. |
| `pegv_31_fourier_full_prefix` | `exact_identity` | `available_exact` | For F_K(theta)=sum_(K<n<4K)mu(n)e(ntheta), S_K(4K-2,H)=integral_0^1 \|F_K(theta)\|^2 sum_(h<=H)cos(2pi h theta)dtheta. | This concerns the full base prefix; varying x remains a bilinear partial-sum problem. |
| `pegv_32_variational_large_sieve` | `literature_guard` | `guard_active` | Lewko-Lewko Lemma 23 controls squared variations of linear interval sums after averaging over primitive characters. | The present target is one quadratic edge coefficient vector against nested threshold features, with no character or separated-frequency average; the theorem does not supply row 15. |
| `pegv_33_gowers_uniformity` | `literature_guard` | `guard_active` | Tao-Teravainen give quantitative U^k uniformity of Mobius with iterated-logarithmic decay for every fixed k. | Even granting a generalized-von-Neumann transfer to the weighted four-linear form, logarithmic decay does not bridge the K^3 to K^(1+epsilon) power gap. |
| `pegv_34_averaged_chowla` | `literature_guard` | `guard_active` | Matomaki-Radziwill-Tao give o(H^k X) after averaging k translation shifts. | The edge-pair lift has varying cutoffs and tail-Gram weights; even a formal K-scale transfer gives only o(K^3) after the 1/K curvature normalization, not K^(1+epsilon). |
| `pegv_35_signed_route_guard` | `proof_guard` | `guard_active` | The energy gate is sufficient and exposes a strong weighted three/four-point correlation; the direct signed pairing (11.22Z.99) remains a potentially weaker route. | Do not silently replace the actual goal by the stronger energy condition. |
| `pegv_36_open_edge_gate` | `open_gate` | `open` | Prove \|C_(alpha,K)\|=O_epsilon(K^(1+epsilon)) on dyadic K for every epsilon>0 and each member of one cofinal positive-alpha sequence, or bypass it by directly bounding the signed planar pairing. | No cited theorem currently supplies this weighted off-diagonal edge correlation. |
| `pegv_37_finite_validation` | `finite_validation` | `validated_finite` | Independent exact finite checks reproduce the pointwise Vaughan identity, directed and symmetrized band formulas, edge-Gram lift, diagonal split, unit-cutoff cancellation, indefinite minor, and Fourier identity. | Finite algebra checks do not prove the all-scale edge gate. |
| `pegv_38_proof_boundary` | `proof_guard` | `guard_active` | The edge-Gram lift, O(K log K) diagonal, Vaughan symmetrization, and theorem-fit exclusions are exact. | No off-diagonal edge cancellation, curvature-energy estimate, full Burnol bound, RH, PF-infinity, Lambda<=0, or Clay-prize result is proved. |

Summary:

- rows: 38
- exact reductions: 28
- proof guards: 5
- literature guards: 3
- open edge gates: 1
