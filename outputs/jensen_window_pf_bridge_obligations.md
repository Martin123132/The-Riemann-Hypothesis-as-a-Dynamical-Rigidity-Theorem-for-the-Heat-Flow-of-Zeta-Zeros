# Jensen-Window PF Bridge Obligations

Date: 2026-07-24

Status: theorem-obligation ledger. This is not a proof of PF-infinity, Jensen
hyperbolicity, Laguerre-Polya membership, RH, or `Lambda <= 0`; it decomposes
the open Jensen-window PF bridge into exact reformulations, finite evidence,
open obligations, and rejection tests.

## Purpose

The target `target_jensen_window_pf_bridge` is concrete but still too broad to
attack as one sentence. This ledger splits it into checkable obligations so a
future proof attempt can be rejected or advanced at the right layer.

Machine-readable ledger:

```text
work/rh_compute/results/jensen_window_pf_bridge_obligations.json
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_bridge_obligations.py
```

Current result:

```text
validated Jensen-window PF bridge obligations: 16 obligations, 0 issues, 3 open obligations
```

## Exact And Finite Layers

```text
jwpf_01_window_pf_jensen_equivalence:
  finite PF-infinity of B^{d,n,0} is the Toeplitz-total-positivity form of
  the real-nonpositive-root condition for one fixed Jensen window

jwpf_02_degree2_signed_hankel_contact:
  degree 2 matches the m=1 signed-Hankel condition exactly

jwpf_03_low_degree_extra_toeplitz_obligations:
  degrees 3 and 4 introduce extra banded Toeplitz obligations

jwpf_04_current_finite_pf_sturm_evidence:
  the earlier Arb/Sturm/PF manifests remain finite evidence; a separate
  real-zero-band/sector theorem now proves every shifted Xi Jensen polynomial
  through degree 361 uniformly on 0<=t<=1/5
```

These rows organize what we already know. None has `would_close_target=true`.
The promoted bounded-degree theorem is recorded in

```text
outputs/jensen_window_pf_newman_zero_slab_degree71_sector_certificate.md
work/rh_compute/results/jensen_window_pf_newman_zero_slab_degree71_sector_certificate.json
python work/rh_compute/scripts/check_jensen_window_pf_newman_zero_slab_degree71_sector_certificate.py
outputs/jensen_window_pf_newman_real_zero_band_degree361_sector_certificate.md
work/rh_compute/results/jensen_window_pf_newman_real_zero_band_degree361_sector_certificate.json
python work/rh_compute/scripts/check_jensen_window_pf_newman_real_zero_band_degree361_sector_certificate.py
```

It closes the actual Xi quartic and quintic layers on the positive Newman
target interval, but its finite cutoff does not prove degree 362, an unbounded
cofinal degree sequence, all-degree Jensen hyperbolicity, or PF-infinity.

## Rejected Antecedent And Open Bridge Obligations

```text
jwpf_05_all_order_shifted_sign_consistency:
  rejected: the actual endpoint sequence has Q_(10,n)(-100)<0 for
  n=0,1,2,3

jwpf_05b_weaker_xi_specific_antecedent:
  identify a weaker Xi/Phi-specific condition that survives those order-ten
  failures, excludes the exact outer-contact prefix whose 4,043 supported
  finite signed-Hankel minors coexist with a nonhyperbolic adjacent quintic,
  and is provably satisfied for every required degree and shift

jwpf_06_sign_regular_to_jensen_pf_conversion:
  legacy id: convert the weaker jwpf_05b structure, without all-shift
  signed-Hankel positivity, into every binomially weighted Jensen-window
  Toeplitz conclusion or directly into Jensen hyperbolicity; at quartic
  contact its low-degree specialization must force `x_5<=U`

jwpf_07_binomial_weight_and_shift_uniformity:
  handle binomial weights binom(d,j) and all shifts n uniformly
```

The central open theorem is `jwpf_06_sign_regular_to_jensen_pf_conversion`.
Its legacy identifier is retained for dependency stability; its admissible
antecedent is now `jwpf_05b_weaker_xi_specific_antecedent`, not the rejected
all-order signed-Hankel hierarchy.
The opposite guard is equally important:

```text
outputs/jensen_window_pf_quartic_outer_contact_normal_form_gate.md
outputs/jensen_window_pf_quartic_outer_contact_length14_survivor_gate.md
outputs/jensen_window_pf_quartic_quintic_polar_contact_lemma.md
```

The normal form supplies an exact negative-discriminant outward collar for
`0<q<15/16`. The survivor lies in that collar, passes all `4,043` supported
finite signed-Hankel minors visible on its prefix, and still has a
nonhyperbolic adjacent quintic. The polar lemma proves that adjacent quintic
hyperbolicity would force the contact threshold `x_5<=U`. A candidate
antecedent must therefore be weaker than the rejected all-shift signed-Hankel
cone on the actual Xi sequence, yet stronger than the survivor precisely in
this degree-coupled direction at unbounded degree. The sector theorem already
excludes this particular low-degree behavior for Xi through degree 361, but it
does not supply the all-degree antecedent required by `jwpf_06`.
The checker permits `would_close_target=true` only on open theorem-obligation
rows, never on finite evidence or countermodel rows.

The theorem machinery audit for this row is:

```text
outputs/jensen_window_pf_theorem_machinery_fit_matrix.md
work/rh_compute/results/jensen_window_pf_theorem_machinery_fit_matrix.json
python work/rh_compute/scripts/check_jensen_window_pf_theorem_machinery_fit_matrix.py
```

Current audit result:

```text
validated Jensen-window PF theorem machinery fit matrix: 11 rows, 0 issues, 0 ready-to-apply rows
```

The structural ansatz workbench for this row is:

```text
outputs/jensen_window_pf_structural_ansatz_matrix.md
work/rh_compute/results/jensen_window_pf_structural_ansatz_matrix.json
python work/rh_compute/scripts/check_jensen_window_pf_structural_ansatz_matrix.py
```

Current ansatz result:

```text
validated Jensen-window PF structural ansatz matrix: 6 ansatz rows, 0 issues, 0 ready-to-apply rows
```

## Downstream Conditional Layer

```text
jwpf_08_jensen_to_laguerre_polya_limit:
  once all-degree/all-shift Jensen hyperbolicity is proved noncircularly,
  document the limiting theorem and the normal-family/growth hypotheses for
  this heat-flow normalization
```

This row is conditional. It cannot be invoked before the all-degree Jensen
theorem is actually proved.

## Rejection Tests

```text
jwpf_09_finite_rectangle_promotion_rejected:
  finite Jensen-window PF/Sturm rectangles cannot be promoted to all-shift
  Jensen hyperbolicity by wording alone

jwpf_10_ordinary_coefficient_pf_route_separated:
  legacy id retained: finite PF evidence for c_k=A_k/k! remains separate and
  nonpromotable, but all-order PF-infinity of c is exactly equivalent to every
  binomially weighted Jensen window being PF-infinity
```

The exact correction is recorded in:

```text
outputs/jensen_window_pf_coefficient_pf_equivalence_gate.md
work/rh_compute/results/jensen_window_pf_coefficient_pf_equivalence_gate.json
python work/rh_compute/scripts/check_jensen_window_pf_coefficient_pf_equivalence_gate.py
```

Thus the coefficient-PF and Jensen-window formulations share one all-order
endpoint target. What remains missing is a noncircular proof of that endpoint
from Xi/Phi-specific or signed-Hankel structure; the current finite evidence on
either side does not supply it.

The same endpoint also has an exact nonlinear Hankel formulation:

```text
jwpf_11_edrei_log_stieltjes_endpoint_equivalence:
  for a_r=p_(r+1)=(-1)^r[z^r]H'/H, the all-order Stieltjes moment property
  is exactly equivalent to H in LP+, coefficient PF-infinity, and every
  Jensen-window PF target
```

The exact source and indexing gate are:

```text
outputs/jensen_window_pf_edrei_stieltjes_equivalence_gate.md
work/rh_compute/results/jensen_window_pf_edrei_stieltjes_equivalence_gate.json
python work/rh_compute/scripts/check_jensen_window_pf_edrei_stieltjes_equivalence_gate.py
```

The 4,205 finite power-Hankel rows remain nonpromotable. A strict proof would
need both determinant columns `s=1,2` at every order, or an equivalent
nonnegative S-fraction or Stieltjes representation of `H_0'/H_0`.

The Stieltjes ratio has an exact Pick-kernel formulation:

```text
jwpf_12_phi_pick_kernel_endpoint_equivalence:
  P_Phi(z)=-Im(F_0'(z)conj(F_0(z)))>0 for every Im(z)>0
  is exactly equivalent to the common LP+/coefficient-PF/Jensen/Stieltjes
  endpoint for the nonexponential Phi transform
```

The explicit polarized double integral, first real-axis wall, mixture guard,
and positive-resolvent alternative are:

```text
outputs/jensen_window_pf_phi_pick_kernel_target.md
work/rh_compute/results/jensen_window_pf_phi_pick_kernel_target.json
python work/rh_compute/scripts/check_jensen_window_pf_phi_pick_kernel_target.py
```

The heat-flow shortcut is separately rejected:

```text
jwpf_13_generic_backward_stieltjes_invariance_rejected:
  the repeated-zero LP+ model (1+beta*z)^2 exits every shifted 2x2
  Stieltjes wall immediately under backward radial heat
```

The exact hierarchy and countermodel are:

```text
outputs/jensen_window_pf_edrei_heat_flow_boundary_gate.md
work/rh_compute/results/jensen_window_pf_edrei_heat_flow_boundary_gate.json
python work/rh_compute/scripts/check_jensen_window_pf_edrei_heat_flow_boundary_gate.py
```

The general boundary flux and collision sensor are:

```text
outputs/jensen_window_pf_edrei_hankel_boundary_flux_gate.md
work/rh_compute/results/jensen_window_pf_edrei_hankel_boundary_flux_gate.json
python work/rh_compute/scripts/check_jensen_window_pf_edrei_hankel_boundary_flux_gate.py
```

They prove that fixed canonical order necessarily escapes at an
infinite-support boundary. An exact orthogonal-polynomial Schur quotient then
proves that a repeated boundary atom forces exponential high-shift
determinant-flux growth. This forces any backward route to supply
Xi/Phi-specific all-order rigidity; the resulting subexponential fixed-size
shifted determinant-flux criterion is exactly equivalent to `Lambda<=0`, and
its Phi-specific estimate remains open.

The exact raw-moment transfer and collision-resolution barrier sharpen what
that estimate must contain:

```text
outputs/jensen_window_pf_edrei_raw_moment_collision_resolution_gate.md
work/rh_compute/results/jensen_window_pf_edrei_raw_moment_collision_resolution_gate.json
python work/rh_compute/scripts/check_jensen_window_pf_edrei_raw_moment_collision_resolution_gate.py
```

At a repeated leading atom the shifted `2x2` determinant is exponentially
small relative to its rank-one summands, with condition number of order
`(beta/q)^s`. Symmetrically split simple factors converge quadratically in
coefficient and compact-open topology but do not reveal their simple-pair
asymptotic until
`s=2*log(beta/epsilon)/log(beta/q)+O(1)`. Therefore, fixed algebraic-accuracy
raw-moment asymptotics cannot close the collision criterion as independent
entrywise estimates. The same crossover survives any fixed larger simple
prefix and hence every fixed collision rank. A Phi-specific route must instead preserve the
determinant cancellation structurally or provide gap-adapted exponentially
small remainders.

The Pick endpoint now also has a published arithmetic-Hankel formulation:

```text
jwpf_14_suzuki_arithmetic_hankel_equivalence:
  on a cofinal omega_n decreasing to zero, prove det(I+/-K_n[t])!=0 for
  every finite t; Suzuki's published criterion also states the terminal
  limit J_n(t;r,r)->0 for every upper-half-plane r
```

The later exact reduction sharpens the cofinal logic:

```text
jwpf_15_suzuki_determinant_only_reduction:
  the all-time determinant family alone implies RH by causal-multiplier
  continuation and shifted-zero accumulation; the terminal premise then
  follows after RH; at fixed omega the determinant property is equivalent
  to reduced innerness and an exact horizontal zero-cancellation rule
```

The exact xi coordinate, Suzuki equivalence, Fredholm-series Hankel contact,
truncation-path spectral reformulation, and shortcut countermodels are:

```text
outputs/jensen_window_pf_xi_pick_suzuki_hankel_bridge.md
work/rh_compute/results/jensen_window_pf_xi_pick_suzuki_hankel_bridge.json
python work/rh_compute/scripts/check_jensen_window_pf_xi_pick_suzuki_hankel_bridge.py
outputs/jensen_window_pf_suzuki_spectral_frontier.md
work/rh_compute/results/jensen_window_pf_suzuki_spectral_frontier.json
python work/rh_compute/scripts/check_jensen_window_pf_suzuki_spectral_frontier.py
outputs/jensen_window_pf_suzuki_determinant_only_reduction.md
work/rh_compute/results/jensen_window_pf_suzuki_determinant_only_reduction.json
python work/rh_compute/scripts/check_jensen_window_pf_suzuki_determinant_only_reduction.py
outputs/jensen_window_pf_suzuki_fixed_omega_phase_diagram.md
work/rh_compute/results/jensen_window_pf_suzuki_fixed_omega_phase_diagram.json
python work/rh_compute/scripts/check_jensen_window_pf_suzuki_fixed_omega_phase_diagram.py
outputs/jensen_window_pf_suzuki_jordan_totient_sign_scout.md
work/rh_compute/results/jensen_window_pf_suzuki_jordan_totient_sign_scout.json
python work/rh_compute/scripts/check_jensen_window_pf_suzuki_jordan_totient_sign_scout.py
outputs/jensen_window_pf_suzuki_cofinal_monotonicity_hierarchy.md
work/rh_compute/results/jensen_window_pf_suzuki_cofinal_monotonicity_hierarchy.json
python work/rh_compute/scripts/check_jensen_window_pf_suzuki_cofinal_monotonicity_hierarchy.py
outputs/jensen_window_pf_suzuki_jordan_error_kernel_reduction.md
work/rh_compute/results/jensen_window_pf_suzuki_jordan_error_kernel_reduction.json
python work/rh_compute/scripts/check_jensen_window_pf_suzuki_jordan_error_kernel_reduction.py
outputs/jensen_window_pf_suzuki_cofinal_l2_hierarchy.md
work/rh_compute/results/jensen_window_pf_suzuki_cofinal_l2_hierarchy.json
python work/rh_compute/scripts/check_jensen_window_pf_suzuki_cofinal_l2_hierarchy.py
outputs/jensen_window_pf_jordan_muntz_causal_energy_bridge.md
work/rh_compute/results/jensen_window_pf_jordan_muntz_causal_energy_bridge.json
python work/rh_compute/scripts/check_jensen_window_pf_jordan_muntz_causal_energy_bridge.py
outputs/jensen_window_pf_jordan_muntz_burnol_hardy_intertwiner.md
work/rh_compute/results/jensen_window_pf_jordan_muntz_burnol_hardy_intertwiner.json
python work/rh_compute/scripts/check_jensen_window_pf_jordan_muntz_burnol_hardy_intertwiner.py
outputs/jensen_window_pf_burnol_cell_energy_tail_obstruction.md
work/rh_compute/results/jensen_window_pf_burnol_cell_energy_tail_obstruction.json
python work/rh_compute/scripts/check_jensen_window_pf_burnol_cell_energy_tail_obstruction.py
outputs/jensen_window_pf_burnol_tail_discrepancy_dyadic_reduction.md
work/rh_compute/results/jensen_window_pf_burnol_tail_discrepancy_dyadic_reduction.json
python work/rh_compute/scripts/check_jensen_window_pf_burnol_tail_discrepancy_dyadic_reduction.py
outputs/jensen_window_pf_weighted_fractional_autocorrelation_gram_bridge.md
work/rh_compute/results/jensen_window_pf_weighted_fractional_autocorrelation_gram_bridge.json
python work/rh_compute/scripts/check_jensen_window_pf_weighted_fractional_autocorrelation_gram_bridge.py
```

Suzuki's local interval for `nu*omega>1` is unconditional, but it is not the
global theorem. Along the continuous compact self-adjoint truncation path,
all-`t` determinant nonvanishing is exactly pointwise strict contractivity
`||K[t]||<1` at every finite `t`. It does not require, and in the desired
Hermite-Biehler case cannot have, one positive gap uniform in `t`.
Hilbert-Schmidt domination and Suzuki's `M_v^2 exp(4vt)<1` high-contour
certificate both have provably finite reach. Total nonnegativity,
real-boundary unimodularity without pole-free contour motion, and finite
`t` grids also do not close the surviving all-time determinant obligation.
For one fixed pair, contractivity is not promoted directly to the terminal
limit; the redundancy is a cofinal-sequence consequence.

At fixed `omega`, the determinant property is equivalent to reduced
meromorphic innerness and to multiplicity-preserving horizontal
zero cancellation. Outside the countable set of exact equal-height
zero-pair displacements, it is simply the shifted zero-free-half-plane
condition. Suzuki's scalar formulation leaves the direct arithmetic target

```text
x^(-1/2)*1_(1,infinity)(x)-h_omega^<1>(x) in L2(1,infinity)
```

for a sequence `omega->0`; eventual one-sign of `h_omega^<1>` is sufficient.
The accompanying 6000-point positive scout is useful reconnaissance but
proves neither condition on the unbounded half-line.

Suzuki's full logarithmic smoothing hierarchy broadens the target. If
`H_(omega,k+1)` is the logarithmic antiderivative of `H_(omega,k)`, then
eventual one-sign behavior at any level `k>=1` gives the same fixed-shift
innerness by Landau's theorem. Composed with the phase diagram,

```text
RH iff some sequence omega_j->0 has eventually one-signed
       H_(omega_j,k_j) for arbitrary selected integers k_j>=1.
```

Every smoothing weight remains signed. The `k=2` route is smoother than the
sampled `k=1` function, but is still an open unbounded arithmetic problem.

Abel summation sharpens that problem again: the Suzuki kernel annihilates
the full positive residue main term of the cumulative Jordan-totient mass.
Every smoothing level is exactly a signed convolution of the remaining
Mobius summatory error. Elementary absolute estimates reach only
`O(x^(1/2-omega)*(1+log x))`, so the missing theorem must provide genuine
cancellation or order, not main-term domination.

There is also an exact energy formulation at every smoothing order. If
`P_(omega,k)` is the logarithmic polynomial determined by the first `k`
central Taylor coefficients of the shifted xi quotient, then

```text
H_(omega,k)(exp(t))-P_(omega,k)(t) in L2(0,infinity)
iff the fixed-shift quotient is inner.
```

Consequently, such an estimate along any sequence `omega_j->0`, at
arbitrary selected orders `k_j`, is RH-equivalent. The Jordan-error
reduction makes the residual an explicit signed Mobius-error convolution.
Its squared-energy estimate is open. Boundary unimodularity, one fixed
shift, and every finite numerical energy cutoff are nonpromotion guards.
The analogous classical Riesz-smoothed totient bounds assume RH and extra
zero-derivative hypotheses, so they do not supply this antecedent.

The unsmoothed Jordan error gives a sharper finite target. Its full
generalized Muntz expansion has natural partial sums `E_(omega,N)`, and
finite Mellin-Plancherel gives an exact weighted all-height norm. The
internally audited causal criterion reduces the surviving obligation to

```text
sup_N ||E_(omega,N)||_H < infinity
```

for each shift in one explicit sequence `omega_j->0`. Equivalently, one
must bound uniformly in `N` the weighted integral of
`|zeta(1/2+omega+i*t)M_N(1/2+omega+i*t)|^2` over every height. The finite
Brownian max-kernel is totally nonnegative, but this proves only energy
nonnegativity. Equal causal and anti-causal boundary norms also show why
boundary Plancherel alone cannot provide the positive-time Hardy estimate.

The exact Hardy intertwiner makes this a published Nyman-style target rather
than only a new Jordan formulation:

```text
(1-2omega)||t^(-omega)f_(2omega,N)||_2
 <=||E_(omega,N)||_H
 <=(1+2omega)||t^(-omega)f_(2omega,N)||_2.
```

Here `f_(2omega,N)` is Burnol's finite natural fractional-part sum.
Burnol's sequence theorem and the Balazard-Saias conditional convergence
estimate source-back the cofinal equivalence. They do not supply the
unconditional uniform bound, which remains precisely the open obligation.

Reciprocal cells now make that obligation completely discrete. With

```text
S_(omega,N)(k)
 =sum_(d<=N)mu(d)d^(-2omega){k/d},

Q_(omega,N)
 =sum_(k>=1)k^(2omega-2)S_(omega,N)(k)^2,
```

the weighted Burnol norms are uniformly bounded in `N` if and only if
`Q_(omega,N)` is. Moreover, the stable prefix `k<=N` forces

```text
|sum_(d<=N)mu(d)d^(-(1+2omega))-1/zeta(1+2omega)|
 =O(N^(-1/2-omega)).
```

Abel summation then makes `Re(s)>1/2+omega` zero-free. This is an exact
necessary obstruction, not an estimate we have proved. The live target is
the uniform `Q_(omega,N)` bound, including its short-multiplicative-interval
cancellation, along one explicit cofinal shift sequence.

The finite-tail reduction separates that target without loss. If
`r_(omega,N)=1/zeta(1+2omega)-A_(omega,N)` and

```text
H_(omega,N)(k)
 =sum_(N<d<=k)mu(d)d^(-2omega)floor(k/d),
```

then

```text
S_(omega,N)(k)-S_(omega,infinity)(k)
 =H_(omega,N)(k)-k*r_(omega,N).
```

Consequently uniform `Q_(omega,N)` is equivalent to all three of:

```text
Q_(omega,infinity)<infinity,
r_(omega,N)=O(N^(-1/2-omega)),
sup_N sum_(j>=0)(2^j*N)^(2omega-2)
 V_(omega,N)(2^j*N)<infinity.
```

Here `V_(omega,N)` is the dyadic square of the omitted-divisor
discrepancy. Partitioning by `q=floor(k/d)` writes it using weighted Mobius
increments on `(k/(q+1),k/q]`. The three-way split is exact; none of its
three estimates is proved by the split.

The same norm has an exact stationary Gram form on logarithmic divisors.
Its zeta spectral density is nonnegative and its diagonal is uniformly
bounded, but all possible growth remains in a signed Mobius off-diagonal.
Thus kernel positivity does not close the uniform upper bound; the missing
arithmetic input may equivalently be sought as a weighted multiplicative
large-sieve or bilinear cancellation estimate.

The autocorrelation asymptotic now splits that kernel as an
Ornstein-Uhlenbeck leading term plus an exponentially smaller pointwise
remainder. The leading Gram is exactly the reciprocal-Mobius tail energy

```text
sum_(k=1)^N [k^(1-alpha)-(k-1)^(1-alpha)]
 |sum_(d=k)^N mu(d)/d|^2.
```

Uniform boundedness for every member of one cofinal sequence
`alpha_j->0` is RH-equivalent. The reduction does not prove that
boundedness and does not close the full Burnol norm: the natural remainder
has a sign-indefinite remainder density, with an explicit negative witness.
Generic logarithmic-spacing Hilbert bounds still cost order
`N^(1-alpha)`. The open handoff is therefore a direct comparison or a
Mobius-specific remainder cancellation estimate.

The positive OU energy has a further exact reduction. Weighted
Hardy/Copson identities make its limiting tail norm equivalent to

```text
sum_(k>=1)M(k)^2/k^(2+alpha).
```

Its cofinal RH criterion is therefore equivalently the dyadic Mertens mean square

```text
sum_(K<=k<2K)M(k)^2
 =O_epsilon(K^(2+epsilon)),
```

or one signed origin-anchored integrated two-point Mobius-correlation
bound. Generic Hardy again costs `N^(1-alpha)`. Published
averaged-Chowla estimates average every base shift, while almost-all
short-interval theorems can discard the exceptional origin slice; neither
proves this cumulative anchored estimate.

The exact weighted-prefix transfer now connects this target to the
`q=floor(k/d)` square function:

```text
outputs/jensen_window_pf_mertens_weighted_prefix_affine_defect_reduction.md
work/rh_compute/results/jensen_window_pf_mertens_weighted_prefix_affine_defect_reduction.json
python work/rh_compute/scripts/check_jensen_window_pf_mertens_weighted_prefix_affine_defect_reduction.py
```

The weighted prefix `A_alpha(k)=sum_(n<=k)mu(n)n^(-alpha)` has the same
reciprocal-zeta Hardy energy as the Mertens coordinate. On the first
post-prefix block,

```text
D_(omega,N)(k)=A_alpha(k)-A_alpha(N)-k*r_(omega,N).
```

A constant-tail countermodel proves that the discrepancy and tail rate
miss one constant anchor mode. Thus a proof through reciprocal cells must
also control `N^((alpha-1)/2)A_alpha(N)` in a summable dyadic theorem.

The complete adjacent-cutoff tail sequence identifies that theorem exactly:

```text
outputs/jensen_window_pf_mertens_anchor_logarithmic_tail_energy_reduction.md
work/rh_compute/results/jensen_window_pf_mertens_anchor_logarithmic_tail_energy_reduction.json
python work/rh_compute/scripts/check_jensen_window_pf_mertens_anchor_logarithmic_tail_energy_reduction.py
```

Weighted Hardy/Copson transfer and the stable-prefix identity give

```text
sum_(N>=1)N^(alpha-2)|A_alpha(N)|^2<infinity
 iff
sum_(N>=1)P_(alpha/2,N)/N<infinity.
```

Thus the anchor is not a separate all-cutoff coordinate: it is the
logarithmically summed stable-prefix tail energy. A critical scalar
countermodel has `sup_N P_N<infinity` but divergent `sum_N P_N/N`, so the
forced pointwise tail rate is one logarithm short. The open obligation is
the Mobius-specific logarithmic gain, not another generic Hardy estimate.

The dyadic cosine-mode reduction makes that gain finite-dimensional at
each scale:

```text
outputs/jensen_window_pf_mertens_dyadic_cosine_mode_bottleneck.md
work/rh_compute/results/jensen_window_pf_mertens_dyadic_cosine_mode_bottleneck.json
python work/rh_compute/scripts/check_jensen_window_pf_mertens_dyadic_cosine_mode_bottleneck.py
```

Exact Neumann cosine Parseval proves every mode
`r>=ceil(sqrt(K))` automatically summable. The open obligation is therefore
the dyadic mean plus `O(sqrt(K))` smooth anchored Mobius modes. Davenport's
arbitrary logarithmic savings retain the power `K^(1-alpha)` and do not
close this series. A square-root additive-twist bound is a sufficient
conditional calibration only; it is not proved.

The local-path transference splits that obligation exactly:

```text
outputs/jensen_window_pf_mertens_local_path_cosine_transference.md
work/rh_compute/results/jensen_window_pf_mertens_local_path_cosine_transference.json
python work/rh_compute/scripts/check_jensen_window_pf_mertens_local_path_cosine_transference.py
```

The mean is the affine defect `sqrt(K)r_K-tau_(K,0)` and every
nonconstant mode is a DCT coefficient of the anchored reciprocal-weighted
local Mertens path. Full-path Abel comparison with ordinary Mertens
increments is exact but alone does not control the low-mode projector. A
non-Mobius scalar construction has zero dyadic endpoint tails and zero
mean modes while its first nonconstant series diverges. Thus endpoint
control and coefficient-envelope or generic large-sieve bounds cannot
discharge the obligation; at this stage the missing input is signed,
Mobius-specific, and presented in two components.

The centered quotient now resolves that projection issue:

```text
outputs/jensen_window_pf_mertens_centered_bridge_vaughan_handoff.md
work/rh_compute/results/jensen_window_pf_mertens_centered_bridge_vaughan_handoff.json
python work/rh_compute/scripts/check_jensen_window_pf_mertens_centered_bridge_vaughan_handoff.py
```

The Abel map preserves constants and is uniformly invertible modulo
constants. Hence the nonconstant gate is exactly equivalent to centered
local Mertens energy with Brownian-bridge kernel `min(i,j)-ij/K`.
Its diagonal is summable, and its off-diagonal correlation has the exact
pre-collapse handoff `C_K=-TI_K+TII_K`. Separate absolute Type I/II
bounds and the published averaged-Chowla scale remain one full power
short. The affine mean defect remains an independent open component.

The affine component can now be localized and merged with that bridge:

```text
outputs/jensen_window_pf_mertens_affine_tent_bridge_handoff.md
work/rh_compute/results/jensen_window_pf_mertens_affine_tent_bridge_handoff.json
python work/rh_compute/scripts/check_jensen_window_pf_mertens_affine_tent_bridge_handoff.py
```

Writing `u_K=sqrt(K)c_(K,0)`, the normalized scale filter
`u_K-(1/2)u_(2K)` is invertible and cancels the mean kernel's infinite
plateau exactly, leaving a compact tent on `(K,4K)`. Its signed square and
the Brownian bridge combine into one local kernel with summable diagonal.
After splitting the support into two standard dyadic intervals, the exact
handoff is `O_(loc,K)=-TI_(loc,K)+TII_(loc,K)`. This is now one arithmetic
obligation rather than two independently presented pieces. Before the
ordinary-Mertens transfer, the affine rank-one term also restores the
bridge's constant direction and gives an equivalent positive compact
tail-lattice square. A second exact Abel isomorphism gives the equivalent
ordinary-Mobius criterion
`sum_K K^(-2-alpha) sum_t |beta_K+sum_(i=t)^(K-1)mu(K+i)|^2<infinity`,
where `beta_K` is a bounded-coefficient anchor supported on `[2K,4K)`.
The sufficient scale calibration is `H_K=O(K^(2+epsilon))` with
`epsilon<alpha`; coefficientwise control gives only `O(K^3)`. None of
the three coordinates supplies a Type I/II gain, and separate absolute
estimates and averaged Chowla remain one full power short.

The ordinary anchored criterion now has a pre-square vector Vaughan
form:

```text
outputs/jensen_window_pf_mertens_ordinary_vector_vaughan_handoff.md
work/rh_compute/results/jensen_window_pf_mertens_ordinary_vector_vaughan_handoff.json
python work/rh_compute/scripts/check_jensen_window_pf_mertens_ordinary_vector_vaughan_handoff.py
```

The suffix energy is `H_K=||B_K||_2^2` for a bounded
`K`-dimensional feature sum. Its explicit Gram diagonal is at most
`(5/4)K^2`. Componentwise finite Vaughan on the two standard intervals
then gives `B_K=-V_(I,K)+V_(II,K)` before squaring. The remaining
obligation is a weighted vector square-function estimate for this signed
difference. Separate Type I and Type II norm bounds are not equivalent,
and coefficientwise control still gives only `H_K=O(K^3)`.

The nested vector geometry is now diagonalized exactly:

```text
outputs/jensen_window_pf_mertens_mixed_boundary_sine_vaughan_handoff.md
work/rh_compute/results/jensen_window_pf_mertens_mixed_boundary_sine_vaughan_handoff.json
python work/rh_compute/scripts/check_jensen_window_pf_mertens_mixed_boundary_sine_vaughan_handoff.py
```

The min-kernel has half-odd sine angles `(2r-1)pi/(2K+1)`. Unfolding
the future anchor makes every spectral coefficient `X_(K,r)` one
ordinary-Mobius test with `O(K^(-1/2))` deterministic coefficients.
The exact target is
`sum_K K^(-alpha)sum_r|X_(K,r)|^2/(2r-1)^2<infinity`, and finite
Vaughan gives `X_(K,r)=-T_(I,K,r)+T_(II,K,r)` before squaring.
Square-root-plus-epsilon cancellation would suffice for
`2epsilon<alpha`; Davenport, Parseval, and scalar all-interval
logarithmic savings remain one power short.

The high half-odd modes are now removed exactly:

```text
outputs/jensen_window_pf_mertens_spectral_anchor_high_mode_reduction.md
work/rh_compute/results/jensen_window_pf_mertens_spectral_anchor_high_mode_reduction.json
python work/rh_compute/scripts/check_jensen_window_pf_mertens_spectral_anchor_high_mode_reduction.py
```

Writing `X_(K,r)=Y_(K,r)+beta_K h_(K,r)`, current-block Parseval and the
endpoint envelope give a `K/R^2` current tail and a
`beta_K^2/(KR)` anchor tail. The unconditional cutoff
`R_(alpha,K)=ceil(K^(1-alpha/2))` bounds the normalized high part by
`2K^(-1)+8K^(-alpha/2)`, which is summable on dyadic `K`. Thus for each
fixed positive `alpha` the bridge obligation is exactly the low-mode
square with `r<=R_(alpha,K)`. This reduction is not uniform at
`alpha=0`, and it proves no low-mode Mobius cancellation.

The surviving low modes now have one exact ordinary-Mobius kernel:

```text
outputs/jensen_window_pf_mertens_truncated_half_odd_kernel_handoff.md
work/rh_compute/results/jensen_window_pf_mertens_truncated_half_odd_kernel_handoff.json
python work/rh_compute/scripts/check_jensen_window_pf_mertens_truncated_half_odd_kernel_handoff.py
```

Its infinite odd completion is the rescaled `min(i,j)` kernel, while
its reduced ordinary-support diagonal is uniformly below `7pi^2/24`.
The endpoint column is an increasing suffix-Abel weight. Collecting the
two finite Vaughan intervals gives
`<u_I,G u_I>+<u_II,G u_II>-2<u_I,G u_II>`, so the bridge obligation is
exactly one weighted signed off-diagonal estimate. Kernel positivity,
bounded trace, and separate Type I/II bounds do not supply that
estimate.

That off-diagonal now has an exact two-axis variation handoff:

```text
outputs/jensen_window_pf_mertens_shift_kernel_variation_handoff.md
work/rh_compute/results/jensen_window_pf_mertens_shift_kernel_variation_handoff.json
python work/rh_compute/scripts/check_jensen_window_pf_mertens_shift_kernel_variation_handoff.py
```

For each shift `h=m-n`, the kernel row has height `O(1/K)` and

```text
Var_(n in I_(K,h))W_(K,h)(n)<6rho_K,
rho_K=pi^2*K/(2K+1)^2.
```

There is a dual `O(1/K)` variation estimate at every fixed base point.
The corresponding Abel formulas expose maximal fixed-origin
two-point correlations and local Mertens intervals. Absolute
summation along either coordinate remains power-short even under
conditional square-root input. The bridge therefore requires
joint cancellation across shifts, base points, Vaughan types, and possibly
dyadic scales.

Those two coordinates now admit one exact planar Abel handoff:

```text
outputs/jensen_window_pf_mertens_planar_abel_handoff.md
work/rh_compute/results/jensen_window_pf_mertens_planar_abel_handoff.json
python work/rh_compute/scripts/check_jensen_window_pf_mertens_planar_abel_handoff.py
```

After zero extension of the triangular support,

```text
O_(alpha,K)=sum_(n,h)S_K(n,h)Delta Wtilde_K(n,h),

sum_(n,h)|Delta Wtilde_K(n,h)|
 <3*pi^2*(1+log(2R))/K.
```

The current interior uses the `L1` norm of the finite odd Dirichlet
kernel; all transition and future pieces telescope by convexity.
Weighted Cauchy-Schwarz reduces the bridge to the sufficient
curvature-energy estimate
`E_(alpha,K)=O_epsilon(K^(1+epsilon))`. That arithmetic estimate
remains open. Its stronger maximal form already contains RH-scale
block-Mertens input, so it is recorded as a conditional calibration
rather than a proof.

The curvature energy now has an exact edge-Gram handoff:

```text
outputs/jensen_window_pf_mertens_planar_edge_gram_vaughan_handoff.md
work/rh_compute/results/jensen_window_pf_mertens_planar_edge_gram_vaughan_handoff.json
python work/rh_compute/scripts/check_jensen_window_pf_mertens_planar_edge_gram_vaughan_handoff.py
```

Represent each triangular-band pair by an edge and each planar prefix
by its threshold indicator. The resulting Gram diagonal satisfies

```text
D_(alpha,K)<(27/2)*pi^2*K*(1+log(2R)).
```

It is therefore already `O_epsilon(K^(1+epsilon))`. The sufficient
energy gate is now exactly equivalent to

```text
|C_(alpha,K)|=O_epsilon(K^(1+epsilon))
```

for the signed off-diagonal edge-pair term. This contains one-vertex collisions
and four-distinct-endpoint pairs; their signed contributions may cancel.
Finite Vaughan decomposition is valid here with `UV<=K` on `n>K`, and
symmetrization retains all four directed terms. The `U=V=1` audit
shows that their large deterministic pieces cancel only after they are
combined. The diagonal is closed, but the signed edge-pair estimate
remains open.

The edge form also has a complete endpoint-incidence expansion and a
strictly weaker direct signed projection:

```text
outputs/jensen_window_pf_mertens_planar_incidence_antidiagonal_handoff.md
work/rh_compute/results/jensen_window_pf_mertens_planar_incidence_antidiagonal_handoff.json
python work/rh_compute/scripts/check_jensen_window_pf_mertens_planar_incidence_antidiagonal_handoff.py
```

The exact endpoint formulas satisfy `C=C_3+C_4`. A bounded-sign
threshold witness has `C_3=210`, `C_4=-210`, and `C=0`, proving that
separate collision and four-distinct estimates can lose decisive
cancellation. Fair random cuts reproduce the strong energy gate
exactly, so they are a bilinear coordinate rather than a gain.

Before absolute mixed curvature, the current-current interior instead
has the exact odd-frequency projection

```text
O_(CC,int,K)
 =(4*pi^2/N^3)sum_(r<=R)sigma_(K,r)Y_(K,r).
```

This isolates the weaker open target
`sum_(r<=R)|Y_(K,r)|^2=O_epsilon(N^6*R^(-1)*K^epsilon)`.
Menon's recent almost-all short-interval results
(`arXiv:2607.15574`) have logarithmic gains and the wrong averaging
quantifiers for this anchored quadratic projection.

The complementary cells now have exact joined coordinates:

```text
outputs/jensen_window_pf_mertens_planar_boundary_flux_anchor_handoff.md
work/rh_compute/results/jensen_window_pf_mertens_planar_boundary_flux_anchor_handoff.json
python work/rh_compute/scripts/check_jensen_window_pf_mertens_planar_boundary_flux_anchor_handoff.py
```

One modewise flux law covers every curvature region. Because the mixed
stencil reaches `m+2=2K`, the top current-interior line already contains
the current-future shoulder `A_CF`. The lossless formula is

```text
O=O_(CC,int)+B_CC-A_CF+beta*gamma+(c/2)(beta^2-Q),
```

and the complement is the joined endpoint-spectral projection
`sum_r Z_(K,r)/q_r^2`. The subtraction `B_CC-A_CF` must precede absolute
values. Orthogonal completion gives the equivalent low-mode energy

```text
L=E_perp+c(beta+gamma/c)^2.
```

This closes the transition/future algebra, not its arithmetic estimate.
The joined projection now has a positive suffix-energy reduction:

```text
outputs/jensen_window_pf_mertens_planar_boundary_suffix_energy_handoff.md
work/rh_compute/results/jensen_window_pf_mertens_planar_boundary_suffix_energy_handoff.json
python work/rh_compute/scripts/check_jensen_window_pf_mertens_planar_boundary_suffix_energy_handoff.py
```

After retaining `B_CC-A_CF`, finite suffix Abel summation gives

```text
R_B=(E_B-D_B)/2+C_delta,
```

where `E_B>=0`, `D_B<5pi^2/12`, and
`|R_B-E_B/2|=O(1+log R)` unconditionally. Thus the direct route no longer
needs a separate `Z` estimate: its complete boundary/future part is
subpower exactly when `E_B` is subpower. Uniform odd-sine control further
gives `E_B<=pi(1/2+2pi)K^(-2)H_K`, where `H_K` is exactly the earlier
affine-tent anchored energy. The live options are the strong signed edge-pair
bound, the lossless low-mode energy bound, or the `Y` projection together
with this anchored `H_K/E_B` weighted Mertens-energy estimate. In
particular, the anchored `E_B` gate introduces no extra logarithmic weight
beyond the earlier `H_K` criterion.

Its exact strength is now classified:

```text
outputs/jensen_window_pf_mertens_planar_boundary_suffix_localization_gate.md
work/rh_compute/results/jensen_window_pf_mertens_planar_boundary_suffix_localization_gate.json
python work/rh_compute/scripts/check_jensen_window_pf_mertens_planar_boundary_suffix_localization_gate.py
```

The odd-sine tail gives `Delta_i>=pi^2/(4N^2)` outside an
`O(K^(alpha/2))` terminal collar. Bounded adjacent Mobius increments recover
that collar, so all-epsilon `E_B` subpower control along one fixed cofinal
positive-alpha sequence is RH-equivalent. At the same time, a bounded
terminal-anchor witness has `K^2 E_B/H_K->0`, ruling out a uniform reverse
norm comparison.

This changes the route decision. Proving `Y_(K,r)` and `E_B` separately is
not a weaker direct shortcut: the `E_B` half is already RH-strength. A
remaining candidate was to retain signed cancellation in
`O_(CC,int)+E_B/2`; the next gate tests that possibility exactly.

That final candidate join also closes exactly:

```text
outputs/jensen_window_pf_mertens_planar_joined_energy_equivalence_gate.md
work/rh_compute/results/jensen_window_pf_mertens_planar_joined_energy_equivalence_gate.json
python work/rh_compute/scripts/check_jensen_window_pf_mertens_planar_joined_energy_equivalence_gate.py
```

For `J_B=O_(CC,int)+E_B/2`,

```text
J_B=L/2+(D_B-D_L)/2-C_delta.
```

Here `L=E_perp+c(beta+gamma/c)^2`, both diagonals are uniformly bounded,
and `C_delta=O(log R)`. Therefore `J_B`, the full signed off-diagonal, and
`L` have equivalent all-epsilon and fixed-alpha weighted dyadic criteria.
The joined proposal returns the old lossless gate; it is not a fourth route.

The corresponding finite float64 scout is:

```text
outputs/jensen_window_pf_mertens_planar_curvature_energy_scout.md
work/rh_compute/results/jensen_window_pf_mertens_planar_curvature_energy_scout.json
python work/rh_compute/scripts/check_jensen_window_pf_mertens_planar_curvature_energy_scout.py
```

Its 28 rows through `K=1024` have max `K*V_K<6.571` and
`E_K/K<0.831`. These observations do not prove a uniform bound or an
asymptotic exponent.

A second guard blocks an attempted transfer from limiting fixed-shift
innerness:

```text
outputs/jensen_window_pf_fixed_shift_inner_reciprocal_boundary_separation_gate.md
work/rh_compute/results/jensen_window_pf_fixed_shift_inner_reciprocal_boundary_separation_gate.json
python work/rh_compute/scripts/check_jensen_window_pf_fixed_shift_inner_reciprocal_boundary_separation_gate.py
```

The exact four-zero model has finite causal `H2` energy and a rational
inner fixed-shift quotient, but infinite reciprocal square energy on the
same shifted line. Thus the logarithmic stable-prefix series cannot be
deduced from limiting Jordan/Burnol innerness by a universal same-shift
norm comparison. The cofinal full-Burnol and direct Mobius routes remain
open.

The alternative noncircular handoff applies Vaughan before correlation
collapse, with bounded shifted-Mobius tests. Exact dyadic reindexing gives

```text
sum_(h<X)sum_(Y<=X-h)C_h(Y)=B_1(X)-TI_X+TII_X.
```

The new live bilinear target is
`-TI_X+TII_X=O_epsilon(X^(2+epsilon))`, preserving signed cancellation
jointly over shift and terminal length. Applying Vaughan after collapse is
endogenous: its test contains `M(n-1)`, and generic Cauchy-Schwarz returns
the target Mertens energy.

The rigorous route counterexample is recorded in:

```text
outputs/jensen_window_pf_endpoint_order10_counterexample.md
work/rh_compute/results/jensen_window_pf_endpoint_order10_counterexample.json
python work/rh_compute/scripts/check_jensen_window_pf_endpoint_order10_counterexample.py
```

## Integration Points

This obligation ledger is tied to:

```text
outputs/jensen_window_pf_bridge_target.md
outputs/signed_hankel_jensen_dependency_graph.md
outputs/jensen_window_pf_theorem_machinery_fit_matrix.md
outputs/jensen_window_pf_structural_ansatz_matrix.md
outputs/jensen_window_pf_coefficient_pf_equivalence_gate.md
outputs/jensen_window_pf_edrei_heat_flow_boundary_gate.md
outputs/jensen_window_pf_phi_pick_kernel_target.md
outputs/jensen_window_pf_xi_pick_suzuki_hankel_bridge.md
outputs/jensen_window_pf_suzuki_spectral_frontier.md
outputs/jensen_window_pf_suzuki_determinant_only_reduction.md
outputs/jensen_window_pf_suzuki_fixed_omega_phase_diagram.md
outputs/jensen_window_pf_suzuki_jordan_totient_sign_scout.md
outputs/jensen_window_pf_suzuki_cofinal_monotonicity_hierarchy.md
outputs/jensen_window_pf_suzuki_jordan_error_kernel_reduction.md
outputs/jensen_window_pf_suzuki_cofinal_l2_hierarchy.md
outputs/jensen_window_pf_weighted_fractional_autocorrelation_gram_bridge.md
outputs/jensen_window_pf_weighted_autocorrelation_ou_tail_energy_reduction.md
work/rh_compute/results/jensen_window_pf_weighted_autocorrelation_ou_tail_energy_reduction.json
python work/rh_compute/scripts/check_jensen_window_pf_weighted_autocorrelation_ou_tail_energy_reduction.py
outputs/jensen_window_pf_ou_mertens_mean_square_reduction.md
work/rh_compute/results/jensen_window_pf_ou_mertens_mean_square_reduction.json
python work/rh_compute/scripts/check_jensen_window_pf_ou_mertens_mean_square_reduction.py
outputs/jensen_window_pf_mertens_weighted_prefix_affine_defect_reduction.md
work/rh_compute/results/jensen_window_pf_mertens_weighted_prefix_affine_defect_reduction.json
python work/rh_compute/scripts/check_jensen_window_pf_mertens_weighted_prefix_affine_defect_reduction.py
outputs/jensen_window_pf_mertens_anchor_logarithmic_tail_energy_reduction.md
work/rh_compute/results/jensen_window_pf_mertens_anchor_logarithmic_tail_energy_reduction.json
python work/rh_compute/scripts/check_jensen_window_pf_mertens_anchor_logarithmic_tail_energy_reduction.py
outputs/jensen_window_pf_mertens_dyadic_cosine_mode_bottleneck.md
work/rh_compute/results/jensen_window_pf_mertens_dyadic_cosine_mode_bottleneck.json
python work/rh_compute/scripts/check_jensen_window_pf_mertens_dyadic_cosine_mode_bottleneck.py
outputs/jensen_window_pf_mertens_local_path_cosine_transference.md
work/rh_compute/results/jensen_window_pf_mertens_local_path_cosine_transference.json
python work/rh_compute/scripts/check_jensen_window_pf_mertens_local_path_cosine_transference.py
outputs/jensen_window_pf_mertens_centered_bridge_vaughan_handoff.md
work/rh_compute/results/jensen_window_pf_mertens_centered_bridge_vaughan_handoff.json
python work/rh_compute/scripts/check_jensen_window_pf_mertens_centered_bridge_vaughan_handoff.py
outputs/jensen_window_pf_mertens_truncated_half_odd_kernel_handoff.md
work/rh_compute/results/jensen_window_pf_mertens_truncated_half_odd_kernel_handoff.json
python work/rh_compute/scripts/check_jensen_window_pf_mertens_truncated_half_odd_kernel_handoff.py
outputs/jensen_window_pf_mertens_shift_kernel_variation_handoff.md
work/rh_compute/results/jensen_window_pf_mertens_shift_kernel_variation_handoff.json
python work/rh_compute/scripts/check_jensen_window_pf_mertens_shift_kernel_variation_handoff.py
outputs/jensen_window_pf_mertens_planar_abel_handoff.md
work/rh_compute/results/jensen_window_pf_mertens_planar_abel_handoff.json
python work/rh_compute/scripts/check_jensen_window_pf_mertens_planar_abel_handoff.py
outputs/jensen_window_pf_mertens_planar_edge_gram_vaughan_handoff.md
work/rh_compute/results/jensen_window_pf_mertens_planar_edge_gram_vaughan_handoff.json
python work/rh_compute/scripts/check_jensen_window_pf_mertens_planar_edge_gram_vaughan_handoff.py
outputs/jensen_window_pf_mertens_planar_incidence_antidiagonal_handoff.md
work/rh_compute/results/jensen_window_pf_mertens_planar_incidence_antidiagonal_handoff.json
python work/rh_compute/scripts/check_jensen_window_pf_mertens_planar_incidence_antidiagonal_handoff.py
outputs/jensen_window_pf_mertens_planar_boundary_flux_anchor_handoff.md
work/rh_compute/results/jensen_window_pf_mertens_planar_boundary_flux_anchor_handoff.json
python work/rh_compute/scripts/check_jensen_window_pf_mertens_planar_boundary_flux_anchor_handoff.py
outputs/jensen_window_pf_mertens_planar_boundary_suffix_energy_handoff.md
work/rh_compute/results/jensen_window_pf_mertens_planar_boundary_suffix_energy_handoff.json
python work/rh_compute/scripts/check_jensen_window_pf_mertens_planar_boundary_suffix_energy_handoff.py
outputs/jensen_window_pf_mertens_planar_boundary_suffix_localization_gate.md
work/rh_compute/results/jensen_window_pf_mertens_planar_boundary_suffix_localization_gate.json
python work/rh_compute/scripts/check_jensen_window_pf_mertens_planar_boundary_suffix_localization_gate.py
outputs/jensen_window_pf_mertens_planar_joined_energy_equivalence_gate.md
work/rh_compute/results/jensen_window_pf_mertens_planar_joined_energy_equivalence_gate.json
python work/rh_compute/scripts/check_jensen_window_pf_mertens_planar_joined_energy_equivalence_gate.py
outputs/jensen_window_pf_mertens_planar_curvature_energy_scout.md
work/rh_compute/results/jensen_window_pf_mertens_planar_curvature_energy_scout.json
python work/rh_compute/scripts/check_jensen_window_pf_mertens_planar_curvature_energy_scout.py
outputs/jensen_window_pf_fixed_shift_inner_reciprocal_boundary_separation_gate.md
work/rh_compute/results/jensen_window_pf_fixed_shift_inner_reciprocal_boundary_separation_gate.json
python work/rh_compute/scripts/check_jensen_window_pf_fixed_shift_inner_reciprocal_boundary_separation_gate.py
outputs/sign_regularity_theorem_fit_matrix.md
work/rh_compute/results/proof_claim_ledger.json
python work/rh_compute/scripts/check_jensen_window_pf_bridge_target.py
python work/rh_compute/scripts/check_signed_hankel_jensen_dependency_graph.py
python work/rh_compute/scripts/check_jensen_window_pf_theorem_machinery_fit_matrix.py
python work/rh_compute/scripts/check_jensen_window_pf_structural_ansatz_matrix.py
python work/rh_compute/scripts/check_jensen_window_pf_coefficient_pf_equivalence_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_edrei_stieltjes_equivalence_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_edrei_heat_flow_boundary_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_phi_pick_kernel_target.py
python work/rh_compute/scripts/check_jensen_window_pf_xi_pick_suzuki_hankel_bridge.py
python work/rh_compute/scripts/check_jensen_window_pf_suzuki_spectral_frontier.py
python work/rh_compute/scripts/check_jensen_window_pf_suzuki_fixed_omega_phase_diagram.py
python work/rh_compute/scripts/check_jensen_window_pf_suzuki_jordan_totient_sign_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_suzuki_cofinal_monotonicity_hierarchy.py
python work/rh_compute/scripts/check_jensen_window_pf_suzuki_jordan_error_kernel_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_suzuki_cofinal_l2_hierarchy.py
python work/rh_compute/scripts/check_sign_regularity_theorem_fit_matrix.py
```

## Boundary

Passing this checker means the Jensen-window bridge has a reproducible
obligation decomposition. It does not prove any open obligation, and it does
not promote finite signed-Hankel, Jensen-window, Sturm, PF, or coefficient-PF
evidence into a proof of `Lambda <= 0`.
