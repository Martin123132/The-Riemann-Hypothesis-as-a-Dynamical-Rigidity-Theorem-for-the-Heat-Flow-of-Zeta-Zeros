# Jensen-Window PF Newman Strict-Laguerre Correlation Target

Date: 2026-07-25

Status: exact strict-Laguerre/Wiener equivalence with a generic-kernel
guard. This is not a proof of RH or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_newman_strict_laguerre_correlation_target.json
python work/rh_compute/scripts/jensen_window_pf_newman_strict_laguerre_correlation_target.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_strict_laguerre_correlation_target.py
```

Current result:

```text
validated Jensen-window PF Newman strict-Laguerre correlation target: 10 rows, 0 issues, 1 strict-Laguerre equivalence, 1 exact correlation identity, 1 Wiener-density equivalence, 1 RH-equivalent density target, 1 exact strict-log-concavity/positive-definiteness countermodel, 2 non-promotion gates, 1 open Xi handoff
```

## Strict Laguerre Target

Set

```text
phi_t(u)=exp(t*u^2)*Phi(u), u in R
H_t(x)=integral_0^infinity phi_t(u)*cos(x*u)du
L_t(x)=H_t'(x)^2-H_t(x)*H_t''(x)
K_(1,t)(v)=integral_R phi_t(s+v)*phi_t(s-v)*s^2 ds
```

The positive-boundary attainment theorem turns the first Laguerre
inequality into an exact endgame:

```text
Lambda<=0 if and only if L_t(x)>0 for every real x and every 0<t<=1/5. The forward implication uses the simple real-zero factorization for t>Lambda; the reverse implication uses the finite multiple zero of H_Lambda forced when Lambda>0 together with the published bound Lambda<=1/5.
```

This is weaker in form than proving every `H_t` is Laguerre-Polya:
only one strict real-axis inequality is requested, but it must hold
for the complete positive-time continuum.

## Correlation Identity

The midpoint/difference change of variables gives exactly

```text
L_t(x)=integral_R K_(1,t)(v)*cos(2*x*v)dv; equivalently, Fourier[K_(1,t)](xi)=L_t(xi/2).
```

There is no missing normalization factor. Wiener's theorem now yields

```text
Translations of K_(1,t) are dense iff L_t has no real zero. Since L_t(0)>0, zero-freeness is equivalent to L_t(x)>0 for all x.
Lambda<=0 if and only if, for every 0<t<=1/5, the translations of K_(1,t) are dense in L1(R).
```

Primary sources: https://arxiv.org/abs/1309.0055 and
https://arxiv.org/abs/1606.05011.

At a hypothetical positive boundary the correlation kernel is already
positive definite, but its transform touches zero at the multiple root.
The missing property is therefore zero-freeness or translate density,
not ordinary nonnegativity.

## Exact Shape Guard

Let

```text
T_a(y)=(a-|y|)_+ = 1_[-a/2,a/2]*1_[-a/2,a/2]
G_sigma(x)=exp(-x^2/(2*sigma^2))
K_(a,sigma)(x)=integral_(-a)^a (a-|y|)*exp(-(x-y)^2/(2*sigma^2))dy
```

For `a=1`, `sigma=2`, conditional differentiation gives

```text
With p_x(y) proportional to T_a(y)G_sigma(x-y), (log K)''=Var_x(Y)/sigma^4-1/sigma^2
Since Y lies in [-a,a], Var_x(Y)<=a^2; for sigma>a, (log K)''<=(a^2-sigma^2)/sigma^4<0.
(log K)''<=-3/16
```

so `K` is smooth, positive, even, and strictly log-concave. Its
Fourier transform is nonnegative, hence `K` is positive definite, but

```text
8*sqrt(2)*sqrt(pi)*exp(-2*xi**2)*sin(xi/2)**2/xi**2
xi=2*pi*n/a for every nonzero integer n
```

The zeros are double, so the translates of `K` are not dense. This
blocks promotion from strict log-concavity plus positive definiteness.
This Gaussian model alone does not reproduce the Xi-specific tail.
The separate theta-tail weighted countermodel closes promotion from
uniform strong and root-variable log-concavity plus theta-type decay;
it still leaves Xi arithmetic or modular correlation structure open.

## Live Handoff

Prove uniformly for every 0<t<=1/5 that Fourier[K_(1,t)] has no real zero, equivalently that translations of K_(1,t) are dense in L1(R). A viable proof must use Xi-specific structure beyond generic strict log-concavity and positive definiteness. The exact higher-shift expansion must be grouped through the theta modular identity before the spectral transform: its ordinary termwise Bessel-transform sum diverges already at zero frequency. Continuum-cell subtraction repairs the endpoint transform and gives a normally convergent coupled matrix at t=0, but each cell block retains a nonzero exp(-5u) tail and cannot be deformed to any t>0. The theta-curvature probability/operator gate now proves that the full positive primitive S_t=exp(tu^2)R is strictly decreasing and convex, that dmu_t=2S_t''du is a probability, and that its transform C_t is strictly positive. Its fixed theta weights satisfy w_1>2799/2800 with uniform C2 omitted-component budgets, isolating a dominant first block on compact frequency; a supplementary C3 bound makes the corresponding J_t' error explicit. The resulting B_0/B_1 disjunction is a rigorous pointwise exclusion test. An exact rational origin estimate and a 160-bit Arb/Taylor certificate with 1900 rational leaf boxes now promote it to (H_t,H_t')!=(0,0) for every |x|<=38 and 0<=t<=1/5. Its raw moment bars still grow polynomially and leave oscillatory phase control necessary for |x|>38. That remaining range requires oscillatory tail control and modular endpoint cancellation inside the corrected high-frequency partition. It reduces multiple contact exactly to J_t(x)=16*x^4*H_t(x) and J_t'(x) vanishing together, while a smooth Neumann-lift guard proves that primitive-transform positivity alone is insufficient. Its ratio G_t=8H_t/C_t has a weighted Doob diffusion with nonnegative instantaneous generator, but a quadratic nodal-loss guard shows that operator positivity and Sturm monotonicity still permit contact; modular endpoint asymptotics also make C_t eventually log-convex, blocking a global contracting-drift shortcut. A separate quadratic backward-heat audit shows that top-time Pick positivity, collision-free speed bounds, and cutoff-weighted energy continuity cannot be bridged through a collision without separation- and cutoff-uniform Xi input. The cancellation-preserving modular partition now gives an exact direct-C1 finite-to-infinite contract. Explicit heat-polynomial envelopes and the reflected-phase minimum prove effective C*N^A*exp(-c*N^(4/3)) derivative tails, replacing the unsupported fixed additive saddle collar by N=max(N_sad,ceil(K*(1+x)^(3/4))). A weighted Cauchy-Schwarz reduction now converts each nonsmooth L1 term to an analytic quadratic integral. The cancellation-stable forward-remainder identity, explicit arithmetic and outer-u tails, and a complete 223-entry Arb matrix now give full m=9 derivative budgets for N=4 through 10. The direct-tail explicit-constant successor now gives cap-free m=9 d0/d1 upper budgets for every N>=2. With N(x)=ceil((1+x)^(3/4)), its adaptive monotonic successor proves both omitted errors below the quarter-gamma scale exp(-pi*x/8)/4 for every x>=245. A 192-bit N=7 retained cover then certifies all 1240 initial boxes on [1/155,1/5]x[38,69], proving that Q_31=[1/155,1/4]x[0,69] has no first-jet contact and has zero winding. Thus Q_1 through Q_31 are rigorous by containment. On the bounded bridge, the ordinary positive-half-line theta series is cheaper: a directed direct-L1 theorem shows that six ordinary terms leave both J and J' errors below 10^-11*exp(-pi*x/8) on 38<=x<=245. A 352-bit shared Taylor certificate then validates all 414 half-unit panels and all 7102 time-frequency leaves in the missing low-time and right-hand strips, with zero subdivisions or unresolved leaves. Combined with Q31, this proves no contact and zero winding on Q_207=[1/1035,1/4]x[0,245]. Thus Q_1 through Q_207 are rigorous by containment. The same ordinary-theta tail formula now also has a cofinal saddle-scale successor: for K_h=ceil(sqrt(x/8+2log(1+x)+h)) and h>=0, retaining N_(F,h)=K_h-1 ordinary terms leaves both first-jet errors below exp(-3h-pi*x/8)/50 on x>=245. For fixed h this count is only sqrt(pi/2) times the Riemann-Siegel saddle count asymptotically, and h can be increased to follow any prescribed positive error tolerance. Thus the earlier x^(3/4) modular count remains valid but is no longer the efficient cofinal error handoff. What remains is retained first-jet separation for this saddle-scale ordinary main, with a terminating analytic margin for x>=245 and 0<t<=1/5. On the critical overlap L>=50, tL<=25, the exact normalizer transfer gives the direct sufficient target T_L[O_h]>2*exp(-6h)/(625*x^(23/2)). It also gives the algebraic identity Z_t=O_h+e_F=J_hat+r_RS and a scaled first-jet comparison below 5600*exp(-3L/4). This does not make the current quantitative targets equivalent: the certified ordinary envelope is O(exp(-3h-23L/4)), while the corrected-Riemann-Siegel envelope is O(exp(-3L/4)), a certified-envelope cost O(exp(5L+3h)). One must therefore prove either the tiny direct ordinary lower bound or the existing corrected target T_L[J_hat]>32000000*exp(-3L/2), without silently transferring the smaller threshold. An exact one-sided score bridge now defines dnu_t=-f_t'/f_t(0), writes F_t=H_t+iY_t=i*f_t(0)*(1-chi_t)/x with Y_t>0, and converts multiple contact to the joint characteristic equations B_t=E_nu[sin(xU)]=0 and B_t'=E_nu[U*cos(xU)]=0. The same score law gives the signed-Hankel coefficients as normalized odd moments and every shifted Jensen window as an explicit positive generalized-Laguerre scale mixture. A uniform/Beta pushforward makes this one positive Abel measure with moments M_k=k!*A_k, Bessel transform 2*H_t(sqrt(-z)), and a closed tail-source flow; its first nontrivial Jensen condition is a lower concentration threshold not supplied by positivity. Full scale support rules out a common interlacer for the complete component family in degree at least two, and opposite-sign Abel minors reject the bare fractional kernel as a variation-diminishing shortcut. The same Abel density now gives an exact radial probability in dimension 2n+2 for every coefficient shift n, with normalized Fourier profile mathcal_F_t^(n)(-|xi|^2)/A_n and the classical dimension walk. A two-Gaussian Newman-flow countermodel retains this entire radial ladder and arbitrarily strong log-concavity while failing every shifted quadratic Jensen test. Therefore own-dimension positive definiteness and dimension walking do not close the coefficient branch. The separate Xi squared-variable log-concavity and reciprocal-defect heat theorems close degrees two and three. The Xi real-zero-band/sector theorem now closes every shifted Jensen layer through degree 361 uniformly on 0<=t<=1/5. The first unproved direct layer is degree 362, but a bounded successor is not enough: an unbounded cofinal sequence or an Xi-specific all-degree input must enter beyond that point. Theta arithmetic and the Xi double-exponential tail remain natural sources for such an input. Its sine observable obeys a closed radial backward-heat equation. On the corrected overlap, restoring the exact lift amplitude makes its first jet uniformly equivalent to T_L[Z_t], with squared constants (9-sqrt(17))/8 and (9+sqrt(17))/8. The raw score phase is rapidly asymptotically flat, and the existing smoothed-triangle countermodel has a score-characteristic double zero, so neither an absolute phase-speed floor nor generic strong log-concavity closes the target. The primary surviving route is an arithmetic C1 separation theorem for this explicit Xi characteristic expression on |x|>38, using the existing corrected Riemann-Siegel partition and a phase-critical-value avoidance theorem on its residual scaled layer, or a coupled first-correlation square. In that partition, |x|<=38 and the dominant ray L>=50, tL>=25 are closed; every fixed tL>=c_*+epsilon is closed only above an existential L_epsilon. The bounded-L band beyond the closed x<=245 core, the finite L_epsilon shoulders, and L>=50 with 0<tL<=c_*+o(1) remain open. The former sharper sufficient subtarget M_t(x)=-L_t'(x)>0 has now been rejected for Xi: Arb certifies M_0(1401016343/100000)<0, and continuity preserves that sign for sufficiently small positive t. The surviving routes are direct zero-freeness of Fourier[K_(1,t)] or corrected C1 double-zero transversality; do not impose global monotonicity of L_t.

The exact route refinement is recorded in
`outputs/jensen_window_pf_newman_theta_bessel_higher_shift_regularization_gate.md`.
The convergent endpoint matrix and its positive-time obstruction are in
`outputs/jensen_window_pf_newman_theta_cell_renormalization_gate.md`.
The positive primitive normalizer and characteristic `C1` contact
reduction are in
`outputs/jensen_window_pf_newman_theta_curvature_probability_operator_gate.md`.
The exact origin collar, diagnostic route map, and rigorous 1,900-box
promotion through `|x|=38` are in
`outputs/jensen_window_pf_newman_theta_compact_transversality_scout.md`
and
`outputs/jensen_window_pf_newman_theta_compact_transversality_interval_certificate.md`.
The cancellation-preserving direct-C1 contract, exact derivative-tail
scale, subordinate fixed-collar stress, and rigorous finite compact
quadratic pilot are in
`outputs/jensen_window_pf_newman_theta_adaptive_modular_c1_remainder_contract.md`,
`outputs/jensen_window_pf_newman_theta_modular_tail_derivative_envelope_gate.md`,
`outputs/jensen_window_pf_newman_theta_modular_tail_derivative_budget_scout.md`,
and
`outputs/jensen_window_pf_newman_theta_modular_tail_arb_quadratic_pilot.md`.
The stable forward-remainder theorem, arbitrary-N direct-tail
compiler, cofinal quarter-gamma tail theorem, six-term finite-bridge
tail and interval theorems, tunable cofinal square-root tail theorem,
exact ordinary-to-corrected-Riemann-Siegel C1 transfer guard,
full compact and outer derivative budgets,
and the completed finite Q31 retained cover are in
`outputs/jensen_window_pf_newman_theta_forward_remainder_tail_gate.md`,
`outputs/jensen_window_pf_newman_theta_switch_defect_explicit_constant_gate.md`,
`outputs/jensen_window_pf_newman_theta_adaptive_gamma_scale_tail_gate.md`,
`outputs/jensen_window_pf_newman_theta_forward_six_term_bridge_tail_gate.md`,
`outputs/jensen_window_pf_newman_theta_forward_six_term_finite_bridge_interval_certificate.md`,
`outputs/jensen_window_pf_newman_theta_forward_adaptive_sqrt_tail_gate.md`,
`outputs/jensen_window_pf_newman_theta_forward_sqrt_to_corrected_rs_C1_transfer_gate.md`,
`outputs/jensen_window_pf_newman_one_sided_phase_moment_bridge_gate.md`,
`outputs/jensen_window_pf_newman_score_abel_radial_dimension_lift_gate.md`,
`outputs/jensen_window_pf_newman_zero_slab_degree71_sector_certificate.md`,
`outputs/jensen_window_pf_newman_real_zero_band_degree361_sector_certificate.md`,
`outputs/jensen_window_pf_newman_theta_stable_remainder_outer_tail_gate.md`,
`outputs/jensen_window_pf_newman_theta_stable_remainder_arb_quadratic_matrix.md`,
`outputs/jensen_window_pf_newman_theta_full_derivative_budget_certificate.md`,
and
`outputs/jensen_window_pf_newman_theta_modular_retained_q31_interval_certificate.md`.
The rejected monotonicity subtarget, its misleading finite diagnostics,
and the rigorous Xi counterexample are in
`outputs/jensen_window_pf_newman_strict_laguerre_monotonicity_scout.md`.
