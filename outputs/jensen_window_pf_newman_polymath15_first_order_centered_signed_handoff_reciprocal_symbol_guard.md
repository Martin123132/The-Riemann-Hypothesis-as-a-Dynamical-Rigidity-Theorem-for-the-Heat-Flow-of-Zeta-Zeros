# Newman Signed-Handoff Reciprocal-Symbol Guard

Date: 2026-07-28

Status: exact continuous leading-symbol and first-jet no-gain guard.
This is not a proof of a discrete B-process theorem, the Newman
boundary, PF-infinity, RH, or a prize-level result.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_signed_handoff_reciprocal_symbol_guard.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_signed_handoff_reciprocal_symbol_guard.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_signed_handoff_reciprocal_symbol_guard.py
```

Primary sources:

```text
https://arxiv.org/abs/1904.12438
https://arxiv.org/abs/2306.05599
```

## Exact Signed Handoff

```text
For w_t(u)=exp[t*log(u)^2/4], the exact zeta-handoff difference is S_N=sum_(n>=1)F_N(n)n^(-s_*), where F_N(u)=w_t(u)*1_(u<=N)-1. Equivalently, S_N=sum_(n<=N)(w_t(n)-1)n^(-s_*)-sum_(n>N)n^(-s_*).
```

## Pi Provenance

The 2*pi in a_omega^2=omega/(2*pi) comes from the standard Poisson phase exp(-2*pi*i*nu*u), and the pi/4 in the global factor is the negative-curvature stationary signature. These are the same completed-zeta and Riemann-Siegel constants used by the Polymath-15 approximation; no new geometric pi is introduced.

## Reciprocal Symbol

```text
Write s_*=sigma-i*omega, a_omega^2=omega/(2*pi), y=nu/a_omega, and u_nu=a_omega/y=a_omega^2/nu. For a positive Poisson mode nu, stationary phase gives the leading term Q(omega)*K_N(nu)*nu^(-conj(s_*)), where Q(omega)=exp(i[omega*(2log(a_omega)-1)-pi/4]) and K_N(nu)=y^(2sigma-1)F_N(u_nu).
On the interior reciprocal side nu>a_omega^2/N, so u_nu<N and F_N(u_nu)=w_t(u_nu)-1. Hence K_N=y^(2sigma-1)(w_t(u_nu)-1)=w_t(nu)*y^Delta*[1-w_t(u_nu)^(-1)], Delta=2sigma-1-t*log(a_omega). This factorization is exact.
After division by the dominant weighted reciprocal saddle, the signed weighted-minus-unweighted symbol is exactly Sigma_N(nu)=1-exp[-t*log(u_nu)^2/4]. It lies strictly between zero and one for t>0 and u_nu!=1. Thus it has neither a sign change nor a zero on the active positive-radius block.
```

## Active Radius

```text
For u_nu=N^r and c=tL with log(N)=L/2+o(1), exp[-t*log(u_nu)^2/4]=N^(-c*r^2/8+o(1)). At r=r_*=125662/155153, Sigma_N=1-N^(-c*r_*^2/8+o(1)). The relative defect therefore tends to zero while the signed symbol tends to one.
At c=2 the weighted reciprocal saddle at radius 2-r_* has N-exponent 8523351684/24072453409, while the subtracted ordinary reciprocal term has exponent 29491/155153. Their exact gap is 3947734561/24072453409=r_*^2/4=0.163993860281979..., so the subtraction is a relative power-small correction. It supplies zero exponent gain against the existing positive deficit d_2.
At c=c_* the relative exponent gap is 1989040967/9549356844=0.208290568620833..., again with the weighted reciprocal term dominant. The conclusion is not an artifact of testing only c=2.
```

```json
{
  "c2_relative_gap": {
    "decimal": "0.163993860281979",
    "exact": "3947734561/24072453409"
  },
  "c2_unweighted_dual_exponent": {
    "decimal": "0.190076891842246",
    "exact": "29491/155153"
  },
  "c2_weighted_dual_exponent": {
    "decimal": "0.354070752124225",
    "exact": "8523351684/24072453409"
  },
  "c_star": {
    "decimal": "2.540223984760009",
    "exact": "4911678521/1933561194"
  },
  "c_two_deficit": {
    "decimal": "0.065088263870695",
    "exact": "3133668399/48144906818"
  },
  "cstar_relative_gap": {
    "decimal": "0.208290568620833",
    "exact": "1989040967/9549356844"
  },
  "cutoff_separation_exponent": {
    "decimal": "0.190076891842246",
    "exact": "29491/155153"
  },
  "endpoint_separation_at_c2": {
    "decimal": "0.815088263870695",
    "exact": "78484697025/96289813636"
  },
  "r_star": {
    "decimal": "0.809923108157754",
    "exact": "125662/155153"
  },
  "r_star_dual": {
    "decimal": "1.190076891842246",
    "exact": "184644/155153"
  }
}
```

## First X Derivative

```text
On a fixed-t, fixed-N chart and for fixed mode nu, put A=log(a_omega), U=log(u_nu)=2A-log(nu), z=log(nu/a_omega), q=tU^2/4, and A_x=omega_x/(2omega). Then K_x/K=2sigma_x*z-(2sigma-1)A_x+tU*A_x*exp(q)/(exp(q)-1). This is the exact first-x derivative of the interior leading amplitude.
For T_nu=Q(omega)K_N(nu)nu^(-conj(s_*)), T_(nu,x)/T_nu=-sigma_x*U-(2sigma-1)A_x+tU*A_x*exp(q)/(exp(q)-1)+i*omega_x*U. The signed factor alone obeys partial_x log(Sigma_N)=tU*A_x/(exp(q)-1), which is power-small on the active block. The first jet therefore inherits the same no-gain verdict.
```

## Endpoint And Cutoff

```text
The reciprocal cutoff is nu=a_omega^2/N, asymptotically radius one. The active dual radius is 2-r_*=184644/155153 and its exact power separation from the cutoff is 1-r_*=29491/155153=0.190076891842246.... The floor seam and adjacent recurrence cannot be identified with this interior block.
The certified first correction has T|d_n|<100 with T asymptotic to pi*N^2, so it is relative O(N^-2). The Polymath-15 endpoint and one adjacent block have scale N^(-1/2-c/8+o(1)); at c=2 this is N^(-3/4+o(1)), separated from the current N^(d_2+o(1)) active-block envelope by the exact exponent 78484697025/96289813636=0.815088263870695.... First derivatives add logarithms only. These terms remain mandatory at the cutoff, but they cannot supply a uniform leading-symbol cancellation at r_*.
```

## Route Decision

```text
The proposed signed reciprocal mechanism is falsified at its first required test. On the active block the exact normalized symbol tends to one, its x derivative is a power-small perturbation, the normalizer is conjugate-locked, and the endpoint/correction scales do not reverse the leading symbol. No N^(-d_2-eta) pairwise gain is produced.
Retire signed primal/tail pairwise cancellation as a route to lower c_*. The reciprocal coordinate may still reorganize a dual exponential sum, but any improvement now requires a new arithmetic estimate inside that dual sum, such as bilinear or additive-energy cancellation. The primary proof programme returns to the direct Xi Abel-phase/contact theorem unless such an input is derived independently.
```

## Boundary

This artifact proves the exact continuous signed reciprocal symbol, its fixed-chart first-x derivative, the active rational exponent audit, and a no-gain verdict for signed pairwise cancellation. It does not prove a discrete endpoint-complete Poisson/B-process remainder theorem or rule out new arithmetic cancellation inside the dual sum. The effective L_epsilon, Abel-scalar gap, inner degree theorem, contact exclusion, final Newman conclusion, PF-infinity theorem, RH proof, and prize-level conclusion remain open.
