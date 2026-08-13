# Six-Moment Morse-Fresnel Endpoint-Composition Retention Gate

Date: 2026-08-02

Status: six exact endpoint-composed value identities, exact roster transfer,
and retained endpoint floors; `0 oscillatory interior bounds`, `0 signed flow
bounds`, and this is not a proof of RH.

## Exact Composition

Retain q=1, j=0,...,5, the fixed roster T={m,...,n}, kappa=1/(2pi i), e_B=e(alpha log B), and all notation of Sections 11.167-11.170.

Set F_j=sum_(r in T)U_(j,r), J_j=kappa sum_(r in T)e(phi_r(alpha/r)) integral c'_(j,r)e(-y^2/2), and let T_j be the twice-integrated outside-roster integral tail.

Define

```text
L_j=delta_(j0)/2-calB_j(1)+kappa sum_(r in T)c_(j,r)(y_1).

U_j=A_j(B)/2+e(-alpha log B)calB_j(B)-kappa sum_(r in T)c_(j,r)(y_B).
```

Direct substitution into finite starred Poisson summation gives, without any
absolute value,

```text
H_j=F_j+L_j+e(alpha log B)U_j+J_j+T_j, 0<=j<=5.
```

The stable endpoint formulas are

```text
With d_j=A_j'(1)=-sigma delta_(j0)+delta_(j1), L_j=delta_(j0)/2-kappa delta_(j0)S_1(alpha)+kappa^2[d_jS_2(alpha)+alpha delta_(j0)S_3(alpha)]+kappa sum_T c_(j,r)(y_1).

For x_B=alpha/B, U_j=A_j(B)/2+kappa A_j(B)S_1(x_B)-kappa^2[A_j'(B)S_2(x_B)+(alpha/B^2)A_j(B)S_3(x_B)]-kappa sum_T c_(j,r)(y_B).
```

For lambda=log B, A_j'(B)=B^(-1)exp(t lambda^2/4-sigma lambda){j lambda^(j-1)+(t lambda/2-sigma)lambda^j}, with the first term absent for j=0.

## Endpoint Retention

All S_k and endpoint c values are real. Thus only the physical half endpoint and kappa^2 terms contribute to Re L_j and Re U_j.

|Re L_0-1/2|<2603/(7200alpha)<1/10000, hence |L_0|>4999/10000.

For every 0<=j<=5, |Re U_j-A_j(B)/2|<15653 A_j(B)/(1800alpha)<A_j(B)/100, hence |U_j|>49A_j(B)/100.

The displayed floors use only alpha>4096, log(B)>1/10, |g'|<=51/100, pi>3, Z_2<=6u/alpha, and Z_3<=10u^2/alpha^2; the physical chart has alpha>2^49.

Thus endpoint composition does not make the hard endpoints small. They must
be retained in the comparison carrier.

## Roster Audit

At a physical endpoint mu, q_r(mu)=alpha/mu-r and the exact Morse identity gives b_(j,r)(y_mu)/y_mu=-A_j(mu)/q_r(mu). Hence c_(j,r)(y_mu)=-A_j(mu)/q_r(mu)-b_(j,r)(0)/y_mu.

S_1(x)+sum_(r in T)1/(x-r)=pi cot(pi x). Substitution gives an equivalent full-cotangent form, but its cotangent pole and the b(0)/y pole cancel only after recombination. The stable digamma/Hurwitz package must be used at saddle crossings.

1/2-kappa pi cot(pi x)=1/[1-e(x)] and 1/2+kappa pi cot(pi x)=1/[1-e(-x)]. These are audit identities, not separately bounded endpoint terms.

For any uniformly nonstationary mode moved across the roster boundary, its exact Morse representation U+Q equals its exact twice-integrated representation: I=U+kappa c(y_1)-kappa e_B c(y_B)+kappa e(phi_0)integral c'e(-y^2/2)=kappa[e(phi)A/q]_1^B-kappa^2[e(phi)(D_rA)/q]_1^B+kappa^2 integral e(phi)C_rA.

Therefore the complete decomposition of H_j is roster-invariant. The separate packages L_j and U_j are bookkeeping-dependent and must not be interpreted without F_j, J_j, and T_j.

## Carrier And Remainder

Define P_j^ec=F_j+L_j+e(alpha log B)U_j and rho_j^osc=J_j+T_j. Then H_j=P_j^ec+rho_j^osc exactly.

The outside tail satisfies |T_j|<2K_jh^2 with K=(6,14,55,336,2738,27936). The grouped interior J_j remains unbounded at h^2 scale.

## Eight-Observation Composition

After the fixed source normalization and common c_xi phase, let y=p+r be the twelve-real moment vector induced respectively by P^ec and rho^osc.

Using Psi'=y^T M_xi y+ell_xi^T y and ell_xi=U_1t_T, one has exactly Psi'=p^T M_xi p+ell_xi^T p+2p^T M_xi r+r^T M_xi r+ell_xi^T r.

The fixed terminal coordinates enter through t_T and hence through ell_xi in this identity. They do not belong inside the scalar Poisson endpoint package and no terminal cancellation is assumed before the observation map.

## Next Target

Insert the explicit carrier p into the eight observations and retain its endpoint, Hermitian, transpose, and terminal terms. Then prove a coefficient-aware bound for the five terms involving r, starting with grouped control of J_0,...,J_5; do not demand an h^2 bound for the retained endpoint carrier itself.

## Pi Provenance

`kappa=1/(2pi i)` is forced by `e(x)=exp(2pi i x)`, while
`pi cot(pi x)` is the symmetric partial-fraction sum over integral Poisson
modes. No circle or fitted geometric constant is introduced.

## Proof Boundary

This proves six exact stable endpoint-composed identities, exact one-mode roster transfer, quantitative retention of the lower order-zero and all six factored upper endpoint packages, and the exact carrier/remainder expansion through the eight-observation flow including terminal coordinates. It does not bound the grouped oscillatory interior, prove an endpoint-composed h^2 remainder, evaluate the retained carrier's signed current, prove a signed flow or Phi_B bound, exclude contact, establish a retained aggregate or Xi theorem, prove Q209, Lambda<=0, PF-infinity, RH, or a prize-level conclusion.

## Reproduce

```text
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_endpoint_composition_retention_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_endpoint_composition_retention_gate.py
```
