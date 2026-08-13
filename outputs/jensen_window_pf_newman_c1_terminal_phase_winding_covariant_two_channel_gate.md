# Newman C1 Terminal-Phase Winding and Covariant Two-Channel Gate

Date: 2026-08-04

Status: complete fixed-N terminal-phase winding and the exact two-channel source-phase matrix proved; joined physical-current inequality open; not a proof of RH.

## Complete-cell phase winding

On q=1, a(L)^2=exp(L)+1/(32L^2), T_0=2pi a(L)^2, Omega=T_0-epsilon, and 0<epsilon<7/(8L^2x), with x=4pi exp(L).

The derivative of a(L)^2 is exp(L)-1/(16L^3)>0 for L>=50. Thus a is continuous, strictly increasing, and unbounded.

```text
Delta phi_N
 =[2pi(2N+1)+epsilon(L_N)-epsilon(L_(N+1))]log N
 >4pi N log N
 >2pi.
```

A continuous real phase whose endpoint increment exceeds 2pi assumes a representative of every residue class modulo 2pi in the cell interior: choose k=floor((phi_left-theta)/(2pi))+1, so phi_left<theta+2pi k<=phi_left+2pi<phi_right. Therefore exp(i phi_N) covers the whole unit circle independently of the half-open endpoint convention.

No fixed proper phase sector contains the terminal phase on a complete fixed-N cell. In particular sin(phi_N) takes both signs and vanishes, so the raw negative anchors in Section 11.213 cannot have a cell-uniform physical projection by a fixed terminal sector.

## Covariant two-channel form

Let h=Re(conjugate(tau_0)mathcal J_H^0), t=tau_0 mathcal J_T^0, and z=nu c_xi. The ideal joined mixed current is R_id(z)=(1/2){|z|^2 h+Re(z^2 t)}.

```text
Q(h,t)=(1/2)[[h+Re(t), -Im(t)],
             [-Im(t), h-Re(t)]].

eigenvalues(Q)=(h-|t|)/2, (h+|t|)/2.
```

For z=r exp(i gamma), R_id=(r^2/2){h+Re(exp(2i gamma)t)}. Negativity for the actual source phase is exactly h+Re(exp(2i gamma)t)<0.

Negativity for every nonzero source phase is equivalent to negative definiteness of Q, hence to the single phase-invariant inequality h+|t|<0. Nonpositivity is equivalent to h+|t|<=0.

Because every complete cell contains phi_N=3pi/2 modulo 2pi, the centre-only phase-uniform margin h_0+|t_0| equals rho(A_H+A_T)>0 there. Static negative centres cannot prove a source-phase-uniform current theorem on a complete cell.

The phase-uniform criterion is stronger than the theorem at one known physical source phase. Its failure for the bare centres does not disprove the complete current: the signed reciprocal residuals, physical corrections, fixed affine row, and moved-tail defect have not been estimated.

## Gate rows

| id | role | state | claim |
|---|---|---|---|
| tpw_01_q1_law | physical phase law | proved | The q=1 cutoff and physical frequency have an exact one-sided epsilon law. |
| tpw_02_monotone | cutoff monotonicity | proved | The saddle cutoff is strictly increasing for L>=50. |
| tpw_03_cells | complete cells | proved | Every sufficiently large integer cutoff has a unique complete L-cell. |
| tpw_04_phase | terminal phase | proved | The physical terminal phase is Omega log N on a fixed-N cell. |
| tpw_05_endpoint | endpoint identity | proved | The exact cell endpoint phase increment is explicit. |
| tpw_06_epsilon | epsilon guard | proved | The one-sided physical correction cannot remove the macroscopic phase increment. |
| tpw_07_turn | full-turn theorem | proved | Every complete fixed-N cell gains more than one full terminal turn. |
| tpw_08_circle | circle coverage | proved | The terminal phase covers every direction on every complete cell. |
| tpw_09_sector | sector obstruction | proved | A fixed proper terminal phase sector cannot be uniform on a complete cell. |
| tpw_10_current | two-channel current | proved | The ideal joined current has one Hermitian scalar and one transpose spin-two scalar. |
| tpw_11_matrix | covariant matrix | proved | The source phase is represented by an exact real symmetric 2x2 quadratic form. |
| tpw_12_spectrum | matrix spectrum | proved | The two invariant eigenvalues are determined by h and |t|. |
| tpw_13_actual | actual-phase criterion | proved | The physical source phase has an exact scalar sign criterion. |
| tpw_14_uniform | phase-uniform criterion | proved | A source-phase-uniform theorem is equivalent to h+|t|<0. |
| tpw_15_centres | centre-only obstruction | proved | The raw anchor centres cannot make the covariant matrix negative on a complete cell. |
| tpw_16_target | covariant current target | open | The live target is a joined phase-correlated Hermitian/transpose current estimate. |

## Next action

Work with the complete physical pair h=Re(conjugate(tau_0)mathcal J_H) and t=tau_0 mathcal J_T, including the correction polynomials, fixed affine row, and moved-tail defect before projection. First derive the exact actual-source-phase scalar h+Re(exp(2i gamma_phys)t) in the contact-box normalization. Then test whether reciprocal stationary recomposition can prove its signed reserve. Use h+|t| only if a genuinely source-phase-uniform theorem is sought; do not estimate the two channels or their anchor disks independently.

## Pi provenance

The factor 2pi in T_0=2pi a^2 and alpha_P=xi/(2pi) comes from the completed-zeta saddle and the Fourier character e(y)=exp(2pi i y). The threshold 2pi is one period of exp(i phi_N). The rational bounds pi>3 and pi<22/7 only certify a deliberately coarse full-turn inequality. No circle, polygon, or fitted plotting constant is inserted into the model.

## Proof boundary

This gate proves complete-cell winding of the physical q=1 terminal phase, failure of any fixed proper terminal phase sector on such a cell, the exact Hermitian/transpose source-phase quadratic matrix, its invariant eigenvalues, and a centre-only obstruction to source-phase-uniform negativity. It proves no signed reciprocal residual estimate, actual-source-phase current inequality, correction-row bound, all-q transport, contact-box exclusion, Q209 shell, cofinal successor, Lambda<=0, RH, or prize-level conclusion.

built Newman C1 terminal-phase winding/covariant two-channel gate: 16 rows, 0 issues, 1 full-cell winding theorem, 1 circle-surjectivity theorem, 1 phase-uniform criterion, 1 centre-only obstruction, 1 live covariant current target
