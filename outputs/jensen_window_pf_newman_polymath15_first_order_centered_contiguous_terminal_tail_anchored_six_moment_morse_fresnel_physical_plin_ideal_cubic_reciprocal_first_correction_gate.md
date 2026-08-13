# Ideal-Cubic Reciprocal First-Correction Gate

Date: 2026-08-02

Status: exact reciprocal cells, dual positive Morse chart, returned leading carrier, and cancellation of the first full-line sequential stationary correction. This is not a proof of a finite-cell h^2 estimate, signed flow, or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_ideal_cubic_reciprocal_first_correction_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_ideal_cubic_reciprocal_first_correction_gate.py
```

## Reciprocal Cells

```text
C_q(alpha_P)={r in Z_(>0): q-1/2<=alpha_P/r<q+1/2}.
C_q(alpha_P)={floor(alpha_P/(q+1/2))+1,...,floor(alpha_P/(q-1/2))}.
The lower derivative boundary is included and the upper is excluded, so adjacent cells are disjoint and every equality alpha_P/r=q+1/2 transfers to C_(q+1).
The cells q>=1 partition every positive mode with alpha_P/r>=1/2; intersection with the fixed roster handles its at-most-one top edge.
```

## Dual Morse Chart

```text
For rho=qr/alpha_P, Phi_q(r)=alpha_P[log(alpha_P/r)-1]+qr=alpha_Plog q+alpha_P[rho-1-log rho].
zeta(rho)=sgn(rho-1)sqrt(2[rho-1-log rho]), eta=sqrt(alpha_P)zeta, so Phi_q=alpha_Plog q+eta^2/2.
dr/deta=[sqrt(alpha_P)/q]J(rho), with J(1)=1, J'(1)=2/3, and J''(1)=-5/18.
The original u chart has F_-=e(-1/8), the reciprocal r chart has F_+=e(+1/8), and F_-F_+=1.
At r_*=alpha_P/q, [qA(q)/sqrt(alpha_P)][sqrt(alpha_P)/q]=A(q), so the sequential leading stationary symbol is exactly e(alpha_Plog q)A(q).
```

## First Correction

```text
mathcal D_qA=q^2A''(q)+2qA'(q)+A(q)/6.
For b_r(y)=[sqrt(alpha_P)/r]A((alpha_P/r)v)J(v), b_(alpha_P/q)''(0)=q mathcal D_qA/alpha_P^(3/2).
The reciprocal leading amplitude after its Jacobian is d_q(eta)=J(rho)A(q/rho)/rho, with d_q''(0)=mathcal D_qA/alpha_P.
[sqrt(alpha_P)/q]b_(alpha_P/q)''(0)=d_q''(0)=mathcal D_qA/alpha_P.
For A(q)=exp(S(log q))P(log q), mathcal D_qA=exp(S)[P''+(2g+1)P'+(g^2+g+t/2+1/6)P](log q), the same removable-saddle operator as Section 11.173.
The positive reciprocal Gaussian contributes -kappa d_q''(0)/2, while the leading transform of the negative-chart c-prime correction contributes +kappa[sqrt(alpha_P)/q]b''(0)/2; their sum is exactly zero.
```

The first full-line sequential stationary correction cancels coefficientwise for an arbitrary twice-differentiable amplitude. For the physical ideal cubic it is exactly the cancellation of the reciprocal correction to mathcal F against the saddle symbol of mathcal J.

## Remaining Boundary

The proof uses the full-line Gaussian moment coefficients at one continuous reciprocal saddle. The exact discrete half-open cell has integer boundary transfers which are not bounded here.

The physical u integral has incomplete Fresnel endpoints. Their tails, the lower and upper endpoint packages, and boundary q=1,B cells must be recomposed before the full-line cancellation becomes a finite theorem.

In the complete-line interior symbol the order-alpha_P^(-1) correction vanishes, so the next formal interior term begins at second correction order. No O(alpha_P^(-2)) finite-cell remainder is claimed until the discrete boundaries and incomplete Fresnel tails are bounded.

## Handoff

Derive an exact finite reciprocal-cell identity with the two endpoint transfers exposed. Subtract the returned physical carrier and the cancelling first correction, then bound the remaining cell defect in a form that is summable across q without absolute family aggregation.

## Proof Boundary

This gate proves the exact half-open reciprocal cells, dual Morse phase and Jacobian, leading carrier return, equality of the two first stationary differential symbols, and their coefficientwise cancellation. It proves no finite reciprocal-cell remainder, incomplete-Fresnel or endpoint-complete h2 bound, q-family cancellation theorem, quadratic residual bound, signed flow estimate, Phi_B upper bound, contact exclusion, retained aggregate or Xi theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion.
