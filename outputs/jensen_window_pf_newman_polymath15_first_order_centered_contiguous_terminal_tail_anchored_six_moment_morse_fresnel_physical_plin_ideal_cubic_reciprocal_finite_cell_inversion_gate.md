# Ideal-Cubic Reciprocal Finite-Cell Inversion Gate

Date: 2026-08-02

Status: exact complete-mode reassembly, half-open cell Poisson formula, carrier Fourier inversion, and global internal-boundary cancellation. This is not a proof of a reciprocal exterior-tail bound, signed flow, or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_ideal_cubic_reciprocal_finite_cell_inversion_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_ideal_cubic_reciprocal_finite_cell_inversion_gate.py
```

## Complete Real-Mode Package

```text
F_A(u)=1_[1,B](u)A(u)e(alpha_Plog u), and I_A(r)=integral_1^B A(u)e(alpha_Plog u-ru)du=Fhat_A(r).
For real r>0, I_A(r)=U_A(r)+kappa e(-r)c_A(y_1(r))-kappa e(alpha_Plog B-rB)c_A(y_B(r))+J_A(r).
For integer r and integer B, the endpoint phases reduce to 1 and e_B, so I_A(r)=U_A(r)+kappa c_A(y_1)-kappa e_Bc_A(y_B)+J_A(r).
```

The real-r extension must use e(-r) and e(alpha_Plog B-rB) at the Morse endpoints. Replacing them off lattice by their integer values 1 and e_B changes the second Poisson transform.

## Exact Reciprocal Cell

```text
a_q=alpha_P/(q+1/2), b_q=alpha_P/(q-1/2), and C_q(alpha_P)={r in Z_(>0):a_q<r<=b_q}.
C_q(alpha_P)={floor(a_q)+1,...,floor(b_q)}.
b_(q+1)=a_q; the shared integer, when present, belongs to C_(q+1).
S_q[A]:=sum_(r in C_q)I_A(r)=tau_q[A]+lim_(R->infinity)sum_(s=-R)^R B_(q,s)[A].
tau_q[A]=(1/2)1_(b_q in Z)I_A(b_q)-(1/2)1_(a_q in Z)I_A(a_q).
B_(q,s)[A]=integral_(a_q)^(b_q)I_A(r)e(sr)dr.
```

## Dual Kernel And Carrier

```text
B_(q,s)[A]=integral_1^B F_A(u)K_q(u-s)du, where K_q(x)=kappa[e(-a_qx)-e(-b_qx)]/x and K_q(0)=b_q-a_q.
The s=q phase has r saddle alpha_P/q, value alpha_Plog q, and curvature q^2/alpha_P.
On [a_q,b_q], alpha_P/r lies in [q-1/2,q+1/2]; hence |s-alpha_P/r|>=1/2 for every integer s!=q.
omega_B(s)=1 for 1<s<B, 1/2 for s=1 or B, and 0 outside [1,B].
PV integral_R I_A(r)e(sr)dr=omega_B(s)A(s)e(alpha_Plog s).
E_q^ext[A]=PV integral_(R setminus [a_q,b_q])I_A(r)e(qr)dr, so B_(q,q)=omega_B(q)A(q)e(alpha_Plog q)-E_q^ext[A].
S_q[A]-omega_B(q)A(q)e(alpha_Plog q)=tau_q[A]+sum_(s!=q)B_(q,s)[A]-E_q^ext[A], with symmetric dual limits.
```

The full-line sequential transform is exact Fourier inversion, so it returns the starred physical carrier with no full-line remainder. The first-correction cancellation in Section 11.179 is its first local coefficient. Finite reciprocal cells differ only through explicit tie, alias, and exterior-tail terms.

## Global Recomposition

```text
The cells 1<=q<=B tile the continuous reciprocal band (alpha_P/(B+1/2),2alpha_P] and its integer modes. For alpha_-<=alpha_P<=alpha_+, every such mode lies in the fixed roster T={m,...,n}, m=max(1,floor(alpha_-/(2B))) and n=ceil(2alpha_+).
Internal tie terms cancel because a_q=b_(q+1); only the two outer reciprocal-band half terms remain.
sum_(q=1)^B S_q[A]=tau_band[A]+sum_(s=1)^B omega_B(s)A(s)e(alpha_Plog s)-sum_(s in Z)^sym E_s^band[A], where E_s^band is the Fourier-inversion exterior of [alpha_P/(B+1/2),2alpha_P].
If M_A=sum_(r in Z)^sym I_A(r) is the complete mode sum and C_band[A]=M_A-sum_(q=1)^B S_q[A], then sum_(s in Z)^sym E_s^band[A]=tau_band[A]+C_band[A].
All internal cell aliases and cellwise exterior tails recompose exactly; the only reciprocal-frequency boundary is the two-sided exterior of the complete physical band.
```

## Remaining Boundary

No uniform estimate is supplied for the noncentral aliases or the global reciprocal-band exterior. Endpoint jumps make the exterior Fourier inversion principal-value at q=1,B, so its two sides must not be bounded independently before endpoint composition.

## Handoff

Recompose the global lower and upper reciprocal-frequency exterior with the nonphysical roster fringes, mathcal T, the hard endpoint halves, mathcal L+e_Bmathcal U, and the terminal recurrence. Then derive a signed bound for that two-sided package, using the alias derivative gap only after internal cell cancellation.

## Pi Provenance

The kernel and inversion factors use the already fixed character e(x)=exp(2pi i x), so kappa=1/(2pi i). The endpoint half weights are the Dirichlet Fourier-inversion values. No circle, polygon, or fitted numerical constant is introduced.

## Proof Boundary

This gate proves the exact off-lattice endpoint phases, integer endpoint reassembly, reciprocal half-cell Poisson formula and tie correction, finite Fourier kernel, central carrier inversion, alias derivative gap, physical-band tiling, internal tie cancellation, and reduction of all cellwise defects to one global reciprocal-frequency exterior package. It proves no bound for that exterior, nonphysical roster fringes, endpoint-complete h2 remainder, q-family signed aggregate, quadratic residual, signed flow, Phi_B, contact exclusion, retained aggregate or Xi theorem, Q209, cofinal descendant theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion.
