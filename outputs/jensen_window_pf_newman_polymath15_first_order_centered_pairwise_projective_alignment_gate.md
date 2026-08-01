# Pairwise Projective-Alignment Gate

Date: 2026-07-28

Status: exact alignment and pole identities with two route guards;
not a proof of the signed Xi occupation bound, contact exclusion,
`Lambda<=0`, PF-infinity, RH, or a Clay-prize result.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_pairwise_projective_alignment_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_pairwise_projective_alignment_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_pairwise_projective_alignment_gate.py
```

## Pi Provenance

The pi in a^2=x/(4*pi)+t/16 and u_x=1/(8*pi*a^2) is inherited from the completed-zeta factor pi^(-s/2)*Gamma(s/2)*zeta(s) and its Riemann-Siegel saddle. It is not fitted from a circle, prime-curvature arc, polygon, or plot.

## Exact Alignment

Fix one N=floor(a) chart on L>=50, 0<tL<=25, and q=2*t*L^2>=1. For 1<=n<m<=N-1, put u_i=log(a/i), k_i=u_i-u_N=log(N/i), c_s=Re(s_*'), b=Im(s_*')<0, and h_i=c_s*k_i-b*u_i*tan(theta_i) whenever X_i!=0.

Here e_i=v_a+Im[d_(i,x)/(1+d_i)], so nu_i=b*u_i+e_i and |e_i|<8449/x^2.

```text
At h_i=H put w_i=c_s*k_i-H and e_i=nu_i-b*u_i. Then h_(i,x)=-b^2*u_i^2-w_i^2+c_(s,x)*k_i-(b_x/b+u_x/u_i)w_i-b*u_i*e_i-e_i*w_i^2/(b*u_i).
For 1<=n<m<=N-1 put delta=u_n-u_m=k_n-k_m>0, S=u_n+u_m, and A_i=b^2*u_i^2+w_i^2. At h_n=h_m=H, D_nm(H):=h_(m,x)-h_(n,x) equals delta*[b^2*S+c_s*(c_s*(k_n+k_m)-2H)-c_(s,x)+c_s*b_x/b+u_x*(c_s*u_N+H)/(u_n*u_m)]+e_n*A_n/(b*u_n)-e_m*A_m/(b*u_m).
```

## Frozen Model

If c_s=c_(s,x)=b_x=u_x=e_n=e_m=0, then D_nm=b^2*(u_n^2-u_m^2)=b^2*delta*S>0. This is the favorable frozen-scale correction-free Sturm orientation.

## Moving-Scale Obstruction

```text
Restoring only the exact moving scale u_x>0 gives D_nm=delta*[b^2*S+u_x*H/(u_n*u_m)].
Put H_mov=-b^2*S*u_n*u_m/u_x<0. Then D_nm(0)=b^2*delta*S>0, D_nm(H_mov)=0, and D_nm(2H_mov)=-b^2*delta*S<0. In that same reduced model h_(i,x)=-b^2*u_i^2-H^2+u_x*H/u_i. It is strictly negative at H=0 and at H=2H_mov<0 for both atoms. Thus strict individual clockwise motion and the actual logarithmic distance order do not fix pairwise crossing orientation once u_x is restored.
```

## Projection Pole

At a finite projective alignment, psi_(m,x)-psi_(n,x)=D_nm(H)/(1+H^2).

```text
At a common projection pole X_n=X_m=0, division-free currents give psi_(i,x)=-nu_i/(b*u_i)=-1-e_i/(b*u_i), and hence psi_(m,x)-psi_(n,x)=e_n/(b*u_n)-e_m/(b*u_m).
The two finite-chart limits H->+infinity and H->-infinity of D_nm(H)/(1+H^2) both equal the exact pole value e_n/(b*u_n)-e_m/(b*u_m).
```

Choose any 0<epsilon<min(8449/x^2,-b*delta,-b*u_m). The two hypothetical assignments (e_n,e_m)=(epsilon,0) and (0,epsilon) preserve nu_m-nu_n>0 and nu_n,nu_m<0. At the pole their relative rates are respectively epsilon/(b*u_n)<0 and -epsilon/(b*u_m)>0. Therefore the currently proved residual bounds and strict absolute angular-rate order do not determine the pole sign. These assignments are source-bound guards, not asserted Xi correction values.

## Legacy Symmetry Boundary

The reconstructed prime-curvature matrix is exactly radial and transpose-symmetric by construction, while every tested raw symmetric kernel is indefinite. That static p<->q symmetry supplies neither the moving-scale term u_x*H/(u_n*u_m) nor an ordering of e_i/u_i. The pictures remain a useful falsification laboratory, but they cannot repair either pairwise obstruction without an exact identity to the Xi contact scalar.

## Aggregate Occupation Identity

```text
On a pole-free local chart let I be the interior atoms with c_i!=0, choose A<min_(i in I)h_i and B>max_(i in I)h_i, and put P_x(s)=sum_(i in I)c_i*1_(A<s<h_i). Then exactly sum_(i in I)d_i=A*sum_(i in I)c_i+integral_A^B P_x(s)ds.
Let D_perp be the sum of d_i over interior atoms with c_i=0. With X_I=sum_(i in I)c_i, the endpoint-complete scalar is mathcal_C_N=D_edge+D_perp+A*X_I+integral_A^B P_x(s)ds, while mathsf_X=C_edge+X_I.
```

The threshold mass P_x(s) is permutation-free. At a pairwise alignment the two order-cell descriptions join without a singular term: h_2(c_1+c_2)-(h_2-h_1)c_1=h_1c_1+h_2c_2=h_1(c_1+c_2)-(h_1-h_2)c_2.

When an atom reaches c_i=0 it is retained in D_perp division-free and the next local projective chart is chosen before continuing. The exact pole-rate formula, not a tangent quotient, is the chart join.

The aggregate formula keeps (C_edge,D_edge)=(c_0+c_N,d_0+d_N) exact. At a cutoff, recompute terminal centering and transport this block through J_a^(adj); never replace the generally complex within-chart J_a^(der) by J_a^(adj), and never assign separate endpoint or terminal signs.

## Route Decision

Retire a uniform pairwise projective Sturm sign as a consequence of the currently proved current bounds. This does not prove that a stronger arithmetic Xi comparison is false. The next admissible target is a one-sided estimate for the signed occupation integral integral P_x(s)ds together with D_edge and D_perp, using actual Xi amplitudes and oscillation before absolute values. Pairwise alignments are harmless permutation joins; projection poles and cutoff edges remain exact separate joins.

The q<1 layer remains a separate multiplicity-compatible parabolic/Hermite or boundary-degree theorem. No pairwise or occupation identity here proves t=0 endpoint simplicity.

## Boundary

This proves the exact finite-alignment relative-current identity, the favorable frozen-scale model, the exact moving-scale sign-reversal guard, the division-free pole rate and its two-sided chart limit, the residual-order insufficiency guard, and a permutation-free local occupation identity. The guards are not actual Xi counterexamples and do not rule out a stronger source-specific theorem. This does not prove a signed Xi occupation bound, an edge-block sign, the Abel-scalar gap, successor count, q<1 closure, finite-height effectivity, contact exclusion, Lambda<=0, PF-infinity, RH, or a Clay-prize conclusion.
