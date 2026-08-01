# Newman First-Order Oriented Successor-Winding Reduction

Date: 2026-07-26

Status: exact oriented successor and one-sided winding reduction.
Three Xi boundary obligations remain open. This is not a proof of
`Lambda<=0`, RH, or the Clay prize.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_oriented_successor_winding_reduction.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_oriented_successor_winding_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_oriented_successor_winding_reduction.py
```

## Orientation

```text
bottom x:0->R; right t:lower->upper; top x:R->0; axis t:upper->lower
ind_(x,t)(H,H_x/ell)=+floor(m/2)>0
Writing coordinates as (t,x) reverses both contact degree and boundary orientation; it must not be mixed with the standard counterclockwise edge order.
```

The standard `(x,t)` orientation matches the boundary order used by
the finite Q208 polygon and by every ray-aligned rectangle.

## Scaled Pruefer Form

```text
W_J=J+i*J_x/ell, ell(t,x)>0
Q_ell=J^2+(J_x/ell)^2
partial_x arg(W_J)=(ell*J*J_xx-J*J_x*ell_x-ell*J_x^2)/(ell^2*Q_ell)
partial_t arg(W_J)=(ell*J*J_xt-J*J_x*ell_t-ell*J_x*J_t)/(ell^2*Q_ell)
At J=0 and J_x!=0, partial_x arg(W_J)=-ell<0.
```

## Signed Crossings

```text
For a closed regular path W=u+iv, wind(W,0) is the oriented intersection count with the positive imaginary ray: sum_(u=0,v>0) sign(-d_s u).
For W_J on x increasing, each zero J=0,J_x>0 contributes -1.
For the same path traversed with x decreasing, each zero J=0,J_x>0 contributes +1.
On x=R with t increasing, a transversal crossing J=0,J_x>0 contributes -sign(J_t).
A vertical tangency is interpreted by oriented intersection number, equivalently by a small regular ray rotation. Corners are assigned once by a half-open edge convention.
```

The vertical connector is part of the integer. Being contact-free
does not by itself make its open-path phase contribution zero.

## Successor Chain

```text
P_j=[t_j,T]x[0,R_j], D_j=[t_(j+1),t_j]x[0,R_j], S_j=[t_(j+1),T]x[R_j,R_(j+1)]
partial P_(j+1)=partial P_j+partial D_j+partial S_j
Since 1/4>1/5>=Lambda, P_0 contains no heat contact and wind(V_H(partial P_0),0)=0.
On S_j, L>=L_j and t>=t_(j+1), so tL>=t_(j+1)L_j=25 and L>=101.
The dominant-saddle theorem gives no contact in S_j, hence wind(V_H(partial S_j),0)=0.
w_(j+1)=w_j+kappa_j, where kappa_j=wind(V_H(partial D_j),0).
kappa_j=sum_(p in D_j) floor(m_p/2) is a nonnegative integer.
```

All internal boundary segments cancel in the stored exact chain audit.
The initial rectangle `P_0` is already contact-free because its lower
time is `1/4>1/5`; Q208 remains a compact calibration, not the logical
source of the `P_0` degree.

## One-Sided Integer Trap

```text
Boundary homotopy gives wind(proxy_j)=kappa_j>=0. Therefore the strictly weaker estimate wind(proxy_j)<1 (equivalently total oriented phase change <2*pi) forces kappa_j=0; exact zero winding need not be proved directly.
For a derivative-compatible regular proxy, kappa_j=N_up(t_j;R_j)-N_up(t_(j+1);R_j)-sum_(t on V_j; J=0,J_x>0) sign(J_t), plus the explicitly assigned finite-shoulder/join intersection count.
```

This is the useful relaxation: a rigorous oriented phase upper bound
below one turn suffices. Computing exact zero winding is unnecessary.

## First-Order Transfer

```text
||(r_[1],r_[1],x/L)||_2^2<50000000000*exp(-5L/2)
J_[1]^2+(J_[1],x/L)^2>50000000000*exp(-5L/2)
||Delta V_J||_2^2/50000000000*exp(-5L/2)<13/500<1
The full first jet, first-order proxy, and adjacent cutoff charts are joined by nonvanishing boundary homotopies. They have the same integer winding; cutoff joins add no independent integer theorem.
```

## Remaining Xi Input

```text
On q=2t_jL^2>=1, prove the first-order boundary margin using the saddle-centered logarithmic moment.
On q<1, prove a multiplicity-compatible boundary margin; its constants may deteriorate as t_j tends to zero.
Certify the finite L<50, oscillatory L_epsilon, core-to-main, and vertical-connector phase cells with exact joins.
For every j, prove the composite successor proxy has wind<1, or equivalently an oriented crossing count <=0. Boundary nonvanishing alone does not imply this inequality.
```

## Scope Guard

```text
F(t,x)=x^2-2t solves F_t=-F_xx.
(t,x)=(0,0) has standard local degree +1.
On [-1,1]x[-1/4,1/4], F+iF_x is boundary-nonzero but has winding +1.
A boundary norm floor, top/right nonvanishing, or cutoff continuity does not determine the winding integer.
```

## Live Handoff

```text
Construct one continuous successor boundary proxy, use signed phase cells to upper-bound its winding by a number below 1, and split only the hard bottom estimates into q>=1, q<1, and finite shoulders. The positive contact degree then forces zero without an exact phase evaluation.
```

This reduction does not prove either bottom margin, the finite
shoulders, the one-sided phase estimate, `Lambda<=0`, or RH.
