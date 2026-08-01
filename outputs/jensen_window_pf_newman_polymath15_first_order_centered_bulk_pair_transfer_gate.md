# Newman First-Order Centered Bulk Pair-Transfer Gate

Date: 2026-07-26

Status: exact reflected/Abel bulk coordinates and a uniform
adjacent locked-pair determinant theorem. The aggregate Xi
cancellation problem remains open; this is not a proof of RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_bulk_pair_transfer_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_bulk_pair_transfer_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_bulk_pair_transfer_gate.py
```

## Reflected Bulk Coordinates

On

```text
L=log(x/(4*pi))>=50, 0<=tL<=25, a^2=x/(4*pi)+t/16, N=floor(a), 1<=n<=N-1
```

put `ell_n=log(n/a)`. The centered pair has the exact
division-free reflected-sharp form

```text
2*(X,A_a)^T=(g_0,G_a)^T+(conj(g_0),conj(G_a))^T+sum_n{(1,-s_*'*ell_n)^T*f_n+(1,-conj(s_*')*ell_n)^T*conj(f_n)}, ell_n=log(n/a)
```

For `s_*'=c+i*b` and `f_n=x_n+i*y_n`, the real
component transfer is

```text
For s_*'=c+i*b and f_n=x_n+i*y_n, V_n=(X_n,A_n)^T=T_(ell_n)(x_n,y_n)^T, T_ell=[[1,0],[-c*ell,b*ell]], det(T_ell)=b*ell
```

Here `b=-1/2-t*C_x/4<0`; for every strict interior
`ell_n<0`, the determinant `b*ell_n` is positive. This is
only a local orientation statement.

Finite Abel summation gives a second exact bulk basis:

```text
F_k=sum_(n<=k)f_n, h_k=log((k+1)/k): M_(1,a)=ell_N*F_N-sum_(k=1)^(N-1)h_k*F_k; A_a=Re(G_a-s_*'*ell_N*F_N+s_*'*sum_(k=1)^(N-1)h_k*F_k)
```

The weights `h_k=log(1+1/k)` are positive and independent
of the saddle chart. The chart dependence is confined to
the terminal/endpoint block already controlled by the adjacent
stability theorem.

## Exact Ratio Lock

Consecutive corrected coefficients are not arbitrary:

```text
R_n=f_(n+1)/f_n=rho_n*exp(i*delta_n)=exp[-s_*h_n+(t/4)*(2log(n)h_n+h_n^2)]*(1+d_(n+1))/(1+d_n)
d_(n+1)-d_n=alpha'*t^2*(h_n^2-2*h_n*alpha_n)/8, alpha_n=alpha-log(n)
```

The imported critical-frame bounds give

```text
|d_n|<2189/x<1/2, |d_(n+1)-d_n|<44*h_n/x, log(rho_n)<-(2/5)h_n, hence 0<rho_n<=exp(-2h_n/5)<1 and 1-rho_n>h_n/4
```

At `L=50`, `x=6.51529311985421213e+22` and the checked upper
bound for `|d_n|` is `3.35978744122717451e-20`.
The complete corrected logarithmic ratio per `h_n` is at most
`-4.12499999999999978e-01<-2/5`.

## Locked-Pair Determinant

Let `R(delta)` be the real rotation matrix. Then

```text
V_n+V_(n+1)=B_n*(Re f_n,Im f_n)^T, B_n=T_(ell_n)+rho_n*T_(ell_(n+1))*R(delta_n)
det(B_n)=b*(ell_n+rho_n^2*ell_(n+1))+rho_n*[b*(ell_n+ell_(n+1))*cos(delta_n)+c*(ell_(n+1)-ell_n)*sin(delta_n)]
```

Writing `b=-B`, `ell_n=-L_0`,
`ell_(n+1)=-M_0`, and `h_n=L_0-M_0`, this becomes

```text
With b=-B, ell_n=-L_0, ell_(n+1)=-M_0, h_n=L_0-M_0>0: det(B_n)=B*(L_0+rho_n^2*M_0)+rho_n*[B*(L_0+M_0)*cos(delta_n)+c*h_n*sin(delta_n)]
```

The phase can be eliminated without estimating it:

```text
C_x>0, D_x>0, B=1/2+t*C_x/4>=1/2, c=t*D_x/4<=1/(4x)<=h_n/16. Therefore det(B_n)>=h_n*[B*(1-rho_n)-c*rho_n]>=h_n^2/16.
```

Thus every adjacent corrected pair is orientation-positive
for every common phase. The determinant also yields

```text
||B_n||_F<3L, so ||(X_n+X_(n+1),A_n+A_(n+1))||_2>h_n^2*|f_n|/(48L).
```

This is a genuine two-frequency first-jet theorem. It holds
on the full displayed first-order domain, hence also on its
`q=2tL^2>=1` sublayer.

## Falsification Gate

Local orientation is not aggregate separation. For generic
component coefficients one has the exact family

```text
For arbitrary c, b!=0, and ell_j!=0, any target w_j=(P_j,Q_j) is realized by f_j=P_j+i*(Q_j/ell_j+c*P_j)/b. The nonzero targets (1,0),(0,1),(-1,-1) sum to zero.
```

Even preserving each internal adjacent ratio is insufficient:

```text
For any three invertible locked-pair matrices B_1,B_2,B_3, set u_j=B_j^(-1)w_j for w_1=(1,0), w_2=(0,1), w_3=(-1,-1). Every pair obeys its internal ratio lock and is nonzero, yet sum_j B_j u_j=0. Cross-block Xi ratios are essential.
```

The latter construction breaks only the ratios *between*
successive pairs. Keeping those relative ratios still leaves a
common absolute phase. In exact matrix notation,

```text
For W_k=B_(2k-1)u_k and u_(k+1)=rho_(2k-1)rho_(2k)*R(delta_(2k-1)+delta_(2k))*u_k, det(W_k,W_(k+1)) is the quadratic form rho_(2k-1)rho_(2k)*u_k^T*B_(2k-1)^T*J*B_(2k+1)*R(delta_(2k-1)+delta_(2k))*u_k
```

and the generic common-phase guard is

```text
If the frame drift and all three ratio phases are zero, positive pair matrices may be diag(p_1,q_1),diag(p_2,q_2). Then det(W_1,W_2)=(p_1*q_2-q_1*p_2)*x*y, which takes both signs under the common phase whenever p_1*q_2!=q_1*p_2.
```

Thus neither another isolated determinant estimate nor
relative ratio phases alone can supply the orientation. The
normalizer/endpoint phase anchor is part of the indispensable
Xi input.

## Xi-Specific Handoff

```text
Use the full consecutive ratio chain together with the exact normalizer/endpoint phase anchor, not isolated pair invertibility or relative phases alone, to prove a quantitative half-plane, signed-area, or one-sided winding bound for the paired block polygon plus the recurrent edge block on the contact value band.
```

A concrete next calculation is to substitute the actual
normalizer phase and recurrent endpoint block into the displayed
turning form, then reduce it in the saddle variable. Its purpose
is to decide whether the anchored stationary-phase block polygon
admits a quantitative half-plane or one-sided winding theorem.
The unanchored sign is now rigorously rejected, so it must not be
rediscovered by a large numerical sweep.

The unpaired last saddle remains composed through

```text
The disjoint adjacent blocks stop before the unpaired edge. The certified recurrence |A_(a,N+1)-A_(a,N)|<5000exp(-7L/4) must be retained for that last-saddle/endpoint block.
```

This artifact does not prove the aggregate contact-band scalar
bound, its signed crossing count, the `q<1` chart, finite phase
cells, contact exclusion, `Lambda<=0`, PF-infinity, RH, or a
Clay-prize conclusion.
