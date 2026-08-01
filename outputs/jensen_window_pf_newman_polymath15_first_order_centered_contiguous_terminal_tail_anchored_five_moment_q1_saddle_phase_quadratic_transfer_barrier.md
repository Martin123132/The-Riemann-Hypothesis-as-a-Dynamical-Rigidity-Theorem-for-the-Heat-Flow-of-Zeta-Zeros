# Physical q=1 Quadratic Transfer Barrier

Date: 2026-08-01

Status: exact algebraic reduction and route guard. This proves no physical
current sign and is not a proof of RH.

## Discriminant Collapse

```text
For M_0=f+ig, M_1=p+iq, and M_2=r+is, 4P_bulk^(0)=-f*r-q^2+2u_x*f*g.

The Section 11.162 discriminant factors identically as D_hat=(f*r+q^2-2u_x*f*g)|M_0|^4=-4P_bulk^(0)|M_0|^4.
```

D_hat is the primary quadratic current multiplied by |M_0|^4, not an independent sixth-degree positivity invariant. At M_0=0, D_hat=0 while P_bulk^(0)=-(Im M_1)^2/4. The discriminant therefore loses the strict exceptional-fibre information carried by P_bulk^(0).

## Exact Frequency Flow

```text
For M_k(xi)=omega_a H_k(xi), partial_xi M_k=-i M_(k+1), 0<=k<=4. Hence, if M_0=f+ig, M_1=p+iq, M_2=r+is, and y=Im M_3, then f_xi=q, g_xi=-p, q_xi=-r, and r_xi=y.

H_k(Omega)-H_k(T_0)=i*epsilon*integral_0^1 H_(k+1)(T_0-theta*epsilon)dtheta.

4 partial_xi P_bulk^(0)=q*r-f*y+2u_x(q*g-f*p).

Since Omega=T_0-epsilon, 4[P_bulk^(0)(Omega)-P_bulk^(0)(T_0)]=-epsilon integral_0^1 K(T_0-theta*epsilon)dtheta, where K=q*r-f*y+2u_x(q*g-f*p).
```

The independent endpoint expansion is

```text
For delta f,delta g,delta q,delta r, the exact endpoint difference is 4 delta P=-delta f*r-f*delta r-delta f*delta r-2q*delta q-(delta q)^2+2u_x[delta f*g+f*delta g+delta f*delta g].
```

## Terminal Scale

```text
At n=N, C_N w_N=-Q_RS(p)exp(W_theta), |Q_RS|=1, |W_theta|<6h, and C_N=(1+d_N)/(1+d_1).
The uniform |d_j|<1/2 gives |C_N|<3, hence A_N=|w_N|>exp(-6h)/3.
Since A_N=A_a exp[B_a u_N+u_N^2/(8L^2)], B_a<1/2, u_N<2h, and L>=50, the exponent is <2h. Therefore A_a>exp(-8h)/3>1/6.
```

## Positive-Mass Barrier

```text
Let K=M+1, B=N-K, and I={n integer: ceil(a/4)<=n<=floor(a/2)}. From K^18<=a and a>72000000000 one has K<a/4 and B>a/2, so I is contained in {1,...,B}.
The block has #I>=a/4-1>a/5. For n in I, log(2)<=u_n=log(a/n)<=log(4) and A_n>=A_a.
Writing d=log(2), A_0^+>a*A_a/5 and A_3^+>a*A_a*d^3/5.

Put A_j^+=sum_(n<=B)u_n^j A_n. Pointwise triangle inequality on the exact line integral gives |P(Omega)-P(T_0)|<=B_abs, where B_abs=epsilon[A_0^+A_3^++A_1^+A_2^++4u_x A_0^+A_1^+]/4.
Because a^2=x/(4pi)+1/(32L^2), x<4pi a^2. Thus epsilon>3/(8L^2x)>3/(32pi L^2a^2).
The positive A_0^+A_3^+ term alone forces B_abs>3A_a^2 log(2)^3/(3200pi L^2)>log(2)^3/(38400pi L^2).
Since h=1/a, a^2>exp(L), exp(L)/L^2 is increasing for L>=50, exp(1)>2, log(2)>69/100, and pi<22/7, B_abs/(h^2/400)>490000000.
```

This lower-bounds the stated positive-mass upper majorant. It does not lower-bound the true signed transfer error. It rules out spending an h^2/400 reserve through this canonical pointwise triangle estimate, not through cancellation-preserving analysis of the signed line integral.

## Route Decision

Do not estimate D_hat as a separate positive invariant and do not transfer P_bulk^(0) from T_0 to Omega by replacing every moment with its positive mass. The first step is algebraically redundant and the second has a certified majorant more than 490000000 times the entire h^2/400 reserve.

Retain K(xi)=q*r-f*Im(M_3)+2u_x[q*g-f*p] with its signs on the whole interval T_0-epsilon<=xi<=T_0. Insert its moments into the endpoint-composed Hermitian/transpose or balanced Mangoldt/Poisson representation before taking absolute values.

```text
A P-only benchmark is |integral_0^1 K(T_0-theta*epsilon)dtheta|<=h^2/(100epsilon), which would give |P(Omega)-P(T_0)|<=h^2/400. The actual proof may instead transfer E+P_bulk^(0)+R_corr jointly and use cancellation between those composed terms.
```

## Pi Provenance

The `pi` in `x<4pi*a^2`, `u_x=h^2/(8pi)`, and `T_0=2pi*a^2`
is inherited from the completed-zeta, Riemann--Siegel, and Poisson
normalizations already traced in Formal Core Section 11.118. No circle,
polygon, or fitted geometric constant is introduced here.

## Boundary

This proves the algebraic collapse of D_hat, exact frequency flow and signed transfer identities, A_a>1/6, a macroscopic positive-mass block, and failure of one explicit absolute-mass transfer majorant at the required scale. It does not prove the true saddle transfer is large, rule out oscillatory or composed transfer, sign P_bulk^(0) at T_0 or Omega, bound E or R_corr, upper-bound Phi_B, establish a signed Type-I/II, Vaughan, or Poisson estimate, exclude contact, transfer the retained model to Xi, prove Q209, Lambda<=0, PF-infinity, RH, or a prize-level result.

## Reproduce

```text
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_five_moment_q1_saddle_phase_quadratic_transfer_barrier.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_five_moment_q1_saddle_phase_quadratic_transfer_barrier.py
```
