# Contiguous Terminal-Tail Growing-Prefix Finite-Height Gate

Date: 2026-07-31

Status: retained first-order `q=1` growing-prefix current theorem. This is not
a proof of the remaining bulk cross-current estimate, an Xi-level current
theorem, `Lambda<=0`, or RH.

## Domain

q=2tL^2=1, L>=50, theta=(1-p)/2 in [0,1], a=N+theta, h=1/a<1/72000000000

Take M>=1, K=M+1, h*K^18<=1, and group the retained endpoint-terminal edge with n=N-1,...,N-M.

Equivalently, `1<=M<=floor(a^(1/18))-1`; this range tends to infinity.

## Exact Carrier Ratio

```text
ell_j=log(1-h*j), k=m+theta, Delta_m=ell_k-ell_theta, chi=alpha-log(a), r=-exp(2*pi*i*p)
L_m=log[(f_(N-m)/f_N)/r^m]=-Delta_m/2-(t/2)chi*Delta_m+(t/4)(ell_k+ell_theta)Delta_m+i[R_m-(pi*t/8)Delta_m]+log(1+d_k)-log(1+d_theta)
R_m=2*pi*h^-2*Delta_m+2*pi*m*(h^-1-theta)+pi*(m^2+4*m*theta)
L_(m,x)=-s_*'*Delta_m+d_(k,x)/(1+d_k)-d_(theta,x)/(1+d_theta)+i*h*m/2
raw_log_ratio-L_m-m*Log(r)=-2*pi*i*[m*N+m*(m+1)/2], N=h^-1-theta in Z
```

The cubic phase limit is

```text
R_m/h -> -(2*pi/3)[(m+theta)^3-theta^3]
```

The inherited source bounds and exact logarithm estimates give

```text
|Delta_m|<2*h*k, |Delta_m+h*m|<h^2*k^2, |ell_k|<2*h*k, |U-k|<h*k^2
d_k-d_theta=alpha'*t^2*(w_k-w_theta)*(w_k+w_theta)/8; the analogous x derivative is the sum of the alpha'' quadratic and (alpha')^2 linear terms
|d_k|<h^2/10
|d_(k,x)|<h^4/100
|W_k|<12*h*k^3
|W_(k,x)|<3*h^2*k^2
|exp(W_k)-1|<24*h*k^3
```

## Four Coordinates

```text
Z=X+iY, V=Z_x/h, U=-ell_k/h, K_m=-Delta_m, gamma=-Im(s_*'); C=X; D=(Re(s_*')/h)K_m X+gamma UY; E=Re V; N=(Re(s_*')_x/h^2)K_mX+(Re(s_*')/h)K_mReV+[gamma_x*(-ell_k)/h^2+gamma/(8*pi)]Y+gamma U ImV
Z_0=-Q*r^m=X_0+iY_0, V_0=-i*k*Z_0/2; (C_0,D_0,E_0,N_0)=(X_0,kY_0/2,kY_0/2,Y_0/(16*pi)-k^2X_0/4)
```

For each carrier,

```text
|C-C_0|<24*h*k^3
|D-D_0|<50*h*k^4
|E-E_0|<18*h*k^4
|N-N_0|<40*h*k^5
```

After adding the certified edge errors and summing `m=1,...,M`,

```text
|eps_C|<124*h*K^4
|eps_D|<300*h*K^5
|eps_E|<268*h*K^5
|eps_N|<940*h*K^6
```

## Current Margin

```text
Delta P=C*eps_N+N*eps_C+eps_C*eps_N-D*eps_E-E*eps_D-eps_D*eps_E
|Delta P|<5020*h*K^7+196960*h^2*K^10
|Delta P|<6000*h*K^7<1/400
```

The nontrivial leading prefix has `P_M<-1/200`; therefore

```text
a^2*J_(tail,M)/S_a^2<-1/200+1/400=-1/400
```

For every admissible M>=1, the retained first-order q=1 contiguous terminal-prefix current satisfies a^2*J_(tail,M)/S_a^2<-1/400.

4^18=68719476736<72000000000, so M=3 is admissible throughout L>=50. At L>=50 the theorem includes M=3, so the canonical prefix through N-3 is finite-height clockwise even though the edge plus N-3 alone has both limiting orientations.

## Boundary

Every pi is inherited from the completed-zeta/Riemann-Siegel source, a^2=x/(4*pi)+t/16, and the published C_0 recurrence; no fitted geometry introduces it.

Join this growing but short terminal block to the remaining n<=N-M-1 Abel/cross-current bulk without termwise absolute values, then use the existing C1 whole-jet residual transfer.

This theorem does not sign the remaining nonterminal bulk, extend the prefix
to all `N-1` carriers, prove the `q>1` grouped-tail theorem, control the omitted
Xi remainder inside one component, prove an Abel gap, winding cap, contact
exclusion, `Q209`, `Lambda<=0`, PF-infinity, RH, or a prize-level conclusion.

## Reproduce

```text
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_growing_prefix_finite_height_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_growing_prefix_finite_height_gate.py
```
