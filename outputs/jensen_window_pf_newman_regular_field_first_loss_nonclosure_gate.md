# Newman Regular-Field First-Loss Nonclosure Gate

Date: 2026-08-03

Status: exact arbitrary-field first-loss countermodel theorem; actual Xi field and degree-uniform Jensen remainder remain open; not a proof of RH.

## Construction

Fix any real `ell`, any `c>0`, and any integer `m>=2`. Set

```text
A=ell^2+2,
B=ell^2-ell+2=(ell-1/2)^2+7/4,
alpha=1-1/A,
beta=1+1/B.
```

Then `A>1`, `B>0`, and `0<alpha<1<beta`. Define

```text
F_(m,ell,c)(s)=(s+c^2)^m(s+alpha*c^2)(s+beta*c^2).
```

Every zero is negative and every coefficient is strictly positive. At the multiplicity-`m` zero `rho=-c^2`, with `F=(s-rho)^m U`, the two simple roots contribute

```text
rho/(rho+alpha*c^2)=A,
rho/(rho+beta*c^2)=-B,
rho U'(rho)/U(rho)=A-B=ell.
```

Thus this family realizes every finite real regular field, not only both signs.

## Exact Heat Boundary

Lift the product to the even real-rooted polynomial

```text
H_(m,ell,c)(x)=F_(m,ell,c)(-x^2).
```

It has multiplicity `m` at both `c` and `-c`. For `tau>=0`,

```text
exp(-tau D_x^2)
 =lim_(N->infinity)[(1-sqrt(tau/N)D_x)(1+sqrt(tau/N)D_x)]^N.
```

Every first-order factor preserves real-rootedness, so `H_tau` is real-rooted for every nonnegative time. Near `c`, the two parabolic limits are

```text
tau^(-m/2) H_tau(c+sqrt(tau)z)/V(c)
 -> exp(-D_z^2)z^m=2^(m/2)He_m(z/sqrt(2)),

sigma^(-m/2) H_(-sigma)(c+sqrt(sigma)z)/V(c)
 -> exp(D_z^2)z^m.
```

The first polynomial has `m` simple real roots. The second has no real root for even `m` and only the root zero for odd `m`, so a nonreal pair appears for every sufficiently small negative time. Heat monotonicity therefore makes zero the exact Newman-style threshold of every model.

## What The Field Does

If `H(x)=(x-c)^m V(x)`, then

```text
B_c=V'(c)/V(c)=(m+4ell)/(2c).
```

Writing `P_m^-=exp(-D_z^2)z^m`, the exact Appell-Hermite recurrence

```text
P_(m+1)^-=z P_m^--2(P_m^-)'
```

shows that every member of the forward cluster has

```text
x_j(tau)=c+sqrt(tau)h_j+2B_c tau+O(tau^(3/2)).
```

So `ell` changes the common drift but cannot change the real-versus-nonreal side of the contact. For the multiplicity-three first-Xi-jet case, `B_c=(3+4ell)/(2c)`.

## Route Decision

This is a genuine nonclosure theorem for the local route. First-loss heat dynamics, evenness, finite genus-zero negative-root geometry, and positive `s`-coefficients do not impose any sign, finite interval, or one-sided bound on `ell`. The earlier two sign examples were not exceptional; every real value occurs inside an exact Newman-boundary model.

The result does not say the actual Xi field is arbitrary. Any closing estimate must use structure absent here: the theta/Fourier source, global Xi zero placement, or zero first-jet boundary winding. The cofinal Jensen route also still requires a remainder uniform in degree.

## Audit

The builder certifies 10 exact rational field models and 15 Hermite multiplicities. It checks 10 exact positive-time evolved polynomials, all with every root real. The independent checker reconstructs separate rational models and the forward/backward local polynomials.

```text
python work/rh_compute/scripts/jensen_window_pf_newman_regular_field_first_loss_nonclosure_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_regular_field_first_loss_nonclosure_gate.py
```

Current result:

```text
validated regular-field first-loss nonclosure gate: 16 rows, 3 sources, 10 arbitrary-field models, 15 multiplicities, 10 positive-time samples, 15 forward Hermite checks, 15 backward nonreal checks, 2 open handoffs, 0 Xi field bounds, 0 degree-uniform bounds
```

No `pi` enters this construction. Its constants are rational functions of the prescribed field; any `pi` in the Xi programme remains tied to the separate Fourier/theta normalization.
