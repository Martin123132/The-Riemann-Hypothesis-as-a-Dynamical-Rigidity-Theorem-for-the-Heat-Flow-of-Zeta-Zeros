# Quartic Outer-Contact Normal-Form Gate

Date: 2026-07-25

Status: exact contact-level reduction; not a proof artifact for
the complete contact-cell sign of `Delta_14`, which remains open.

```text
work/rh_compute/results/jensen_window_pf_quartic_outer_contact_normal_form_gate.json
python work/rh_compute/scripts/jensen_window_pf_quartic_outer_contact_normal_form_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_quartic_outer_contact_normal_form_gate.py
```

## Root Chart

Put `delta=1-a` and `q=C/(4 delta^2)`. On the left outer branch,

```text
0<delta<1, 0<q<1,
p=4*delta**2*q - 3*delta**2 + 2*delta + 1,
Disc_(b,c)=-16*delta**2*(q - 1).
```

The right outer branch has `delta=-r`. The exact factorization of
`x_3-x_2` contains a quadratic polynomial convex in `q`; both
endpoint values are negative for `0<r<1`. Hence that branch has
`x_3<x_2` and is excluded whenever the required contraction
ordering is imposed.

## Scale Removal

The contact contractions have the exact form

```text
1-x_j=delta^2 e_j,  j=2,3,4,5,
x_5-U=delta^2(T-s).
```

Thus the outward region is `s<T`. The first two order-three gaps
have different contact scales:

```text
G_1=delta^6 h_1,
G_2=delta^4 h_2,
h_1=8*q**3/(2*delta**2*q - 3*delta**2 + 3)**3.
```

The second normalized gap splits exactly into an outer-threshold
slack plus a positive remainder:

```text
h_2=x_4^2*e_3*(T-s)+8*delta**2*q**3/(27*(2*delta**2*q - 2*delta**2 + delta + 1)**3)
```

## Adjacent Quintic Phase

Put `epsilon=T-s=(x_5-U)/delta^2`. The adjacent normalized
quintic has the exact discriminant factorization

```text
Disc(P_5)=positive_prefactor*epsilon^2*R(epsilon),
R(epsilon)=R_0+R_1*epsilon+R_2*epsilon^2,
R_0=256*delta**6*q**4*(16*q - 15)*(2*delta**2*q - 2*delta**2 + delta + 1)**2,
R_1=-96*delta**3*(delta - 1)**3*(15*q**2 - 40*q + 24)*(2*delta**2*q - 2*delta**2 + delta + 1)*(4*delta**2*q - 3*delta**2 + 2*delta + 1)**2,
R_2=9*(delta - 1)**6*(4*delta**2*q - 3*delta**2 + 2*delta + 1)**4.
```

The prefactor is strictly positive on the open contact chart,
and

```text
Disc_epsilon(R)=-147456*delta**6*(delta - 1)**6*(q - 6)**2*(q - 1)**3*(2*delta**2*q - 2*delta**2 + delta + 1)**2*(4*delta**2*q - 3*delta**2 + 2*delta + 1)**4>0.
```

For `0<q<15/16`, `R_0<0<R_2`. Hence `R` has exactly one
positive root `epsilon_+`, and

```text
0<epsilon<epsilon_+ => Disc(P_5)<0.
```

This is an exact negative-discriminant collar of strict outward
contacts. It proves adjacent degree-five nonhyperbolicity there,
not that an Xi contact enters the collar.

## First Extension

At `x_6`, the scaled-defect wall and order-four cap compete:

```text
W_6=delta^4 w_6,
C_6=delta^2 c_6,
Y_6=W_6/C_6=delta^2 w_6/c_6,
0<G_3<min(W_6,C_6).
```

A complete contact certificate must therefore split the exact
regimes `Y_6<=1` and `Y_6>=1`; the fixed-contact theorem lies in
the former regime.

For every later step, normalized gaps and defects obey

```text
h_(k-3)=y_k*h_(k-4)^2/(x_(k-2)^3*h_(k-5))
e_k=(e_(k-1)^2-h_(k-3))/(x_(k-1)^2*e_(k-2))
x_k=1-delta^2 e_k.
```

Finally,

```text
Delta_14=delta^4 delta14_hat.
```

This removes the degenerating contact scale. It does not prove
that `delta14_hat<0` over every admissible contact-tail cell. The
complete outer-contact obstruction, Xi threshold, PF-infinity,
`Lambda<=0`, RH, and Clay-prize problem remain open.
