# Jensen-Window PF Scaled Double-Zero Boundary-Layer Lemma

Date: 2026-07-11

Status: exact scaled double-zero boundary layer and finite-threshold
exhaustion theorem. This is not a proof of PF-infinity, RH, or
`Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_scaled_double_zero_boundary_layer_lemma.json
python work/rh_compute/scripts/jensen_window_pf_scaled_double_zero_boundary_layer_lemma.py
python work/rh_compute/scripts/check_jensen_window_pf_scaled_double_zero_boundary_layer_lemma.py
```

Current result:

```text
validated Jensen-window PF scaled double-zero boundary-layer lemma: 10 rows, 0 issues, 3 exact scaling identities, 1 heat PDE, 1 double-zero transversality, 1 universal boundary layer, 1 D^(-3/2) gap law, 1 external-field D^(-2) collision law, 1 exact toy family, 1 threshold-exhaustion theorem, 1 open uniform handoff
```

## First Degree Correction

For every fixed shift and compact `z` set,

```text
J_(D,n,lambda)(z/D)=sum_(j=0)^D ((D)_j/D^j)*A_(n+j)(lambda)*z^j/j!
(D)_j/D^j=1-j*(j-1)/(2*D)+O_j(D^(-2))
J_(D,n,lambda)(z/D)=F_(n,lambda)(z)-z^2*F_(n,lambda)''(z)/(2*D)+O_K(D^(-2))
J_D(z/D)=F-z^2*F''/(2*D)+(z^3*F'''/3+z^4*F''''/8)/D^2+O_K(D^(-3))
```

The `O_K(D^(-2))` remainder follows by splitting the entire series
on a larger compact disk and applying the falling-factorial expansion
to the coefficient core.

## Double-Zero Layer

The coefficient flow is

```text
partial_lambda F_n=(4*n+2)*partial_z F_n+4*z*partial_z^2 F_n
F(rho)=F'(rho)=0 => partial_lambda F(rho)=4*rho*F''(rho)
```

Let `rho<0` be a nondegenerate double zero at `lambda_*`. Under

```text
lambda=lambda_*+tau/D, z=rho+eta/sqrt(D),
```

the scaled Jensen polynomial has the universal limit

```text
If lambda=lambda_*+tau/D and z=rho+eta/sqrt(D), then 2*D*J_D(z/D)/F''_*(rho)->eta^2+8*rho*tau-rho^2.
```

Consequently

```text
eta_+/- -> +/-sqrt(rho^2-8*rho*tau); at tau=0, w_+/-=rho/D+/-|rho|/D^(3/2)+O(D^(-2)).
lambda_D=lambda_*+rho/(8*D)-[rho*(6*n+1)/64+rho^2*F'''_*(rho)/(48*F''_*(rho))]/D^2+o(D^(-2)).
z_D=rho+rho*(1-2*n)/(4*D)+o(D^(-1)) in the scaled Jensen variable.
```

Because `rho<0`, the finite-degree collision lies on the bad side of
the entire-function boundary. This is why every fixed degree can remain
strict at the limiting boundary.
Writing `F_*(z)=(z-rho)^2 U(z)` gives
`F'''_*(rho)/F''_*(rho)=3 U'(rho)/U(rho)`. Thus the leading layer is
universal, while the global root external field first enters at order
`D^(-2)`.

## Exact Model Check

The solvable family

```text
F_lambda(z)=12*lambda**2 + 4*lambda + z**2 + z*(12*lambda + 2) + 1
partial_lambda F=2*F'+4*z*F''
J_D(z/D)=12*lambda**2 + 12*lambda*z + 4*lambda + z**2 + 2*z + 1 - z**2/D
```

has a double zero `rho=-1` at `lambda=0`. Its near collision threshold is

```text
(sqrt(2)*sqrt(D - 1) - sqrt(2*D + 1))/(6*sqrt(2*D + 1))
=-h/8 + h**2/64 - 5*h**3/256 + O(h**4)
```

and its exact unscaled boundary gap is

```text
2/(sqrt(D)*(D - 1))
```

matching the universal `-1/(8D)` shift and `2/D^(3/2)` gap.

## Newman Threshold Exhaustion

Let `S_D` be the parameter set where degree `D` is hyperbolic and set

```text
Theta_D=inf{ell:[ell,infinity) is contained in the degree-D hyperbolicity set}
Theta_D<=Theta_(D+1)<=Lambda
sup_D Theta_D=Lambda
```

Polar descent gives the nesting. Closedness and Jensen's theorem force
the supremum to equal the full Newman boundary. A single endpoint test
does not bound `Theta_D`; one must control the full forward parameter ray.

The remaining proof obligation is now quantitative: a degree- and
height-uniform bound on the scaled collision layer and its remainders.
