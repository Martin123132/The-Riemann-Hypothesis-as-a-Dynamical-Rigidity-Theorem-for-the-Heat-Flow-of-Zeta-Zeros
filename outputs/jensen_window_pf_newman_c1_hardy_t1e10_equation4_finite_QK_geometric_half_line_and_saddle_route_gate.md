# Finite Q_K geometric half-line compression and saddle route

Date: 2026-08-27

Status: exact finite compression certified; saddle-ray deformation open

Section 11.484 writes every explicit A33 label as one convergent half-line
integral.  Because the physical roster is finite, those integrals can be
summed without any interchange theorem.  Put

```text
s=1/2+it,
lambda=pi(1+i)/sqrt(2),
K_0=exp(3pi*t/4+i3pi/8),
M=(B-A)/2+1=2481423.                                (GH1)
```

For one odd label,

```text
S_alpha=K_0 integral_0^infinity
        x^(-s)exp(-pi*x^2-lambda*alpha*x)dx.         (GH2)
```

Therefore, exactly,

```text
S_W=K_0 integral_0^infinity x^(-s)exp(-pi*x^2)G_W(x)dx,

G_W(x)=sum_(j=0)^(M-1)exp[-lambda(A+2j)x]
      =exp(-lambda*A*x)
       (1-exp(-2lambda*M*x))/(1-exp(-2lambda*x)),    (GH3)

G_W(0)=M.                                           (GH4)
```

For real `x>0`, `|exp(-2lambda*x)|<1`, so the denominator in (GH3) is
nonzero.  The apparent endpoint singularity is removable by (GH4).  With
`Hardy_t[X]=2Re[exp(i theta(t))X]`, the corrected physical roster is

```text
H(t)Q_K=Hardy_t[S_W].                               (GH5)
```

Thus 2,481,423 ill-conditioned direct Kummer calls have been replaced by one
exact scalar integral and a stable finite geometric quotient.  At the
surrogate `t=5.0`, roster `1..7`, direct summation of the
individual Kummer labels and independent quadrature of (GH3) differ by only
`4.1740077599530154e-31`.

The one-label exponent in (GH2) is

```text
F_alpha(x)=-pi*x^2-lambda*alpha*x-s Log(x).          (GH6)
```

Its saddles solve

```text
2pi*x^2+lambda*alpha*x+s=0,
D_alpha=lambda^2 alpha^2-8pi*s
       =-4pi+i*pi*(pi*alpha^2-8t).                  (GH7)
```

The transition is therefore exactly

```text
alpha_0=sqrt(8t/pi).
```

At `t=10^10`, Arb gives

```text
alpha_0=[159576.9121605730711759784239737527473903434524659738630663703318682632 +/- 2.97e-65],
159577-alpha_0=[0.08783942692882402157602624725260965654753402613693362967 +/- 1.88e-57]. (GH8)
```

Hence `159577` is the first odd label above the saddle transition; `159575`
is below it.  Arb also places both roots of (GH7) in the open third quadrant
at the previous odd label, the first two roster labels, and the upper roster
endpoint.

There is a revealing algebraic candidate contour.  On the ray
`x=exp(-3i*pi/4)y`, using `Log(x)=log(y)-3i*pi/4`, the whole large prefactor
cancels pointwise:

```text
K_0 x^(-s)exp(-pi*x^2-lambda*alpha*x)dx
 =y^(-s)exp[i*pi*(alpha*y-y^2)]dy.                  (GH9)
```

The real stationary equation of the right-hand side is

```text
2pi*y^2-pi*alpha*y+t=0,                             (GH10)
```

with the same threshold `alpha_0`.  Equation (GH9) explains the source's
transition geometry and removes the exponentially large normalization at
the integrand level.  It is not yet permission to rotate the integral:
`exp(-pi*x^2)` grows in intervening Stokes sectors, and the branch at zero
must be connected with the correct multiplier.  A valid next proof must
derive that connection by a finite steepest-descent contour or an exact
parabolic-cylinder connection formula, including every large-arc and branch
contribution, before using (GH9) as an integral identity.

Pi provenance: all occurrences come from the RSI Gaussian/sine kernel,
quarter-turn contour, and gamma/Kummer normalization.  The transition value
`sqrt(8t/pi)` is derived from the discriminant (GH7), not inserted or fitted.

Proof boundary: exact finite geometric compression, removable endpoint,
surrogate quadrature check, exact saddle polynomial, and actual threshold
arithmetic only.  The saddle-ray integrand identity is algebraic, but its
contour deformation and Stokes multiplier remain open.  No actual-height
quadrature, `B_W`, `Q_K`, `D_K`, `Delta_KU`, `J_Z`, non-A, all-height,
`Lambda<=0`, PF-infinity, RH, or prize-level enclosure is proved.
