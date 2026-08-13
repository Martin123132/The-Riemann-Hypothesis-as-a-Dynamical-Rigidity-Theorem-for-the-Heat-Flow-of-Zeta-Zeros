# Equation-(10) global theta-current reduction

Date: 2026-08-10

Status: exact algebraic and roster reduction validated; not a proof of the required uniform estimate

The paper's exact paired collection in equation (10) contains both

```text
alpha_E cos(pi*j*alpha_E*x/2)
```

and

```text
i*j sin(pi*j*alpha_E*x/2).
```

Before equation (51), the paper explicitly drops the second term and applies
the saddle approximation (50).  Equation (62), equations (126)--(127), and
the source Gaussian-sum path inherit that non-exact step.

For `u=pi*alpha_E*j*x/2`, the complete pair has the exact identities

```text
(alpha_E+j)e^(iu)+(alpha_E-j)e^(-iu)
  =2[alpha_E cos(u)+i*j sin(u)],

P(alpha_E,j;x)
  =4/(i*pi*x) d/d(alpha_E)
     [e^(i*pi*(alpha_E^2+j^2)*x/4) cos(pi*alpha_E*j*x/2)].
```

Thus the omitted sine term is not an unrelated correction: it completes the
retained cosine term into a pivot derivative current.

The source checkpoint proves that its 36 blocks partition the contiguous odd
roster

```text
159577, 159579, ..., 5122421
```

with no gap or overlap.  It contains 2,481,423 alpha terms:
101 direct terms followed by
7,034 paired collections representing
2,481,322 terms.  The two phase branches require
14,068 Gaussian phase sums.

Define the finite incomplete theta sum

```text
Theta_[A,B](x,y)
 =sum_(A<=alpha<=B, alpha odd)
    exp(i*pi*alpha^2*x/4+i*pi*alpha*y/2).
```

Then, exactly,

```text
sum alpha exp(i*pi*alpha^2*x/4)
 =2/(i*pi) d/dy Theta_[A,B](x,y)|_(y=0).
```

This is the first global coordinate that preserves every alpha and both
members of every equation-(10) pair before asymptotic approximation.  It is
therefore a viable route for the signed source-aligned error.

The finite atlas explains why this matters: its independent block triangle is
[0.13303034632772588889577778961552984056331811542139107007600333460949171898645177828910589990288151071169985270 +/- 1.56e-95] to
[0.15526911168332570866524206560505587017936704396427745771141840173422307665265110948460749888620594547545527539 +/- 1.56e-95], while only
[0.046676278888853190883810738547243994065873870564711966046631841291001587350584384053636955042345017900784087929 +/- 1.06e-94] to
[0.049885739432795757895306957463493408724662413234293278453121755392441111230634118209931142229329333425497311454 +/- 1.23e-94] survives in the signed sum.

There is an important guard.  A continuous pivot derivative sampled on a
step-two lattice is not automatically an endpoint telescope.  The next
theorem must apply an endpoint-complete finite Poisson, Abel, or
Euler--Maclaurin transform to the incomplete theta derivative, carry the
transition and global `REM` terms, and identify its dual saddle range with the
classical `T_upper` sum.  Only then can the source saddle truncation and
Gaussian-evaluator defects be inserted.

Proof boundary: exact symbolic identities, paper/source audit, and finite
roster coverage only.  No finite-Poisson remainder, source-error bound,
height-uniform Hardy theorem, `Lambda<=0`, RH, or prize-level conclusion is
proved.
