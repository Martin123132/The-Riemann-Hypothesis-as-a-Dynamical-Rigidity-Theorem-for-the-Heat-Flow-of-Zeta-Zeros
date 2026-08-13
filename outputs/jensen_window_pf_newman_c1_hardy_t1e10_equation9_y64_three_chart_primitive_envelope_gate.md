# Exact y=64 three-chart primitives and envelopes

Date: 2026-08-11

Status: exact affine primitive and explicit saved-height chart envelopes
validated; not a proof of the remaining z integral or grouped mode sum

For fixed `(z,m)`, put

```text
k=x(z)/beta,        p=d_m+beta*tanh(z/beta),
y_*=p/k,            c=sigma/C,
E(y)=exp(i[k y^2/2-py]),
H0(L,U)=integral_L^U E(y)dy.                            (Y3E1)
```

Direct differentiation gives the exact affine primitive

```text
K(L,U)=integral_L^U(1+cy)E(y)dy
      =(1+c y_*)H0(L,U)-i(c/k)[E(U)-E(L)].              (Y3E2)
```

Using `pi*C^2=8beta^3` and `sigma=C/(2beta^2)`, (Y3E2) becomes

```text
C sigma K
 =2m sigma H0/x+2[E(U)-E(L)]/(i pi x).                 (Y3E3)
```

Thus (Y3E3) is exactly the source's paired Fresnel and endpoint current.  In
particular

```text
K(0,y_B)=K(0,64)+K(64,y_B),                            (Y3E4)
```

and the artificial `E(64)` endpoint cancels algebraically.  The split creates
no new source term and no double counting.

On the nonstationary chart `g_m(z)<=-delta`, one integration by parts uses
`phi'(64)>=delta` and `phi''=k>0`.  Evaluating the affine derivative and
curvature integrals before taking uniform bounds gives

```text
|K(64,y_B)|
 <= [32.76633261811716084385647121951513797828056082199418986287651812436831629855674161148090106291960692 +/- 1.08e-98] <33. (Y3E5)
```

On the endpoint/Fresnel chart, the normalized lower endpoint remains inside
`|q_0|<=[4.016744664136184053250960463946153597949110319858742833308741957427711151919819424997063124302504227 +/- 1.45e-99]`.  Completing the square and
using `|int_Q^infinity exp(iq^2/2)dq|<=2/Q` outside that compact interval
gives

```text
|K(64,y_B)|
 <= [561.2125999917513959820295904311901564902255551031971748985080208861439167441547275769823466478556500 +/- 3.47e-97] <600. (Y3E6)
```

This is an envelope only; the incomplete Fresnel function remains exact in
the working representation.

On the ordinary Morse chart, define the grouped main by replacing only `H0`
in (Y3E2) with the full Gaussian

```text
H_full=exp(-i p^2/(2k))sqrt(2pi/k)exp(i*pi/4),          (Y3E7)
```

while retaining the exact endpoint current in (Y3E2).  Then

```text
|K-K_Morse,grouped|
 <= [33.06232863894525889234696334879823972428052282455538145150293241944331943216434390199765530220176018 +/- 1.08e-98] <34. (Y3E8)
```

and the grouped main itself is below

```text
[165.7507167102790739554723056825340089120697918759980032263907285218626059961651795383323096097814814 +/- 4.64e-98] <200.    (Y3E9)
```

The constants in (Y3E5)--(Y3E9) are per-mode absolute envelopes.  They are
not summed over 84 modes and are not claimed to be small enough for the final
source error.  Their purpose is to close the local chart analysis without
losing the exact endpoint/Fresnel pairing; the next estimate must recover
cancellation in the remaining `z` integral and mode sum.

Pi provenance: every pi in (Y3E3) and (Y3E7) comes from the Kummer quadratic
phase, Fourier--Poisson character, and standard Fresnel Gaussian.  No fitted
geometric normalization is introduced.

Proof boundary: exact source-normalized primitives and explicit per-mode
saved-height envelopes only.  No completed z integral, grouped 84-mode sum,
join to all lower-interior modes, `T_upper` theorem, `Lambda<=0`, RH, or
prize-level conclusion is proved.
