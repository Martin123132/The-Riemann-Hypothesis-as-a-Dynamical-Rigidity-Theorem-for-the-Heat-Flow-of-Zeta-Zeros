# Exact finite-endpoint tail decomposition

Date: 2026-08-13

Status: exact reassembly reduction proved; not a proof of the quantitative endpoint-tail bounds

Let

```text
a_m=m sqrt(2/x),
E_D=exp(i*pi*q_D^2/2),
F_D=integral_0^(q_D) exp(i*pi*q^2/2)dq.
```

Since `F(+infinity)=(1+i)/2` and
`F(-infinity)=-(1+i)/2`, every finite current splits exactly as

```text
P_m=(E_B-E_A)/(i*pi)+a_m(F_B-F_A)
   =P_bulk+P_A+P_B,                                    (ET1)

P_bulk=a_m(1+i),
P_A=-E_A/(i*pi)-a_m[F_A+(1+i)/2],
P_B=+E_B/(i*pi)-a_m[(1+i)/2-F_B].                      (ET2)
```

Section 11.349 closes `P_bulk` by an exact Gamma integral.  Equations
(ET1)--(ET2) isolate the only remaining finite-endpoint objects without
freezing either at the outer saddle.

There is a mandatory summation guard.  Quadratic completion gives

```text
-m^2/x+q_D^2/2=xD^2/4-mD.                             (ET3)
```

Both source endpoints are odd, so for integer `m`

```text
(-1)^m exp(-i*pi*mD)=1.                               (ET4)
```

The bare endpoint exponential is therefore independent of `m`.  It is not
an independently summable infinite Poisson series.  The lower and upper
terms in (ET2) must remain paired with their Fresnel tails and then be
reassembled with the symmetric zero, negative, and outer-positive modes plus
the endpoint half-current before absolute values are taken.  This explains
why the frozen endpoint-current model failed while the exact bulk succeeded.

After parity, the endpoint phase is

```text
Phi_D(x)=pi D^2 x/4+(t/2)log((1-x)/x),
Phi_D'(x)=pi D^2/4-t/[2x(1-x)].                        (ET5)
```

The lower `A` phase is the characteristic fold already covered locally by
the Airy/logistic atlas.  The upper `B` phase is the nonstationary endpoint
piece already bounded on the transition chart.  The next quantitative task
is to perform their complete symmetric roster reassembly, not to assign
termwise endpoint errors.

Pi provenance: all occurrences derive from the equation-(9) quadratic
Kummer phase and integer Fourier-Poisson character.  No geometric fit is
used.

Proof boundary: exact per-mode decomposition, endpoint phase, and summation
guard only.  No complete symmetric endpoint-tail bound, ordinary/fold
reassembly, complete `T_upper` theorem, height-uniform theorem, or
`Lambda<=0` theorem is proved.  No claim of PF-infinity, RH, or a prize-level
conclusion is made.
