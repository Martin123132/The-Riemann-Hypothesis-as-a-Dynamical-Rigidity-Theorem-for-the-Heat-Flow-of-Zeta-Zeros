# Characteristic canonical-core interval quadrature

Date: 2026-08-10

Status: finite characteristic-core quadrature validated; not a proof of the
outer normal-saddle join or the complete Hardy transformation

For the grouped canonical phase from Section 11.314, define

```text
I_M(Z,Y)=sum_m integral_(-Z)^Z integral_0^Y
 exp(i[z^3/3-lambda*z-(z+d_m)y+y^2/(4beta)]) dy dz.     (CQ1)
```

Two exact one-dimensional representations are used.  Integrating `z` on the
full real line first gives

```text
I_M(infinity,Y)=2pi integral_0^Y Ai(-lambda-y)
 exp(i*y^2/(4beta))D_M(y)dy.                            (CQ2)
```

Integrating `y` first on `[0,Y]`, complete the square:

```text
H_beta(p;Y)=exp(-i*beta*p^2)sqrt(2beta)
 integral_(-sqrt(2beta)p)^((Y-2beta*p)/sqrt(2beta))
 exp(i*q^2/2)dq.                                        (CQ3)
```

Then `I_M(Z,Y)` is the `z` integral of
`exp(i[z^3/3-lambda*z]) sum_m H_beta(z+d_m;Y)`.

The detuning lattice gives strict branch margins

```text
sqrt(lambda)-d_max=[0.0058655409515117871899773863921807083 +/- 3.37e-30] >0.0058,
sqrt(lambda)+d_min=[0.032867370546508168938542494450478896 +/- 3.37e-30] >0.032.
```

Thus neither leading Airy branch has a normal stationary point on
`0<=y<=64`.  This is an arithmetic lattice effect; the dangerous continuous
threshold lies just outside the final integer mode.

At `Z=4`, `Y=64`, 32 deterministic interval panels give

```text
I_M(4,64)=[-46.273430907931427642709548334683336 +/- 1.02e-12]
          +i*[1.4147856469615762362957843469812907 +/- 1.02e-12].             (CQ4)
```

The exact finite-`t` logit phase remains quadratic in `y`.  A degree-17 odd
Taylor form for `tanh(z/beta)` removes interval dependency cancellation; a
Cauchy estimate on `|w|=1` bounds its contribution to the full rectangle by

```text
[3.9987949799942440767110646728760747e-35 +/- 2.33e-68] <5e-35.                 (CQ5)
```

After restoring the exact amplitude, the finite-`t` rectangle divided by the
constant `sqrt(2)/pi` is

```text
[-46.273728557590014905179857409197297 +/- 9.85e-13]
+i*[1.4147878359740294320447177389450181 +/- 9.85e-13].                  (CQ6)
```

Its direct correction to (CQ4) has magnitude and relative size

```text
[0.00029765770783340073839338263406162355 +/- 2.02e-12] <0.00031,
[6.4295376875348422659241795842910427e-6 +/- 4.36e-14] <6.6e-6.       (CQ7)
```

This cancellation-preserving comparison is much sharper than integrating the
pointwise phase-error majorant absolutely.

The full-`z` value from (CQ2) is

```text
[-43.332933185562355910492316316787781 +/- 1.05e-23]
+i*[1.4281900292033471330446424551391455 +/- 1.05e-23].                     (CQ8)
```

Consequently the canonical `|z|>4`, `0<=y<=64` contribution has relative
size

```text
[0.063517115097500926000702623852936513 +/- 2.36e-14] <0.064.         (CQ9)
```

It is controlled and explicitly measured, but it is not negligible.  The
next theorem must place (CQ9) on a legal steepest-descent contour and join the
`y>64` component to every displaced alpha saddle.

Pi provenance: `pi` in (CQ2) is the standard Airy Fourier inversion factor;
`beta`, `d_m`, and the physical `sqrt(2)/pi` factor retain the Kummer and
Fourier--Poisson normalization traced in Sections 11.308 and 11.314.  No
geometric fit is introduced.

Proof boundary: rigorous finite quadrature and exact-transform algebra at
`t=10^10` only.  No height-uniform canonical bound, `y>64` saddle partition,
complete `T_upper` assembly, `Lambda<=0`, RH, or prize-level conclusion is
proved.
