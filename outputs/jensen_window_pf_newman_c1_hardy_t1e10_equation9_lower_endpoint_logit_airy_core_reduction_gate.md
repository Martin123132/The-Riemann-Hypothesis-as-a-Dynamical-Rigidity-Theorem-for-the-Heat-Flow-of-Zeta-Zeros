# Exact lower-endpoint logit fold and Airy core

Date: 2026-08-10

Status: exact fold phase and certified Airy core validated; not a proof of the grouped fold remainder

Restrict the joint equation-(9) phase to the lower alpha endpoint `alpha=C`
and introduce the logit coordinate

```text
u=log((1-x)/x),       x=1/(1+exp(u)),
eta=pi*C^2/(8t).
```

After removing the constant `pi*C^2/8-pi*m*C`, the phase is exactly

```text
Phi_C(u)-Phi_C(0)=t[u/2-eta*tanh(u/2)].                 (LA1)
```

This is independent of the Poisson mode except for the removed constant.  Its
stationary equation is

```text
sech^2(u/2)=1/eta.                                     (LA2)
```

At the fold `eta=1`,

```text
t[u/2-tanh(u/2)]=t[u^3/24-u^5/240+...].                (LA3)
```

Put `mu=eta-1` and scale

```text
u=2z/(t*eta)^(1/3),
lambda=mu*t^(2/3)*eta^(-1/3).                          (LA4)
```

The linear and cubic phase becomes the canonical Airy polynomial exactly:

```text
t[-mu*u/2+eta*u^3/24]=z^3/3-lambda*z.                 (LA5)
```

A global real derivative bound `|tanh^(5)|<=512` gives

```text
|R_phase(z)|
 <=(64/15)t^(-2/3)eta^(-2/3)|z|^5.                    (LA6)
```

For `t=10^10`, `C=159577`, interval arithmetic certifies

```text
eta    =[1.0000011009042588332499966510739060892901771786178344435167606786115783877181726187489164531305393951710723053 +/- 4.73e-110],
mu     =[1.1009042588332499966510739060892901771786178344435167606786115783877181726187489164531305393951710723053436773e-6 +/- 3.63e-111],
lambda =[5.1099430394918356283059191674846647509663170203654407926790704695202619930481405568500939226158509559393019513 +/- 1.69e-104].
```

The two exact logit stationary points reproduce the prior portcullis roots:

```text
N_-=[39852.391386231238961843136786826147844269890017700142822067081499723395630016783577515625960639418285931660383 +/- 1.35e-103],
N_+=[39936.108613768761038156863213173852155730109982299857177932918500276604369983216422484374039360581714068339617 +/- 1.35e-103].
```

On the natural core `|z|<=4`, the conservative phase
remainder satisfies

```text
|R_phase|<=[0.00094128618812900810795259480279158435017965493452879447449818214531710186157684518019028494708880597560803135400 +/- 5.59e-113]<0.001. (LA7)
```

The transformed Kummer measure is also exact and regular:

```text
[x(1-x)]^(-1/4)|dx|
 =2^(-3/2)cosh(u/2)^(-3/2)du.                          (LA8)
```

There is a crucial grouping guard.  Since `C` is odd,

```text
(-1)^m exp(-i*pi*m*C)=1                                (LA9)
```

for every integer mode.  The lower endpoint current is therefore
mode-independent after parity and cannot be summed separately.  Equations
(LA1)--(LA8) must be applied inside the exact grouped boundary/Fresnel current,
where that apparent divergence cancels.

The canonical fold phase is no longer the unknown.  The next obligation is
to transport the full grouped amplitude through the Airy core, bound its
variation and the `|z|>4` contour tails explicitly, and join the result to the
ordinary interior and nonstationary ranges.

Proof boundary: exact lower-boundary fold algebra and a finite certified Airy
phase core at `t=10^10`.  No grouped fold-uniform remainder, complete
`T_upper` assembly, source-aligned height-uniform error, `Lambda<=0`, RH, or
prize-level conclusion is proved.
