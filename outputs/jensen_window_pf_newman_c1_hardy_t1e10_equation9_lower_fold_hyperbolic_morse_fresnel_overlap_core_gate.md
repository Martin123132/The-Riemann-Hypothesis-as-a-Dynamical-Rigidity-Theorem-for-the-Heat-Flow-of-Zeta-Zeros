# Characteristic hyperbolic Morse-Fresnel overlap core

Date: 2026-08-12

Status: exact hyperbolic fold form and a compact three-height overlap core
certified; this is not a proof of the complete Airy/logistic integral join

Keep the paired finite-Poisson current from Section 11.337 and set

```text
u=sqrt(t/2)s,       Q=sqrt(pi)q,
Y=y/sqrt(2 beta),  beta^3=pi C^2/8.                    (HF1)
```

The universal logistic phase and exact Fresnel square then combine to the
exact hyperbolic phase

```text
Phi-Phi_center=(Q^2-Q_0^2-u^2)/2.                     (HF2)
```

There is no finite-`t` phase remainder in (HF2).  On the same physical alpha
slice as the Airy chart,

```text
Q=Q_C(u)+sqrt(2x(u))Y.                                (HF3)
```

After extracting the common `sqrt(2)/pi` factor, the exact amplitude is

```text
A_t(u,Y)=r^(1/4)v^(-1/4)(dv/ds)sqrt(2x)
 [1+(Q_C+sqrt(2x)Y)sqrt(rx/t)].                       (HF4)
```

For mode `39894`, the canonical Airy-Fresnel phase under
`z=-d+u/sqrt(2 beta)` is

```text
Phi_0-Phi_0,center
 =ell u+u^2/(2C)+u^3/(3C sqrt(pi))-uY+Y^2/2,
ell=2(t-tau)/(C sqrt(pi)).                             (HF5)
```

Thus the exact-minus-canonical phase is

```text
R_t(u,Y)=R_boundary(u)
 +[sqrt(2x)Q_C+u]Y+[x-1/2]Y^2.                        (HF6)
```

The exact mode carrier and (HF5) carrier differ by less than
`[1.927664248648459925938916988109521117005123713366041604455528786067000151054148308610432238711012051e-12 +/- 2.78e-112]` at all three certified
heights.

## Certified core

The interval cover `|v-1|<=0.00284` contains `|u|<=200` at

```text
t=tau-pi/16,       t=tau,       t=tau+pi/16,
0<=Y<=64/sqrt(2 beta).
```

Across all three heights, independent interval panels prove

```text
|R_boundary| < [0.01536731587133653482373777987057258695825489721298228790507948813802307596406206049027532877662302396 +/- 1.01e-102],
|sqrt(2x)Q_C+u| < 0.09515299883509335645733350327458814621461868910046177916228771209716796875000000000000000000000000000,
|x-1/2| < 0.0007121265119786404640400689913803944364190101623535156250000000000000000000000000000000000000000000000,
|R_t(u,Y)| < [0.08171152709230543202643484554114436721180975375124409349935328080261027164117908397556391400506615203 +/- 2.38e-102],
|A_t(u,Y)-1| < [0.0009537211312279485575704770189319514074805086597975165938466469197985668669856127522644649348364241925 +/- 5.25e-104].   (HF7)
```

The logarithmic quotient
`[v-1-log(v)]/(v-1)^2` was evaluated by a degree-18 convergent
series with an explicit geometric tail.  This avoids interval cancellation
at `v=1`.

## Interpretation

The Airy fold and ordinary-Morse region are not two unrelated phases.  They
are the same exact quadratic form `(Q^2-u^2)/2` cut by a lower boundary whose
tangent is nearly the characteristic line `Q=-u`.  Straightening that
boundary produces the cubic Airy term in (HF5).  The retained Fresnel current
is therefore the correct transition object on both sides of the event.

## Boundary

This gate proves an exact coordinate identity and compact pointwise phase and
amplitude bounds.  It does not integrate the exact-minus-canonical core,
control `|u|>200`, complete the one-mode Airy/logistic remainder, propagate
through 399 events, establish complete `T_upper`, or prove `Lambda<=0`, RH,
or a prize-level conclusion.
