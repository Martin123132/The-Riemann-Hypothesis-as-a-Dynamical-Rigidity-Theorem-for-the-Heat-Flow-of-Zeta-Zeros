# Exact half-Kummer reflection and one-branch reduction

Date: 2026-08-13

Status: exact-lemma certificate; not a proof of the quantitative half-domain remainder

For an odd source index `alpha`, write the equation-(9) integrand as

```text
K_alpha(x)=alpha exp(i*pi*alpha^2*x/4)
 exp[i(t/2)log((1-x)/x)]/[x(1-x)]^(1/4).
```

Odd squares satisfy `alpha^2=1 (mod 8)`.  Therefore

```text
K_alpha(1-x)=exp(i*pi/4) conjugate(K_alpha(x)),         (HR1)

exp(-i*pi/8) integral_0^1 K_alpha(x)dx
 =2 Re[exp(-i*pi/8) integral_0^(1/2)K_alpha(x)dx].     (HR2)
```

Equation (HR2) is termwise exact, so it remains exact after the finite odd
source-roster sum and its endpoint half-weights.  It also survives symmetric
finite Poisson summation and the already certified interchange:

```text
Q_main=2 Re[exp(-i*pi/8) Q_half],

Q_half=integral_0^(1/2) W_t(x)
 {H_x+lim_(M->infinity)sum_(m=-M)^M I_m(x)}dx.       (HR3)
```

This changes the stationary ledger decisively.  A positive Poisson mode has

```text
x_m=2*pi*m^2/(t+2*pi*m^2),
alpha_m=2m+t/(pi*m).                                   (HR4)
```

Since `x_m` is strictly increasing,

```text
x_m<1/2 iff m<sqrt(t/(2*pi)).                          (HR5)
```

At `t=10^10`,

```text
sqrt(t/(2*pi))=39894.228040143267793994605993438186847585863116493465766592582967065792589930184,
39894<sqrt(t/(2*pi))<39895.                            (HR6)
```

Consequently the half-integral has only one positive stationary branch:

```text
1..621:          B-endpoint side (alpha_m>B),
622..39852:      lower interior,
39853..39894:    A-endpoint side (alpha_m<A),
m>=39895:        no interior x saddle on 0<x<=1/2.     (HR7)
```

For `m<=0`, `partial_alpha Phi=pi(alpha*x/2-m)>0`, so no joint
saddle exists either.  The formerly listed reflected transition
`39895..39936` and upper interior `39937..2560588` live in `x>1/2`.
They are already represented exactly by the conjugate in (HR2), and must not
be counted as a second stationary source family in (HR3).

The continuous involution explains the geometry:

```text
m*=t/(2*pi*m),  alpha_(m*)=alpha_m,  x_(m*)=1-x_m.     (HR8)
```

It need not preserve integers because (HR1)--(HR3), not rounded reciprocal
pairing, performs the exact discrete reassembly.

The quantitative problem is now smaller.  On the half-domain there are only
two endpoint transitions: the `B` crossing `621|622` and the combined `A`
endpoint/`x=1/2` crossing `39894|39895`.  The latter is the characteristic
fold already represented by the certified atlas.  The former is
noncharacteristic and should be treated as one bulk-plus-Fresnel crossing,
together with the nonpositive and outer-positive nonstationary completion.

Pi provenance: every `pi` in (HR1)--(HR8) is inherited from the equation-(9)
Kummer phase, the odd-square character, or the Fourier-Poisson phase.  No
geometric fit is introduced.

Proof boundary: exact source reflection, half-domain Poisson representation,
and stationary classification only.  No half-domain nonstationary bound,
`B` crossing theorem, quantitative `A` splice, complete `T_upper`,
height-uniform theorem, `Lambda<=0`, PF-infinity, RH, or prize-level
conclusion is proved.
