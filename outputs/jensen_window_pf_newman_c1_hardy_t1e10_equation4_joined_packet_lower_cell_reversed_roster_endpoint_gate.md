# Reversed-roster endpoint certificate for the lower finite cell

Date: 2026-08-27

Status: rigorous production lower-cell complex enclosure certified; ordinary
Gamma-subtracted join remains open.

Put `s=1/2+it`, `p=t/(2*pi)`, `L=621.5`, and reverse the odd finite roster as

```text
D_W(y)=sum_(r=1)^M exp[2*pi*i(q_+-r)y],
q_+=2561211.5,       M=2481423.                     (LC1)
```

For the `r`th label absorb `y^(-it)` into the phase.  Its derivative is

```text
phi_r'(y)=-2*pi*d_r(y),
d_r(y)=p/y+y-q_++r.                                 (LC2)
```

On `0<y<=L`, `d_r` decreases and the nearest label satisfies

```text
d_1(L)=[230.67967651384986136385152879102227610049461715927188523481993397904479849613982 +/- 4.85e-78],
|phi_1'(L)|=[1449.4031541367613813287260204170388508006059462702672448196113428978683551005089 +/- 4.64e-77].
```

Thus no lower-cell label has a stationary point.  At the half-integer endpoint,

```text
exp[2*pi*i(q_+-r)L]=exp(2*pi*i*q_+*L)(-1)^r.        (LC3)
```

Starting with `h_0(y)=y^(-1/2)`, define

```text
h_(n+1,r)=d/dy[h_(n,r)/d_r].                        (LC4)
```

Repeated integration by parts is exact for the finite roster.  Every zero-end
boundary term vanishes, while every `L`-end term is a finite linear combination
of

```text
sum_(r=1)^M (-1)^r/(r+a)^k,
a=p/L+L-q_+=[229.67967651384986136385152879102227610049461715927188523481993397904479849613982 +/- 4.85e-78].   (LC5)
```

These finite alternating sums are evaluated without iterating the roster by
the shifted eta/Hurwitz identity

```text
sum_(r=1)^M (-1)^r/(r+a)^k
 =-eta_k(a+1)-eta_k(a+1+M),                         (LC6)

eta_k(z)=2^(-k)[zeta(k,z/2)-zeta(k,(z+1)/2)],       k>1,
eta_1(z)=[psi((z+1)/2)-psi(z/2)]/2.
```

The order-18 endpoint sum is

```text
[3.7858004113710693147087482562564388530657859837305712373852139310276415033137470e-6 +/- 2.79e-86]
+ i [1.3336666061723477872466968106445732537183169685396993387709478448318633284642851e-5 +/- 4.62e-85].
```

For the discarded integral, split at `y=300`.  Below the split,
`d_r>=d_1>=p/(2y)` gives a termwise exact-monomial integral.  Above it,
4096 monotone slabs use

```text
sum_(r=1)^M d_r^(-k)
 <=min[M d_1^(-k), d_1^(-k)+d_1^(1-k)/(k-1)].       (LC7)
```

The certified remainder is

```text
small-y: [1.8772031095951948680344038303794069914666396894854731393349866206257210162595984e-139 +/- 4.04e-219]
upper:   [2.2999900816639904578846322633575209547735068587333424457584597439764680969354272e-16 +/- 3.66e-96]
total:   [2.2999900816639904578846322633575209547735068587333424457584597439764680969354272e-16 +/- 3.66e-96].
```

Therefore the complete lower finite cell is

```text
V_L = [3.7858004113710693147087482562564388530657859837305712373852139310276415033137470e-6 +/- 2.30e-16]
    + i [1.3336666061723477872466968106445732537183169685396993387709478448318633284642851e-5 +/- 2.30e-16].           (LC8)
```

This is a complex ball retained for later addition.  No absolute value of
`V_L` is substituted into the joined packet.

Pi provenance: every `pi` in (LC1)--(LC8) comes from the inherited quadratic
Riemann-Siegel/Fresnel phase, the exact Fourier spacing, or `p=t/(2*pi)`.
No geometric fit or inserted circle constant is used.

Proof boundary: exact reversed-roster identity, nonstationary lower-cell phase,
finite alternating-Hurwitz endpoint expansion, and rigorous production complex
lower-cell enclosure only.  This gate does not yet enclose the ordinary
Gamma-subtracted packet, the complete joined packet, `J_Z`, or `D_K`, and it
does not prove a non-A, all-height, `Lambda<=0`, PF-infinity, RH, or prize-level
conclusion.
