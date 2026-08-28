# Finite branch Fresnel connection and mode roster

Date: 2026-08-27

Status: exact connection and mode roster certified; joined amplitude bound open

Write

```text
s=1/2+it,
c_alpha=pi*alpha*(1+i)/sqrt(2),
d=exp(-pi*t/4-i*pi/8),
K_0=exp(3pi*t/4+i3pi/8),                            (FC1)

I_-(alpha)=integral_0^infinity x^(-s)exp(-pi*x^2-c_alpha*x)dx,
I_+(alpha)=integral_0^infinity x^(-s)exp(-pi*x^2+c_alpha*x)dx,
J(alpha)=integral_0^infinity x^(-conj(s))
         exp(-pi*x^2+i*c_alpha*x)dx.                (FC2)
```

Scaling DLMF 12.5.1 to (FC2) and applying the minus-sign version of the exact
connection formula DLMF 12.2.19 gives, for every positive odd `alpha`,

```text
I_-(alpha)=-i exp(-pi*t)I_+(alpha)+kappa(t)J(alpha), (FC3)

kappa(t)=i sqrt(2pi)(2pi)^(it)exp(-pi*t/2)
         /Gamma(1/2+it),                            (FC4)

|kappa(t)|=sqrt(1+exp(-2pi*t)).                     (FC5)
```

The odd-label phase `exp(i*pi*alpha^2/4)=exp(i*pi/4)` is what makes (FC4)
independent of the label.  Three altered `(t,alpha)` rows directly verify
(FC3), (FC5), and the prefactor identity

```text
K_0[-i exp(-pi*t)] = d.                             (FC6)
```

Consequently the connection formula is exactly the finite-label split found
in Section 11.484:

```text
S_alpha=K_0 I_-(alpha)
       =d I_+(alpha)+K_0 kappa(t)J(alpha),

d I_+(alpha)=S_alpha-T_alpha,
K_0 kappa(t)J(alpha)=T_alpha.                       (FC7)
```

Unlike the rejected third-quadrant rotation of `I_-`, the `I_+` contour can
be rotated inside `0<=arg(x)<pi/4`, where the Gaussian supplies quadratic
decay.  Taking the endpoint as a Fresnel/Abel limit and integrating the tail
by parts gives the exact oscillatory representation

```text
d I_+(alpha)=integral_0^infinity[Fresnel]
             y^(-s)exp(i*pi*(alpha*y-y^2))dy.       (FC8)
```

The phase derivative is eventually `-2pi*y`, so (FC8) is a convergent
improper oscillatory integral; near zero the amplitude is `y^(-1/2)`.
Because the roster is finite, summing (FC8) is exact:

```text
B_W=integral_0^infinity[Fresnel]
    y^(-s)exp(-i*pi*y^2)D_W(y)dy,                   (FC9)

D_W(y)=sum_(j=0)^(M-1)exp(i*pi*(A+2j)y)
      =exp(i*pi*(A+M-1)y)sin(pi*M*y)/sin(pi*y),     (FC10)
```

with removable values supplied by the finite sum at integer `y`.  Equations
(FC3)--(FC10) own the branch multiplier and Stokes complement that were open
in Section 11.485.

The one-label stationary equation in (FC8) is

```text
2pi*y^2-pi*alpha*y+t=0,
y_+/-=(alpha+/-sqrt(alpha^2-8t/pi))/4.              (FC11)
```

At `t=10^10`, Arb proves

```text
y_-(B)=[621.5560036102371463889412039565834444235567629462757775115581974443559 +/- 3.84e-68],
y_-(A)=[39852.39138623123896184313678682614784426989001770014282206708149972340 +/- 4.38e-66],
y_+(A)=[39936.10861376876103815686321317385215573010998229985717793291850027660 +/- 4.38e-66],
sqrt(t/(2pi))=[39894.22804014326779399460599343818684758586311649346576659258296706579 +/- 2.60e-66]. (FC12)
```

The induced integer roster is therefore

```text
ordinary lower branch: 622..39852, 39231 modes,
transition block:       39853..39936, 84 modes,
transition split:       39853..39894 and 39895..39936,
                        42 modes on each side.       (FC13)
```

This independently recovers the existing physical ownership indices:
the extended P/Gamma roster is `622..39936`, the paired A roster is
`39853..39936`, the classical target ends at `39894`, and the extra Gamma
block is `39895..39936`.  The agreement certifies the route geometry and
index ownership, not a new amplitude estimate.  The next step is to express
the already-owned `G+A_transition` carriers on the common Fresnel integral
(FC9), subtract them before norms, and enclose only the joined remainder.

Pi provenance: DLMF's `pi` is the standard gamma/parabolic-cylinder
normalization; all project-specific occurrences descend from the RSI
Gaussian, quarter-turn, and Riemann-Siegel phase.  No fitted constant or
circle construction is introduced.

Proof boundary: exact parabolic-cylinder connection, branch/Stokes ownership,
finite Fresnel/Dirichlet representation, and actual saddle-mode roster only.
No common-integral carrier subtraction, joined amplitude bound, actual-height
`B_W`, `Q_K`, `D_K`, `J_Z`, non-A, all-height, `Lambda<=0`, PF-infinity, RH,
or prize-level enclosure is proved.
