# Beta^-4 scaled Airy-branch amplitude envelope

Date: 2026-08-13
Status: rigorous saved-height envelope; not a proof of the ordinary splice

For branch sign `sigma in {+1,-1}`, set

```text
W_sigma(X)=Ai(-X)-i*sigma*Bi(-X),
W_sigma,X=d_X W_sigma,
xi=2X^(3/2)/3,

S_sigma(X,y)=e^(-i*sigma[xi-pi/4])
  [U W_sigma+V W_sigma,X]/2.                           (AE1)
```

By the exact Airy--Hankel connection in Section 11.378, `S_sigma` is the
complete beta-minus-four branch amplitude after extracting only the cubic
carrier `exp(i*sigma[xi-pi/4])`.  There is no asymptotic replacement in
(AE1).

On all 399 event cells and `0<=y<=64`, Section 11.379 gives

```text
[9.11373501846528081509557975247001911858917670407189115712808383608025308676795678667665257453894321481239e-5 +/- 3.71e-109] <=X<=[179.782621284236435878329113335532092186189432792103078919475031389492390462386608929729662870931897000187 +/- 1.06e-102],
|lambda|<116.                                           (AE2)
```

Direct rational bounds for (11.377.6) give

```text
|U-1|<1.6e-5,       |V|<5.38e-4,
|U_y|<7.1e-8,       |V_y|<3.05e-6.                    (AE3)
```

For negative Airy argument, the DLMF modulus definitions identify

```text
|W_sigma|=M(-X),       |W_sigma,X|=N(-X).              (AE4)
```

A deterministic Arb atlas uses 8,000 equal `log X` panels below `X=1` and
120,000 equal cubic-phase panels above it.  For the derivative modulus on
`X>=1`, DLMF 9.8.21 and its signed next-term remainder give

```text
N(-X)^2<=sqrt(X)/pi [1+7/(32X^3)].                    (AE5)
```

The right side has no interior maximum and its upper endpoint dominates at
the saved `X_max`.  Combining this with the direct low-`X` atlas gives,
without a numerical monotonicity assumption,

```text
|W_sigma|<[0.710032507114397415382001099715125747025012969970703125000000000000000000000000000000000000000000000000000 +/- 3.57e-8]<0.711,
|W_sigma,X|<[2.06591285917196858596799449919172193408804238934648758868715076704740206572831757172853809804182207996152 +/- 5.80e-105]<2.07.        (AE6)
```

The derivative after extracting the cubic carrier uses

```text
rho_sigma=W_sigma,X-i*sigma*sqrt(X)W_sigma.            (AE7)
```

The same interval atlas encloses the derivative after extracting the cubic
carrier,

```text
|rho_sigma|<[0.511777876876294612884521484375000000000000000000000000000000000000000000000000000000000000000000000000000 +/- 1.71e-6]<0.512. (AE8)
```

Since `W_sigma,XX=-X W_sigma`, differentiation of (AE1) is exact:

```text
S_sigma,y=e^(-i*sigma[xi-pi/4])/2
 {U_y W_sigma+V_y W_sigma,X
   +(U-i*sigma*sqrt(X)V)rho_sigma}.                   (AE9)
```

Combining (AE3), (AE5), and (AE7) proves uniformly

```text
|S_sigma|<[0.355576905294536463575169301042928996710157733826201884238193135616978756783054745755623797103437198318877 +/- 1.79e-8]<0.356,
|S_sigma,y|<[0.257739643802453347333714717057624321211842046796196278033065294538575616472879843484031488397341955831502 +/- 8.58e-7]<0.258. (AE10)
```

Pi provenance: the extracted `pi/4` is the forced Hankel phase, and the
`pi` inside `beta` and the event-cell endpoints comes from the Kummer/Fourier
normalization.  The Airy modulus definitions are recorded from `https://dlmf.nist.gov/9.8`;
the finite inequalities in (AE5)--(AE10) combine its signed asymptotic
remainder with direct Arb enclosures.

This gate controls the exact canonical branch amplitude on the compact fold
strip.  It does not yet compare that amplitude with the exact finite-height
logistic/Gamma carrier, integrate the opposite branch over every owned
corridor, close `Q_K-T` or `T_upper`, prove a height-uniform theorem,
`Lambda<=0`, PF-infinity, RH, or a prize-level result.
