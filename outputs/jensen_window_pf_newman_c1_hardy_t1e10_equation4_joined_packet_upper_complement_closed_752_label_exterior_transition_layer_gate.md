# Closed 752-label exterior transition layer

Date: 2026-08-27

Status: the exact dual-endpoint regrouping, its rigorous transition-packet
complex ball, and both eight-width far-launch tests are certified; the two far
integrals remain open; this is not a proof of the complete joined packet or RH.

The failed unsplit eight-round integration-by-parts route begins too close to
the endpoint saddles.  Replace that split by

```text
C=L+3/16=621.6875=621+11/16,
N=752=47*16.
```

At `L=621.5` the Fourier step is `-1`.  At `C` it is the primitive order-16
root `exp(2*pi*i*11/16)`, because `gcd(11,16)=1`.  Thus the same roster
vanishes exactly at both endpoints:

```text
sum_(j=0)^751 (-1)^j=0,
sum_(j=0)^751 exp(2*pi*i*j*C)=0.                 (ET1)
```

All 230 stationary labels lie in `[L,C]`; their extreme lower roots are
`[621.50016499467458464125532308222573458314325504102629961178993300508841771420086 +/- 5.92e-79]` and
`[621.55576081192507870131278288344926918563290127646748069839952373328749026286757 +/- 3.61e-78]`.  The packet

```text
G_752=sum_(q=q_+)^(q_++751)
      integral_L^C y^(-s)exp(-i*pi*y^2+2*pi*i*q*y)dy             (ET2)
```

is evaluated with one grouped geometric amplitude.  At 384 bits on 192 Arb
panels,

```text
G_752 = [0.00018532739869309836058169750757857242633607806086105562725497482925875668962334382 +/- 1.14e-55]
      + i [0.00027936104831354108310086432421683262755421947744964903790671465960471682153248850 +/- 1.14e-55],

|G_752| = [0.00033524444815864031122611645089865484790157372321229048443712932692246354379252809 +/- 1.58e-55].
```

For `Q(y)=y+p/y`, `p=t/(2*pi)`, define

```text
chi=delta/sqrt(abs(Q'(y))),
delta=q-Q(y).
```

The comparison `chi>8` is made without decimal threshold fitting by proving
`delta^2-64*abs(Q'(y))>0`.  For the finite continuation at `C`, the certified
scaled detuning is
`[8.4545838510485978847962962402297009495878377055946537648332498503320650292980522 +/- 3.89e-80]` and its squared
margin is `[30794.367607571050821427114358660767055821490722091384880204603751969531300339486 +/- 3.88e-76]`.
The preceding 1/16-grid point `L+2/16` fails, with margin
`[-182229.84045883826978416735980525167500694551419313066533962088748284717701073297 +/- 2.59e-75]`.

For the remote tail beginning at `q_++752` from `L`, the scaled detuning is
`[8.1380564649104406365297930127439823569359429535964720970136877720054769506660166 +/- 3.16e-80]` and its squared margin is
`[9177.8389118786377246123189318906665773212333964537792994049962355565303759151349 +/- 3.60e-77]`.  The preceding multiple
of 16, `736`, fails, with margin
`[-7280.4114396781667117444321467966205874629388544495203730807658771140360722083907 +/- 4.11e-77]`.  Monotonicity of
the positive detuning on the lower saddle branch makes these the first
1/16-grid split and the least multiple-of-16 roster meeting the threshold.

The exact upper-complement decomposition is therefore

```text
T_upper=G_752+F_C+I_L,                              (ET3)

F_C=sum_(q=q_+)^(q_++751) integral_C^a (...),
I_L=sum_(q>=q_++752)       integral_L^a (...).
```

Both retained pieces are genuinely far nonstationary candidates for the
already-certified endpoint-Hurwitz algebra and repeated integration by parts.
They are not bounded by this gate.

Pi provenance: every `pi` comes from the inherited Fresnel phase, exact
Fourier characters, or `p=t/(2*pi)`.  The split `3/16` comes from the first
rational Fourier grid point meeting the stated squared-detuning inequality;
`752` is the least multiple of its exact root order 16 meeting the remote
inequality.  No circle, polygon, fitted constant, or visual pattern supplies
`pi` here.

Proof boundary: exact dual-endpoint cancellation, the exact decomposition,
stationary-root census, rigorous `G_752` complex enclosure, and the two
eight-width launch inequalities only.  No quantitative enclosure of `F_C`,
`I_L`, the lower complementary tail, complete ordinary or joined packet,
`J_Z`, or `D_K` is proved, nor is any non-A, all-height, `Lambda<=0`,
PF-infinity, RH, or prize-level conclusion.
