# Eight-round closure of the far upper-complement remainders

Date: 2026-08-27

Status: both far upper remainders and the complete upper complement at the
production height are rigorously enclosed as complex balls; the lower
complement and complete joined packet remain open; this is not a proof of RH.

Section 11.502 supplies the exact decomposition

```text
T_upper=G_752+F_C+I_L,                              (FR1)
```

where `C=621+11/16`, `F_C` contains 752 labels on `[C,a]`, and `I_L`
starts at `q_++752` on `[L,a]`.  Both launches exceed eight quadratic widths
and the full interval lies below `sqrt(p)`, so `Q` decreases throughout.

The endpoint sums are evaluated before any norm.  At `C`, the finite roster
uses the primitive order-16 root with numerator 11 and 47 cycles; at `a` it
uses 376 alternating cycles.  The remote endpoints at `L` and `a` use the
alternating infinite Lerch/digamma formulas.  Every finite formula through
power 16 overlaps its direct 752-term Arb sum, and every infinite row obeys
the shifted Lerch recurrence, including exact cancellation of the `k=1`
Hurwitz poles.

Apply eight exact integration-by-parts rounds using

```text
h_0=y^(-1/2),       h_(n+1)=d/dy[h_n/(q-y-p/y)].    (FR2)
```

On 8192 fourth-power endpoint-clustered slabs, the retained complex results
are

```text
F_C endpoint partial = [3.87555794344493853928813960640579152234872720532556666090747967244141497e-6 +/- 4.48e-78]
                     + i [-1.38592282533207387796164545570979502234037469127378212911190503787055996e-6 +/- 2.33e-78],
F_C remainder <= [1.74059857308995518791986489757469035948411640896219199320366540758648383e-17 +/- 4.08e-89],
F_C ball = [3.87555794344493853928813960640579152234872720532556666090747967244141497e-6 +/- 1.75e-17]
         + i [-1.38592282533207387796164545570979502234037469127378212911190503787055996e-6 +/- 1.75e-17];

I_L endpoint partial = [-1.61172374629028454582842978335087643309643634449050430664513319699000654e-6 +/- 1.80e-78]
                     + i [-5.90075204188353427680533261072081862451883296946478993356891817790393831e-6 +/- 6.25e-78],
I_L remainder <= [2.97311462255729774535247807324865831673383605381719226191011596912730600e-17 +/- 7.55e-90],
I_L ball = [-1.61172374629028454582842978335087643309643634449050430664513319699000654e-6 +/- 2.98e-17]
         + i [-5.90075204188353427680533261072081862451883296946478993356891817790393831e-6 +/- 2.98e-17].
```

Joining these balls to the certified transition packet before taking a norm
gives the complete upper complement

```text
T_upper = [0.000187591232890253014575157217401627341425330351721890689609237175734208098 +/- 4.74e-17]
        + i [0.000272074373446325474946097346150402013907360269788910465844033836388942323 +/- 4.74e-17],

|T_upper| = [0.000330476830267260523244164958366910454790055240437141037546098232269287109 +/- 6.58e-17].               (FR3)
```

The independent checker changes to nine recurrence rounds, 10000
fifth-power-clustered slabs, direct finite endpoint sums, and explicit
alternating eta values.  Both changed-order family balls and their complete
upper join overlap the production enclosures.  A separate altered-height
positive-denominator direct quadrature overlaps its fifth-order recurrence
ball.

Pi provenance: every `pi` comes from the inherited Fresnel/Fourier phase,
`phi'=2*pi*D`, exact root-of-unity characters, or `p=t/(2*pi)`.  The endpoint
`C` and count 752 are inherited from the exact order-16/eight-width gate.  No
circle, polygon, fitted constant, or visual pattern supplies `pi`.

Proof boundary: rigorous endpoint tables through power 16, eight-round
complex enclosures of `F_C` and `I_L`, and the complete production-height
upper-complement ball only.  No lower complementary-tail, complete ordinary
or joined packet, `J_Z`, `D_K`, non-A, all-height, `Lambda<=0`, PF-infinity,
RH, or prize-level conclusion is proved.
