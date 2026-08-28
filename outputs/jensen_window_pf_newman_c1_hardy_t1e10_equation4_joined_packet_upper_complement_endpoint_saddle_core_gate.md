# Grouped upper-complement endpoint-saddle core

Date: 2026-08-27

Status: the complete 230-label stationary core is rigorously enclosed on a
rational endpoint interval; its nonstationary continuation remains open.

The complementary half-lattice gate leaves exactly the upper labels

```text
q=q_+,q_++1,...,q_++229,       q_+=2561211.5.
```

All of their lower saddle roots lie in the rational interval

```text
L=621.5 <= y <= B=621.5625=L+1/16.
```

The extreme roots are `[621.50016499467458464125532308222573458314325504102629961178993300508841771420086 +/- 5.92e-79]` and
`[621.55576081192507870131278288344926918563290127646748069839952373328749026286757 +/- 3.61e-78]`.  The first-label root remains
below `B` by `[0.0067391880749212986872171165507308143670987235325193016004762667125097371324336003 +/- 2.17e-83]`, while at `B`
the nearest phase is already nonstationary by
`[27.756031598419002372836641297342834233954904672518891224228313374799248240762753 +/- 3.23e-79]`.

Because the roster length is even and `L` is a half integer, its grouped
Dirichlet amplitude has the exact endpoint zero

```text
sum_(j=0)^229 exp(2*pi*i*j*L)=sum_(j=0)^229 (-1)^j=0.
```

With `x=y-L`, the builder evaluates the complete finite block through the
single stable geometric amplitude

```text
(1-exp(2*pi*i*230*x))/(1+exp(2*pi*i*x))
```

and the endpoint-centered logarithmic phase.  Thirty-two Arb panels at 384
bits give

```text
P_core = [0.00035610961020398829219946871323619826184984813157583078802006037010285327589948581 +/- 1.03e-56]
       + i [0.00025364593807039957779411465513483891276542545784342859557571554590754766316996976 +/- 1.03e-56],

|P_core| = [0.00043720740659239690000111216470913579459817770338562357342509108297817709461923998 +/- 1.43e-56].
```

This yields the exact split `T_upper=P_core+R_upper`.  Every phase retained in
`R_upper` is nonstationary on its own interval: the 230 labels continue from
`B`, and labels from `q_++230` upward start at `L`.

Pi provenance: every `pi` comes from the inherited Fresnel phase, exact
half-integer Fourier spacing, or `p=t/(2*pi)`.  The rational split is
`L+1/16`; no fitted geometric constant is used.

Proof boundary: exact 230-label endpoint cancellation and rational split,
production saddle census, and rigorous grouped complex enclosure of `P_core`
only.  No quantitative enclosure of `R_upper`, the lower complementary tail,
complete ordinary or joined packet, `J_Z`, or `D_K` is proved, nor is any
non-A, all-height, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion.
