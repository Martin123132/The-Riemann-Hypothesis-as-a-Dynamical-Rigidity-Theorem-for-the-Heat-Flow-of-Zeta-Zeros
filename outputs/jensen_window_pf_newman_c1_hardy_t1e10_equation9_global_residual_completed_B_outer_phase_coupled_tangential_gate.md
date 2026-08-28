# Phase-coupled tangential bound for the complete analytic outer B block

Date: 2026-08-13

Status: saved-height, cutoff-uniform and Abel-uniform analytic outer-block
bound; not a proof of the complete completed-B estimate

The three preceding gates bound, for every `m>=B=5122421`, the paired first
three rational currents, the exact two-sign Fresnel remainder, and the exact
Fresnel-to-boundary dictionary through two derivatives.  Let their joined
Kummer-weighted amplitude be `U(x)`, and put

```text
A(x)=eta_delta(x) U(x),
h(x)=Phi_B'(x)=pi B^2/4-t/[2x(1-x)].                 (OP1)
```

The exterior domain is `[delta,x_low] union [x_high,1/2]`.  Two exact
tangential integrations by parts give boundary majorants

```text
|A/h|,
|A'|/|h|^2+|A||h'|/|h|^3,                           (OP2)
```

and integral remainder

```text
integral { |A''|/|h|^2+3|A'||h'|/|h|^3
           +|A||h''|/|h|^3+3|A||h'|^2/|h|^4 } dx. (OP3)
```

The amplitude and phase are not replaced by unrelated global suprema.
Instead, (OP3) is enclosed on `45` deterministic Arb slabs:
`16` across the endpoint cutoff, `8`
on the left exterior, and `21` on the right exterior.
The exterior slabs are dyadic in normalized distance from the B saddle, so
the growth of the outer amplitude remains coupled to the much faster growth
of `|h|`.

All actual exterior boundary terms are retained.  The cutoff has
`A(delta)=A'(delta)=0`; the remaining three boundaries contribute

```text
sum |A/h| <= [1.853131996579961802983005733396561662989932684026019694869423779833932521281748856375774283892589714e-7 +/- 3.49e-102],
sum second-boundary <= [3.176298426940954484070060334754488309584155666210105492651602055332633780216535489124139847149765601e-11 +/- 2.14e-105].         (OP4)
```

The complete integral remainder is

```text
[1.093392166959339003678436259879049835623893098138848381265324004965112908475938709224102184479070136e-10 +/- 1.38e-104].                         (OP5)
```

After the exact equation-(9) half-domain projection factor,

```text
E_outer,m>=B <= [6.565486022379221416510266647544862923073521082626919637883972959976130875503541450674936869045802908e-10 +/- 1.25e-104]
              < 8e-10.                 (OP6)
```

This proves the complete analytic `m>=B` contribution to the released
completed-B exterior current is negligible at the saved height, uniformly
in every finite cutoff and Abel weight bounded by one.  No process or mode
enumeration is involved.  The same derivative majorants are summable and
independent of the regulator, so dominated convergence also gives this
outer block its own Abel limit without separating any unsummed endpoint
piece.

The finite block `m<B`, including its exact crossing completion, is still
open and must be compressed before norms.  The common-regulator joined
remainder is also open.  Thus this does not prove a complete exterior-B
bound, the `1.4058e-4` joined target, `T_upper`, a height-uniform theorem,
`Lambda<=0`, PF-infinity, RH, or a prize-level conclusion.

Pi provenance: `pi` in (OP1)--(OP6) is inherited from the equation-(9)
Fresnel phase and exact physical projection.  No fitted geometric constant
is introduced.
