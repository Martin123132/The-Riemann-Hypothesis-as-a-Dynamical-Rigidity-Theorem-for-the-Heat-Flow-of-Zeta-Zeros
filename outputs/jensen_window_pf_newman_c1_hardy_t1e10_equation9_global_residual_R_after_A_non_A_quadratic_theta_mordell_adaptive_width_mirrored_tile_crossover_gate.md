# Mirrored adaptive-width Mordell tile crossover

Date: 2026-08-25

Status: certified local width-lattice crossover; finite cover remains open

## Exact mirrored map

The right child weight from Section 11.470 is independent of `x`.  The left
map has the matching cancellation

```text
(AW1) p_(2,L)=3*31916-3/2-A/2-s.
```

After the integer curvature shift and conjugation, the left normalized child
is

```text
(AW2) a_(3,L)=1-[(3/2+A/2+s)x-31916]/(3x-1),
      t_(3,L)=x/[2(3x-1)]-1.
```

Its parent main, both joined Mordell endpoint currents, the physical main,
and both physical endpoint currents are retained in the same original-current
orientation as Section 11.469.

## Local width search

The fixed `s` half-width is `10^-23`.  For integer `M>=2`, the right search
uses `[2/5+10^-23,2/5+M*10^-23]` and the left search uses
`[2/5-M*10^-23,2/5-10^-23]`.  The declared complete-source radius target is
`0.01`.  Every row is a rigorous source
box; a failing row means only that this enclosure does not certify the target.

| side | units | source radius | terminal propagation | target met |
|---|---:|---:|---:|---|
| right | `2` | `0.003665095791802741587162017822265625` | `0.002723075573033298477781949387122040207032` | `True` |
| right | `3` | `0.007330240558076183375602052905151140294038` | `0.00408467526704316304109187996118635055609` | `True` |
| right | `4` | `0.01098609041946474636219921450219771941192` | `0.005446316370188182716394198479292754200287` | `False` |
| left | `2` | `0.002987630454299505366327904809509163897019` | `0.00272307296849107682498725147013374225935` | `True` |
| left | `3` | `0.0059752980232588015496730804443359375` | `0.004084671360112302612721535410855722147971` | `True` |
| left | `4` | `0.008955457247793674468994140625` | `0.005446311160790287021760125441005584434606` | `True` |
| left | `5` | `0.011950860658544115722179412841796875` | `0.006807992370525038292039532450417027575895` | `False` |

The largest passing right and left lattice widths are respectively
`3` and
`4` units.  Their immediate next integer
candidates are explicitly certified but fail the declared radius target.
This is lattice maximality under the current wall-anchored Taylor enclosure,
not continuous or method-independent maximality.

## First positive floor-wall anchor

At `m_3=1`, both sides have `t=1/(2*496283)`.  The exact normalized anchors
have minimal phase periods

```text
right: 1488849=3*496283,
left:  2977698=6*496283.
```

Both periods exceed the current length `496283`, so residue-period compression
has no complete cycle there.  Conversely, the next Mordell step contracts the
active index to one.  The missing bridge is therefore a cancellation-preserving
small-`t` endpoint expansion or grouped third step, not further fixed-wall
Taylor tiling.

## Validation and boundary

The production search uses 320 bits, `Y=10`, and tolerance `1e-55`.  The
independent checker rederives the left map and Jacobian, verifies both wall
anchor periods, and replays each passing/failing boundary pair at 384 bits,
`Y=11`, and tolerance `1e-65`.

Exact mirrored left child and endpoint formulas, exact local Jacobians and first-wall anchor periods, and discrete maximal passing/failing complete-source tiles on the declared 1e-23 integer lattice for the diagnostic radius target 0.01 only. No continuous maximal-width theorem, finite recursive cover, cancellation-preserving third Mordell endpoint expansion, physical quadrature, non-A bound, joined R_after_A, R_Dir, Q_K-T, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.
