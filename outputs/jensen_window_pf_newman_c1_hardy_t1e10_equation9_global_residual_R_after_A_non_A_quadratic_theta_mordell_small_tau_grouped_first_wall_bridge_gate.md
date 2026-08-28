# Grouped small-t first-floor-wall bridge

Date: 2026-08-25

Status: certified local `m3=0|1` wall crossing; finite recursive cover remains open

## Exact grouped jump

For the affine current `W_n(p;a,t)`, the third Mordell floor changes from
`m=0` to `m=1` at `2*n*t=1`, with `n=496283` and

```text
t*=1/992566.
```

The new child term and the change of the upper joined endpoint must be kept
together.  The two upper Mordell arguments satisfy `u_+(0)=u_+(1)+1`.
Applying the exact unit recurrence simultaneously to `h` and `h_z` reduces
the endpoint coefficient to `(p'+1)/t`, exactly the coefficient of the new
child term.  Their phase difference is

```text
2*n=992566,
```

an even integer in the `exp(i*pi*phase)` convention.  Hence on both physical
sides

```text
Delta upper endpoint = - Delta third main.
```

At the two wall anchors the child affine weights are respectively `1/3` and
`-1/3`.  This proves that the floor switch is removable only after the main
and endpoint currents are grouped; it does not bound either piece separately.

## Common finite-current enclosure

Instead of evaluating the ill-conditioned small-`t` Mordell pieces, both
floor branches are represented by the analytic finite current

```text
W_n(p;a,t)=sum_(k=0)^n (p+k) exp(2*pi*i*(a*k+t*k^2)).
```

Four base moments are computed rigorously at each exact wall anchor.  The
production evaluator resets the exact quadratic phase every
`16` terms, so rectangular complex-ball wrapping cannot
accumulate across the full sum.  For `rho=p-p0`, `alpha=a-a0`, and
`beta=t-t0`, the common enclosure is

```text
B0 + 2*pi*i*(alpha*B1+beta*B2) + R2,

|R2| <= 2*pi^2*(delta_a^2*S2
                 +2*delta_a*delta_t*S3
                 +delta_t^2*S4).
```

This is one enclosure of the full current across `m=0|1`, rather than a hull
of separately bounded Mordell components.

## Certified source boxes

Both physical boxes use x half-width `1/1000000000000000000000000` and s half-width
`1/1000000000000000000000000` around their exact first-wall anchors.  All earlier Mordell
mains and joined endpoints are retained through the physical source.

| side | child deviation | source radius | target |
|---|---:|---:|---:|
| right | `0.0002502150705942783907517745767279393476201` | `0.010659942519851028919219970703125` | `<0.06` |
| left | `0.002002960256450679946138571096980740549043` | `0.0577669341000728378543449537119158776477` | `<0.06` |

The independent checker uses 384 bits, block size 13, `Y=11`, and tolerance
`1e-65`; it recomputes every base moment with different reset boundaries,
rederives the exact jump algebra, and replays both complete source boxes.

## Boundary

The exact m3=0|1 third-main/upper-endpoint jump cancellation and one endpoint-complete physical source box of x half-width 1e-24 and s half-width 1e-24 around each first positive third-floor wall only. The common finite-current Taylor enclosure avoids separately estimating small-t Mordell pieces. No cover between the normalization wall and these first floor walls, continuous adaptive-width theorem, subsequent floor wall, physical non-A quadrature, non-A bound, joined R_after_A, R_Dir, Q_K-T, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.
