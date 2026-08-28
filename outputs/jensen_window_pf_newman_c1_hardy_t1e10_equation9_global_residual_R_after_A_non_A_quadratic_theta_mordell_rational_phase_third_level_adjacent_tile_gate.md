# Rational-phase terminal and adjacent third-level tile

Date: 2026-08-25

Status: certified local recursive-tile lemma; finite-domain tiling remains open

## Generic rational-phase terminal

For rational `(a_0,tau_0)`, the phase `e(a_0 k+tau_0 k^2)` has exact period
`P` when

```text
(RP1) 2*tau_0*P is an integer,
      a_0*P+tau_0*P^2 is an integer.
```

Each power moment is grouped exactly by residues modulo `P`.  If
`rho=p-p_0`, `alpha=a-a_0`, and `beta=tau-tau_0`, the common evaluator is

```text
(RP2) B_j=M_j+rho F_j,

(RP3) W=B_0+2*pi*i(alpha B_1+beta B_2)+R_2,

(RP4) |R_2| <= 2*pi^2(delta_a^2 S_2
                       +2 delta_a delta_tau S_3
                       +delta_tau^2 S_4).
```

This reproduces the old alternating period-two terminal and both conjugate
period-three terminals.  The root-of-unity coefficients are stored exactly;
`pi` enters only through `e(u)=exp(2*pi*i*u)` and Taylor's exponential
remainder.

## Dependency cancellation

On the right branch the second affine weight simplifies before interval
evaluation:

```text
(RP5) p_2=2*tau_1*31916-a_1=A/2+s+1-2*31916.
```

Thus `p_2` is independent of `x`.  On the certified adjacent tile this reduces
the inherited interval radius by a factor
`159580800000000000000000000060000000000000000000001/800000000000000000000060000000000000000000001`.

## Next branch walls

After integer-curvature normalization and conjugation, the two open
period-three children have

```text
(RP6) t_R=(5x-2)/(2(1-2x)),
      t_L=(2-5x)/(2(3x-1)).
```

For child length `n=496283`, their floor walls are

```text
(RP7) x_R(q)=(2n+q)/(5n+2q),
      x_L(q)=(2n+q)/(5n+3q).
```

The first positive right wall is `992567/2481417`,
at distance `1/12407085` from `2/5`.
The adjacent `10^-23` tile is therefore strictly inside the `m_3=0` branch.
A uniform march at this width would require more than
`20000000000000000000000/2481417` tile units merely
to reach that first wall, so fixed-width enumeration is not a viable cover
strategy.

## Adjacent source tile

The certified box is

```text
x in [2/5+10^-23,2/5+2*10^-23],
s in [1/3-10^-23,1/3+10^-23].
```

Its rational-phase terminal has period `3` and deviation
below `0.0001217796417834267504436834839154357723601`.  The second recursive main,
its two joined Mordell endpoint currents, the first physical main, and both
physical endpoint currents are all retained.  The complete source radius is
below `0.003665095791802741587162017822265625` and its source box overlaps
the previously certified wall tile at their shared boundary.

Formally taking one more Mordell step would have child index zero, but its main
coefficient is already bounded only by
`3.19062159152189843618594266546176e+32`.
That `t^(-3/2)` inflation would require cancellation against the small-`t`
endpoint terms.  The rational-phase terminal avoids discarding that
cancellation; it does not claim a useful third Mordell decomposition.

## Scale audit

| width unit | delta p | delta a | delta tau | child deviation | quadratic remainder | p-radius reduction |
|---:|---:|---:|---:|---:|---:|---:|
| `1e-20` | `1/100000000000000000000` | `23936810000000000000000003/2999999999999999999400000000000000000000` | `5/19999999999999999996` | `0.1254796850040465483466789464728208258748` | `0.003703746967587377356134759054384630871937` | `159580800000000000000000060000000000000000001/800000000000000000060000000000000000001` |
| `1e-21` | `1/1000000000000000000000` | `239368100000000000000000003/299999999999999999994000000000000000000000` | `5/199999999999999999996` | `0.01221463127332179064710171445540254353546` | `0.00003703746967587377301924650430109409171564` | `15958080000000000000000000600000000000000000001/80000000000000000000600000000000000000001` |
| `1e-22` | `1/10000000000000000000000` | `2393681000000000000000000003/29999999999999999999940000000000000000000000` | `5/1999999999999999999996` | `0.001218129755061350336098335844781104242429` | `0.000000370374696758737745015541619961196850852` | `1595808000000000000000000006000000000000000000001/8000000000000000000006000000000000000000001` |
| `1e-23` | `1/100000000000000000000000` | `23936810000000000000000000003/2999999999999999999999400000000000000000000000` | `5/19999999999999999999996` | `0.0001217796417834267504436834839154357723601` | `0.000000003703746967587377714853212216580824467371` | `159580800000000000000000000060000000000000000000001/800000000000000000000060000000000000000000001` |
| `1e-24` | `1/1000000000000000000000000` | `239368100000000000000000000003/299999999999999999999994000000000000000000000000` | `5/199999999999999999999996` | `0.00001217763084111559083034907746689867735768` | `3.70374696758737768900381807429870962764e-11` | `15958080000000000000000000000600000000000000000000001/80000000000000000000000600000000000000000000001` |
| `1e-26` | `1/100000000000000000000000000` | `23936810000000000000000000000003/2999999999999999999999999400000000000000000000000000` | `5/19999999999999999999999996` | `0.0000001217759417402061204792167170563033096187` | `3.703746967587377820270272703075074048151e-15` | `159580800000000000000000000000060000000000000000000000001/800000000000000000000000060000000000000000000000001` |

## Proof boundary

One exact rational-phase Taylor terminal theorem, exact next-wall formulas for the two local period-three children, and one endpoint-complete adjacent right tile of width unit 1e-23 only. No finite recursive cover, useful third Mordell decomposition at small tau, physical quadrature, non-A bound, joined R_after_A, R_Dir, Q_K-T, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.
