# Right-side affine Mordell wall and period-three recursive box

Date: 2026-08-24

Status: certified right-side local interval lemma; the left side remains open

## Status

Certified at `t=10^10` on

```text
x in [2/5,2/5+10^-23],
s in [1/3-10^-23,1/3+10^-23].
```

This is one side of the interior normalization wall. It is not a complete
two-sided neighbourhood, a physical quadrature, a non-A bound, or an RH-level
result.

## Exact wall partition

After the first conjugation and integer shifts,

```text
(RW1) a=31916/x-(A+2s)/2-1,
      tau=1/(2x)-1,
      W_original=conjugate(W_992568(31916;a,tau)).
```

At `x=2/5`, `tau=1/4` and
`floor(2*992568*tau)=496284`. For every `x>2/5` in the box,
`tau<1/4` and the floor is `496283`. The closed source box is therefore the
union of the wall line and the open-right branch, not one silently frozen
floor branch.

## Period-three child enclosure

Both branches have child parameters near

```text
(p_0,a_0,tau_0)=(95747/6,1/3,-1).
```

Since `exp(2*pi*i*(k/3-k^2))=omega^k`, where
`omega=exp(2*pi*i/3)`, every base moment is represented exactly as
`u+v*omega` by splitting `k` modulo three. If `rho=p-p_0`,
`alpha=a-a_0`, and `beta=tau+1`, then

```text
(RW2) W=B_0(rho)+2*pi*i[alpha B_1(rho)+beta B_2(rho)]+R_2,
      B_j(rho)=M_j+rho F_j.
```

Only the quadratic phase remainder is bounded absolutely:

```text
(RW3) |R_2|<=2*pi^2[delta_a^2 S_2
                     +2 delta_a delta_tau S_3
                     +delta_tau^2 S_4].
```

For the open-right branch, the total period-three child deviation is below
`0.0000608889426035426342622043638375117780015`.

## Endpoint-complete recursion

Each fixed-floor branch uses the affine Mordell identity with its own child
length, main coefficient, and both joined `h,h_z` endpoint currents. The
open branch is bounded on its closed analytic parameter hull but is asserted
only for `x>2/5`; the wall line supplies the excluded endpoint. Their complete
normalized-parent balls overlap and are joined by interval hull:

```text
(RW4) T_right = hull(T_m=496283, T_m=496284),
      SourceBox=P_1*conjugate(T_right)+Endpoint_1.
```

The complete physical source radius is below
`0.0036649473840952850878238677978515625` and contains the exact source current
at `(x,s)=(2/5,1/3)`.

## Scale audit

| half-width | max delta p | max delta a | max delta tau | child deviation | two-level propagation |
|---:|---:|---:|---:|---:|---:|
| `1e-20` | `7979040000000000000000001/4000000000000000000100000000000000000000` | `7979040000000000000000001/2000000000000000000050000000000000000000` | `5/39999999999999999996` | `0.06181395513255632384597149098226509522647` | `1.382202056345179563834335567662492394447` |
| `1e-21` | `79790400000000000000000001/400000000000000000001000000000000000000000` | `79790400000000000000000001/200000000000000000000500000000000000000000` | `5/399999999999999999996` | `0.006098061051182326712527537182495507295243` | `0.1363567904138750785936196052716695703566` |
| `1e-22` | `797904000000000000000000001/40000000000000000000010000000000000000000000` | `797904000000000000000000001/20000000000000000000005000000000000000000000` | `5/3999999999999999999996` | `0.000608972760497499552082012463927185308421` | `0.01361704488918107774919619146203331183642` |
| `1e-23` | `7979040000000000000000000001/4000000000000000000000100000000000000000000000` | `7979040000000000000000000001/2000000000000000000000050000000000000000000000` | `5/39999999999999999999996` | `0.0000608889426035426342622043638375117780015` | `0.001361518147396043650759756005186318361666` |
| `1e-24` | `79790400000000000000000000001/400000000000000000000001000000000000000000000000` | `79790400000000000000000000001/200000000000000000000000500000000000000000000000` | `5/399999999999999999999996` | `0.00000608881092589218972886948993972211496839` | `0.0001361499513243836763463218142433674984204` |
| `1e-26` | `7979040000000000000000000000001/4000000000000000000000000100000000000000000000000000` | `7979040000000000000000000000001/2000000000000000000000000050000000000000000000000000` | `5/39999999999999999999999996` | `0.00000006088801759101361702205499346410300098853` | `0.000001361497463487094360462472016071178160246` |
| `1e-28` | `797904000000000000000000000000001/40000000000000000000000000010000000000000000000000000000` | `797904000000000000000000000000001/20000000000000000000000000005000000000000000000000000000` | `5/3999999999999999999999999996` | `6.08880166743345367015921765538972662557e-10` | `0.00000001361497442989527011045479296070734309687` |

## Validation

- Production uses 320 bits, Mordell cutoff 10, and one numerical worker.
- The checker derives the branch floors independently, checks every residue
  power sum against direct witnesses, verifies both base period-three currents,
  and replays both endpoint-complete branches at 384 bits and cutoff 11.
- Every `pi` comes from `e(u)=exp(2*pi*i*u)` or the inherited exact Mordell
  transform; no circle or polygon constant is inserted.

## Proof boundary

One right-sided width-1e-23 physical box at the x=2/5 normalization wall, its exact split between child indices 496283 and 496284, two endpoint-complete affine Mordell branch boxes, and their complete source-current hull only. No left-side branch, full two-sided interior neighbourhood, recursive tiling, physical quadrature, non-A bound, joined R_after_A, R_Dir, Q_K-T, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.
