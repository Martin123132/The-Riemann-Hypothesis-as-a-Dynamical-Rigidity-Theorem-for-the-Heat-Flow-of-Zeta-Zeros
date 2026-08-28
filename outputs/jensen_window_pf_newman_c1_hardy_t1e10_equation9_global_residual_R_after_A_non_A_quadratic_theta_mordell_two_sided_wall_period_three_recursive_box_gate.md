# Two-sided affine Mordell wall and period-three recursive box

Date: 2026-08-24

Status: certified two-sided local interval lemma; recursive tiling remains open

## Exact left normalization

For `x<2/5`, the right coordinate has `tau_R>1/4`. Apply the exact parity
half-shift and conjugation:

```text
(TW1) a_L=1/2-a_R=3/2+A/2+s-31916/x,
      tau_L=1/2-tau_R=3/2-1/(2x).
```

The new conjugation cancels the first conjugation, so the normalized left
parent is already in the original transformed-current orientation. At the
wall `(a_L,tau_L)=(1/3,1/4)`.

## Left child and orientation witness

The open left child has index `496283` and tends to

```text
(p_0,a_0,tau_0)=(47873/3,2/3,-1).
```

It uses the conjugate period-three root `omega^2`. Exact residue pairs and the
same cancellation-preserving Taylor form give

```text
(TW2) W=B_0(rho)+2*pi*i[alpha B_1(rho)+beta B_2(rho)]+R_2,

(TW3) |R_2|<=2*pi^2[delta_a^2 S_2
                     +2 delta_a delta_tau S_3
                     +delta_tau^2 S_4].
```

At width `10^-23`, the left child deviation is below
`0.0000608889227922001195769663439705254859291`.

A separate `496284` left-wall branch is completed with its own Mordell main
and joined endpoints. It overlaps the conjugate of the independently
certified right-wall branch, proving the orientation match before hulling.

## Two-sided source box

The completed transformed branches satisfy

```text
(TW4) T_two=hull(T_left_open,T_left_wall,
                 conjugate(T_right_wall),conjugate(T_right_open)),

(TW5) SourceBox=P_1*T_two+Endpoint_1.
```

The certified physical box is

```text
x in [2/5-10^-23,2/5+10^-23],
s in [1/3-10^-23,1/3+10^-23].
```

Its complete source radius is below `0.005268418157356791198253631591796875`
and it contains the exact centre source current.

## Scale audit

| half-width | max delta p | max delta a | max delta tau | left child deviation | two-level propagation |
|---:|---:|---:|---:|---:|---:|
| `1e-20` | `7979039999999999999999999/3999999999999999999900000000000000000000` | `7979039999999999999999999/1999999999999999999950000000000000000000` | `5/39999999999999999994` | `0.06181393496016751587518456290126778185368` | `1.382201605276853317860741299227811396122` |
| `1e-21` | `79790399999999999999999999/399999999999999999999000000000000000000000` | `79790399999999999999999999/199999999999999999999500000000000000000000` | `5/399999999999999999994` | `0.00609805906647013871296758580342611821834` | `0.1363567460343613757522973628510953858495` |
| `1e-22` | `797903999999999999999999999/39999999999999999999990000000000000000000000` | `797903999999999999999999999/19999999999999999999995000000000000000000000` | `5/3999999999999999999994` | `0.0006089725623515477572963994212784655246651` | `0.01361704045850290097086077167887196992524` |
| `1e-23` | `7979039999999999999999999999/3999999999999999999999900000000000000000000000` | `7979039999999999999999999999/1999999999999999999999950000000000000000000000` | `5/39999999999999999999994` | `0.0000608889227922001195769663439705254859291` | `0.001361517704400957507263059120816706126789` |
| `1e-24` | `79790399999999999999999999999/399999999999999999999999000000000000000000000000` | `79790399999999999999999999999/199999999999999999999999500000000000000000000000` | `5/399999999999999999999994` | `0.000006088808944790465003146610889483980599834` | `0.0001361499070256024261291183385935710248305` |
| `1e-26` | `7979039999999999999999999999999/3999999999999999999999999900000000000000000000000000` | `7979039999999999999999999999999/1999999999999999999999999950000000000000000000000000` | `5/39999999999999999999999994` | `0.00000006088799778003214949530940867591932708081` | `0.000001361497020500081686091541077754385469234` |
| `1e-28` | `797903999999999999999999999999999/39999999999999999999999999990000000000000000000000000000` | `797903999999999999999999999999999/19999999999999999999999999995000000000000000000000000000` | `5/3999999999999999999999999994` | `6.088799686335342734405182722669681166794e-10` | `0.00000001361497000002522423192216676152491161567` |

## Validation

- Production uses 320 bits, Mordell cutoff 10, and one numerical worker.
- The independent checker derives the double-conjugation map and both floor
  branches, verifies the `omega^2` residue pairs against direct witnesses, and
  replays the four completed parent branches and source assembly at 384 bits.
- Every `pi` comes from `e(u)=exp(2*pi*i*u)` or the inherited exact Mordell
  transform; no circle or polygon constant is inserted.

## Proof boundary

One two-sided width-1e-23 physical box at the x=2/5 normalization wall, the exact left parity/double-conjugation map, both 496283 open child branches, two 496284 wall orientation representations, and their endpoint-complete source hull only. No recursive tiling beyond this local box, physical quadrature, non-A bound, joined R_after_A, R_Dir, Q_K-T, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.
