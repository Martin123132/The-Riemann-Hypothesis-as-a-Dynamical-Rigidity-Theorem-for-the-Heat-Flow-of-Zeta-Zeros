# Order-Eleven Curvature Bridge Target

Date: 2026-07-22

Status: the former continuum target is discharged and the resulting
fixed-order-eleven composition is independently checked. This is not a proof
of PF-infinity, RH, or `Lambda <= 0`.

## Coordinate

```text
X(t)=9*B(t)-z(t-1)+2*z(t)-z(t+1)
y(t)=2*z(t)-w(t)+log(1-exp(-X(t)))
Q_(10,n)=A_(n+9)^10*exp(y(n+9))
E_n=log(Q_(10,n)*Q_(10,n+2)/Q_(10,n+1)^2)=10*log(x_k)+Y_k, Y_k=y(k-1)-2*y(k)+y(k+1), k=n+10
```

## Transfer

Under the inherited order-ten `z` curvature theorem,

```text
min(X_j,X_j^(1))>1/j, j>=1252
phi(min(X_j,X_j^(1)))<j
|Y_k-Y_k^(1)|<37/k^2 for every integer k>=1253
exact scaled transfer=3.6000240211299395138963759577008330105802313148259E+1<37
```

The power envelope has eighteen exact rational rows.

## Endpoint Budget

```text
y_1''(t)<=6000/t^2 for every real t>=1252
y_1''(t)<=6000/t^2 => Y_k^(1)<=6000*[-log(1-1/k^2)]<6001/k^2, k>=1253
[z_1''(t)<=4200/t^2 on t>=1251 and y_1''(t)<=6000/t^2 on t>=1252] => Y_k<6001/k^2+37/k^2=6038/k^2<6100/k^2, k>=1253
-10*log(x_k)>=10*d_k>=2510/(250*(2*k+1)), k>=320
6100/k^2<2510/(250*(2*k+1)), k>=1253
[z_1''(t)<=4200/t^2 on t>=1251 and y_1''(t)<=6000/t^2 on t>=1252] => Q_(11,n)(-100)>0 for every n>=1243
Q_(11,n)(-100)>0 for every 0<=n<=1242
[Q_(10,n)(lambda)>0 for every n>=4 and -100<=lambda<=0 and Q_(11,n)(-100)>0 for every n>=4] implies Q_(11,n)(lambda)>0 for every n>=4 and -100<=lambda<=0
```

## Resolution

The finite endpoint prefix through `n=1242`, the delayed order-ten heat ray,
and the shifted order-eleven heat theorem are rigorous inputs. The hash-bound
four-range composition now proves `y_1''(t)<=6000/t^2` for every real
`t>=1252`. Composing that theorem with the exact transfer and finite prefix
proves `Q_(11,n)(-100)>0` for every integer `n>=0`; delayed heat descent plus
the four-row lambda-zero complement then proves `Q_(11,n)(0)>0` for every
integer `n>=0`. The initial-minor transfer gives arbitrary-column signs
through order eleven. Order twelve, PF-infinity, RH, and `Lambda<=0` remain
open.

## Current Certificate Route

The active lower-bridge construction uses exact H0-H23 anchors on a step-two
lattice, H24 Taylor propagation to half-grid H0-H16 jets, and common-variable
Taylor models of degrees `(16,15,14)`. The canonical source is complete and
independently validated at `2233/2233` rows through `t=5708`. The canonical
segment cache validates all `279/279` segments and 17,792 quarter blocks
over `1252<=t<=5700`. Hash-bound pilots pass at
`t=1252`, the near/broad transition `t=1500`, the canonical far regime
`t=2200`, and the stress cell at `t=3000`.

The independent complete-prefix checker and final composer certify
`y_1''(t)<=6000/t^2` on the full finite lower bridge. A separate outward-rounded
dimensionless certificate proves the stronger `t^2*y_1''(t)<2000` on every
saddle mode `u>=20`. The finite saddle certificate independently covers all
17,999 mode blocks on `2001/1000<=u<=20`, with largest scaled upper
`482.7782145186298121`. The final piece was the compact real-t handoff

```text
5700 <= t <= V'(2001/1000)
V'(2001/1000) <= 38019.6621635193385678801171533426331770390296897595492274576
```

which is now complete and independently checked.

## Compact Resolution

The direct localized interval route remains retired because dependency
wrapping destroys the shifted cancellation. Exact H0-H23 stencil pilots now
identify source enclosure quality, especially quadrature panel count, as the
repairable obstruction:

```text
t=20004 ordinary p896/b180/n80 exact unit stencil: 17626.4157734123 (fails)
t=20004 p1536/b320/n128 exact unit stencil:          42.9019020606 (passes)
p896/b200/n96 step-four lattice, representative maxima:
  t=6000   41.5096789031
  t=10000  42.2081422734
  t=15004  43.2479096753
  t=20004  73.8573237086
  t=25004  787.0094323583
  t=30004  11582.2001559632 (fails)
  t=35004  192022.526214650 (fails)
  t=37800  1293682.15867585 (fails)
```

The upper-region profile ladder now closes the three failing centers on the
same step-four geometry. The cheapest passing tested contract is p896/b200/n108:

```text
p896/b200/n108 step-four lattice:
  t=30004  106.6387983799
  t=35004  645.4712927743
  t=37800  1949.4593338882
p896/b200/n104 at t=37800: 10601.1209490943 (fails)
```

Both true endpoint collars are also rigorous pilots. The lower collar beginning
at `t=5700` has maximum scaled upper `41.4262719267`; the upper collar through
`t=38020`, including the special H24 tile `38028..38029`, has maximum
`2122.8651669101`.

The active full-source contract is therefore a step-four adaptive lattice with
p896/b200/n96 on anchors `5692..25016` and p896/b200/n108 on
`25020..38028`. Its append-only exact H0-H23 cache currently validates
all `8085/8085` rows through `t=38028`. The final cache SHA-256 is
`63dead138ec03af796ffd37a53307ac21895c31feec960b819df7d3013e2d5bb`.
All 4,832 lower-profile rows and all 3,253 upper-profile rows are present.
Exact rebuilds match at transition rows 4831 and 4832, namely `t=25016` under
p896/b200/n96 and `t=25020` under p896/b200/n108.

The strengthened independent checker audits the outer source contract,
manifest parameters, cache hash, lower/upper row counts, next task, and live
exact rebuilds. The manifest now distinguishes the historical inner quadrature
kernel `contract_id` from the authoritative adaptive-lattice
`source_contract_id`; the inherited `order10` and `step8` tokens do not describe
the outer order-eleven step-four geometry.

The complete source is propagated across every real interval between anchors.
The canonical block driver and independent checker validate all `2020/2020`
segments and 129,280 quarter blocks over `5700..38020`. The final cache SHA-256
is `70062a5bf5c906ef894a559db1cc8ebab3379d69ffcaf7f6f88561986b7ae9a2`.
Its maximum scaled upper is `2122.8651669101`, with minimum margin
`3877.1348330899`; both occur in final segment 2019, block 63, at anchor
`152079/4`. Canonical segment 1207 still exactly matches its independent pilot.

The compact promoter and global four-source composer now pass on the live
artifacts. Their fail-closed gates still reject a 2019-segment partial cache,
a missing compact theorem, and hidden full-kernel, heat-forward, or RH claims.
The endpoint and lambda-zero checkers additionally reject higher-order leakage
and loss of the fixed-order nonpromotion countermodel.

The implemented post-source contract is recorded in
`outputs/jensen_window_pf_compound_order11_compact_block_handoff.md`. Named
pilots at the lower endpoint, profile transition, and upper endpoint pass all
192 quarter blocks. Five deliberate corruptions of ordering, profile, source
hash, propagation distance, and quarter adjacency are rejected.

## Compact Source Reproduction

```powershell
python -u work/rh_compute/scripts/jensen_window_pf_compound_order11_compact_adaptive_point_h0_h23_cache.py --workers 4 --runtime-seconds 9000
python -u work/rh_compute/scripts/check_jensen_window_pf_compound_order11_compact_adaptive_point_h0_h23_cache.py --require-complete
python -u work/rh_compute/scripts/check_jensen_window_pf_compound_order11_compact_adaptive_point_h0_h23_cache.py --skip-live-rebuilds --rebuild-index 4831 --rebuild-index 4832
```

Compact segment reproduction and resume commands are in the handoff document.

## Proof Boundary

The resolved chain establishes contiguous and arbitrary-column signed-Hankel
positivity at lambda zero through fixed order eleven. It does not establish
order twelve or any all-order limit, so it is not PF-infinity, the missing
Jensen bridge, `Lambda<=0`, RH, or a Clay-prize proof.
