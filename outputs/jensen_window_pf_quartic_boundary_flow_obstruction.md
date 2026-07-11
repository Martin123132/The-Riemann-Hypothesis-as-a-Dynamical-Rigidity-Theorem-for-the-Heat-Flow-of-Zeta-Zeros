# Jensen-Window PF Quartic Boundary-Flow Obstruction

Date: 2026-07-10

Status: exact quartic boundary-flow countermodel. This blocks promotion
of the cubic cone but is not a counterexample to the zeta trajectory and
is not a proof or disproof of RH.

```text
work/rh_compute/results/jensen_window_pf_quartic_boundary_flow_obstruction.json
python work/rh_compute/scripts/jensen_window_pf_quartic_boundary_flow_obstruction.py
python work/rh_compute/scripts/check_jensen_window_pf_quartic_boundary_flow_obstruction.py
```

Current result:

```text
validated Jensen-window PF quartic boundary-flow obstruction: 10 rows, 0 issues, 4 exact quartic identities, 1 hyperbolic boundary point, 4 positive ratio margins, 3 strict cubic margins, 1 negative quartic derivative, 1 blocked promotion, 1 open coupled-invariant handoff
```

## Quartic Frontier

After positive rescaling, a shifted quartic window is

```text
J_4=1+4*w+6*x*w^2+4*x^2*y*w^3+x^3*y^2*z*w^4.
Disc(J_4)=256*x^6*y^2*Q(x,y,z)
```

The derivative cubic satisfies

```text
Disc(partial_w J_4)=-6912*x^6*y^2*F(y,z)
F(s,t)=s^2*t^2-6*s*t+4*s+4*t-3
```

The exact heat vector for `(x,y,z)` makes `Q'` linear in the next
contraction `u`.

## Exact Boundary Countermodel

Take `k=1` and

```text
(1+13*w/20)^2*(1+21*w/50)*(1+57*w/25)
x=48901/60000
y=2147067000/2391307801
z=u=104061328/112253067
```

The displayed factorization gives four negative roots with one double
root, so `Q=0` is a genuine hyperbolic quartic boundary point. The
contractions are nondecreasing and satisfy every k=1 pointwise wall.
All neighboring cubic tests are strict:

```text
F(x,y)=-260340783523/71739234030000<0
F(y,z)=-457533698448089/268431634803275667<0
F(z,u)=-242320485975382436128527365651/158778927046920986937948622307121<0
```

Nevertheless, direct substitution into the heat vector gives

```text
Q=0, Q'/r_1=-13108711376416987159336748097/20606742971316325673502124987495<0.
```

Because `Disc(J_4)=256*x^6*y^2*Q`, the discriminant immediately becomes
negative. Thus the ratio cone plus all cubic inequalities is not a
quartic invariant.

## Consequence

Degree four needs a new coupled condition involving the quartic frontier
and the next contraction. This abstract guard says nothing adverse about
the actual zeta coefficient trajectory; it only rejects an insufficient
proof shortcut. PF-infinity, RH, and `Lambda <= 0` remain open.
