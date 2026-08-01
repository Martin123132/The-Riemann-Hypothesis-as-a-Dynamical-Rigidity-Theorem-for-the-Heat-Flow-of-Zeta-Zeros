# Newman Crossing Slope-Gap Reduction

Date: 2026-07-25

Status: exact reduction, impossibility gate, and finite diagnostic note; not a proof of the Xi separation theorem, `Lambda<=0`, or RH.

## Real Component Jet

Write the corrected complex main as

```text
E=sum_j e_j,                    e_j=a_j exp(i theta_j),
e_j'=(u_j+i v_j)e_j,           a_j>0,
X=Re(E),                        U=Re(E').
```

Set

```text
p_j=cos(theta_j),
q_j=u_j cos(theta_j)-v_j sin(theta_j).
```

Then `X=sum a_j p_j` and `U=sum a_j q_j`. For `L>0`,

```text
Q_L=X^2+(U/L)^2
   =a^T (p p^T+q q^T/L^2) a                         (1)
   =sum_j a_j^2(p_j^2+q_j^2/L^2)
    +2 sum_(j<k) a_j a_k(p_j p_k+q_j q_k/L^2).
```

This is the requested exact diagonal/off-diagonal expansion. Its matrix has
rank at most two. On the crossing hyperplane `p.a=0`, it collapses to

```text
Q_L=(q.a)^2/L^2.                                    (2)
```

For at least three components, the kernel of (2) inside `p_perp` has
dimension at least `m-2`. Therefore no argument using only a positive lower
eigenvalue of this unrestricted crossing form can work.

## Slope-Gap Identity

Put `c_j=a_j p_j`, `d_j=a_j q_j`, and, when `c_j!=0`,
`h_j=d_j/c_j=u_j-v_j tan(theta_j)`. Split the indices into
`P={j:c_j>0}`, `N={j:c_j<0}`, and `Z={j:c_j=0}`, and define

```text
M_+=sum_P c_j,                 M_-=sum_N (-c_j),
M=(M_++M_-)/2,                X=M_+-M_-,
h_+=(sum_P c_j h_j)/M_+,
h_-=(sum_N (-c_j) h_j)/M_-,
D_0=sum_Z d_j.
```

Direct substitution gives the robust near-crossing identity

```text
U=M(h_+-h_-)+(X/2)(h_++h_-)+D_0.                  (3)
```

At an exact crossing, `M_+=M_-=M`, so

```text
U=M(h_+-h_-)+D_0,
Q_L=(M(h_+-h_-)+D_0)^2/L^2.                       (4)
```

Equation (4) is the surviving reduction. It uses the crossing constraint
instead of asking an indefinite full form to become positive.

A strong sufficient condition is immediate. If `D_0=0` and all positive-side
slopes exceed all negative-side slopes by at least `gamma>0`, or conversely,
then

```text
|U|>=M gamma,                 Q_L>=M^2 gamma^2/L^2. (5)
```

The actual weighted-mean gap in (4) is weaker than this pairwise ordering and
is the better theorem target.

## Half-Plane Form

Let `g_j=(p_j,q_j/L)`. If some unit vector `eta` and `kappa>0` obey

```text
eta.g_j>=kappa  for every j,
```

then positive amplitudes give

```text
sqrt(Q_L)=||sum_j a_j g_j||>=kappa sum_j a_j.      (6)
```

This is exactly the component-level version of the rotating half-plane cone
needed on a successor strip. It is sufficient, not asserted for Xi.

## Centered Arithmetic Frame

At `t=0`, the Dirichlet carriers have
`v_n=beta_0'+log(n)/2`. For any real center `mu`, put
`v_mu=beta_0'+mu/2` and
`R_(0,mu)=e_0'-i v_mu e_0` for the endpoint pseudo-component. Then

```text
E'=i v_mu E
   +(i/2) sum_(n=1)^N (log(n)-mu)e_n+R_(0,mu).    (7)
```

At `X=0`, with `E=iY`,

```text
U=-v_mu Y
  -(1/2) Im sum_(n=1)^N (log(n)-mu)e_n
  +Re R_(0,mu).                                   (8)
```

This identifies the missing input as a centered logarithmic trigonometric
moment separation with the endpoint retained, rather than as a generic
frequency-ordering claim.

## Exact Generic Obstruction

The four-carrier thought experiment

```text
E_*(x)=3 exp(i*pi/3-4ix)+3 exp(2i*pi/3-3ix)
      +7 exp(-i*pi/3-2ix)+7 exp(4i*pi/3-ix)       (9)
```

has positive amplitudes, speeds `-4<-3<-2<-1`, and no component with zero
real part at `x=0`. Nevertheless,

```text
E_*(0)=-4i sqrt(3),       E_*'(0)=-5i,
Re(E_*''(0))=-21,         M=5,
h_+=h_-=-sqrt(3)/5.                              (10)
```

Thus the slope gap can vanish exactly under all of those generic assumptions.
The new arithmetic theorem must use the actual zeta phases, amplitudes, and
endpoint relation.

## Corrected-Crossing Diagnostics

These are reproducible `t=0` route-shaping diagnostics, not interval
certificates and not evidence for the `L>=50` theorem.

| label | N | crossing | h_+-h_- | pairwise separation | common half-plane | scaled cancellation |
|---|---:|---:|---:|---|---|---:|
| ordinary_negative | 6 | 517.2201245916915343451981788255318076677 | -0.59717593741384642765405841263372255 | False | False | 12.7114793148485356900585646120 |
| ordinary_positive | 6 | 519.7486269586060360224240422895051497992 | 0.46465122614549325639932939409349042 | False | False | 17.6069535681214689479909286949 |
| lehmer_left | 33 | 14010.12594382089318324866358493301988117 | -0.030869418555740143033175118301043354 | False | False | 464.759238023827015611671363910 |
| lehmer_right | 33 | 14010.20090109869997081492153406456328151 | 0.031400855619841948649429859759450416 | False | False | 454.759131470830693542715372267 |

The coarse/fine precision comparison passed with maximum root delta
`0.0` and maximum
relative quantity delta
`0.0`.

## Xi Handoff

At an exact corrected-main crossing, the live target

```text
Q_L>8000000 exp(-3L/2)
```

is equivalent to

```text
|M(h_+-h_-)+D_0|
  >2000 sqrt(2) L exp(-3L/4),       L>=50.         (11)
```

This is a sharper statement of the missing theorem, not its proof. A viable
proof may bound the centered log moment in (8), prove the weighted slope gap
(11), or establish an equivalent Xi-specific half-plane cone on newly
entering high-frequency strips.

There is an important quantifier guard. If (11) is imposed uniformly at one
fixed `x` for every positive `t` tending to zero, continuity promotes it to
a nonzero endpoint first jet and therefore adds a high-zero simplicity
burden. The existing positive-boundary delta-localization gate already
proves that such a uniform floor is sufficient but not logically necessary
for `Lambda<=0`. The multiplicity-compatible alternative is to use (11), or
a weaker entry certificate, only where new cells enter and then transport
their descendants with a time-adaptive relative estimate. The finite
diagnostics, exact countermodel, Q208 theorem, and one finite successor do
not prove either route, an all-`j` successor, `Lambda<=0`, RH, PF-infinity,
or the Clay prize.

Machine-audited files:

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_crossing_slope_gap_reduction.json
work/rh_compute/scripts/jensen_window_pf_newman_polymath15_crossing_slope_gap_reduction.py
work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_crossing_slope_gap_reduction.py
```
