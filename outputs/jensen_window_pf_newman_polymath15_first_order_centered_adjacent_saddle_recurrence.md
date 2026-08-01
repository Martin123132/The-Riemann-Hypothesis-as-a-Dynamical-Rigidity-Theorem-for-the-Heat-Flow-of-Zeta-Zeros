# Jensen-Window PF Newman Polymath-15 First-Order Centered Adjacent-Saddle Recurrence

Date: 2026-07-26

Status: exact adjacent-saddle identities and selected high-height
diagnostics. The uniform differentiated-ratio theorem remains open;
this is not a proof of `Lambda<=0` or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_adjacent_saddle_recurrence.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_adjacent_saddle_recurrence.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_adjacent_saddle_recurrence.py
```

## Exact Adjacent Recurrence

Put `n=N+1`, `p=1-2(a-N)`,
`r=(n-a)/a=(p+1)/(2a)`, and
`R(p)=exp(pi*i*(p^2/2+p+3/8))`. Differentiating the exact
`C_0` recurrence gives

```text
R'''=(-3*pi^2*(p+1)-i*pi^3*(p+1)^3)*R; R''''=(pi^4*(p+1)^4-6i*pi^3*(p+1)^2-3*pi^2)*R

J_a=H_a(p+2)+H_a(p)=R(p)*[1-r/2-(2*pi*i/3)*a^2*r^3], r=(N+1-a)/a

J_(a,x)/R=(1+r)/(16*pi*a^2)-pi*a^2*r^4/3+i*(-r/2+r^2/2+r^3/12)
```

If `H_a=F+F'''/(12pi^2a)` and `g_N=-kappa_NH_a(p)`,
then `kappa_(N+1)=-kappa_N` and

```text
g_(N+1)-g_N=kappa_N*J_a.
```

No subtraction of two long finite sums occurs.

## Centered Scalar Jump

The moment changes by
`M_(1,a;N+1)-M_(1,a;N)=log((N+1)/a)f_(N+1)`.
Consequently

```text
Delta A_a=Re[-s_*'*ell*f_(N+1)+kappa_N*(J_(a,x)+mu_a*J_a)]

Q_N=f_(N+1)+kappa_N*J_a=E_[1],N+1-E_[1],N; Delta A_a=Re(Q_(N,x)-lambda_a*Q_N-e_(N+1)*d_(N+1,x))

Delta A_a=Delta U-u_a*Delta X+v_a*Delta Y-Re(e_(N+1)*d_(N+1,x))
```

The complex half itself is chart dependent: its adjacent mismatch
can have a large imaginary part. The exact real trace and its real
derivative are the relevant quantities; replacing them by a complex
absolute-value estimate would throw away the cancellation.

## Leading Cancellation

For the reflected positive sharp block, set `epsilon=1/a` and
`w=N+1-a`. The exact adjacent log-ratio has first coefficient

```text
-(1/2-i*pi*t/8)w-i*pi*t*w/8-i*(2*pi/3)w^3
 =-w/2-i*(2*pi/3)w^3,
```

which is exactly the first coefficient of
`log(1-w/(2a)-i*(2pi/3)w^3/a)`. The heat shift cancels too.
The finite correction starts only at `a^-2`. Therefore

```text
For epsilon=1/a and w=N+1-a, [epsilon]log(main_sharp/endpoint_C1)=0
main_sharp/endpoint_C1=1+O(a^-2); after centered x-differentiation Delta A_a=O(a^-7/2)=O(exp(-7L/4))
```

This explains why the separate finite and endpoint pieces live at
`exp(-3L/4)` while their centered real combination is expected at
`exp(-7L/4)`. A uniform remainder constant is still required.

## Selected High-Height Audit

Independent 110- and 155-digit runs agree with maximum relative drift `0.0`.

| N | theta | c=tL | L | finite e3 | endpoint e3 | Delta A e7 |
|---:|---:|---:|---:|---:|---:|---:|
| 100000000000 | 0 | 0 | 50.65687205 | -0.191341716 | 0.191341716 | 1.781994112 |
| 100000000000 | 0 | 1 | 50.65687205 | -0.191341716 | 0.191341716 | 1.780505559 |
| 100000000000 | 0 | 25 | 50.65687205 | -0.191341716 | 0.191341716 | 1.744780292 |
| 100000000000 | 0.5 | 0 | 50.65687205 | 0.230969883 | -0.230969883 | 0.06315348448 |
| 100000000000 | 0.5 | 1 | 50.65687205 | 0.230969883 | -0.230969883 | 0.06332573523 |
| 100000000000 | 0.5 | 25 | 50.65687205 | 0.230969883 | -0.230969883 | 0.06745975311 |
| 100000000000 | 0.9 | 0 | 50.65687205 | -0.0161958709 | 0.0161958709 | -0.01037345218 |
| 100000000000 | 0.9 | 1 | 50.65687205 | -0.0161958709 | 0.0161958709 | -0.01042998459 |
| 100000000000 | 0.9 | 25 | 50.65687205 | -0.0161958709 | 0.0161958709 | -0.0117867625 |
| 10000000000000 | 0 | 0 | 59.86721242 | -0.191341716 | 0.191341716 | 1.781994112 |
| 10000000000000 | 0 | 1 | 59.86721242 | -0.191341716 | 0.191341716 | 1.780734567 |
| 10000000000000 | 0 | 25 | 59.86721242 | -0.191341716 | 0.191341716 | 1.750505495 |
| 10000000000000 | 0.5 | 0 | 59.86721242 | 0.230969883 | -0.230969883 | 0.06315348448 |
| 10000000000000 | 0.5 | 1 | 59.86721242 | 0.230969883 | -0.230969883 | 0.06329923511 |
| 10000000000000 | 0.5 | 25 | 59.86721242 | 0.230969883 | -0.230969883 | 0.06679725024 |
| 10000000000000 | 0.9 | 0 | 59.86721242 | -0.0161958709 | 0.0161958709 | -0.01037345218 |
| 10000000000000 | 0.9 | 1 | 59.86721242 | -0.0161958709 | 0.0161958709 | -0.0104212873 |
| 10000000000000 | 0.9 | 25 | 59.86721242 | -0.0161958709 | 0.0161958709 | -0.01156933014 |

Here `finite e3` and `endpoint e3` multiply the real pieces by
`exp(3L/4+cL/16)`, while `Delta A e7` multiplies their sum by
`exp(7L/4+cL/16)`. Across these selected rows,

```text
max edge-scaled magnitude = 0.230969883127807635880681556221 < 0.24
max scalar-jump magnitude = 1.78199411194874888375743432551 < 2
```

These are diagnostics at eighteen prescribed points, not an
interval proof over `theta`, `L`, or `t`.

## Route Decision

The edge contribution cannot be assigned a sign by bounding the
last centered finite term and endpoint separately. Their
`exp(-3L/4)` pieces are equal-and-opposite in the proof-facing
real coordinate. The adjacent recurrence must be differentiated
before absolute values are taken.

The next technical theorem is explicit:

```text
Prove |rho-1|<=C0/a^2 and |rho_x|<=C1/a^3
uniformly for L>=50, 0<tL<=25, 0<=N+1-a<=1.
Deduce |A_(N+1)-A_N|<=C_A*exp(-7L/4)
with any explicit constant small at the exp(-5L/4) contact scale.
```

Once that chart-stability lemma is certified, the remaining
RH-level object is the chart-invariant bulk scalar lower bound and
its oriented positive-crossing count. This artifact supplies
neither theorem and does not prove contact exclusion, `Lambda<=0`,
RH, PF-infinity, or a Clay-prize conclusion.
