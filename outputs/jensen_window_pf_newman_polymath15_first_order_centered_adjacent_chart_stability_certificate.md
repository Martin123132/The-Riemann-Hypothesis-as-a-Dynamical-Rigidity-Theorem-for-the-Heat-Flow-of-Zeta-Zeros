# Jensen-Window PF Newman Polymath-15 First-Order Centered Adjacent-Chart Stability Certificate

Date: 2026-07-26

Status: uniform adjacent-chart theorem on `L>=50`,
`0<=tL<=25`. This closes one endpoint bookkeeping problem;
it is not the bulk scalar lower bound and not a proof of RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_adjacent_chart_stability_certificate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_adjacent_chart_stability_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_adjacent_chart_stability_certificate.py
```

## Scaled Domain

Put

```text
L>=50, 0<=tL<=25, epsilon=1/a, 0<=w=N+1-a<=1, epsilon<=exp(-25)
h=1/2-i*pi*t/8, n=N+1.
```

The 256-bit interval endpoint is `[1.388794386496402059466176374608685691039976038020505558354779996088137e-11 +/- 1.87e-81]`;
in particular `epsilon<1.4e-11` and all logarithm disks below
have radius less than `4e-11`.

## Exact Corrected Ratio

Reflect the entering negative sharp to the positive sharp, which
does not change its real trace. Dividing by the adjacent
`C_0+C_1/a` recurrence block gives

```text
rho=exp(Omega)*(1+d_+)/B, B=1-epsilon*w/2-i*(2*pi/3)*epsilon*w^3
Omega=h*r0+R_M-delta/2+i*pi*t*r1/8+t*(r1-delta)^2/4-i*pi*Psi
```

The heat/logarithm combination simplifies exactly:

```text
-h*delta+i*pi*t*(r1-delta)/8
 =-delta/2+i*pi*t*r1/8.
```

There is therefore no uncancelled `t*delta` term. The scaled
blocks are

```text
delta=log(1+epsilon*w); Psi=epsilon^-2*(2delta-2epsilon*w+epsilon^2*w^2); r0=epsilon^2/(4*pi*i)+epsilon^2/(2*pi*i-epsilon^2); r1=epsilon^2/[2*(2*pi*i+h*epsilon^2)]+epsilon^2/[2*pi*i+(h-1)*epsilon^2]+log(1+h*epsilon^2/(2*pi*i))/2
```

These forms are analytic at `epsilon=0`; bounding an apparent
`1/epsilon` derivative term by term before this rewrite would
miss the cancellation.

## Value Bound

Write `B=1+epsilon*B1`,
`B1=-w/2-i*(2pi/3)w^3`. The exact first coefficients of
`-delta/2-i*pi*Psi` and `log B` agree. Taylor's bound

```text
|log(1+z)-z+z^2/2|<=|z|^3/[3(1-|z|)]
```

gives the geometric coefficient budget

```text
1/4+pi/2+|B1|^2/2+three explicit remainders
 = [5.186239300556402808485452271825124637990800286934980516319957486169029 +/- 3.89e-70] < 5.2.
```

The remaining scaled terms satisfy

```text
|h*r0|/epsilon^2 < |h|*3/(4pi),
|R_M|/epsilon^2 < |h|^2/(2pi),
|i*pi*t*r1/8|/epsilon^2 < pi/16,
|t*(r1-delta)^2/4|/epsilon^2 < 0.13,
|log(1+d_+)|/epsilon^2 < 0.12.
```

Their interval sum is

```text
[5.806753732218342642502703806471324781277737981561575927296266048261212 +/- 7.31e-71] < 8,
|Z|=|log rho|<8*epsilon^2.
```

## Derivative Bound

Differentiation is performed after the first coefficient has
cancelled. For `w` the coefficient derivative costs

```text
1/2+2pi+|B1||B1'|+0.4
 = [24.78144804675178895332146594409987779440883054090474023846761447598611 +/- 4.46e-69] < 30,
|Z_w|<30*epsilon^2.
```

For `epsilon`, direct differentiation of the displayed scaled
rational forms gives the termwise budget

```text
geometric <14, r0 <1, R_M <1, heat-r1 <1,
square <1, log(1+d_+) <1;
interval sum = [14.80000000000000000000000000000000000000000000000000000000000000000000 +/- 3e-73] < 15.
|Z_epsilon|<100*epsilon.
```

The saved constant leaves a factor-four margin. Since

```text
w_x=-epsilon/(8*pi), epsilon_x=-epsilon^3/(8*pi)
```

one obtains

```text
|Z_x|<[1.193662073244473391153922866778970142587972861528969144535929510183526 +/- 1.28e-70]*epsilon^3
<2*epsilon^3,
|rho-1|<11*epsilon^2, |rho_x|<3*epsilon^3.
```

## Centered Endpoint Rate

The endpoint prefactor has a second exact cancellation. With
`chi=alpha-log(a)`,

```text
mu_a=K_x/K+i*(alpha-log(a))/2+i*t*alpha'*(alpha-log(a))/4=O(epsilon^2)
```

and the `-pi/8` in `K_x/K` cancels the `pi/8` from
`i*chi/2`. The scaled alpha remainder gives
`|mu_a|<epsilon^2`. The explicit recurrence polynomial gives
`|J_(a,x)/J_a|<[0.5000000000397790913777049427839907511384308974120272886353038688349288 +/- 2.73e-71]epsilon<0.51epsilon`; together with
`|lambda_a|<3epsilon^2`,

```text
|lambda_a+mu_a+J_(a,x)/J_a|<2*epsilon
```

## Chart-Stability Theorem

The inherited recurrence endpoint bound, enlarged from 50 to
100 to include `B`, and the ratio estimates imply

```text
|Delta X|<1100*exp(-5L/4), |Delta U|<2500*exp(-7L/4)
```

The exact projection identity is

```text
Delta A_a=Delta U-u_a*Delta X+v_a*Delta Y
             -Re(e_(N+1)d_(N+1,x)).
```

The entering derivative correction costs at most
`1500exp(-7L/4)`. The imported `u_a`, `v_a`, and complex-main
bounds put both frame terms below `5exp(-7L/4)`. Hence

```text
raw constant = 4005.000000000000000000000000000000000000000000000000000000000000000000
|A_(a,N+1)-A_(a,N)|<5000*exp(-7L/4)
```

At `L=50` this is smaller than the `exp(-5L/4)` scale by
`[6.943971932482010297330881873043428455199880190102527791773899980440684e-8 +/- 6.94e-79] < 1e-7`, and the
ratio decreases thereafter.

This proves adjacent chart stability only. It removes the
last-saddle/endpoint bookkeeping obstruction but supplies no
lower bound or sign for the chart-invariant bulk centered
moment, no positive-crossing count, no q<1 theorem, and no
contact exclusion, `Lambda<=0`, RH, PF-infinity, or
Clay-prize conclusion.
