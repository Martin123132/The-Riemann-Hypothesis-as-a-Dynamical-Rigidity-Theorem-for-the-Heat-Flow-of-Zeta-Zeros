# Exact fixed-selector height transport of the lower fold

Date: 2026-08-11

Status: exact height-transport identity and complete 84-mode saved-height
derivative validated; not a proof of a nonzero height-cell enclosure

While the odd source endpoint selector `C` is fixed,

```text
beta^3=pi*C^2/8,
lambda(t;C)=[pi*C^2/8-t]/beta,
d lambda/dt=-1/beta.                                  (HT1)
```

For the finite transform `G_Y(d)` of Section 11.324,

```text
partial_lambda G_Y
 =Ai(-lambda-Y)w_Y-Ai(-lambda)
  +G_Y'/(2beta)+i*d*G_Y,                              (HT2)

dG_Y/dt=-(partial_lambda G_Y)/beta.                    (HT3)
```

Equation (HT2) follows from
`partial_lambda Ai(-lambda-y)=partial_y Ai(-lambda-y)` and one exact
integration by parts.  Both finite endpoint values are retained.

The current saved height lies in the fixed-selector cell

```text
[9999760347.647943710077056053348795787925647462426307094224362837236653 +/- 4.21e-61]<t<
[10000011009.04258833249996651073906089290177178617834443516760678611578 +/- 4.07e-60],                       (HT4)
```

at distance

```text
[11009.04258833249996651073906089290177178617834443516760678611578384866 +/- 2.14e-61]            (HT5)
```

below its upper selector fold.  All 84 values of `G_64'`,
`partial_lambda G_64`, and `dG_64/dt` were evaluated as complex balls.  The
grouped canonical height derivative is

```text
d/dt [2pi sum_m G_64(d_m)]
 =[0.1305127144667297220085814791523352302470003090059549893087442709631377 +/- 1.33e-25]
  +i*[-0.0001213396458241485134591364143593757347215943672136758824701598851589646 +/- 1.33e-25], (HT6)
```

with magnitude

```text
[0.1305127708723700739560575864728139845147250996108012813981227205416919 +/- 1.33e-25].        (HT7)
```

The independently assembled grouped chain-rule residual contains zero in
both components.  This gives an exact local transport law, but a point
derivative is not yet a nonzero-radius height theorem.  The next obligation
is a second-derivative or interval-ODE bound on a selector-stable height cell,
followed by subdivision up to the selector boundary.

Proof boundary: exact fixed-selector height identity and complete derivative
at `t=10^10` only.  No nonzero-radius height enclosure, selector-transition
join, complete `T_upper`, `Lambda<=0`, RH, or prize-level conclusion is
proved.
