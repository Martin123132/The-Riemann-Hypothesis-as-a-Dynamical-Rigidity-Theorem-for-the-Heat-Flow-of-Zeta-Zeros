# Ray-Aligned Energy Endpoint-Margin Gate

Date: 2026-07-28

Status: exact endpoint margins and whole-collar conditioning
audit. This is not a proof of the descendant theorem,
Lambda<=0, or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_ray_aligned_energy_endpoint_margin_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_ray_aligned_energy_endpoint_margin_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_ray_aligned_energy_endpoint_margin_gate.py
```

Current result:

```text
validated ray-aligned energy endpoint-margin gate: 18 rows, 10 compact edge records, 1 left raw margin, 1 right normalized margin, 1 exact contact floor, 1 open full-collar bulk budget, 2 nonpromotion guards, 0 pointwise contact exclusions
```

## Active Stages

The quantitative compact endpoint source applies to

```text
j>=25
L_j=101+j>=126
t_(j+1)=25/L_j and t_j=25/(L_j-1)
W_j=R_j-38=4*pi*exp(L_j)-38
```

t_j=25/(100+j)<=1/5 exactly when j>=25; earlier finite stages remain separate base/transition bookkeeping and are not used in this cofinal endpoint estimate

Pi provenance remains unchanged: `pi` is the ordinary circle
constant inherited from `L=log(x/(4*pi))`.

## Left Endpoint

At `x=38`, all ten stored compact edge boxes use the
derivative branch and tile `0<=t<=1/5`. They prove

```text
|J_(1,t)'|/B_1>7/4
|J_t'-J_(1,t)'|<B_1
B_1(t,38)>=4*38^3/2800
|J_t'(38)|>3*38^3/2800=20577/350
```

The last line is a margin for the full exact heat-flow
characteristic derivative, not only its first theta block.

Because

```text
J_t'(x)=64x^3H_t(x)+16x^4H_t'(x)
|J_t'(38)|<=sqrt((64*38^3)^2+(16*38^4/s_pf(t,38))^2)*sqrt(E_pf(t,38))
```

we obtain

```text
E_pf(t,38)>M_38^2/[4096*38^6+256*38^8/s_pf(t,38)^2], M_38=20577/350
E_pf(t,38)>=M_38^2/[4096*38^6+256*38^8*(log(38/(4*pi))^2+L_j/50)]
```

No unknown jet norm was used as a divisor.

## Right Endpoint

At `x=R_j`, the complete dominant-ray theorem gives

```text
tL_j>=25 and L_j>=126
C=(H_x^2-H H_xx)/A_t^2>3L_j^2/40
Z=H/A_t, g=H_x/A_t, h_2=H_xx/A_t, C=g^2-Zh_2
```

The main-sum, remainder, and normalizer derivative bounds
compose to

```text
|H_xx|/(A_t L_j^2)<13/20
s_pf^2 L_j^2=q/(1+q)>=6300/6301, q>=6300
C/L_j^2<=y^2/(s_pf^2L_j^2)+(13/20)y
at y=99/1000 the upper side is 51906789/700000000=3/40-593211/700000000
sqrt(E_pf(t,R_j))>(99/1000)A_t(R_j)
```

The factor `A_t(R_j)` is retained. It is not normalized away.

## Exact Contact Floor

The prior energy gate now gives the fully sourced implication

```text
integral_38^R_j B_pf dx>=[m_(L,j)+(99/1000)A_t(R_j)]^2/W_j
```

A sufficient strict theorem would be

```text
D_j(t)<=U_j(t), with U_j derived independently from Xi source bounds and D_j=d/dt integral E_pf+2[J_pf(R_j)-J_pf(38)]
U_j(t)<2[m_(L,j)+(99/1000)A_t(R_j)]^2/W_j
```

That upper budget is not currently proved.

## Conditioning

The left-only certified floor has the exact asymptotic

```text
m_(L,j)^2/W_j~[9/(231853260800*pi)]*exp(-L_j)/L_j
```

The currently certified raw floor loses one exponential radius factor. This is a conditioning statement about the available sufficient bound, not a claim that the actual endpoint energy has this exact size.

| stage j | L_j | bottom t | certified left norm | log10(left-only floor) |
|---:|---:|---:|---:|---:|
| 25 | 126 | 0.1984126984126984 | 9.093316658558458e-07 | -67.902869954066773 |
| 100 | 201 | 0.12437810945273632 | 7.6868810626580444e-07 | -100.62090040381037 |
| 1000 | 1101 | 0.022706630336058128 | 3.6542330980843115e-07 | -492.13184194570636 |
| 10000 | 10101 | 0.0024750024750024749 | 1.2360567975861653e-07 | -4401.7236947141373 |

The table is a deterministic conditioning diagnostic, not a
finite-stage proof or counterexample.

## Interior Numerator Gap

For fixed `t`, the dominant theorem begins at

```text
x_dom(t)=4*pi*exp(25/t)
x_dom(t_(j+1))=R_j
x_dom(t_j)=R_(j-1)
[x_dom(t),R_j]
[38,x_dom(t))
```

The compact certificate stops at x=38 and the dominant-ray certificate starts at x_dom(t). Neither supplies a full-collar upper bound U_j for D_j. Pointwise noncontact in other bands does not itself bound the integrated H_x and H_xx bulk.

## Route Decision

Substituting D_j=2 integral B_pf is an identity, not the strict source-level upper estimate required by the test.

The right margin is relative to A_t(R_j); A_t cannot be dropped or replaced by one without a proved comparison.

A finite evaluation of rho_j cannot prove the all-stage bulk budget.

Retain both endpoint margins, but do not launch a whole-collar bulk computation until a source-level interior majorant exists. The preferred refinement is a localized energy criterion on bounded logarithmic or phase cells with independently certified endpoint margins; otherwise return to the contact-normal C1 arithmetic target.

Open handoff:

Derive an exact localization of the energy/current identity on bounded cells and inventory which cell endpoints already have quantitative first-jet margins from compact, oscillatory-zeta, dominant-saddle, or phase-cell sources. Determine whether the cell boundary fluxes telescope without requiring the same open interior theorem.

## Boundary

This artifact proves the compact left derivative margin, its parabolic-frequency energy conversion, the dominant-ray normalizer-relative right energy margin, the instantiated two-endpoint contact floor, its cofinal conditioning, and the moving dominant-tail coverage geometry. It does not prove a full-collar Xi bulk upper bound, the strict no-contact budget, a localized telescoping theorem, Q209, the cofinal descendant theorem, Lambda<=0, PF-infinity, RH, or a Clay-prize conclusion.
