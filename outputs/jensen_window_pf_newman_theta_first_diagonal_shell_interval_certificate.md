# Newman Theta First Diagonal-Shell Interval Certificate

Date: 2026-07-24

Status: rigorous first-stage interval theorem, not a proof of RH.
Stages `j>=2`, `Lambda <= 0`, and RH remain open.

```text
work/rh_compute/results/jensen_window_pf_newman_theta_first_diagonal_shell_interval_certificate.json
python work/rh_compute/scripts/jensen_window_pf_newman_theta_first_diagonal_shell_interval_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_theta_first_diagonal_shell_interval_certificate.py
```

## Certified Edge

```text
At t=1/5, (H_t(x),H_t'(x))!=(0,0) for 38<=x<=39
Evenness gives the same contact exclusion for -39<=x<=-38
The prior compact theorem covers |x|<=38 at t=1/5, so (H_(1/5),H_(1/5)')!=(0,0) for every |x|<=39
```

The replay uses the same 160-bit Arb retained integrals, analytic
`10^-800` tail, and bivariate Taylor machinery as the compact
certificate, specialized to the exact time `t=1/5`.

```text
certified boxes=5
subdivisions=0
unresolved boxes=0
minimum normalized disjunction ratio=[1.2104189126201015266000831045663248320677132494 +/- 4.12e-47]
```

| t | x interval | branch | certified ratio lower |
|:---|:---|:---|:---|
| 1/5 | [38,191/5] | derivative | [1.6186252035268730561634382935386137715289823125 +/- 3.88e-47] |
| 1/5 | [191/5,192/5] | derivative | [1.5156167058277156552974231782844860698851352612 +/- 4.00e-48] |
| 1/5 | [192/5,193/5] | derivative | [1.4128873681814139883981530629800038348565640621 +/- 1.26e-47] |
| 1/5 | [193/5,194/5] | derivative | [1.3109896095964587064351245549108035042053820658 +/- 9.02e-49] |
| 1/5 | [194/5,39] | derivative | [1.2104189126201015266000831045663248320677132494 +/- 4.12e-47] |

## First Stage

```text
The published bound Lambda<=1/5 implies that H_t has only simple zeros for every t>1/5
Q_1=[1/5,1/4]x[0,39] contains no common zero of H_t and H_t'
[1/5,1/4]x[-39,39] contains no common zero of H_t and H_t'
Z=H+iH_x is nonzero on partial Q_1 and wind(Z(partial Q_1),0)=0
```

Thus the first member of the independent diagonal exhaustion is
now a theorem, not a diagnostic.

## Proof Boundary

```text
The next diagonal stage Q_2=[1/10,1/4]x[0,40] requires bottom-edge separation at t=1/10, right-edge separation at x=40 for 1/10<=t<=1/5, and zero first-jet winding; none is proved by this first-stage certificate
```

validated Newman theta first diagonal-shell interval certificate: 8 rows, 0 issues, 5 certified edge boxes, 0 subdivisions, 0 unresolved, minimum ratio >6/5, 1 first-stage no-contact theorem, 1 zero-winding composition, 1 open second-stage handoff
