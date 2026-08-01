# Newman Theta Fourth Diagonal-Shell Two-Block Interval Certificate

Date: 2026-07-24

Status: rigorous fourth-stage interval theorem, not a proof of RH.
Stages `j>=5`, `Lambda <= 0`, and RH remain open.

```text
work/rh_compute/results/jensen_window_pf_newman_theta_fourth_diagonal_shell_two_block_interval_certificate.json
python work/rh_compute/scripts/jensen_window_pf_newman_theta_fourth_diagonal_shell_two_block_interval_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_theta_fourth_diagonal_shell_two_block_interval_certificate.py
```

## Inherited Two-Block Theorem

```text
M_k<=2*pi^2*k!*sum_(n>=3)n^4*exp(-pi*n^2)/(4*pi*n^2-9)^(k+1), k=0,1; Using pi>3 and pi^2<10, consecutive majorants have ratio <4*exp(-21)<1/2; hence M_0<40*81*exp(-27)/99<10^-10 and M_1<40*81*exp(-27)/99^2<10^-12; |J_t-J_(<=2,t)|<16*x^4*10^-10 and |J_t'-J_(<=2,t)'|<64*x^3*10^-10+16*x^4*10^-12
```

## Certified Slab

```text
S_4=[1/20,1/5]x[38,42]
(J_t(x),J_t'(x))!=(0,0), equivalently (H_t(x),H_t'(x))!=(0,0), for every (t,x) in [1/20,1/5]x[38,42]
initial boxes=120
certified boxes=120
subdivisions=0
unresolved boxes=0
minimum normalized disjunction ratio=[88928.172125791174079409729388765596692868343484099 +/- 1.92e-46]
```

## Sign Transition

```text
105 boxes prove J_t<0; the 15 boxes covering [1/20,1/5]x[83/2,42] use strict J_t' separation
```

| t interval | x interval | value ratio | derivative ratio |
|:---|:---|:---|:---|
| [1/20,3/50] | [83/2,42] | 0 | [92225.672908130126234020600478693702834135171549686 +/- 4.82e-46] |
| [3/50,7/100] | [83/2,42] | 0 | [91991.110134614203391416746074524658196714453845931 +/- 4.72e-46] |
| [7/100,2/25] | [83/2,42] | 0 | [91756.395624056041290396218356835589421109641659060 +/- 2.15e-46] |
| [2/25,9/100] | [83/2,42] | 0 | [91521.529974211179985564756552583661522375751000444 +/- 3.25e-46] |
| [9/100,1/10] | [83/2,42] | 0 | [91286.513356679302366188625725682044696660118288768 +/- 4.71e-46] |
| [1/10,11/100] | [83/2,42] | 0 | [91051.346569184919949675725672634448833693231223005 +/- 4.74e-46] |
| [11/100,3/25] | [83/2,42] | 0 | [90816.030153936288495598641968781994936043487740997 +/- 2.20e-46] |
| [3/25,13/100] | [83/2,42] | 0 | [90580.564539825259229600245083566954536492423022038 +/- 1.15e-46] |
| [13/100,7/50] | [83/2,42] | 0 | [90344.950241507181090394504954465128052521617625957 +/- 1.66e-46] |
| [7/50,3/20] | [83/2,42] | 0 | [90109.187717200489757440435527287510978362702013603 +/- 1.48e-46] |
| [3/20,4/25] | [83/2,42] | 0 | [89873.277738406763974471802658512082649757942646987 +/- 4.55e-46] |
| [4/25,17/100] | [83/2,42] | 0 | [89637.220678909792203289479603909152933479427350000 +/- 1.76e-46] |
| [17/100,9/50] | [83/2,42] | 0 | [89401.017026695757001131005150767938868579867372504 +/- 3.45e-46] |
| [9/50,19/100] | [83/2,42] | 0 | [89164.667127992715914965917209128407085852078507883 +/- 2.02e-46] |
| [19/100,1/5] | [83/2,42] | 0 | [88928.172125791174079409729388765596692868343484099 +/- 1.92e-46] |

## Fourth Stage

```text
Q_4=[1/20,1/4]x[0,42] contains no common zero of H_t and H_t'
[1/20,1/4]x[-42,42] contains no common zero of H_t and H_t'
Z=H+iH_x is nonzero on partial Q_4 and wind(Z(partial Q_4),0)=0
```

## Proof Boundary

```text
The next stage Q_5=[1/25,1/4]x[0,43] requires a new certificate on [1/25,1/5]x[38,43]; it is not proved here
```

validated Newman theta fourth diagonal-shell two-block interval certificate: 8 rows, 0 issues, 120 certified slab boxes, 0 subdivisions, 0 unresolved, 105 negative-value boxes, 15 derivative-only boxes, minimum disjunction ratio >88000, 1 fourth-stage no-contact theorem, 1 zero-winding composition, 1 open fifth-stage handoff
