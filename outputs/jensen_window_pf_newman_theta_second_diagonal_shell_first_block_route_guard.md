# Newman Theta Second Diagonal-Shell First-Block Route Guard

Date: 2026-07-24

Status: rigorous pointwise countercertificate for the unchanged
first-block sufficient bars, not a proof of RH. `Q_2`,
`Lambda <= 0`, and RH remain open.

```text
work/rh_compute/results/jensen_window_pf_newman_theta_second_diagonal_shell_first_block_route_guard.json
python work/rh_compute/scripts/jensen_window_pf_newman_theta_second_diagonal_shell_first_block_route_guard.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_theta_second_diagonal_shell_first_block_route_guard.py
```

## Exact Point Audit

The ratios use lower bounds to certify a branch and upper bounds
to certify failure of both branches.

| t | x | value ratio interval | derivative ratio interval | result |
|:---|:---|:---|:---|:---|
| 1/10 | 198/5 | [[0.84288431151307944016596013511775480697036568114 +/- 3.32e-48], [0.84288431151307944016596013511775480925091245629 +/- 3.83e-48]] | [[1.0046112346529529272333085056054455079295048144 +/- 2.54e-47], [1.0046112346529529272333085056054455303568951491 +/- 1.26e-47]] | raw_disjunction_certified |
| 1/10 | 397/10 | [[0.82261911920733491888296363287239076485236138220 +/- 3.14e-49], [0.82261911920733491888296363287239076718284969070 +/- 4.98e-48]] | [[0.95590997016615601654304873817611947824794321844 +/- 2.83e-48], [0.95590997016615601654304873817611950120956087898 +/- 3.25e-48]] | raw_disjunction_rigorously_false |
| 1/10 | 199/5 | [[0.80318335372701128293114742598269425866604741786 +/- 2.82e-48], [0.80318335372701128293114742598269426104757076836 +/- 4.36e-49]] | [[0.90798453357514738914518254006276905161040868858 +/- 3.69e-48], [0.90798453357514738914518254006276907511885791843 +/- 1.76e-48]] | raw_disjunction_rigorously_false |
| 1/10 | 399/10 | [[0.78455559814390148074027213381356984907421318167 +/- 7.24e-49], [0.78455559814390148074027213381356985150788909835 +/- 4.34e-48]] | [[0.86086595463177080360978346970996374923832204287 +/- 2.83e-48], [0.86086595463177080360978346970996377330650271754 +/- 2.01e-48]] | raw_disjunction_rigorously_false |
| 1/10 | 40 | [[0.76671444081779264215769884881565070452061105590 +/- 4.88e-48], [0.76671444081779264215769884881565070700758143004 +/- 1.93e-48]] | [[0.81458256996893635990149200005932916081446713324 +/- 4.21e-48], [0.81458256996893635990149200005932918545558035719 +/- 8.45e-49]] | raw_disjunction_rigorously_false |

```text
At t=1/10 and each x in {397/10,199/5,399/10,40}, rigorous Arb upper bounds give |J_1|/B_0<1 and |J_1'|/B_1<1
The unchanged first-block sufficient disjunction cannot certify the complete bottom edge of Q_2
maximum failed upper ratio=[0.95590997016615601654304873817611950120956087898 +/- 3.93e-48]
```

## Proof Boundary

```text
Failure of a sufficient bound does not imply H=H_x=0, a multiple zero, Lambda>0, or failure of RH
Q_2 requires a sharper tail decomposition, additional theta blocks, or a different direct boundary/winding certificate
```

validated Newman theta second diagonal-shell first-block route guard: 6 rows, 0 issues, 5 exact point audits, 4 strict two-branch failures, maximum failed ratio <24/25, 1 second-stage route guard, 1 open replacement handoff
