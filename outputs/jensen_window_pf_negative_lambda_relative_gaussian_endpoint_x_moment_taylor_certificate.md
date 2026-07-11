# Jensen-Window PF Negative-Lambda Relative-Gaussian Endpoint X-Moment Taylor Certificate

Date: 2026-07-09

Status: endpoint x-moment Taylor certificate. This is not a proof
of the full compact interval, a finite-grid interval certificate, a uniform
collar theorem, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_moment_taylor_certificate`.

Proof boundary: this artifact certifies the finite `n<=30` cancellation-reduced
value and derivative cores on the first panel `0<=y<=1` for `T=10000`,
`F_21`. It does not compose the separate `n>30` tail here and does not
certify the five panels with `y>=1`.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_moment_taylor_certificate.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_moment_taylor_certificate.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_moment_taylor_certificate.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian endpoint x-moment Taylor certificate: 7 rows, 0 issues, 65 exact-moment rows, 1 certified first panel, 5 open compact panels, 0 ready-to-apply rows
```

## Certified Setup

```text
T: 10000
index: F_21
alpha: 20.5
y interval: 0<=y<=1
x interval: 0<=x<=0.01
Phi finite terms: 30
subtracted polynomial M: 20
Taylor degree: 64
Arb precision bits: 8192
exact transformed-moment rows: 65
```

Exact transformed moment:

```text
mu_k=T^(-k/2)*lower_gamma(alpha+1+k/2,1)/Gamma(alpha+1)
```

## Complex-Disk Remainder

```text
For |z|<=R<pi/8, Re(exp(4z))>=exp(-4R)*cos(4R)>0; hence |Phi_30(z)| is bounded termwise by the recorded positive majorant.
disk radius R: 0.2
4R<pi/2 certified: True
Re(exp(4z)) lower: 0.313050504004480024148649034681440220229137988792389098912584
finite Phi_30 disk majorant upper: 95.0575638764905004952154227153911180127486223883434865867288
x/R: 0.05000000000000000000000000000000000000000
value integrated Cauchy radius upper: 9.808948189167392704429706266299725940443101491142294878574993990827353E-104
derivative integrated Cauchy radius upper: 3.190489463634446679651346564512226647996756169224178018399655940174371E-102
```

## First-Panel Enclosures

```text
value Taylor-moment integral: [-2.7778525801082570082696403014722026321372004582856996606427526747473232005120313E-86 +/- 1.82E-166]
derivative Taylor-moment integral: [-5.8332978789969920509137104127278740012007559993199535094730146881310238695894843E-85 +/- 4.60E-165]
value total ball: [-2.7778525801082570E-86 +/- 1.81E-103]
derivative total ball: [-5.8332978789969921E-85 +/- 8.10E-102]
value certified negative: True
derivative certified negative: True
value absolute/cap ratio upper: 4.09589998776490160168595418401036799336466150619608535450185E-47
derivative absolute/cap ratio upper: 4.09576479916002087721501830403351576818148465399064522751207E-47
both channels below caps: True
```

Selected coefficient/moment rows:

```text
k=0: moment=[1.61555602394641035019843832303118644468321656671618684684532E-21 +/- 2.94E-81]; value contribution=[+/- 4.61E-2485]; derivative contribution=0
k=1: moment=[1.57725335167772937773895526899416325205170947283247086561860E-23 +/- 3.86E-83]; value contribution=[5.30130219235593332678286374305186517151665473799095135931004E-1323 +/- 6.68E-1385]; derivative contribution=[2.65065109617796666339143187152593258575832736899547567965502E-1323 +/- 3.35E-1385]
k=2: moment=[1.54072154286455442629772161099652967305438211420263405513836E-25 +/- 2.19E-85]; value contribution=[-3.12349483531434022547314586123861692334621155552383750813262E-1321 +/- 3.76E-1381]; derivative contribution=[-3.12349483531434022547314586123861692334621155552383750813262E-1321 +/- 3.76E-1381]
k=20: moment=[1.08718970898723579281404852234354056721160838535616438741172E-61 +/- 5.61E-122]; value contribution=[-5.00903009146369919626391117630062438940379227243420372480039E-1302 +/- 5.17E-1363]; derivative contribution=[-5.00903009146369919626391117630062438940379227243420372480039E-1301 +/- 5.17E-1362]
k=40: moment=[8.19111547267127599948280571775420794900313832353343520156028E-102 +/- 4.43E-162]; value contribution=[-3.95910786268676422492156513050344655620975523804610887962470E-1290 +/- 3.40E-1350]; derivative contribution=[-7.91821572537352844984313026100689311241951047609221775924941E-1289 +/- 3.21E-1349]
k=41: moment=[8.09133881874827338221713547937623030877196092567730705388426E-104 +/- 4.72E-164]; value contribution=[1.13569199506113623169109596561650858610974262249875718431753E-1289 +/- 1.97E-1349]; derivative contribution=[2.32816858987532927496674672951384260152497237612245222785093E-1288 +/- 2.47E-1348]
k=42: moment=[7.99396239602527294884941394505984100728249972468372934439058E-106 +/- 4.34E-166]; value contribution=[-2.77977929690390669652166712350282810282388397383349153736522E-86 +/- 3.16E-146]; derivative contribution=[-5.83753652349820406269550095935593901593015634505033222846697E-85 +/- 1.37E-145]
k=48: moment=[7.45558708319041618083938053550055723773952165645713245819027E-118 +/- 3.22E-178]; value contribution=[8.79519934822957315842193838141632346945230157208885965658478E-96 +/- 5.32E-157]; derivative contribution=[2.11084784357509755802126521153991763266855237730132631758035E-94 +/- 2.93E-154]
k=56: moment=[6.84122131112447221568522290857994408717006980547935177729690E-134 +/- 2.68E-194]; value contribution=[5.19467524141066694473614398361972132320384896320197969152114E-109 +/- 2.00E-170]; derivative contribution=[1.45450906759498674452612031541352197049707770969655431362592E-107 +/- 8.56E-168]
k=64: moment=[6.32036710315282208670952221708957098646619555258371939939812E-150 +/- 2.80E-210]; value contribution=[-2.57263602911077335248908452085282427340241175552410330212490E-121 +/- 1.81E-181]; derivative contribution=[-8.23243529315447472796507046672903767488771761767713056679969E-120 +/- 4.24E-180]
```

## Remaining Work

```text
Compose the already separate n>30 Phi-tail source without changing normalization.
Certify the five compact panels 1<=y<=200 by analytic-domain/Taylor or Bernstein bounds.
Lift the worst-row compact certificate to all required rows and then to the full collar.
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_panel_route_matrix.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_panel_route_matrix.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_parity_repair_matrix.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md
outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md
outputs/signed_hankel_jensen_dependency_graph.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.md
```

Summary:

The endpoint x-moment Taylor certificate closes the finite n<=30 first-panel integration problem for the T=10000, F_21 worst row. A degree-64 Arb Taylor model is integrated by exact transformed Gamma moments, and a rigorous |z|<=0.2 Cauchy majorant bounds the true-function tail. Both first-panel channels are certified negative and far below their global error caps. This retires the endpoint branch obstruction for the finite core only; the separate n>30 tail composition, five y>=1 panels, all-row coverage, and uniform collar theorem remain open.
