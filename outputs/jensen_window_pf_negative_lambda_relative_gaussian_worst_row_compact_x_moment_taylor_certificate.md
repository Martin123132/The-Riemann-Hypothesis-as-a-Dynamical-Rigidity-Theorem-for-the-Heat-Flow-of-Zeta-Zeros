# Jensen-Window PF Negative-Lambda Relative-Gaussian Worst-Row Compact X-Moment Taylor Certificate

Date: 2026-07-09

Status: worst-row compact x-moment Taylor certificate. This is not a proof
of a complete worst-row expectation, an all-row finite-grid certificate,
a uniform collar theorem, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_x_moment_taylor_certificate`.

Proof boundary: this artifact certifies the finite `n<=30` cancellation-reduced
value and derivative cores on all of `0<=y<=200` for `T=10000`, `F_21`.
The separate `n>30` Phi tail and `y>200` far tail are not composed here.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_x_moment_taylor_certificate.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_x_moment_taylor_certificate.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_x_moment_taylor_certificate.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian worst-row compact x-moment Taylor certificate: 7 rows, 0 issues, 129 exact-moment rows, 1 certified compact interval, 0 open compact panels, 0 ready-to-apply rows
```

## Certified Setup

```text
T: 10000
index: F_21
alpha: 20.5
y interval: 0<=y<=200
x interval: 0<=x<=sqrt(0.02)
x right endpoint: 0.1414213562373095048801688724209698078569671875376948073176679737990732
disk radius: 0.38
Taylor degree: 128
exact transformed-moment rows: 129
Arb precision bits: 8192
```

Exact compact moment:

```text
mu_k(200)=T^(-k/2)*lower_gamma(alpha+1+k/2,200)/Gamma(alpha+1)
```

## Complex-Disk Remainder

```text
For |z|<=R<pi/8, Re(exp(4z))>=exp(-4R)*cos(4R)>0; hence |Phi_30(z)| is bounded by the recorded positive termwise majorant.
4R upper: 1.520000000000000000000000000000000000000
pi/2 lower: 1.570796326794896619231321691639751442099
4R<pi/2 certified: True
Re(exp(4z)) lower: 0.01110498340884990488555969851233503491016313452332698847251506357441587
finite Phi_30 disk majorant upper: 1768508.446734926415433097204233572413431274458023462036057859663537778
x_max/R upper: 0.3721614637823934338951812432130784417288610198360389666254420363133507
value integrated Cauchy radius upper: 2.655198459050798383387346946068294311858052394623814612402050072867752E-49
derivative integrated Cauchy radius upper: 1.720472565601475716457020246375098967791947816209975921130278988840058E-47
value tail-radius/cap ratio upper: 3.915048413229948731882778565524718055955968130406233194894092243554826E-10
derivative tail-radius/cap ratio upper: 1.208004651619589888088535471069610017484967372017378180314715491575033E-9
both tail radii below caps: True
```

## Compact Enclosures

```text
value Taylor-moment integral: [-6.5833869236101933120086975942979724158233913100544908541959385799240163069244594E-34 +/- 1.23E-114]
derivative Taylor-moment integral: [-1.3805838741523015640267420526051885567628123046821874080802380609833339891072295E-32 +/- 4.29E-112]
value compact total ball: [-6.58338692361019E-34 +/- 5.97E-49]
derivative compact total ball: [-1.38058387415230E-32 +/- 3.29E-47]
value certified negative: True
derivative certified negative: True
finite-core compact interval certified: True
panel interpolation needed for finite core: False
```

Selected coefficient/moment rows:

```text
k=0: moment=[0.999999999999999999999999999999999999999999999999999999999979 +/- 3.82E-61]; value contribution=[+/- 2.86E-2464]; derivative contribution=0
k=1: moment=[0.0460993168396280314412204760522750604907656687203771316352278 +/- 2.85E-62]; value contribution=[1.54944295517315770638031069808913441075107761439301714156304E-1301 +/- 2.11E-1361]; derivative contribution=[7.74721477586578853190155349044567205375538807196508570781521E-1302 +/- 5.29E-1363]
k=2: moment=[0.00214999999999999999999999999999999999999999999999999999999959 +/- 4.67E-63]; value contribution=[-4.35868111731608048871646620366619436736464693661805108920604E-1299 +/- 2.07E-1360]; derivative contribution=[-4.35868111731608048871646620366619436736464693661805108920604E-1299 +/- 2.07E-1360]
k=20: moment=[1.32721227041473608398437499999999999999999999999997765013906E-26 +/- 3.79E-86]; value contribution=[-6.11489066288183665000513295239539457031930664464454793596955E-1267 +/- 2.34E-1327]; derivative contribution=[-6.11489066288183665000513295239539457031930664464454793596955E-1266 +/- 2.34E-1326]
k=40: moment=[4.69962860638360176510883519049111008644103760821300057971935E-51 +/- 4.21E-111]; value contribution=[-2.27152658625298118051765399630560010177323039153067685181970E-1239 +/- 4.27E-1299]; derivative contribution=[-4.54305317250596236103530799261120020354646078306135370363941E-1238 +/- 1.47E-1298]
k=41: moment=[3.01841883879028318884205657134782361369292791463922312248700E-52 +/- 4.69E-112]; value contribution=[4.23662163919390050007254078614444767010501096305131603791368E-1238 +/- 3.38E-1298]; derivative contribution=[8.68507436034749602514870861159611772371527247425519787772304E-1237 +/- 2.92E-1297]
k=42: moment=[1.95034587164919473252016660405381068587302672412561225142482E-53 +/- 3.52E-113]; value contribution=[-6.78203224787260481759133359940609477593988520207382502315849E-34 +/- 4.15E-94]; derivative contribution=[-1.42422677205324701169418005587527990294737589243550325486328E-32 +/- 2.04E-92]
k=64: moment=[5.28582329569339531311654321901731037518237435064812359044936E-79 +/- 2.22E-139]; value contribution=[-2.15153633200047295023466782013125503417431972029244770410507E-50 +/- 2.99E-110]; derivative contribution=[-6.88491626240151344075093702442001610935782310493583265313621E-49 +/- 4.44E-109]
k=80: moment=[5.85194730834193648846313944573388618018756117966711150383207E-97 +/- 1.18E-157]; value contribution=[-1.43458632041096045135886379072762958716867610770037760461457E-61 +/- 3.48E-121]; derivative contribution=[-5.73834528164384180543545516291051834867470443080151041845827E-60 +/- 3.91E-120]
k=96: moment=[1.85543724313476833512941805830018092717428784151073179613707E-114 +/- 3.10E-174]; value contribution=[-3.11270895006911792525436333909727047341756202453381313594679E-73 +/- 2.26E-133]; derivative contribution=[-1.49410029603317660412209440276668982724042977177623030525446E-71 +/- 2.84E-132]
k=112: moment=[1.49044597402595328258726257535198528261228326038108712299543E-131 +/- 2.35E-191]; value contribution=[5.65365747337953488963435361956829143580164435071900976963027E-83 +/- 1.04E-143]; derivative contribution=[3.16604818509253953819523802695824320404892083640264547099295E-81 +/- 6.21E-142]
k=128: moment=[2.75299452860557002904962075055915296362630262547921497178869E-148 +/- 2.56E-208]; value contribution=[4.15861764717418026652593633688187819120340142809227349514372E-93 +/- 4.74E-153]; derivative contribution=[2.66151529419147537057659925560440204237017691397905503689198E-91 +/- 3.83E-151]
```

## Remaining Splice

```text
Compose the separately certified n>30 Phi/Phi'/Phi(0) tail with matching normalization.
Compose the separately certified y>200 far-tail bound.
Then extend the compact/tail certificate from this worst row to every required grid row and the collar.
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_moment_taylor_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_moment_taylor_certificate.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_arb_chebyshev_interpolant_moment_scout.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_interval_integration_scout.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_plus_tail_budget_certificate.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md
outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md
outputs/signed_hankel_jensen_dependency_graph.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.md
```

Summary:

The worst-row compact x-moment Taylor certificate replaces the six-panel interpolation obligation for the finite n<=30 core. For T=10000 and F_21, the full range 0<=y<=200 maps inside |x|<=sqrt(0.02)<0.38; a degree-128 Arb Taylor model is integrated by exact transformed Gamma moments, and the true-function Cauchy remainder radii consume less than 1.3e-9 of their source caps. Both compact integrals are certified negative. This closes only the finite-core compact source for one row; n>30 tail composition, y>200 far-tail composition, all-row coverage, and the uniform collar theorem remain open.
