# Newman Theta Second Diagonal-Shell Two-Block Interval Certificate

Date: 2026-07-24

Status: rigorous second-stage interval theorem, not a proof of RH.
Stages `j>=3`, `Lambda <= 0`, and RH remain open.

```text
work/rh_compute/results/jensen_window_pf_newman_theta_second_diagonal_shell_two_block_interval_certificate.json
python work/rh_compute/scripts/jensen_window_pf_newman_theta_second_diagonal_shell_two_block_interval_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_theta_second_diagonal_shell_two_block_interval_certificate.py
```

## Two-Block Replacement

```text
Phi(u)=sum_(n>=1)Phi_n(u), Phi_n=(2*pi^2*n^4*exp(9u)-3*pi*n^2*exp(5u))*exp(-pi*n^2*exp(4u))
H_(<=2,t)(x)=integral_0^infinity exp(tu^2)*(Phi_1(u)+Phi_2(u))*cos(xu)du
M_k<=2*pi^2*k!*sum_(n>=3)n^4*exp(-pi*n^2)/(4*pi*n^2-9)^(k+1), k=0,1
Using pi>3 and pi^2<10, consecutive majorants have ratio <4*exp(-21)<1/2; hence M_0<40*81*exp(-27)/99<10^-10 and M_1<40*81*exp(-27)/99^2<10^-12
|J_t-J_(<=2,t)|<16*x^4*10^-10 and |J_t'-J_(<=2,t)'|<64*x^3*10^-10+16*x^4*10^-12
```

The failed first-block curvature-moment bars are not reused. The
actual oscillatory `n=2` contribution is retained, and only the
tiny positive `n>=3` kernel is bounded absolutely.

## Certified Slab

```text
S_2=[1/10,1/5]x[38,40]
initial boxes=25
certified boxes=25
subdivisions=0
unresolved boxes=0
minimum normalized full-kernel ratio=[255984.06834982681328334815749670913675813946043815 +/- 3.71e-45]
minimum normalized value-sign ratio=[39655.801881622599811288517094311291881391746984960 +/- 4.88e-46]
```

| t interval | x interval | branch | ratio lower |
|:---|:---|:---|:---|
| [1/10,3/25] | [38,192/5] | derivative | [437435.33533820617481167366250643985038063802758334 +/- 4.54e-45] |
| [1/10,3/25] | [192/5,194/5] | derivative | [390437.26950638183865761986047072209109613890855896 +/- 2.54e-45] |
| [1/10,3/25] | [194/5,196/5] | derivative | [344362.54223934115702736246275157453643442261933366 +/- 3.67e-45] |
| [1/10,3/25] | [196/5,198/5] | derivative | [300008.20416487294254858103057122162384615169094633 +/- 3.16e-46] |
| [1/10,3/25] | [198/5,40] | derivative | [257991.05426816579916017197049491533630710476060787 +/- 4.33e-46] |
| [3/25,7/50] | [38,192/5] | derivative | [437227.37178744291647209654547862448435569792379576 +/- 1.61e-45] |
| [3/25,7/50] | [192/5,194/5] | derivative | [390121.08629957902629413458701361825360899436710853 +/- 4.68e-45] |
| [3/25,7/50] | [194/5,196/5] | derivative | [343963.78382618895937614995414390717906016854993334 +/- 2.61e-46] |
| [3/25,7/50] | [196/5,198/5] | derivative | [299549.73340720750123034471638408812501215187850309 +/- 7.80e-46] |
| [3/25,7/50] | [198/5,40] | derivative | [257492.50548072165736724110189633586383120271306951 +/- 4.71e-45] |
| [7/50,4/25] | [38,192/5] | derivative | [437015.63741805701832775132582336473741193436646414 +/- 1.78e-45] |
| [7/50,4/25] | [192/5,194/5] | derivative | [389801.56175140957869814212762676770987923853645035 +/- 7.99e-46] |
| [7/50,4/25] | [194/5,196/5] | derivative | [343562.10324437538236281470365073001857265955391586 +/- 2.41e-45] |
| [7/50,4/25] | [196/5,198/5] | derivative | [299088.74446114994499962178870749199488915894339861 +/- 3.65e-45] |
| [7/50,4/25] | [198/5,40] | derivative | [256991.82215673151705410997087821795829943521623272 +/- 1.14e-45] |
| [4/25,9/50] | [38,192/5] | derivative | [436800.13480895675793008040937384309586843935145229 +/- 4.51e-45] |
| [4/25,9/50] | [192/5,194/5] | derivative | [389478.69867223067636511145045363743470439743434321 +/- 5.40e-46] |
| [4/25,9/50] | [194/5,196/5] | derivative | [343157.50399654407482299196019977045321985523907137 +/- 1.49e-45] |
| [4/25,9/50] | [196/5,198/5] | derivative | [298625.24097746628909795027269990137173613946171316 +/- 4.15e-45] |
| [4/25,9/50] | [198/5,40] | derivative | [256489.00832261943421618042056697727627711345593867 +/- 3.54e-45] |
| [9/50,1/5] | [38,192/5] | derivative | [436580.86550679886940896583276590143479912098324146 +/- 2.57e-45] |
| [9/50,1/5] | [192/5,194/5] | derivative | [389152.49988997139743292910900000580616064351235184 +/- 1.05e-45] |
| [9/50,1/5] | [194/5,196/5] | derivative | [342749.98967240757829058653075568733665298308035078 +/- 4.51e-45] |
| [9/50,1/5] | [196/5,198/5] | derivative | [298159.22675894242611013481146334648125079382902799 +/- 9.78e-46] |
| [9/50,1/5] | [198/5,40] | derivative | [255984.06834982681328334815749670913675813946043815 +/- 3.71e-45] |

## Second Stage

```text
J_t(x)<0, equivalently H_t(x)<0, for every (t,x) in [1/10,1/5]x[38,40]
Q_2=[1/10,1/4]x[0,40] contains no common zero of H_t and H_t'
[1/10,1/4]x[-40,40] contains no common zero of H_t and H_t'
Z=H+iH_x is nonzero on partial Q_2 and wind(Z(partial Q_2),0)=0
```

## Proof Boundary

```text
The next stage Q_3=[1/15,1/4]x[0,41] requires a new certificate on [1/15,1/5]x[38,41]; it is not proved here
```

validated Newman theta second diagonal-shell two-block interval certificate: 10 rows, 0 issues, 25 certified slab boxes, 0 subdivisions, 0 unresolved, 25 negative-value boxes, two analytic n>=3 moment bounds, minimum value ratio >39000, 1 second-stage no-contact theorem, 1 zero-winding composition, 1 open third-stage handoff
