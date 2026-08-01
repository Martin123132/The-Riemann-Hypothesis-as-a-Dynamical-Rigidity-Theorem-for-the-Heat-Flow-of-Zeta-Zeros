# Newman Theta Two-Block Shell Frontier Scout

Date: 2026-07-24

Status: finite interval diagnostic, not a proof of RH or `Lambda <= 0`.

This scout measures where the fixed two-block first-jet method continues to certify finite diagonal shells. Passing finitely many rows cannot establish a cofinal theorem.

```text
work/rh_compute/results/jensen_window_pf_newman_theta_two_block_shell_frontier_scout/frontier.jsonl
python work/rh_compute/scripts/jensen_window_pf_newman_theta_two_block_shell_frontier_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_theta_two_block_shell_frontier_scout.py
```

## Fixed Configuration

```text
{"cpu_consecutive_park_samples": 2, "cpu_park_threshold_percent": 85.0, "cpu_sample_every_evaluations": 25, "cpu_sample_interval_seconds": 0.2, "endpoint_serialization_digits": 50, "max_depth": 6, "max_evaluations_per_shell": 20000, "max_shell_seconds": 180.0, "precision_bits": 160, "retained_theta_blocks": [1, 2], "schema_version": 1, "start_j": 5, "tail_moment_one_upper": "1/1000000000000", "tail_moment_zero_upper": "1/10000000000", "time_step": "1/100", "x_step": "1/2"}
```

## Shell Rows

| j | status | initial | evaluated | subdiv | negative | value-separated | derivative-only | minimum ratio | seconds |
|---:|:---|---:|---:|---:|---:|---:|---:|:---|---:|
| 5 | certified | 160 | 160 | 0 | 112 | 140 | 20 | [39673.192546178463746125623164631306575741582437849 +/- 1.94e-46] | 8.923 |
| 6 | certified | 204 | 204 | 0 | 119 | 182 | 22 | [8990.0945486841265568648977962700810179803242025827 +/- 4.62e-47] | 9.904 |
| 7 | certified | 252 | 252 | 0 | 126 | 229 | 23 | [8588.6078859585630749141901202499830514490872451622 +/- 3.61e-47] | 12.421 |
| 8 | certified | 288 | 288 | 0 | 126 | 264 | 24 | [7709.0988208757710009748745200655328551536715495393 +/- 1.77e-47] | 14.036 |
| 9 | certified | 324 | 324 | 0 | 126 | 300 | 24 | [7707.2311134837930375111394906100643420545539351209 +/- 1.56e-47] | 15.800 |
| 10 | certified | 360 | 360 | 0 | 126 | 336 | 24 | [7705.7340864985282124668348550123357088920847118185 +/- 2.32e-47] | 16.886 |
| 11 | certified | 418 | 418 | 0 | 133 | 393 | 25 | [7704.5073578926565039291892558946492767370649831468 +/- 2.93e-47] | 20.054 |
| 12 | certified | 456 | 456 | 0 | 133 | 411 | 45 | [7703.4837848084382958820478304214983513879268872209 +/- 2.61e-47] | 21.176 |
| 13 | certified | 494 | 494 | 0 | 158 | 436 | 58 | [4079.1354819477682457322095592966849104606692033803 +/- 4.28e-48] | 24.540 |
| 14 | certified | 532 | 532 | 0 | 196 | 474 | 58 | [1353.0325730058331130321848955853350392121810597572 +/- 1.13e-47] | 25.557 |
| 15 | certified | 570 | 570 | 0 | 234 | 512 | 58 | [743.08066903553526895579677979320634644396035474167 +/- 3.96e-48] | 26.554 |
| 16 | certified | 608 | 608 | 0 | 272 | 550 | 58 | [642.09171845808649269087196281130867667093087753804 +/- 4.74e-48] | 28.122 |
| 17 | certified | 646 | 646 | 0 | 310 | 588 | 58 | [642.05980411163639831970924670057157501149739321113 +/- 1.36e-48] | 32.745 |
| 18 | certified | 684 | 684 | 0 | 348 | 626 | 58 | [642.03142349272283376108199868495689610355979909638 +/- 2.24e-48] | 32.320 |
| 19 | certified | 722 | 722 | 0 | 386 | 664 | 58 | [642.00602079505582500773437136307463381006394715781 +/- 4.61e-48] | 33.567 |
| 20 | certified | 760 | 760 | 0 | 424 | 702 | 58 | [641.98315055128624290884112046315136387206459207014 +/- 4.82e-48] | 36.667 |
| 21 | certified | 840 | 840 | 0 | 487 | 780 | 60 | [440.73929632072420073612621628637015591812774672174 +/- 5.05e-49] | 37.016 |
| 22 | certified | 880 | 880 | 0 | 486 | 779 | 101 | [180.81950539939247606594251398047442200850003730262 +/- 3.10e-48] | 42.475 |
| 23 | certified | 920 | 934 | 7 | 486 | 779 | 148 | [1.3709314869015039163669568274457936340118136691538 +/- 2.92e-50] | 46.754 |
| 24 | certified | 960 | 1054 | 47 | 486 | 849 | 158 | [1.4014758640469298295401483252276762213618642071908 +/- 2.80e-50] | 52.437 |
| 25 | certified | 1000 | 1174 | 87 | 486 | 928 | 159 | [1.4295761528089343528433856115959048061127350709239 +/- 4.73e-50] | 55.942 |
| 26 | certified | 1040 | 1294 | 127 | 486 | 1008 | 159 | [1.4555135704549989002644257290115198126741933165371 +/- 2.14e-50] | 65.270 |
| 27 | certified | 1080 | 1414 | 167 | 486 | 1054 | 193 | [1.4795293613046688330210544969084168294027622207015 +/- 4.90e-50] | 71.457 |
| 28 | certified | 1120 | 1534 | 207 | 486 | 1054 | 273 | [1.5018283424750188382046630784389549877905209751392 +/- 4.45e-50] | 74.286 |
| 29 | certified | 1160 | 1654 | 247 | 486 | 1054 | 353 | [1.5225895319228132027263760861688558063479439490436 +/- 2.78e-50] | 82.298 |
| 30 | certified | 1200 | 2006 | 403 | 642 | 1210 | 393 | [1.5419658611909556051403134870261010993034534046432 +/- 4.34e-50] | 98.408 |

## Proof Boundary

The rows are rigorous finite interval diagnostics for the stored boxes, but this artifact is not promoted as a theorem for any new `Q_j`. A cofinal result needs an analytic `j`-dependent separation theorem or another uniform mechanism.

validated Newman theta two-block shell frontier scout: 26 shell rows, 0 issues, j=5..30, 26 certified shells, first finite frontier=none
