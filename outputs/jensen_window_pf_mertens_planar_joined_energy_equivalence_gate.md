# Jensen-Window PF Mertens Planar Joined-Energy Equivalence Gate

Date: 2026-07-24

Status: exact equivalence between the proposed joined
current-interior/boundary coordinate and the inherited lossless low-mode
energy, with a summable coordinate defect. This is a route guard, not an
estimate and not a proof of RH, PF-infinity, `Lambda<=0`, or a Clay-prize
result.

```text
work/rh_compute/results/
jensen_window_pf_mertens_planar_joined_energy_equivalence_gate.json
python work/rh_compute/scripts/
jensen_window_pf_mertens_planar_joined_energy_equivalence_gate.py
python work/rh_compute/scripts/
check_jensen_window_pf_mertens_planar_joined_energy_equivalence_gate.py
```

## The Full Off-Diagonal Already Knows The Positive Energy

Retain the notation of Corollaries 11.22Z.8-Z.10. Define

```text
L
 :=sum_(r=1)^R (x_r+beta e_r)^2/q_r^2
  =E_perp+c(beta+gamma/c)^2,

D_L
 :=sum_(r=1)^R D_(C,r)/q_r^2+cQ.                    (PJE.1)
```

The three-block compression from Corollary 11.22Z.8 is exactly

```text
O_(alpha,K)=(L-D_L)/2.                              (PJE.2)
```

The diagonal estimate from Corollary 11.22Z.3 gives

```text
0<=D_L<7pi^2/24,

|O_(alpha,K)-L/2|=D_L/2<7pi^2/48.                  (PJE.3)
```

Thus the complete signed off-diagonal is already the positive low-mode
energy up to a uniformly bounded diagonal.

## The Proposed Join Returns The Same Coordinate

Corollary 11.22Z.9 gives

```text
O_(alpha,K)=O_(CC,int)+R_B,

R_B=(E_B-D_B)/2+C_delta.                            (PJE.4)
```

Define the candidate joined quantity

```text
J_B:=O_(CC,int)+E_B/2.                              (PJE.5)
```

Eliminating `O_(alpha,K)` between (PJE.2) and (PJE.4) gives the exact
identity

```text
J_B
 =L/2+(D_B-D_L)/2-C_delta.                          (PJE.6)
```

The endpoint estimates already proved in Corollary 11.22Z.9 give

```text
|C_delta|<=pi H_R^odd+pi^2/2,
0<=D_B<5pi^2/12.
```

Together with (PJE.3),

```text
|J_B-L/2|
 <=pi H_R^odd+41pi^2/48,

|J_B-O_(alpha,K)|
 <=pi H_R^odd+17pi^2/24.                            (PJE.7)
```

For polynomially bounded `R`, every right side is `O(1+log K)`.
Consequently,

```text
J_B=O_epsilon(K^epsilon) for every epsilon>0
 iff
O_(alpha,K)=O_epsilon(K^epsilon) for every epsilon>0
 iff
L=O_epsilon(K^epsilon) for every epsilon>0.          (PJE.8)
```

For each fixed `alpha>0`, the stronger weighted statement is also exact:

```text
sum_(K dyadic)K^(-alpha)|J_B|<infinity
 iff
sum_(K dyadic)K^(-alpha)|O_(alpha,K)|<infinity
 iff
sum_(K dyadic)K^(-alpha)L<infinity.                  (PJE.9)
```

On one fixed cofinal positive-alpha sequence, the final member is the
prior lossless RH criterion.

## Route Decision

Keeping `O_(CC,int)` signed with `E_B/2` does not create a new weaker
interface. It returns `L/2` up to the same summable endpoint and diagonal
coordinate changes already isolated by the planar algebra. The options are
therefore honest:

```text
prove the RH-strength E_B anchored theorem;
prove the lossless L theorem;
prove the strong signed nested-edge theorem; or
find new structure outside this equivalent planar join.
```

The `Y_(K,r)` projection remains a useful exact coordinate, but pairing its
current-interior value with `E_B/2` cannot lower the theorem strength.

## Finite Mobius Diagnostic

At `alpha=1/2`, direct recomputation gives:

| K | R | L | D_L | O | E_B | O_CC,int | J_B | J_B-L/2 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 16 | 8 | 0.086714164 | 1.276405173 | -0.594845505 | 0.110319464 | -0.001713988 | 0.053445744 | 0.010088662 |
| 32 | 14 | 0.074555912 | 1.268291221 | -0.596867655 | 0.083247421 | -0.000825629 | 0.040798082 | 0.003520126 |
| 64 | 23 | 0.129088701 | 1.278746522 | -0.574828911 | 0.135323930 | 0.000052069 | 0.067714034 | 0.003169684 |
| 128 | 39 | 0.143665376 | 1.315481547 | -0.585908085 | 0.147238007 | 0.000053254 | 0.073672258 | 0.001839570 |
| 256 | 64 | 0.112963133 | 1.324398989 | -0.605717928 | 0.115185299 | -0.000166902 | 0.057425747 | 0.000944181 |
| 512 | 108 | 0.109596358 | 1.323092731 | -0.606748187 | 0.111068291 | -0.000042452 | 0.055491694 | 0.000693515 |
| 1024 | 182 | 0.096999647 | 1.325065222 | -0.614032788 | 0.097890613 | 0.000020860 | 0.048966167 | 0.000466344 |


The final column decreases on this finite grid because the two coordinate
systems are already close there. This is diagnostic only. The exact
identity (PJE.6), not the table, proves the route equivalence.

No estimate for `L`, `E_B`, or `Y_(K,r)`, no RH, no PF-infinity, no
`Lambda<=0`, and no Clay-prize result is proved.

## Claim Ledger

| ID | Role | Status | Statement | Boundary |
|---|---|---|---|---|
| `pjeg_01_inherited_data` | `exact_definition` | `available_exact` | Retain O_(CC,int), E_B, D_B, C_delta, x_r, e_r, beta, Q, c, and D_(C,r) from Corollaries 11.22Z.8-Z.10. | No arithmetic estimate is assumed. |
| `pjeg_02_lossless_energy` | `exact_definition` | `available_exact` | Set L=sum_(r<=R)(x_r+beta e_r)^2/q_r^2=E_perp+c(beta+gamma/c)^2. | L is nonnegative and is the inherited low-mode energy. |
| `pjeg_03_lossless_diagonal` | `exact_definition` | `available_exact` | Set D_L=sum_(r<=R)D_(C,r)/q_r^2+cQ. | D_L is the current/future diagonal in the same mode coordinates. |
| `pjeg_04_block_compression` | `exact_identity` | `available_exact` | O_(alpha,K)=(L-D_L)/2. | This is the exact three-block compression, not an estimate. |
| `pjeg_05_lossless_diagonal_bound` | `exact_bound` | `available_exact` | 0<=D_L<7pi^2/24. | The bound is uniform in alpha, K, and R. |
| `pjeg_06_offdiagonal_lossless_comparison` | `exact_bound` | `available_exact` | |O_(alpha,K)-L/2|=D_L/2<7pi^2/48. | The full signed off-diagonal and positive energy differ boundedly. |
| `pjeg_07_suffix_remainder` | `exact_identity` | `available_exact` | R_B=(E_B-D_B)/2+C_delta and O_(alpha,K)=O_(CC,int)+R_B. | The boundary shoulder remains included in R_B. |
| `pjeg_08_defect_bound` | `exact_bound` | `available_exact` | |C_delta|<=pi H_R^odd+pi^2/2. | The two terms are the shoulder and anchor endpoint corrections. |
| `pjeg_09_joined_candidate` | `exact_definition` | `available_exact` | J_B:=O_(CC,int)+E_B/2. | This is the proposed signed current-interior/boundary join. |
| `pjeg_10_joined_identity` | `exact_identity` | `available_exact` | J_B=L/2+(D_B-D_L)/2-C_delta. | No Cauchy-Schwarz or absolute cell split is used. |
| `pjeg_11_joined_error_bound` | `exact_bound` | `available_exact` | |J_B-L/2|<=pi H_R^odd+41pi^2/48. | Use D_B<5pi^2/12 and D_L<7pi^2/24. |
| `pjeg_12_joined_offdiagonal_bound` | `exact_bound` | `available_exact` | |J_B-O_(alpha,K)|<=pi H_R^odd+17pi^2/24. | This is the same summable coordinate defect as in PBSE.15. |
| `pjeg_13_pointwise_subpower_equivalence` | `exact_equivalence` | `available_exact` | For polynomially bounded R, J_B, O_(alpha,K), and L are all-epsilon subpower simultaneously. | The logarithmic defects are absorbed for every positive epsilon. |
| `pjeg_14_weighted_dyadic_equivalence` | `exact_equivalence` | `available_exact` | For every fixed alpha>0, sum K^(-alpha)L is finite iff sum K^(-alpha)|J_B| is finite, iff sum K^(-alpha)|O_(alpha,K)| is finite. | The K^(-alpha) logarithmic and bounded defects are summable. |
| `pjeg_15_cofinal_composition` | `theorem_composition` | `available_exact` | On one fixed cofinal positive-alpha sequence, the joined J_B criterion is the prior lossless RH criterion. | This is a route identification, not a proof of its antecedent. |
| `pjeg_16_route_collapse` | `route_guard` | `guard_active` | Retaining O_(CC,int) with E_B/2 does not create a new weaker theorem interface; it returns L/2 up to a summable coordinate defect. | Separate current-interior/E_B work may not be advertised as easier. |
| `pjeg_17_finite_mobius_audit` | `finite_validation` | `validated_finite` | At alpha=1/2 and K=16,...,1024, direct Mobius recomputation verifies the joined identity to float64 roundoff. | Finite values prove no asymptotic estimate. |
| `pjeg_18_open_lossless_gate` | `theorem_target` | `open` | Prove dyadic summability of K^(-alpha)L for every member of one fixed cofinal positive-alpha sequence. | This remains the RH-equivalent low-mode arithmetic obligation. |
| `pjeg_19_next_route` | `route_guard` | `guard_active` | Future work must attack L, E_B, or the strong signed edge form honestly at RH strength, or find structure outside this algebraically equivalent join. | Positivity and coordinate changes alone cannot supply cancellation. |
| `pjeg_20_proof_boundary` | `proof_guard` | `guard_active` | The joined identity, diagonal bounds, subpower equivalence, and route collapse are exact. | No L, E_B, or Y estimate, RH, PF-infinity, Lambda<=0, or Clay-prize result is proved. |

Summary:

- rows: 20
- exact/theorem-composed results: 16
- finite validations: 1
- active route/proof guards: 3
- open lossless energy gates: 1
