# Newman C1 Reciprocal-Stationary Disk Geometry Gate

Date: 2026-08-04

Status: exact anchored reciprocal decomposition and disk geometry proved; separated termwise absolute budgets obstructed; signed complete-join estimate open; not a proof of RH.

## Exact combined residuals

```text
W_H:=R_H=A_H+mathcal J_H^0,
W_T:=R_T+Q_T=A_T+mathcal J_T^0.
```

The exact sufficient condition |R_T+Q_T|<=A_T-delta is weaker than the previous separated condition |R_T|+|Q_T|<=A_T-delta and still implies |mathcal J_T^0|>=delta.

For each reciprocal cell,

```text
2M_q-S_q=M_q-tau_q+E_q^(ext)-Lambda_q^(alias).
```

Therefore W_H and W_T are A_H and A_T plus the complete tie-invariant reciprocal joins. Internal aliases, exterior pieces, and tie halves must recompose before a modulus is taken.

## Tangent-disk geometry

```text
|A+J|^2-A^2=|J|^2+2A Re(J).

|A+J|<=A-delta
iff |J|^2+2A Re(J)<=-2A delta+delta^2.
```

Any nonzero J in the closed tangent disk |A+J|<=A satisfies Re(J)<=-|J|^2/(2A)<0. A magnitude-only estimate for J cannot establish this directional condition.

## Absolute-budget obstruction

Let I_N={q in Z:ceil(N/64)<=q<=floor(N/32)}. On L>=50, N>4160 and #I_N>=N/65; every such q is an interior carrier with full starred weight.

Since sigma<1/2+1/(8L), log(a)<(L+1)/2, and (L+1)/(16L)<=51/800<log(10/9), exp(S(log q))>9/(10 sqrt(q))>9/(2 sqrt(N)) on I_N.

For q in I_N, u_q-u_N=log(N/q)>=log 32>3 and u_q+u_N>=u_q-u_N. Thus |P_H^0(log q)|>27/4. The imaginary cubic part of P_T^0 has the same floor and is orthogonal to its real affine part, so |P_T^0(log q)|>27/4.

```text
B_H>A_H,   B_T>A_T,                  L>=50, q=1.
```

At L=50 the 41-term positive Taylor lower sum for exp(25/2) gives (243/520)exp(25/2)>125396 while U(50)<42531. The logarithmic derivative of U is below 4/(L+1)<=4/51<1/4, so the strict budget gap increases for all L>=50.

Therefore B_H>A_H and B_T>A_T throughout the q=1 chart. A proof that first applies the triangle inequality separately to carrier atoms, and then adds any nonnegative variation budget, cannot close either anchor disk. This does not lower-bound |mathscr M_N| or disprove signed cancellation.

## Gate rows

| id | role | state | claim |
|---|---|---|---|
| rsd_01_hermitian | combined residual | proved | The Hermitian disk residual is A_H plus the complete ideal join. |
| rsd_02_transpose | combined residual | proved | The transpose correction belongs inside one combined residual. |
| rsd_03_improvement | optimal handoff | proved | The combined transpose disk strictly improves the separated triangle handoff. |
| rsd_04_cell | reciprocal cell | proved | Each finite reciprocal cell has an exact one-carrier join. |
| rsd_05_global | global recomposition | proved | Internal ties and aliases recompose into one global exterior. |
| rsd_06_anchor_stationary | anchored decomposition | proved | Both disk residuals have a complete tie-invariant reciprocal decomposition. |
| rsd_07_disk | disk geometry | proved | The tangent anchor disk has an exact signed quadratic equation. |
| rsd_08_delta | quantitative disk | proved | The contracted disk gives the exact delta margin. |
| rsd_09_sign | sign necessity | proved | Every nonzero point in the anchor disk has negative raw real part. |
| rsd_10_block | carrier block | proved | A macroscopic interior carrier block has a rational counting floor. |
| rsd_11_amplitude | amplitude floor | proved | The physical q=1 weight has an explicit lower floor on that block. |
| rsd_12_hermitian_poly | Hermitian floor | proved | The Hermitian ideal cubic is uniformly nonzero on the test block. |
| rsd_13_transpose_poly | transpose floor | proved | The transpose ideal cubic has the same modulus floor. |
| rsd_14_budget | absolute carrier budget | proved | Both fully termwise carrier budgets grow at least like exp(L/4). |
| rsd_15_anchor_upper | anchor upper bound | proved | Both anchors have a polynomial q=1 upper envelope. |
| rsd_16_comparison | uniform mismatch | proved | The termwise budgets already exceed the anchors for every L>=50. |
| rsd_17_guard | nonpromotion guard | proved | Separated termwise absolute estimates cannot prove either disk inclusion. |
| rsd_18_target | signed stationary target | open | The live target is a signed complete-join estimate in the tangent half-plane, preferably after physical phase composition. |

## Next action

Use the exact global M-tau+E decomposition to seek a signed estimate for the complete Hermitian or transpose join, with aliases, exterior, tie halves, and terminal phase composed before moduli. The raw disk route requires a negative-real tangent-half-plane estimate; alternatively derive a phase-weighted two-channel inequality that proves the physical current directly. Do not return to separated carrier/variation triangle budgets.

## Proof boundary

This gate proves exact combined anchor residuals, an improved transpose handoff, tangent-disk geometry, and a q=1 source-specific obstruction to separated termwise absolute budgets. It proves no signed reciprocal join, disk inclusion, physical-current sign, all-q transport, contact-box exclusion, Q209 shell, cofinal successor, Lambda<=0, RH, or prize-level conclusion.

built Newman C1 reciprocal-stationary disk geometry gate: 18 rows, 0 issues, 2 combined residuals, 2 optimal disk identities, 2 absolute-budget obstructions, 1 live signed stationary target
