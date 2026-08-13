# Newman C1 Fixed-Cell Low-N Determinant Evaluator Gate

Date: 2026-08-05

Status: finite evaluator contract numerically cross-validated; physical-scale compression and signed determinant open; not a proof of RH.

## Contract

For P(lambda)=sum_(k=0)^5 c_k lambda^k, tabulate mathscr M_N, the p=q, |p-q|=1, and |p-q|>=2 pieces of mathscr F_N, mathscr N_N, and mathscr C_N on lambda^0,...,lambda^6. Degree six is retained because the fixed-cell tangent lift i(lambda-log a)P raises degree by one.

J_q=[max(1,q-1/2),min(N,q+1/2)] and C_p={r>0:p-1/2<=alpha_P/r<p+1/2}. Every roster mode is assigned once. The block sum is checked against both the direct finite Dirichlet kernel and H_N{f_P(1)-f_P(N)}+V_N[P].

The common-cutoff paired remote functional is evaluated without truncating an infinite series: sum_(k>n)1/k^2 and sum_(k>n)cos(2pi k u)/k^2 are replaced by zeta(2) minus finite sums and the periodic Bernoulli B_2 kernel after two exact integrations by parts.

With S and P fixed and z=nu c_xi carrying exp(-i xi log a), partial_xi Re{z L[P]}=Re{z L[i(lambda-log a)P]} on a fixed roster-and-reciprocal-cell chart.

For b=omega_E+p, d_0=f^(0)-p, d_1=f^(1), d_2=f^(2), the evaluator computes Q(b-d)-Q(b) and its diagonal-first primitive/current telescope. A centered xi finite difference independently checks the lifted current while every mode label remains fixed.

## Audit Summary

- Fixtures: `3`.
- Total finite modes: `184`.
- Largest block/direct discrepancy: `5.113e-12`.
- Largest block/variation discrepancy: `3.038e-12`.
- Largest outer near/remote discrepancy: `3.337e-12`.
- Largest determinant-current finite-difference discrepancy: `6.564e-10`.

The fixtures use `N=4,5,7`, retain `a^2-1<alpha_P<=a^2`, and stay at least `0.1` in alpha from every relevant half-open transfer. The stored basis includes the exact ideal relative cubics and a synthetic four-row determinant projection.

## Scope

The fixtures retain the physical scale relation a^2-1<alpha_P<=a^2, the exact ideal cubics, the completed Fourier normalization, and all endpoint halves, but use small N and synthetic four-row polynomial coefficients. They validate the evaluator contract; they are not low-height Xi data and carry no physical sign conclusion.

## Next Action

Extract the four actual physical coefficient rows P_X and the genuine edge vector at each calibrated N about 7.2e10. Replace only the direct block engine by a tie-complete compressed diagonal/adjacent plus rail-compressed far evaluator, while preserving this basis, lift, near/remote, and determinant API as regression oracles.

## Pi Provenance

Every pi is forced by e(x)=exp(2pi i x), alpha_P=xi/(2pi), kappa=1/(2pi i), and the exact Fourier/Bernoulli identities. No geometric or fitted pi is introduced.

## Proof Boundary

This gate proves an executable low-N contract and cross-representation numerical identities to recorded tolerances. It is floating-point quadrature, not interval arithmetic. It proves no physical-N block value or sign, no retained determinant bound, no complete-current inequality, no all-q transport, contact exclusion, Lambda<=0, PF-infinity, RH, or prize-level conclusion.

built fixed-cell low-N determinant evaluator gate: 20 rows, 3 fixtures, 184 finite modes, 3 independent finite-band representations, exact finite paired-remote kernel, primitive/current telescopes and xi finite differences validated, 0 physical signed bounds
