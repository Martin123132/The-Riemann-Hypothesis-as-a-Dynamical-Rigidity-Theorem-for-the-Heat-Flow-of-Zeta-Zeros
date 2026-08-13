# Newman C1 Nonterminal Determinant-Completion Gate

Date: 2026-08-05

Status: the proposed E_FAN+E_quad determinant recombination is falsified as an exact identity; one omitted retained/nonterminal linear package is restored exactly; the arithmetic sign remains open; not a proof of RH.

## Exact Residual Split

```text
R_raw=2p^TJ dot(a+e)+2dot(p)^TJ(a+e)+2(a+e)^TJ dot(a+e)+t_E^T dot(a+e).
R_raw=R_nr+R_quad+R_a*+R_aa+R_aff.
R_nr=2p^TJ dot(e)+2dot(p)^TJe=partial_xi(2p^TJe).
R_quad=2e^TJ dot(e).
R_a*=2(p+e)^TJ dot(a)+2(dot(p)+dot(e))^TJa.
R_aa=2a^TJ dot(a).
R_aff=t_E^T[dot(a)+dot(e)].
```

The residual sector represented by terminal mixed, pure terminal, edge affine, and nonterminal quadratic terms differs from R_raw by exactly R_nr.

## Full-Support Polynomial Owner

With alpha_p=alpha_N, beta_p=beta_N-beta_E, P_nr(lambda)=F_(alpha_p)(lambda)+i(lambda-log a)F_(beta_p)(lambda), and P_E(lambda)=i(lambda-log a)F_(beta_E)(lambda), one has P_lin^[N]=P_nr+P_E.

Using mathscr O_N=mathcal T_N+mathscr R_N, Re{nu c_xi mathscr O_N[P_lin^[N]]}=Re{nu c_xi mathcal T_N[P_nr]}+Re{nu c_xi mathscr R_N[P_nr]}+R_aff. The middle term is R_nr and is neither terminal mixed nor edge affine.

Define E_nr:=2R_nr/|nu|^2. The eight-package ledger omitted E_nr; E_quad^cur contains only R_quad and cannot absorb the retained/nonterminal linear cross.

## Determinant Completion

```text
R_nr+R_quad=2{(p+e)^TJ[dot(p)+dot(e)]-p^TJdot(p)}=partial_xi{(p+e)^TJ(p+e)-p^TJp}.
Define Q_nt=(2/|nu|^2){(p+e)^TJ(p+e)-p^TJp}. Then partial_xi Q_nt=E_nr+E_quad^cur exactly.
```

The proposed completion E_FAN+E_quad^cur is not a determinant-current identity. E_FAN belongs to the terminal/edge-near channel; the exact determinant partner of E_quad^cur is E_nr.

The independent exact-rational specialization gives missing R_nr = `-3410569/60843510`, so the proposed and completed nonterminal values are `-120401707/334639305` and `-278319673/669278610`. This is an exact independent-vector algebraic specialization that disproves the proposed universal bookkeeping identity; it is not asserted to be a physical Xi state.

## Corrected Pointwise Ledger

```text
E_tot^corr=E_F+E_R+E_term+E_Delta+E_aff+E_aa+E_move^cur+E_nr+E_quad^cur.
E_main^corr:=E_FAN+E_nr+E_quad^cur.
```

The bounded six-package fraction remains 289339/858000. Therefore adverse-witness negativity necessarily requires E_main^corr<-(551501/858000)rho A_T, while E_main^corr<-(1130179/858000)rho A_T is triangle-safe sufficient.

The constants do not change because E_nr is restored to the unbounded signed main, not assigned an invented absolute allowance.

## Next Action

Decompose the exact completed nonterminal determinant current E_nr+E_quad^cur through the retained finite, near, and common-cutoff paired-remote functionals. Keep E_FAN separate, and test whether the completed nonterminal current has an exact endpoint transfer or a source-correlated sign before taking moduli.

## Pi Provenance

No new pi is introduced. This correction uses only the symmetric observation matrix J, linearity, and the already normalized full-support functionals.

## Proof Boundary

This gate proves a bookkeeping defect and its exact algebraic repair. It proves no bound or sign for E_nr, E_quad, their completed determinant current, E_FAN, or the complete current; no all-q transport, contact exclusion, Lambda<=0, PF-infinity, RH, or prize-level conclusion.

built Newman C1 nonterminal determinant-completion gate: 17 rows, 0 issues, 8 symbolic audits, 5 exact-rational audits, 7 source audits, 1 missing linear package identified, 1 corrected determinant completion, 1 open arithmetic row
