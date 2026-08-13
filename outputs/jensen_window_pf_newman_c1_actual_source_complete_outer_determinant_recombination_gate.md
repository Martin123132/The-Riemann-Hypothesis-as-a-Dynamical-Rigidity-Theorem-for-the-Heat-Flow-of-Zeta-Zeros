# Newman C1 Complete Outer-Determinant Recombination Gate

Date: 2026-08-05

Status: eight residual and affine packages recombined into one geometric determinant current; moving current remains separate; signed arithmetic open; not a proof of RH.

## Observation Owners

p_X=Re{nu c_xi mathscr M_N[P_X]}, a_X=Re{nu c_xi mathcal T_N[P_X]}, e_X=Re{nu c_xi mathscr R_N[P_X]}, and o_X=a_X+e_X=Re{nu c_xi mathscr O_N[P_X]}.

Because mathscr O_N=mathcal T_N+mathscr R_N, one has a+e=o exactly before any modulus.

## Complete Determinant

```text
R_outer=2p^TJdot(o)+2dot(p)^TJo+2o^TJdot(o)=partial_xi{(p+o)^TJ(p+o)-p^TJp}.
R_nr+R_quad+R_(a*)+R_aa=R_outer.
With b=omega_E+p and t_E=2Jomega_E, R_geom=R_outer+R_aff=partial_xi{(b+o)^TJ(b+o)-b^TJb}.
Define Q_geom=(2/|nu|^2){(omega_E+p+o)^TJ(omega_E+p+o)-(omega_E+p)^TJ(omega_E+p)} and E_geom^cur=partial_xi Q_geom.
E_geom^cur=E_F+E_R+E_term+E_Delta+E_aa+E_nr+E_quad^cur+E_aff. Hence E_tot^corr=E_geom^cur+E_move^cur.
```

The independent exact-rational audit gives outer `-24175003120052/2749383478947105`, affine `140378607632/3710369067405`, and geometric total `15969109027052/549876695789421`. This exact rational-vector audit checks the identities independently of the builder's symbolic names; it is not asserted to be a physical Xi state.

## Near And Paired Remote

```text
Q_near=(2/|nu|^2){Q(omega_E+p+n)-Q(omega_E+p)}.
Q_remote=(2/|nu|^2){Q(omega_E+p+n+c)-Q(omega_E+p+n)}. It contains the near/remote cross and remote self-determinant, so the existing linear remote budget cannot be reused for it without a new proof.
Q_geom=Q_near+Q_remote exactly for mathscr O_N=mathscr N_N+mathscr C_N, with mathscr C_N defined only as a common-cutoff paired-remote observation.
```

The near-first order assigns the near/remote cross to the remote increment, matching the certified rule that the finite mean-one near kernel is composed first and the paired remote tail is attached afterward. Reversing the order changes only cross ownership, not the total.

## Finite-Cell Form

Since mathscr O_N=-mathscr D_N, write d_X=Re{nu c_xi mathscr D_N[P_X]}; then R_geom=partial_xi{Q(b-d)-Q(b)}.

Because mathscr D_N=mathscr F_N-mathscr M_N, b-d=omega_E+2p-f, where f_X=Re{nu c_xi mathscr F_N[P_X]}. Thus the completed geometric current has a purely finite reciprocal-band representation.

With d_X=Re{nu c_xi mathscr D_N[P_X]} and mathscr O_N=-mathscr D_N, Q_geom=(2/|nu|^2){Q(omega_E+p-d)-Q(omega_E+p)}. Since mathscr D_N=-mathscr M_N+H_N{f_P(1)-f_P(N)}+V_N[P], this contains no conditional infinite series.

At either reciprocal boundary, transfer the complete mode and its lift between mathscr F_N and mathscr O_N. The full observation and Q_geom are tie invariant; differentiating m_N or n_N as floor functions is forbidden.

## Corrected Signed Criterion

Sigma_full=-2rho A_TB+E_geom^cur+E_move^cur. On a roster-interior continuous arc carrying the certified values B>49/100 and B<-49/100, the intermediate-value theorem supplies a calibrated point B=-49/100. At that point, |E_move^cur|<rho A_T/100 makes E_geom^cur<-(97/100)rho A_T necessary for negativity, while E_geom^cur<-(99/100)rho A_T is sufficient independently of the moving-current sign.

These 97/100 and 99/100 constants replace the prior signed-main bookkeeping only for the newly completed geometric package. The previously proved component budgets remain valid fallbacks but are not spent in this joined determinant formulation.

## Next Action

Choose between two equivalent arithmetic charts. In the finite chart, insert the Dirichlet-cell formula for mathscr D_N and derive the reciprocal-stationary decomposition of the determinant Q(omega_E+p-d)-Q(omega_E+p). In the near/remote chart, compose the mean-one finite near kernel with the retained carrier before bounding the paired-remote cross and self-determinant. Test both charts on the physical coefficient rows before selecting Type-I/II, Vaughan, adjacent-recurrence, or reciprocal-pair machinery.

## Pi Provenance

No new pi is introduced. The determinant recombination uses only J and linearity. Pi in the finite and paired-remote representations is inherited from e(x)=exp(2pi i x), kappa=1/(2pi i), and the finite Dirichlet kernel.

## Proof Boundary

This gate proves exact package recombination, edge completion, canonical near/paired-remote and finite-cell determinant forms, tie ownership, and geometric-current threshold arithmetic. It proves no sign or bound for E_geom, Q_near, Q_remote, the finite-cell determinant, or the complete current; no all-q transport, contact exclusion, Lambda<=0, PF-infinity, RH, or prize-level conclusion.

built Newman C1 complete outer-determinant recombination gate: 18 rows, 0 issues, 9 symbolic audits, 6 exact-rational audits, 11 source audits, 8 packages collapsed, 2 pointwise packages remain, 1 finite-cell determinant form, 1 open arithmetic row
