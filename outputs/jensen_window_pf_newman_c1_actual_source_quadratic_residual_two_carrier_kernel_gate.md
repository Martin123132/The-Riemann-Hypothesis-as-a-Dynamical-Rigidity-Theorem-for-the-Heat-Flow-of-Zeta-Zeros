# Newman C1 Quadratic-Residual Two-Carrier Kernel Gate

Date: 2026-08-05

Status: exact residual-residual kernel closure; this is not a proof of the joint signed arithmetic estimate or RH.

## Functional And Kernels

Work on one fixed reciprocal-roster cell and the fixed physical coefficient chart. The terminal conditional mathcal T_N block has already been removed, and mathscr R_N is retained as one linear nonterminal functional.

Put x_lambda=lambda-log a and mathcal L_xi[P]=(nu c_xi/|nu|)mathscr R_N[P]. For X in {V,mathcal N,A,Q}, U_X=mathcal L_xi[P_X], e_X=|nu|Re(U_X), and partial_xi mathcal L_xi[P]=mathcal L_xi[i x_lambda P].

```text
K_T(lambda,mu)=C_lambda mathcalN_mu+mathcalN_lambda C_mu-A_lambda Q_mu-Q_lambda A_mu=(R_mu-R_lambda)(C_lambda Q_mu-C_mu Q_lambda)+(R_(lambda,x)+R_(mu,x))C_lambda C_mu.

C_lambda Q_mu-C_mu Q_lambda=C_lambda C_mu[(R_mu-R_lambda)+(delta_mu-delta_lambda)]+C_lambda D_mu-C_muD_lambda.

K_H(lambda,mu)=C_lambda conjugate(mathcalN_mu)+mathcalN_lambda conjugate(C_mu)-A_lambda conjugate(Q_mu)-Q_lambda conjugate(A_mu)=(conjugate(R_mu)-R_lambda)[C_lambda conjugate(Q_mu)-Q_lambda conjugate(C_mu)]+[R_(lambda,x)+conjugate(R_(mu,x))]C_lambda conjugate(C_mu).
```

## Primitive And Current

```text
Q_quad=(1/2)Re{(mathcal L_xi tensor conjugate(mathcal L_xi))[K_H]+(mathcal L_xi tensor mathcal L_xi)[K_T]}.
E_quad^cur=(1/2)Re{(mathcal L_xi tensor conjugate(mathcal L_xi))[i(lambda-mu)K_H]+(mathcal L_xi tensor mathcal L_xi)[i(lambda+mu-2log a)K_T]}.
```

The factor 1/2 is forced by 2Re(z)Re(w)=Re(z conjugate(w)+zw). It is not an adjustable convention.

The Hermitian current multiplier i(lambda-mu) vanishes identically on lambda=mu. The transpose multiplier becomes 2i(lambda-log a), so its diagonal remains unless a separate source identity kills it.

## Exact Kernel Match

```text
For C=1, D=delta=0, R=i*x/2, R_x=-i*u_x/2, one has K_H^0=-(x_lambda+x_mu)^2/4 and K_T^0=-(x_mu-x_lambda)^2/4-i*u_x.
```

Because x_lambda=-u_lambda on a retained physical atom, K_H^0/2=-(u_lambda+u_mu)^2/8 and K_T^0/2=-(u_lambda-u_mu)^2/8-i*u_x/2, exactly the kernels already certified in Section 11.160.

Hermitian products use one weight and one conjugate weight and therefore carry only phase differences. Transpose products use two unconjugated weights and carry only phase sums. The residual-residual determinant creates no third phase family.

The independent finite audit used 3 exact rational complex carriers; both functional differences and the Hermitian current diagonal were `0`, `0`, and `0`, while the transpose diagonal was nonzero.

## Arithmetic Handoff

Keep mathscr R_N[P]=mathcal U_N[P]-mathcal T_1[P]-mathcal U_1[P]+mathcal I_N[P] composed until K_H and K_T have been formed. Expanding its tensor square first creates sixteen coordinate terms and destroys the certified tie-invariant cancellations; those terms are not new ownership packages.

The pointwise route must insert the displayed E_quad^cur kernel into E_main^cur=E_FAN+E_quad^cur before any Hermitian/transpose arithmetic decomposition. A separate absolute |E_quad| allowance is neither supplied nor required by this closure.

Prove or falsify, on every adverse actual-source witness cell, the joint signed inequality E_FAN+E_quad^cur<-(1130179/858000)rho A_T, using the two existing phase families and preserving the grouped mathscr R_N tensor square; alternatively return to the exact endpoint determinant transfer.

No new pi is introduced. The u_x=h^2/(8pi) term in the ideal transpose kernel is inherited from the completed-zeta/Riemann--Siegel saddle normalization already audited in Section 11.160.

This proves exact kernel closure, factorization, normalization, phase multipliers, diagonal classification, and ideal-kernel agreement. It proves no signed Type-I/II, Vaughan, reciprocal-pair, or endpoint determinant estimate; no complete-current sign, contact exclusion, Lambda<=0, PF-infinity, RH, or prize-level conclusion.

built Newman C1 quadratic residual two-carrier kernel gate: 16 rows, 0 issues, 8 symbolic audits, 4 finite-functional audits, 4 source audits, 2 exact kernel families, 1 open signed theorem row
