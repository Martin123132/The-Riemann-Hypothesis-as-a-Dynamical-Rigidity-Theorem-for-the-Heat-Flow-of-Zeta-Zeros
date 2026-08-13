# Full-Support Terminal/Nonterminal Relative-Lift Gate

Date: 2026-08-02

Status: exact endpoint-resolved remainder, tie-safe completed lifts, and degree-five mixed-current compression proved; signed lower/interior estimate open; not a proof of RH.

## Outer Remainder

Write T_mu[P]=e(alpha_Plog mu)kappa A_P(mu)S_1^[N](alpha_P/mu) and U_mu[P]=-e(alpha_Plog mu)kappa^2{A_P'(mu)S_2^[N](alpha_P/mu)+(alpha_P/mu^2)A_P(mu)S_3^[N](alpha_P/mu)}.

Write I_N[P]=kappa^2 sum_(r notin mathcal R_N)^sym integral_1^N e(phi_r)C_rA_Pdu; this series is absolute only after the joined twofold integration by parts.

The exact outer functional is mathscr O_N[P]=T_N[P]+mathscr R_N[P], where mathscr R_N[P]=U_N[P]-T_1[P]-U_1[P]+I_N[P].

For X in {V,mathcal N,A,Q}, rho_X=Re{nu c_xi mathscr R_N[P_X]} and dot rho_X=Re{nu c_xi mathscr R_N[i(lambda-log a)P_X]} on a fixed roster cell.

The complete linear outer functional splits without componentwise budgets as mathscr O_N[P_lin]=T_N[P_lin]+mathscr R_N[P_lin].

## Relative Lifts

T_N[i(lambda-log a)P]=-iu_NT_N[P].

D_Hmathscr R_N[P]:=mathscr R_N[i(lambda-log a)P]+iu_Nmathscr R_N[P]=mathscr R_N[i(lambda-log N)P].

D_Tmathscr R_N[P]:=mathscr R_N[i(lambda-log a)P]-iu_Nmathscr R_N[P]=mathscr R_N[i(lambda-log a-u_N)P].

On a retained carrier atom with lift -iu_q, D_H=-i(u_q-u_N) and D_T=-i(u_q+u_N); the q=N Hermitian relative lift vanishes exactly.

At a reciprocal tie, transfer the complete mode I_P and its lift from mathscr R_N to the retained carrier. Their sum Y and both dot Y+iu_NY and dot Y-iu_NY are invariant; the complement is never differentiated alone.

## Mixed Current

Let Y_X be the retained carrier observation plus mathscr R_N[P_X], and let A_X=tau_0P_X(log N), tau_0=kappa e(alpha_Plog N)exp(S(log N))S_(1,N). Then the complete terminal/nonterminal current uses D_HY_X=dot Y_X+iu_NY_X and D_TY_X=dot Y_X-iu_NY_X.

Re H_(a*)=Re{A_Vconjugate(D_HY_N)+A_Nconjugate(D_HY_V)-A_Aconjugate(D_HY_Q)-A_Qconjugate(D_HY_A)}.

T_(a*)=A_VD_TY_N+A_ND_TY_V-A_AD_TY_Q-A_QD_TY_A.

The actual mixed real current is (1/2)Re{|nu|^2H_(a*)+nu^2c_xi^2T_(a*)}; the fixed affine row and moved-tail defect remain separate.

## Polynomial Compression

B_(H,N)(lambda)=conjugate(C_N)mathcal N(lambda)+conjugate(mathcal N_N)C(lambda)-conjugate(A_N)Q(lambda)-conjugate(Q_N)A(lambda).

B_(T,N)(lambda)=C_Nmathcal N(lambda)+mathcal N_NC(lambda)-A_NQ(lambda)-Q_NA(lambda).

The outer-remainder cross satisfies Re H_(a rho)=Re{conjugate(tau_0)mathscr R_N[i(lambda-log N)B_(H,N)(lambda)]}.

The outer-remainder cross satisfies T_(a rho)=tau_0mathscr R_N[i(lambda-log a-u_N)B_(T,N)(lambda)].

Both B_(H,N) and B_(T,N) have degree at most four; each relative lift has degree at most five, so the source-complete mixed cross remains inside the certified six-moment closure.

## Endpoint Audit

For P_H=i(lambda-log N)P, U_N[P_H]=-i kappa^2 e(alpha_Plog N)exp(S(log N))N^(-1)P(log N)S_2^[N](alpha_P/N); its terminal S_1 and S_3 values vanish.

For P_T=i(lambda-log a-u_N)P, U_N[P_T]=-i kappa^2 e(alpha_Plog N)exp(S(log N)){N^(-1)[P-2u_N(P'+gP)]S_2-2alpha_Pu_NN^(-2)P S_3} at lambda=log N.

The lower conditional piece is -T_1[P_H]=i kappa(log N)P(0)S_1^[N](alpha_P).

The lower conditional piece is -T_1[P_T]=i kappa(log a+u_N)P(0)S_1^[N](alpha_P).

B_(T,N)(lambda)=[R(lambda)-R_N]{C_NQ(lambda)-Q_NC(lambda)}+[R_x(lambda)+R_(N,x)]C_NC(lambda).

At lambda=log N, B_(T,N)=2R_(N,x)C_N^2, while the Hermitian relative factor lambda-log N is zero.

The retained q=N cross is therefore Hermitian-null and transpose-suppressed; it does not cancel the lower endpoint or interior algebraically.

At mu=N, every relative-lift S_2/S_3 coefficient contains N^(-1), u_N, or alpha_Pu_N/N^2; the established boxes give N^(-1)<=(3/2)h and alpha_Pu_N/N^2<(9/2)h.

At mu=1, 0<S_1^[N](alpha_P)<5/N<=(15/2)h, but the relative lifts also contain log N or log a+u_N. The lower endpoint must remain grouped with the interior.

## Handoff

Terminal higher-endpoint and retained q=N cross terms keep explicit h mechanisms. The unresolved completed channel is the lower endpoint plus absolute interior, with both relative phase families still joined.

Expand the two degree-five mathscr R_N functionals into their lower U_1 and absolute-interior pieces, combine them with the finite q<N carrier sums and correction perturbation, and prove or falsify a joined Abel/Mangoldt/Vaughan or reciprocal-pair estimate before taking moduli.

## Pi Provenance

The constants `kappa=1/(2pi i)` and `alpha_P=xi/(2pi)` are inherited from the fixed Fourier character. Relative lifting introduces only logarithmic differences and sums. No circle, polygon, plotted symmetry, or fitted constant is introduced.

## Proof Boundary

This gate proves the exact four-piece outer remainder, its real observation projection, fixed-cell relative lifts, half-open tie-transfer invariance of the completed observation, exact Hermitian and transpose mixed-current compression, two degree-five polynomial functionals, terminal higher-endpoint reductions, lower conditional traces, and the transpose carrier-bilinear factorization. It proves no signed lower-endpoint/interior estimate, complete outer-complement or quadratic residual bound, signed full current, `Phi_B` theorem, contact exclusion, retained aggregate or Xi theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion.
