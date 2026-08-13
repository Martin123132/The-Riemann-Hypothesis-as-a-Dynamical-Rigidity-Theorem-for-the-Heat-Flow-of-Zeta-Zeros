# Full-Support Outer Endpoint Kernel Gate

Date: 2026-08-02

Status: exact pole-free outer-complement and endpoint-linear kernel proved; signed full-support estimate open; not a proof of RH.

This is not a proof of RH. It replaces the split lower and upper outer tails by one stable endpoint package and then carries that package through the physical observation rows.

## Reciprocal Band

m_N=floor(alpha_P/(N+1/2))+1, n_N=floor(2alpha_P), and R_N={m_N,...,n_N}.

For x(u)=alpha_P/u, a_x=x-m_N+1>1/4 and b_x=n_N+1-x>alpha_P on 1<=u<=N.

The complement uses the same half-open reciprocal band as the finite-cell theorem; at a moving integer tie its transfer must be joined to the band before differentiation.

## Pole-Free Kernels

S_1^[N](x)=psi(b_x)-psi(a_x)=PV sum_(r notin R_N)1/(x-r).

S_2^[N](x)=zeta(2,a_x)+zeta(2,b_x)=sum_(r notin R_N)1/(x-r)^2.

S_3^[N](x)=zeta(3,a_x)-zeta(3,b_x)=sum_(r notin R_N)1/(x-r)^3.

At integer x in R_N, S_1^[N](x)=H_(n_N-x)-H_(x-m_N); no cotangent pole is separated.

On the physical segment S_2^[N]<21, |S_3^[N]|<73, and Z_4^[N]<278. These bound kernels, not the complete outer functional.

## Outer Recomposition

I_P(r)=kappa[e(phi_r)A_P/q_r]_1^N-kappa^2[e(phi_r)(D_rA_P)/q_r]_1^N+kappa^2 integral_1^N e(phi_r)C_rA_P du.

C_rA=A''/q_r^2-3A'q_r'/q_r^3-Aq_r''/q_r^3+3A(q_r')^2/q_r^4.

B_mu^[N][P]=e(alpha_Plog mu){kappa A_P(mu)S_1^[N](alpha_P/mu)-kappa^2[A_P'(mu)S_2^[N](alpha_P/mu)+(alpha_P/mu^2)A_P(mu)S_3^[N](alpha_P/mu)]}.

O_N[P]=B_N^[N][P]-B_1^[N][P]+kappa^2 sum_(r notin R_N)^sym integral_1^N e(phi_r)C_rA_P du.

Only S_1 is principal-value conditional; the S_2, S_3, and interior C_r terms are absolutely summable after the joined twofold integration by parts.

## Physical Pullthrough

The complete outer linear correction is Re[nu c_xi O_N[P_lin^[N]]], with no six-moment componentwise split.

P_lin^[N](log N)=F_(alpha_N)(log N)-i u_N F_(beta_N)(log N).

For beta_E=(mathcal N_E,V_E,-Q_E,-A_E), P_E(lambda)=i(lambda-log a)F_(beta_E)(lambda).

B_N^[N][P_E] has exact factor i exp(S(log N)) times {-kappa u_NF S_1-kappa^2 N^(-1)[F-u_N(F'+gF)]S_2+kappa^2 alpha_Pu_NN^(-2)F S_3}.

The conditional S_1 edge coefficient has the exact factor u_N<2h; the remaining edge coefficients contain either N^(-1)<=(3/2)h or alpha_Pu_N/N^2<9h/2.

The F_(alpha_N)(log N) coefficient survives without u_N. Thus terminal centering suppresses the edge-only trace but does not algebraically cancel the complete lower-boundary package.

## Quadratic Channels

R_out=(1/2)Re[|nu|^2 H_out+nu^2 c_xi^2 T_out], where H_out uses O_X conjugate(O_dotY) and T_out uses O_X O_dotY with signs +,+,-,-.

The common phase cancels exactly from H_out and appears squared in T_out.

The retained moving-tail defect remains 2Delta(omega_C)^TJdot(omega_B)+2(omega_E+omega_B+omega_C)^TJdot(omega_C); it is not an outer-mode remainder and is added exactly once.

## Handoff

The genuine edge no longer carries an unsuppressed conditional endpoint trace: its `S_1` coefficient contains `u_N`, and its remaining endpoint coefficients contain `1/N` or `alpha_P u_N/N^2`. The complete alpha-row term still survives. The next theorem must therefore estimate that term jointly with the returned reciprocal-band carrier, Hermitian phase-difference, transpose phase-sum, correction, and moved-tail kernels.

## Pi Provenance

`e(x)=exp(2pi i x)` fixes `kappa=1/(2pi i)`. The cotangent formula is used only as an audit of the symmetric integer partial fraction; the stable formula uses digamma and Hurwitz functions with positive arguments. No circle, polygon, or fitted visual constant enters.

## Proof Boundary

This gate proves a stable pole-free outer-mode representation, exact endpoint reassembly, the one-polynomial physical pullthrough, terminal-centering suppression of the edge-only trace, and exact Hermitian/transpose decomposition of the four quadratic outer pairings. It proves no complete outer-complement estimate, quadratic residual bound, signed full current, `Phi_B` theorem, contact exclusion, retained aggregate or Xi theorem, `Q209`, cofinal descendant theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion.
