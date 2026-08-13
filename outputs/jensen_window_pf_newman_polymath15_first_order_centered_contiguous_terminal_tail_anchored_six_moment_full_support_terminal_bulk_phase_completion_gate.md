# Full-Support Terminal-Bulk Phase Completion Gate

Date: 2026-08-02

Status: exact terminal-bulk phase split, observation coefficient completion, and pure conditional null proved; signed completed estimate open; not a proof of RH.

## Terminal-Bulk Split

For q<N put d_q=-iu_q mathcal N_qw_q, so mathcal N_(<N,xi)=sum_(q<N)Re(d_q).

Im(w_N)Re(d_q)={Im(w_Nd_q)+Im(w_Nconjugate(d_q))}/2.

The live kernel is sum_(q<N)Re{K_H(N,q)w_Nconjugate(w_q)+K_T(N,q)w_Nw_q}, with K_H=S_(1,N)u_qconjugate(mathcal N_q)/(4pi) and K_T=-S_(1,N)u_qmathcal N_q/(4pi).

K_T(N,q)=-conjugate(K_H(N,q)); this relates coefficients but does not cancel the two distinct phase families.

Under w_q=omega_a A_q exp(-i xi u_q), w_Nconjugate(w_q)=|omega_a|^2A_NA_q exp(i xi(u_q-u_N)) is Hermitian phase-difference, while w_Nw_q=omega_a^2A_NA_q exp(-i xi(u_N+u_q)) is transpose phase-sum.

## Ideal Phase Chart

For mathcal N_q^0=-u_q^2/4-i u_(N,x)/2, K_H^0=-S_(1,N)u_q^3/(16pi)+iS_(1,N)u_qu_(N,x)/(8pi) and K_T^0=S_(1,N)u_q^3/(16pi)+iS_(1,N)u_qu_(N,x)/(8pi).

The ideal real parts are opposite and the imaginary parts agree; neither coefficient has a pointwise sign.

At u_q=1, u_(N,x)=0, and w_N=w_q=i, the joined pair equals -S_(1,N)/(8pi), not zero.

## Coefficient Completion

For C(epsilon)=2omega_p^TJdot(epsilon)+2dot(omega_p)^TJepsilon+2epsilon^TJdot(epsilon), split epsilon=a+rho into the terminal conditional observation and every nonterminal residual observation.

C(a+rho)=C(rho)+2(omega_p+rho)^TJdot a+2(dot(omega_p)+dot rho)^TJa+2a^TJdot a.

Put omega_*=omega_p+rho and dot(omega_*)=dot(omega_p)+dot rho. The terminal-dependent quadratic channel is Delta_a=2omega_*^TJdot a+2dot(omega_*)^TJa+2a^TJdot a.

Delta_a=V_*dot a_N+N_*dot a_V-A_*dot a_Q-Q_*dot a_A+dot V_*a_N+dot N_*a_V-dot A_*a_Q-dot Q_*a_A+a_Vdot a_N+a_Ndot a_V-a_Adot a_Q-a_Qdot a_A.

The mixed coefficient of a_V is the complete nonterminal derivative dot N_*=mathcal N_(p,xi)+rho_(N,xi); the pure a-by-dot a term remains separate.

The affine term t_T^Tdot(epsilon) splits linearly as t_T^Tdot rho+t_T^Tdot a and is not absorbed into C or counted twice.

## Pure Conditional Current

For the source-inclusive conditional amplitudes A_X=tau_NP_X(log N), tau_N=kappa S_(1,N)w_N, the frequency lift is A_(dot X)=-iu_NA_X.

H_aa=A_Vconjugate(A_(dot N))+A_Nconjugate(A_(dot V))-A_Aconjugate(A_(dot Q))-A_Qconjugate(A_(dot A)) has Re(H_aa)=0 exactly.

T_aa=-2iu_N(A_VA_N-A_AA_Q).

The carrier identities A_N=R_NC_N and mathcal N_N=R_NQ_N+R_(N,x)C_N give C_Nmathcal N_N-A_NQ_N=R_(N,x)C_N^2 exactly.

Hence T_aa=-2iu_N tau_N^2 R_(N,x)C_N^2, and the actual pure-terminal real current is R_aa=(1/2)Re(T_aa).

The physical terminal boxes give u_N<2h, |R_(N,x)|<h^2/20, |C_N|<3/2, and |kappa|^2=1/(4pi^2)<1/36 because pi>3.

The source-inclusive pure conditional quadratic obeys |R_aa|<h^3 S_(1,N)^2 |w_N|^2/160.

This is h^3 times a squared logarithmic kernel, whereas the mixed terminal-bulk channel is linear in S_(1,N) and contains the complete nonterminal derivative observation; no cancellation follows from scale alone.

## Handoff

Pure conditional self-interaction cannot cancel the live terminal-bulk channel: its Hermitian real part vanishes and its transpose remainder is h^3S_(1,N)^2 scale.

Expand rho into the lower S_1 endpoint, the S_2/S_3 endpoints, and the absolutely convergent interior; then form the completed Hermitian and transpose kernels for dot N_*a_V before taking moduli.

## Pi Provenance

The factor `1/(4pi)` is the polarization half of the fixed conditional factor `1/(2pi)`, itself forced by `kappa=1/(2pi i)` and `e(x)=exp(2pi i x)`. The estimate `|kappa|^2<1/36` uses only `pi>3`. No circle, polygon, plotted symmetry, or fitted constant is introduced.

## Proof Boundary

This gate proves the exact Hermitian phase-difference and transpose phase-sum decomposition of the live terminal-bulk quadrature, their anti-conjugate coefficient relation and explicit noncancellation witness, exact observation-level coefficient completion, exact vanishing of the pure terminal Hermitian real current, collapse of the pure terminal transpose current through `C_N mathcal N_N-A_NQ_N=R_(N,x)C_N^2`, and its `h^3 S_(1,N)^2` envelope. It proves no signed completed terminal/nonterminal estimate, complete outer-complement or quadratic residual bound, signed full current, `Phi_B` theorem, contact exclusion, retained aggregate or Xi theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion.
