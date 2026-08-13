# Physical P_lin And Morse-Amplitude Gate

Date: 2026-08-02

Status: exact physical coefficient factorization and Morse `c'` amplitude;
`0 grouped interior bounds`, `0 signed flow bounds`, and this is not a proof
of RH.

## Terminal-Centered Carrier

Put x=lambda-log(a) and u_N=log(a/N), so log(N)=log(a)-u_N. Writing s_*'=c+ib and s_*''=c_x+ib_x gives the exact identities R=-s_*'x-c*u_N and R_x=-s_*''x-c_x*u_N+i*b*u_(N,x).

Put delta=chi_N-i*b*u_N. Then Q=(R+delta)C+D and N=R*Q+R_x*C. The existing physical estimate is |delta|<4h^2, but no estimate is used in this algebraic gate.

In the centered variable, C=c_0+c_1*x+c_2*x^2 and D=d_0+d_1*x+d_2*x^2, where c_0=1+rho_1*log(a)+rho_2*log(a)^2, c_1=rho_1+2rho_2*log(a), c_2=rho_2, d_0=rho_(1,x)*log(a)+rho_(2,x)*log(a)^2, d_1=rho_(1,x)+2rho_(2,x)*log(a), and d_2=rho_(2,x).

## Two Correction Channels

For any real f=(f_V,f_N,f_A,f_Q), let F_f=f_V*C+f_N*N+f_A*A+f_Q*Q. Use alpha=(N_(p,xi),V_(p,xi),-Q_(p,xi),-A_(p,xi)) and beta=(N_p+A_(T,x),X_T+V_p,-Q_p-X_(T,x),-A_T-A_p). Then P_lin=F_alpha+i*x*F_beta.

For R=r_0+r_1*x and R_x=t_0+t_1*x, every F_f factors exactly as F_f=C*U_f+D*V_f. Here V_f=f_N*R+f_Q and U_f=f_N*[R(R+delta)+R_x]+f_Q(R+delta)+f_A*R+f_V. Thus four physical carrier polynomials reduce to the two correction channels C and D.

Writing U_f=u_0+u_1*x+u_2*x^2 and V_f=v_0+v_1*x, put g_0=r_0(r_0+delta)+t_0, g_1=r_1(2r_0+delta)+t_1, g_2=r_1^2. Then u_0=f_Ng_0+f_Q(r_0+delta)+f_Ar_0+f_V, u_1=f_Ng_1+(f_Q+f_A)r_1, u_2=f_Ng_2, v_0=f_Nr_0+f_Q, and v_1=f_Nr_1.

Let gamma=(u_(alpha,0),u_(alpha,1)+i*u_(beta,0),u_(alpha,2)+i*u_(beta,1),i*u_(beta,2)) and eta=(v_(alpha,0),v_(alpha,1)+i*v_(beta,0),i*v_(beta,1)). Then P_lin=C*sum_(j=0)^3gamma_j*x^j+D*sum_(j=0)^2eta_j*x^j.

## Six Complex Coefficients

```text
p_0=c_0*gamma_0+d_0*eta_0
p_1=c_0*gamma_1+c_1*gamma_0+d_0*eta_1+d_1*eta_0
p_2=c_0*gamma_2+c_1*gamma_1+c_2*gamma_0+d_0*eta_2+d_1*eta_1+d_2*eta_0
p_3=c_0*gamma_3+c_1*gamma_2+c_2*gamma_1+d_1*eta_2+d_2*eta_1
p_4=c_1*gamma_3+c_2*gamma_2+d_2*eta_2
p_5=c_2*gamma_3
```

The top coefficient is exactly p_5=i*rho_2*(X_T+V_p)*(s_*')^2. Hence the degree-five channel vanishes on the total-value fibre X_T+V_p=0; this is a degree reduction, not a remainder bound.

## Correction-Free Audit

At C=1, D=0, s_*'=-i/2, s_*''=0, delta=0, and b=-1/2, one has R=Q=i*x/2 and N=-x^2/4-i*u_(N,x)/2. The induced observation identities are A_p=Q_p, A_(p,xi)=Q_(p,xi), and V_(p,xi)=2A_p.

```text
P_lin=N_(p,xi)-i*A_p*u_(N,x)+{i[N_p+A_(T,x)-Q_(p,xi)]+(X_T+V_p)u_(N,x)/2}x+[A_p+A_T+X_(T,x)]x^2/2-i(X_T+V_p)x^3/4.
```

```text
On X_T+V_p=0 and A_T+A_p=0, the correction-free leading polynomial is quadratic: P_lin=N_(p,xi)-i*A_p*u_(N,x)+i[N_p+A_(T,x)-Q_(p,xi)]x+X_(T,x)x^2/2.
```

The actual correction channels can restore degrees four and five. The ideal reduction is a scale guide, not an exact deletion of C-1 or D.

## Exact Morse Interior

For u_0=alpha/r, lambda=log(u_0v), S(lambda)=t*lambda^2/4-sigma*lambda, z=sgn(v-1)sqrt(2[v-1-log v]), and J=dv/dz, put W_P(v)=exp(S(lambda))*P(lambda)*J(v). Then b_(P,r)(y)=sqrt(alpha)W_P(v)/r and y=sqrt(alpha)z.

For z!=0, c_(P,r)(y)=[W_P(v)-W_P(1)]/(r*z).

For z!=0, c'_(P,r)(y)={z*J*W_P'(v)-W_P(v)+W_P(1)}/[r*sqrt(alpha)*z^2]. With g=t*lambda/2-sigma, this equals {exp(S)[zJ^2(P'+gP)/v+(zJJ'-J)P]+exp(S(lambda_0))P(lambda_0)}/[r*sqrt(alpha)*z^2], where lambda_0=log(alpha/r).

The removable saddle value is c'_(P,r)(0)=exp(S(lambda_0))/[2r*sqrt(alpha)] times {P''+(2g_0+1)P'+[g_0^2+g_0+t/2+1/6]P}(lambda_0), with g_0=t*lambda_0/2-sigma.

At fixed xi the same physical P_lin is used for every Poisson mode. Mode dependence enters only through lambda_0=log(alpha/r), the Morse endpoints, and the common phase. Thus the saddle audit is one second-order differential operator sampled on the logarithmic mode lattice, not six unrelated moment amplitudes.

Neither P_lin(lambda_0)=0 nor the displayed saddle operator vanishes identically. Any zero or discrete-difference gain must be proved from the physical coefficients and mode coupling.

## Route Decision

Use x=lambda-log(a) and the C/D factorization as the proof-facing coefficient chart. Do not expand powers of log(N) and log(a) separately, and do not return to six componentwise Morse estimates.

## Next Target

Insert the certified formulas for rho_1,rho_2 and their x derivatives, bound the six centered coefficients on the full physical chart, and apply the saddle differential operator to P_lin. Then test Abel or reciprocal pairing on that one mode sequence before taking absolute values; keep the four quadratic observation pairings separate.

## Proof Boundary

This proves the exact physical terminal-centered carrier laws, the C/D two-channel factorization, all six centered P_lin coefficients, the total-value factor in the degree-five coefficient, the correction-free cubic and contact-quadratic reductions, and exact off-saddle and saddle Morse c-prime formulas. It proves no coefficient norm on the physical chart, grouped c-prime sum, discrete cancellation theorem, quadratic residual bound, h^2 flow error, Phi_B upper bound, contact exclusion, retained aggregate or Xi theorem, Q209, Lambda<=0, PF-infinity, RH, or a prize-level conclusion.

## Reproduce

```text
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_amplitude_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_amplitude_gate.py
```
