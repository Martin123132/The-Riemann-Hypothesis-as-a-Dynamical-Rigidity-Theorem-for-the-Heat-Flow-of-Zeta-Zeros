# Full-Support Terminal Conditional Quadrature Gate

Date: 2026-08-02

Status: exact terminal quadrature, carrier-half audit, and recurrence defect proved; signed bulk-terminal estimate open; not a proof of RH.

## Physical Carrier Row

For u_q=log(a/q), R_q=s_*'u_q-cu_N and R_(q,x)=s_*''u_q-c_xu_N+ib u_(N,x). Then Q_q=(R_q+delta)C_q+D_q and mathcal N_q=R_qQ_q+R_(q,x)C_q.

With w_q=nu c_xi exp(S(log q))e(alpha_Plog q), mathcal N_(p,xi)=sum_(q=1)^N Re{-iu_q mathcal N_qw_q} after endpoint composition.

Before endpoint composition the reciprocal return and hard endpoint each contribute one half of Re{-iu_Nmathcal N_Nw_N}; together they give the full q=N atom exactly once.

The full q=N lifted mathcal N atom is smaller than (17/5)h^3|w_N|; each reciprocal/hard half is smaller than (17/10)h^3|w_N|.

## Terminal Kernel

For x_N=alpha_P/N, a_*=x_N-floor(alpha_P/(N+1/2)) and b_*=floor(2alpha_P)+1-x_N satisfy 1/4<a_*<2 and b_*>alpha_P(2-1/N)>a_*.

log(b_*/a_*)<S_(1,N)<log(b_*/a_*)+a_*^(-1)-b_*^(-1).

log(((N^2-1)(2-1/N))/2)<S_(1,N)<log(8(N+1)^2+4)+4.

The terminal S_1 kernel is positive and logarithmic in N; it is not an O(1) or O(h) factor.

## Conditional Boundary

If S_(1,N)=S_1^[N](alpha_P/N), the terminal conditional linear trace is L_(N,1)=[S_(1,N)/(2pi)]Im{w_NP_lin^[N](log N)}.

L_(N,1)=[S_(1,N)/(2pi)]mathcal N_(p,xi)Im(w_N)+E_row.

With R_beta the l1 norm of beta_N, |P_lin^[N](log N)-mathcal N_(p,xi)|<(8/5)hR_alpha+(16/5)hR_beta.

|E_row|<[S_(1,N)/(2pi)]h|w_N|{(8/5)R_alpha+(16/5)R_beta}.

Writing mathcal N_(<N,xi)=sum_(q<N)Re{-iu_qmathcal N_qw_q}, L_(N,1)=[S_(1,N)/(2pi)]Im(w_N)mathcal N_(<N,xi)+E_term.

|E_term|<[S_(1,N)/(2pi)]{h|w_N|[(8/5)R_alpha+(16/5)R_beta]+(17/5)h^3|w_N|^2}.

## Recurrence Audit

For Q_RS(p), r(p), and R(p) of Section 11.154, Q_RS(p)r(p)^m=(-1)^m conjugate(R(p-2m-2)).

C_0(p)-Q_RS(p)sum_(m=0)^M r(p)^m=(-1)^(M+1)C_0(p-2M-2)+2i sum_(m=0)^M(-1)^m Im R(p-2m-2).

At p=0 and M=0 the explicit quadrature-defect summand is 2i sin(3pi/8), so it is nonzero.

The known recurrence telescopes the real terminal trace, whereas the conditional S_1 boundary kernel sees the imaginary quadrature because kappa=-i/(2pi).

At the physical anchor, C_Nw_N=-Q_RS(p)exp(W_theta) and |exp(W_theta)/C_N-1|<9h.

The terminal quadrature obeys Im w_N(p=0)>5/8 and Im w_N(p=1)<-1/4.

## Handoff

The first unsuppressed boundary object is [S_(1,N)/(2pi)]Im(w_N)mathcal N_(<N,xi), not a local q=N cancellation.

Estimate this one-versus-many bulk-terminal correlation jointly with the remaining S_2/S_3 and interior outer terms, H_out, T_out, correction perturbation, and moved-tail defect.

The next calculation must estimate the displayed one-versus-many bulk-terminal correlation before adding the remaining endpoint, interior, Hermitian, transpose, correction, and moved-tail channels. The real terminal recurrence cannot be substituted for the required quadrature identity.

## Pi Provenance

The factor `1/(2pi)` is forced by `kappa=1/(2pi i)` under the fixed Fourier character `e(x)=exp(2pi i x)`. The phases `Q_RS` and `R` are inherited from the completed-zeta/Riemann--Siegel endpoint recurrence. No new geometric `pi` is introduced.

## Proof Boundary

This gate proves the exact all-carrier source expansion of `mathcal N_(p,xi)`, reassembly and h-cubed size of the two q=N derivative halves, positivity and logarithmic scale of the terminal `S_1` kernel, the exact conditional quadrature projection, a row-relative bulk-terminal reduction, the complex conjugate recurrence relation, its explicit imaginary defect, and two terminal-factor sign witnesses. It proves no signed bulk-terminal correlation estimate, complete outer-complement or quadratic residual bound, signed current, contact exclusion, retained aggregate or Xi theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion.
