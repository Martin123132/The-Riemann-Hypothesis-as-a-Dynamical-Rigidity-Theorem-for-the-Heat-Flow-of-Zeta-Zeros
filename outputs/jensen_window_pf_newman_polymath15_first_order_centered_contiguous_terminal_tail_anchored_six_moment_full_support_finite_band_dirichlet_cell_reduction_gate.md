# Finite-Band Dirichlet-Cell Reduction Gate

Date: 2026-08-02

Status: exact finite Dirichlet kernel, cell variation, relative-channel and tie reductions proved; signed physical variation estimate open; not a proof of RH.

## Band Kernel

Put m_N=floor(alpha_P/(N+1/2))+1 and n_N=floor(2alpha_P). The reciprocal band is the finite integer interval {m_N,...,n_N}, with m_N>=1 on the physical chart.

Define K_N(u)=sum_(r=m_N)^n_N e(-ru). Then mathscr F_N[P]=integral_1^N f_P(u)K_N(u)du, where f_P(u)=e(alpha_Plog u)exp(S(log u))P(log u).

For u notin Z, K_N(u)=e(-(m_N+n_N)u/2)sin(pi(n_N-m_N+1)u)/sin(pi u); at u in Z its removable value is n_N-m_N+1.

Because all frequencies are integral, K_N is 1-periodic, K_N(-v)=conjugate(K_N(v)), and integral_(-1/2)^(1/2)K_N(v)dv=0 since m_N>=1.

H_N:=integral_0^(1/2)K_N(v)dv=2kappa sum_(m_N<=r<=n_N, r odd)1/r is purely imaginary, and integral_(-1/2)^0K_N(v)dv=-H_N.

## Cell Reduction

Partition [1,N] into [1,3/2], the centered cells [q-1/2,q+1/2] for 2<=q<=N-1, and [N-1/2,N]. Periodicity lets every cell use the same K_N(v).

Define V_N[P] as the lower half-cell integral of [f_P(1+v)-f_P(1)]K_N(v), the sum of centered integrals of [f_P(q+v)-f_P(q)]K_N(v), and the upper half-cell integral of [f_P(N+v)-f_P(N)]K_N(v).

The zero cell mean and half-kernel identity give mathscr F_N[P]=H_N{f_P(1)-f_P(N)}+V_N[P].

Therefore mathscr D_N[P]=-mathscr M_N[P]+H_N{f_P(1)-f_P(N)}+V_N[P]. This finite formula is exactly tau_band-E_band with no conditional exterior sum left.

The returned-carrier join is mathscr M_N[P]-mathscr D_N[P]=2mathscr M_N[P]-H_N{f_P(1)-f_P(N)}-V_N[P].

On an interior cell, V_(N,q)[P]=integral_0^(1/2){Re(K_N)Delta_even f_P+i Im(K_N)Delta_odd f_P}dv, where Delta_even=f(q+v)+f(q-v)-2f(q) and Delta_odd=f(q+v)-f(q-v).

## Relative Channels

For P_H=i(lambda-log N)B_(H,N), f_(P_H)(N)=0 and mathcal T_N[P_H]=0. Hence mathscr R_N[P_H]=mathscr M_N[P_H]-H_Nf_(P_H)(1)-V_N[P_H], and the carrier-only join is 2mathscr M_N[P_H]-H_Nf_(P_H)(1)-V_N[P_H].

For P_T=i(lambda-log a-u_N)B_(T,N), mathscr R_N[P_T]=mathscr M_N[P_T]-H_N{f_(P_T)(1)-f_(P_T)(N)}-V_N[P_T]+2iu_Nmathcal T_N[B_(T,N)]. The carrier-only join adds a second mathscr M_N[P_T].

The transpose endpoint satisfies f_(P_T)(N)=-2iu_N f_(B_(T,N))(N). It must remain with the H_N endpoint coefficient, the conditional terminal term, and the pure terminal transpose block.

Insert the intact ideal cubics P_H^0=-i(x+u_N)(x-u_N)^2/4 and P_T^0=-i(x-u_N)(x+u_N)^2/4+u_(N,x)(x-u_N) into mathscr M_N,H_Nf(1),V_N before attaching Delta_H,Delta_T.

## Tie Transport

When 2alpha_P crosses an integer r_*, n_N gains r_* and mathscr F_N,mathscr D_N jump by +I_P(r_*); the outer remainder jumps by -I_P(r_*).

When alpha_P/(N+1/2) crosses an integer r_*, m_N changes from r_* to r_*+1 and mathscr F_N,mathscr D_N jump by -I_P(r_*); the outer remainder jumps by +I_P(r_*).

On a fixed roster cell, if P and S are held fixed under the auxiliary frequency, f_(i(lambda-c)P)=(partial_xi-ic)f_P and mathscr D_N[i(lambda-c)P]=(partial_xi-ic)mathscr D_N[P]. At a tie the complete transferred mode and lift must be included.

Differentiating the sine-quotient formula while allowing m_N or n_N to move creates boundary distributions. The half-open jump laws above, not an ordinary derivative of floor functions, are the exact global rule.

## Reciprocal Blocks

Let J_q be the nearest-integer physical u-cell, truncated to [1,N], and let C_p={r in Z_(>0):p-1/2<=alpha_P/r<p+1/2}. Then mathcal R_N is the disjoint union of C_p for 1<=p<=N.

Define B_(q,p)[P]=sum_(r in C_p)integral_(J_q)f_P(u)e(-ru)du. The finite band sum is exactly mathscr F_N[P]=sum_(q=1)^N sum_(p=1)^N B_(q,p)[P].

Writing v=alpha_P/r in the p-cell, the phase gradient on J_q is alpha_P/u-r=alpha_P(v-u)/(uv). The diagonal p=q contains the reciprocal stationary point u=v; p=1,N are the two incomplete endpoint diagonals.

For |p-q|>=2, |alpha_P/u-r|>=alpha_P(|p-q|-1)/{(p+1/2)(q+1/2)} throughout the block. These separated blocks are uniformly nonstationary before any modulus is taken.

The p=q+1 and p=q-1 blocks meet the diagonal cells at a shared half-integer boundary and have no positive uniform gradient gap. They must be composed in opposite orientations with the half-open tie transfer and, at the terminal edge, the complex adjacent recurrence.

For u=q+w, phi_r(q+w)-phi_r(q)=(alpha_P/q-r)w+alpha_P{log(1+w/q)-w/q}. This exact form keeps the linear reciprocal offset and negative logarithmic curvature together.

## Route Guards

For the exact test source f(u)=e(ku) with integer k in [m_N,n_N], orthogonality gives mathscr F_N=N-1=mathscr M_N and mathscr D_N=0.

For f(u)=e(ku) with integer k outside [m_N,n_N], including k=0, mathscr F_N=0 while mathscr M_N=N-1, so mathscr D_N=-(N-1) and mathscr M_N-mathscr D_N=2(N-1).

The finite-cell defect has no coefficient-blind smallness or sign theorem. Its behavior is controlled by where the source frequency lies relative to the reciprocal band; the logarithmic physical phase must be retained.

For f_P=e(alpha_Plog u)exp(S(log u))P(log u), f_P'(u)=u^(-1)e(alpha_Plog u)exp(S){P'+[g+2pi i alpha_P]P}. The alpha_P term is large, so V_N cannot be bounded as a slowly varying cell error before reciprocal stationary cancellation.

Return to the reciprocal r-cells inside the finite V_N formula: isolate the stationary r near alpha_P/q, pair the nonstationary modes globally, and estimate the ideal cubic carrier and variation together. Do not take an L1 norm of K_N or a cellwise derivative majorant.

## Handoff

Replace the conditional global exterior by the finite Dirichlet-cell functional. The live Hermitian and transpose joins are now two starred carriers minus one endpoint harmonic and one source-coupled cell variation, plus the explicit transpose terminal term.

Derive a reciprocal-stationary decomposition of V_N[P_H^0] and V_N[P_T^0] that preserves cancellation with 2mathscr M_N. Bound or sign that ideal pair before attaching Delta_H,Delta_T and the remaining current packages.

## Pi Provenance

The sine quotient is only the finite geometric sum for e(x)=exp(2pi i x); kappa=1/(2pi i). No geometric input is added.

## Proof Boundary

No signed physical defect, h^2/800 reserve, completed current, or RH-level conclusion is proved. This gate proves no signed finite-cell defect or completed-current bound, `Phi_B` theorem, contact exclusion, retained aggregate or Xi theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion.
