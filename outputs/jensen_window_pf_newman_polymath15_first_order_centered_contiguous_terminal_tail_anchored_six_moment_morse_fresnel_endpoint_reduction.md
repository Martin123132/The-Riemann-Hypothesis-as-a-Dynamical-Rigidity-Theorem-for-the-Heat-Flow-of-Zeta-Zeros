# Six-Moment Logarithmic Morse-Fresnel Endpoint Reduction

Date: 2026-08-01

Status: exact endpoint and nonstationary-tail reduction with explicit
constant-one residual formulas, `0 evaluated physical remainder constants`,
`0 signed flow bounds`, and this is not a proof of RH.

Primary analytic source used only for comparison and motivation:

```text
https://arxiv.org/abs/1205.0090
```

## Domain And Source Boundary

Use the physical q=1, L>=50, fixed-N chart, B=N-M-1, alpha=xi/(2*pi), and alpha_-<=alpha<=alpha_+ corresponding to T_0-epsilon<=xi<=T_0. The six amplitudes are A_j(u)=(log u)^j exp[t(log u)^2/4-sigma log u], 0<=j<=5.

Vandehey Section 11 motivates an incomplete-Fresnel endpoint chart and explicitly notes that it is computer-evaluable near transition, but its displayed remainder and Corollary 1.4 retain implicit constants. The formulas below are instead exact direct specializations to the logarithmic phase; no sourced big-O term is promoted into a numerical bound.

## Exact Logarithmic Morse Coordinate

For r>0 put u_0=alpha/r, v=u/u_0, h(v)=v-1-log v, and phi_r(u)=alpha log u-ru. Then phi_r(u)=phi_r(u_0)-alpha h(v), where phi_r(u_0)=alpha[log(alpha/r)-1].

Define z(v)=sgn(v-1)sqrt(2h(v)), z(1)=0. Since h>=0, h'= (v-1)/v, z'=(v-1)/(vz)>0, and z tends to -infinity and +infinity at the two ends, z is a smooth increasing bijection (0,infinity)->R.

With J(v)=dv/dz=vz/(v-1) and J(1)=1, put y=sqrt(alpha)z. Then du/dy=[sqrt(alpha)/r]J(v) and phi_r(u)=phi_r(u_0)-y^2/2 exactly.

For delta=v-1, z=delta-delta^2/3+7delta^3/36+O(delta^4) and J=1+2delta/3-5delta^2/36+O(delta^3). These expansions supply the removable values at v=1; they are not asymptotic approximations in the exact transformed integral.

The transformed endpoints are y_1=sqrt(alpha)z(r/alpha) and y_B=sqrt(alpha)z(Br/alpha). The lower saddle crossing r=alpha has y_1=0; the upper crossing r=alpha/B has y_B=0.

## Fixed Transition Roster

Set m=max(1,floor(alpha_-/(2B))), n=ceil(2alpha_+), and T={m,...,n}. This roster is fixed on the whole auxiliary segment. It contains every mode whose saddle can enter [1,B], as well as both endpoint-transition collars.

## Incomplete-Fresnel Transition

Let v(y) be the inverse Morse coordinate and define b_(j,r)(y)=[sqrt(alpha)/r]A_j(u_0v(y))J(v(y)). Then b_(j,r)(0)=[sqrt(alpha)/r]A_j(u_0).

For every r in the fixed roster, I_(j,r)=e(phi_r(u_0))*integral_(y_1)^(y_B) b_(j,r)(y)e(-y^2/2)dy. This identity remains valid when u_0 is just outside [1,B].

Define F(y_1,y_B)=integral_(y_1)^(y_B)e(-y^2/2)dy and U_(j,r)=b_(j,r)(0)e(phi_r(u_0))F(y_1,y_B). Then I_(j,r)=U_(j,r)+Q_(j,r), where Q is the same integral with b(y)-b(0).

F(-infinity,infinity)=e(-1/8)=exp(-i*pi/4), while each half-line integral ending at zero is e(-1/8)/2. Thus the usual full stationary main and both starred half-mains are limits of one continuous function, not separate cutoff rules.

Put c(y)=[b(y)-b(0)]/y for y!=0 and c(0)=b'(0). Exact integration by parts gives |Q_(j,r)|<={|c(y_1)|+|c(y_B)|+integral_(y_1)^(y_B)|c'(y)|dy}/(2*pi). The constant is exactly 1/(2*pi), with no hidden big-O constant.

The identities c(y)=integral_0^1 b'(sy)ds and c'(y)=integral_0^1 s b''(sy)ds make the residual bound regular at y=0. Moreover c_(j,r)(0)=r^(-1){jA_(j-1)(u_0)+[2/3-sigma+(t/2)log u_0]A_j(u_0)}.

Because the roster is fixed and z, J, F, and c use removable charts, U+Q is smooth through y_1=0 and y_B=0. The two sharp jump laws of Section 11.167 are replaced by one exact incomplete-Fresnel transition; no cutoff absence is assumed.

## Reassembled Nonstationary Tail

For r outside T, q_r(u)=alpha/u-r has no zero on [1,B]. The low tail lies below alpha_-/(2B), and the high tail lies above 2alpha_+, so both are uniformly separated from every stationary point on the full xi segment.

Let kappa=1/(2*pi*i), D_rA=(A/q_r)', and C_rA=((D_rA)/q_r)'. Then exactly I_r=kappa[e(phi_r)A/q_r]_1^B-kappa^2[e(phi_r)(D_rA)/q_r]_1^B+kappa^2 integral_1^B e(phi_r)C_rA du.

Writing q'= -alpha/u^2 and q''=2alpha/u^3, C_rA=A''/q^2-3A'q'/q^3-Aq''/q^3+3A(q')^2/q^4.

For x in (m-1,n+1), a_x=x-m+1 and b_x=n+1-x, the symmetric outside-roster sums are S_1(x)=psi(b_x)-psi(a_x), S_2(x)=zeta(2,a_x)+zeta(2,b_x), and S_3(x)=zeta(3,a_x)-zeta(3,b_x). Equivalently S_1=pi*cot(pi*x)-sum_(r=m)^n 1/(x-r), with removable values at every integer in T.

At an integer endpoint mu in {1,B}, x=alpha/mu and e(-rmu)=1. The complete outside-roster endpoint functional is B_j(mu)=e(alpha log mu){kappa A_j(mu)S_1(x)-kappa^2[A_j'(mu)S_2(x)+(alpha/mu^2)A_j(mu)S_3(x)]}. The summed tail equals B_j(B)-B_j(1)+kappa^2 sum_(r outside T)integral e C_rA_j.

For k>=2 put Z_k(x)=zeta(k,a_x)+zeta(k,b_x). The remaining tail has the explicit bound 1/(4*pi^2) integral_1^B {|A_j''|Z_2+[3|A_j'||q'|+|A_j||q''|]Z_3+3|A_j||q'|^2Z_4}du, where x=alpha/u.

The only conditional tail term is S_1, which is evaluated in the same symmetric sense as full Poisson summation. After that endpoint term is reassembled, S_2, S_3, and the integral remainder are absolutely convergent. No logarithmic truncation loss or artificial dyadic endpoint remains.

## Six-Value Handoff

Apply the exact representation independently to j=0,...,5. Section 11.167 already proves that the eight physical flow observations are fixed rows of these six values, so no derivative of the incomplete-Fresnel main is needed and no seventh moment is introduced.

Use the fixed-roster Morse-Fresnel plus twice-integrated tail decomposition as the primary endpoint route. Do not use the inexplicit Vandehey remainder as a numerical certificate, and do not introduce a hard dyadic partition unless a later explicit comparison proves it smaller after all boundary terms are rejoined.

The new quantitative objects are explicit: six transition variation sums and six Hurwitz-zeta tail integrals. Their physical interval evaluation, Hermitian pairing, and endpoint-terminal composition remain open.

## Pi Provenance

The factor 2*pi comes from e(x)=exp(2*pi*i*x); the Fresnel limit e(-1/8)=exp(-i*pi/4) is the negative-curvature signature. The cotangent identity is the symmetric partial-fraction expansion of pi*cot(pi*x). No geometric circle or fitted pi is introduced.

## Next Target

Build interval majorants for the six transition quantities |c(y_1)|+|c(y_B)|+integral|c'|, sum them over the fixed roster without splitting the Hermitian pairs, and evaluate the explicit Hurwitz-zeta tail envelope. Then compose both physical endpoint functionals with the terminal recurrence and insert the six value transforms into the eight observation rows.

Prove Psi(T_0)-epsilon*integral_0^1 Psi'(T_0-theta*epsilon)dtheta<=h^2/400, preferably h^2/800.

## Proof Boundary

This proves an exact fixed-roster Morse-Fresnel representation for every transition mode, an exact twice-integrated reassembly of the infinite nonstationary tail, stable digamma/Hurwitz-zeta endpoint sums, and explicit constant-one residual majorants. It does not evaluate the summed transition variation norms on the physical chart, bound the imaginary Hermitian coefficient, prove a signed offset-sum or signed flow estimate, upper-bound Phi_B, exclude contact, prove a retained aggregate or Xi theorem, Q209, a cofinal descendant theorem, Lambda<=0, PF-infinity, RH, or a prize-level conclusion.

## Reproduce

```text
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_endpoint_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_endpoint_reduction.py
```
