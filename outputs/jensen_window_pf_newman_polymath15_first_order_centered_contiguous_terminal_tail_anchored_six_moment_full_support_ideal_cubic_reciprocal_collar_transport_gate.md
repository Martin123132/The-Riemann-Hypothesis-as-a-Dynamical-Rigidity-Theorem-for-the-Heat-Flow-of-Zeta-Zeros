# Ideal-Cubic Reciprocal Collar Transport Gate

Date: 2026-08-02

Status: exact reciprocal collar, rail transport, corner phases, and finite-moment boundary identities proved; signed tie-complete ideal estimate open; not a proof of RH.

## Collar Decomposition

For r in C_p define I_(q,r)[P]=integral_(J_q)f_P(u)e(-ru)du and mathcal B_(q,p)[P]=sum_(r in C_p)I_(q,r)[P].

For 0<=L<=N-1 define mathcal C_(N,L)[P]=sum_(|p-q|<=L)mathcal B_(q,p)[P] and mathcal E_(N,L)[P]=sum_(|p-q|>=L+1)mathcal B_(q,p)[P]. Then mathscr F_N[P]=mathcal C_(N,L)[P]+mathcal E_(N,L)[P].

The bookmarked diagonal-plus-adjacent object is mathcal C_(N,1)=sum_q mathcal B_(q,q)+sum_(q<N)mathcal B_(q,q+1)+sum_(q>1)mathcal B_(q,q-1).

For fixed p, the far physical cells are two contiguous rails: [1,p-L-1/2] when p>=L+2 and [p+L+1/2,N] when p<=N-L-1. Hence mathcal E_(N,L) is a sum of at most two oriented integrals per frequency cell, not a collection of independently bounded unit cells.

## Rail Transport

When alpha_P/r crosses p+1/2 upward, the mode r moves from C_p to C_(p+1). Its collar jump is 1_(p+L+1<=N)I_(p+L+1,r)-1_(p-L>=1)I_(p-L,r).

For L=1 the exact jump is 1_(p+2<=N)I_(p+2,r)-1_(p-1>=1)I_(p-1,r). The diagonal and two adjacent blocks therefore do not form an internal-tie-invariant object.

The far complement has the opposite rail jump, so mathcal C_(N,L)+mathcal E_(N,L)=mathscr F_N remains invariant under every internal reassignment.

For the exact test source f(u)=u e(ru), every full-cell atom I_(q,r) equals q. At an interior tie the collar jump is 2L+1, proving that no proper fixed-width collar is coefficientwise tie invariant.

At the two outer reciprocal-band crossings the complete mode, not only its near collar, is transferred to or from the outer remainder. The collar and both far rails must travel together.

On a fixed roster-and-cell chart, mathcal C_(N,L)[i(lambda-c)P]=(partial_xi-ic)mathcal C_(N,L)[P] and likewise for mathcal E_(N,L). At a tie the two rail atoms and their relative lifts are transferred explicitly.

## Exact Local Phases

For r in C_q put v=alpha_P/r and u=v(1+y). Then phi_r(u)=alpha_P(log v-1)+alpha_P{log(1+y)-y}. This is the exact negative Morse chart; no quadratic truncation is made.

For p=q+1 let c=q+1/2, u=c-a, v=c+b, and d=(a+b)/(c+b). Then phi_r(u)=alpha_P(log v-1)+alpha_P{log(1-d)+d}, with phi_r'(u)>=0 and equality only at the shared corner a=b=0.

For p=q-1 let c=q-1/2, u=c+a, v=c-b, and d=(a+b)/(c-b). Then phi_r(u)=alpha_P(log v-1)+alpha_P{log(1+d)-d}, with phi_r'(u)<=0 and equality only at the shared corner.

Both corner corrections are strictly negative away from the shared boundary: d/d d[log(1-d)+d]=-d/(1-d) and d/d d[log(1+d)-d]=-d/(1+d). Their cubic terms have opposite orientations and cannot be replaced by one unsigned half-Gaussian.

## Far Rails

On either far rail, v=alpha_P/r lies in the p-cell and the physical variable lies outside [p-L-1/2,p+L+1/2]. Hence |v-u|>=L and q_r(u)=alpha_P/u-r has no zero for every L>=1.

For q_r=alpha_P/u-r, integral_a^b A e(phi_r)du=kappa[Ae(phi_r)/q_r]_a^b-kappa integral_a^b e(phi_r){A'/q_r-Aq_r'/q_r^2}du on each oriented far interval.

The unit-cell boundary terms telescope before estimation because each far side is one contiguous interval. Only u=1, u=N, and the two collar rails p-L-1/2 and p+L+1/2 remain for each frequency cell.

The rail-compressed identity prevents an artificial O(N^2) boundary budget, but no signed or absolute aggregate bound for the remaining p,r sum is proved here.

## Finite First Correction

The inherited full-line sequential coefficient remains -kappa mathcal D_qA/(2alpha_P)+kappa mathcal D_qA/(2alpha_P)=0 for arbitrary twice-differentiable A.

On a finite negative chart, integral_a^b y^2e(-y^2/2)dy=kappa integral_a^b e(-y^2/2)dy-kappa[y e(-y^2/2)]_a^b.

On a finite positive chart, integral_c^d eta^2e(+eta^2/2)deta=-kappa integral_c^d e(+eta^2/2)deta+kappa[eta e(+eta^2/2)]_c^d.

The full-line cancellation uses vanishing regularized boundary terms. Finite diagonal and corner charts retain the four displayed endpoint terms; they may cancel only after the oriented collar rails, outer endpoints, and terminal recurrence are composed.

## Ideal Cubics

Retain P_H^0=-i(x+u_N)(x-u_N)^2/4 and P_T^0=-i(x-u_N)(x+u_N)^2/4+u_(N,x)(x-u_N) inside every block before a modulus is taken.

At x=-u_q, P_H^0=i(u_q-u_N)(u_q+u_N)^2/4 and P_T^0=i(u_q+u_N)(u_q-u_N)^2/4-u_(N,x)(u_q+u_N).

The exact polynomial relation is P_T^0(x)=-P_H^0(-x)+u_(N,x)(x-u_N). The physical interval is not reflection-symmetric, so this identity alone gives no cancellation or sign.

At q=N, P_H^0(log N)=0 and P_T^0(log N)=-2u_Nu_(N,x). The transpose endpoint must remain joined to +2iu_N mathcal T_N[B_(T,N)^0] and the complex terminal recurrence.

Since mathscr F_N=mathcal C_(N,1)+mathcal E_(N,1), the ideal carrier-only joins are mathcal J_H^0=2mathscr M_N[P_H^0]-mathcal C_(N,1)[P_H^0]-mathcal E_(N,1)[P_H^0] and mathcal J_T^0=2mathscr M_N[P_T^0]-mathcal C_(N,1)[P_T^0]-mathcal E_(N,1)[P_T^0]+2iu_Nmathcal T_N[B_(T,N)^0]. The collar rail jump cancels only against the far rail jump inside these complete joins.

## Handoff

The diagonal-plus-adjacent collar is useful only on a fixed cell chart; it is not a globally closed observation. Preserve its two rail transfers and compose the rail-compressed far functional before seeking a uniform sign.

Insert P_H^0 and P_T^0 into the exact diagonal and corner charts, carry the four finite-moment boundary terms into the two rail interfaces, join 2mathscr M_N and the transpose terminal survivor, and then prove or falsify a signed aggregate for the resulting tie-complete collar-plus-far functional.

## Pi Provenance

The phases use only e(x)=exp(2pi i x) and kappa=1/(2pi i). The finite moment identities follow by differentiating e(+-y^2/2); no geometric constant is inserted.

## Proof Boundary

No signed ideal-cubic collar, far aggregate, finite first-correction remainder, h^2 reserve, completed current, or RH-level theorem is proved. This gate proves no signed physical collar or completed-current bound, `Phi_B` theorem, contact exclusion, retained aggregate or Xi theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion.
