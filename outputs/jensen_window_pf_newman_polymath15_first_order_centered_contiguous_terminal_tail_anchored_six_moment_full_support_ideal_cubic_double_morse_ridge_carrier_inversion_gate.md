# Double-Morse Ridge/Carrier Inversion Gate

Date: 2026-08-03

Status: exact carrier-coefficient and finite-diagonal obstruction gate; signed ridge-completed estimate open; not a proof of RH.

## Exact Carrier Coefficient

With Lambda_q^alias=sum_(s!=q)^sym B_(q,s), finite Poisson inversion gives S_q=tau_q+B_(q,q)+Lambda_q^alias and B_(q,q)=M_q-E_q^ext. Hence 2M_q-S_q=M_q-tau_q+E_q^ext-Lambda_q^alias.

For the complete physical band, mathscr F_N=tau_band+mathscr M_N-E_band. Therefore 2mathscr M_N-mathscr F_N=mathscr M_N-tau_band+E_band: exactly one starred carrier remains.

Since mathscr D_N=mathscr F_N-mathscr M_N=tau_band-E_band, the same join is mathscr M_N-mathscr D_N. This is an identity, not an exterior-smallness estimate.

For every starred atom, 2omega_N(q)-omega_N(q)=omega_N(q). Thus the surviving carrier has the original endpoint half-weights, not a new unstarred normalization.

The internal lower/upper tie pair is -1/2+1/2=0. The global formula retains only tau_band at the two outer reciprocal boundaries.

## Hyperbolic Ridge Decomposition

The full diagonal rectangle a_q<=r<=b_q, 1<=u<=N maps to a curvilinear domain Omega_q in (s,d). Because eta+y_U(eta) is strictly increasing, every fixed-s fiber is an interval.

The fibers meeting d=0 are indexed by S_q^Delta=[sqrt(2)eta_-,sqrt(2)eta_+]. On this interval their actual limits satisfy d_-(s)<=0<=d_+(s), including one-sided endpoint cases q=1,N.

K_q(s)=integral_(d_-(s))^(d_+(s))e(sd)dd=kappa{e(sd_+)-e(sd_-)}/s for s!=0, with K_q(0)=d_+(0)-d_-(0).

On S_q^Delta write G_q(s,d)=G_q(s,0)+dH_q(s,d). Define Z_q=e(alpha_Plog q)integral G_q(s,0)K_q(s)ds and V_q=e(alpha_Plog q)integral integral dH_q(s,d)e(sd)dd ds.

Let W_q be the contribution of the remaining s-wings. Then B_(q,q)=Z_q+V_q+W_q and M_q=Z_q+V_q+W_q+E_q^ext.

A direct identification Z_q=M_q is equivalent to the additional identity V_q+W_q+E_q^ext=0. Neither the hyperbolic chart nor finite Fourier inversion supplies that identity.

## Exact Obstruction Witness

Take alpha_P=2, q=2, B>=3, and the piecewise-BV test source F_A(u)=A(u)e(alpha_Plog u)=1_[3/2,5/2](u). Full inversion at the interior point u=2 returns M_2=1.

For this source I_A(r)e(2r)=integral_(-1/2)^(1/2)e(-rt)dt=sin(pi r)/(pi r).

The reciprocal cell is [a_2,b_2]=[4/5,4/3], of exact width 8/15, so B_(2,2)=integral_(4/5)^(4/3)sin(pi r)/(pi r)dr={Si(4pi/3)-Si(4pi/5)}/pi.

Since |sin(pi r)/(pi r)|<=1, |B_(2,2)|<=8/15<1=M_2. Therefore |E_2^ext|=|M_2-B_(2,2)|>=7/15.

High-precision audit values are B_(2,2)=-0.019407578480917249745586237162 and E_2^ext=1.01940757848091724974558623716.

The witness is nonphysical and rejects a coefficient-blind finite-diagonal-equals-carrier identity. It also shows that any ridge-equals-carrier claim needs a separate residual-plus-wings-plus-exterior cancellation; it does not disprove such a special physical cancellation or RH.

## Ideal Joins

The ideal joins are J_H^0=mathscr M_N[P_H^0]-tau_band[P_H^0]+E_band[P_H^0] and J_T^0=mathscr M_N[P_T^0]-tau_band[P_T^0]+E_band[P_T^0]+2iu_N mathcal T_N^[N][B_(T,N)^0]. The transpose terminal block remains explicit.

## Route Decision

Retire direct finite-ridge-to-carrier matching. The exact target is the one-carrier global package mathscr M_N-tau_band+E_band, with the transpose terminal survivor in its own channel.

Before applying a 1/s modulus, pull E_band through the certified full-support exterior/complement identity and join its negative-, low-positive-, and high-positive-frequency pieces to the roster fringes, endpoint halves, and terminal recurrence. A useful estimate must act on that complete non-tautological package.

Prove or rigorously obstruct a signed bound for the complete ideal Hermitian and transpose joins after this exterior composition. Only then return to the central and outer hyperbolic residuals and attach Delta_H, Delta_T.

## Pi Provenance

The normalized sinc and sine-integral witness follows from the fixed character e(x)=exp(2pi i x). The constant pi is inherited from that Fourier convention; no circle, polygon, fitted constant, or plotted symmetry is inserted.

## Proof Boundary

No signed ridge-completed ideal bound, exterior/complement bound, completed-current theorem, Phi_B theorem, contact exclusion, retained aggregate or Xi theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved. This gate is not a proof of RH.
