# Ideal-Cubic Symmetric Outer-Pairing Gate

Date: 2026-08-03

Status: exact symmetric outer decomposition and one-sided-divergence obstruction proved; signed carrier-plus-near estimate open; not a proof of RH.

## Exterior Pullback

Full-support Poisson summation gives mathscr F_N[P]+mathscr O_N[P]=mathscr M_N[P], where mathscr O_N is the common symmetric sum over all integer modes outside mathcal R_N.

Because mathscr F_N=tau_band+mathscr M_N-E_band, one has mathscr O_N=E_band-tau_band=-mathscr D_N exactly.

The ideal Hermitian join is mathcal J_H^0=mathscr M_N[P_H^0]+mathscr O_N[P_H^0].

The ideal transpose join is mathcal J_T^0=mathscr M_N[P_T^0]+mathscr O_N[P_T^0]+2iu_N mathcal T_N^[N][B_(T,N)^0].

These are recombinations of the same Fourier involution. They do not imply that either complete join is small or has a sign.

## Canonical Symmetric Split

Put m_N=floor(alpha_P/(N+1/2))+1 and n_N=floor(2alpha_P), with mathcal R_N={m_N,...,n_N}. For every R>=n_N, the finite symmetric outer sum has one canonical disjoint split.

sum_(-R<=r<=R,r notin mathcal R_N)I_P(r)=sum_(r=1)^(m_N-1)I_P(r)+I_P(0)+sum_(k=1)^n_N I_P(-k)+sum_(k=n_N+1)^R{I_P(-k)+I_P(k)}.

Define mathscr N_N[P]=sum_(r=1)^(m_N-1)I_P(r)+I_P(0)+sum_(k=1)^n_N I_P(-k). It is finite and equals integral_1^N f_P(u)K_near,N(u)du.

K_near,N(u)=sum_(j=-(m_N-1))^n_N e(ju). Its mean on every integer unit cell is exactly 1 because its only zero Fourier mode has coefficient 1.

The remaining remote object is the paired series mathscr C_N[P]=sum_(k=n_N+1)^infinity{I_P(-k)+I_P(k)}. The common symmetric cutoff, not a rearrangement of two one-sided series, defines the passage to this limit.

## Absolutely Convergent Remote Pairs

For every positive integer k, I_P(-k)+I_P(k)=2integral_1^N f_P(u)cos(2pi ku)du.

If f_P is C^2 on [1,N], integer endpoints make both sine boundary values zero, and I_P(-k)+I_P(k)={2/(2pi k)^2}{f_P'(N)-f_P'(1)-integral_1^N f_P''(u)cos(2pi ku)du}.

Hence |mathscr C_N[P]|<={|f_P'(N)-f_P'(1)|+||f_P''||_1}/{2pi^2} sum_(k>n_N)k^(-2)<={|f_P'(N)-f_P'(1)|+||f_P''||_1}/{2pi^2 n_N}.

This proves absolute convergence of the paired remote tail. It is not an h^2 estimate: alpha_P occurs in f_P' and f_P'', so the displayed derivative numerator must still be analyzed in the physical scaling.

## One-Sided Obstruction

For Delta_P=f_P(N)-f_P(1), one integration by parts at integer endpoints gives I_P(k)=-kappa Delta_P/k+O(k^(-2)) and I_P(-k)=kappa Delta_P/k+O(k^(-2)), with kappa=1/(2pi i).

If Delta_P is nonzero, the positive and negative one-sided tails have opposite nonzero harmonic leading terms. Neither tail converges by itself; only their common symmetric pairing cancels the 1/k term.

For x=log(u/a), u_N=log(a/N), S(0)=0, and P_H^0=-i(x+u_N)(x-u_N)^2/4, one has P_H^0(log N)=0 and P_H^0(0)=i log N{log(a^2/N)}^2/4.

Thus Delta_H=-i c_H, where c_H=log N{log(a^2/N)}^2/4>0 in the physical regime a>=N>=2.

Consequently I_H(k)=c_H/(2pi k)+O(k^(-2)) and I_H(-k)=-c_H/(2pi k)+O(k^(-2)). The two separate outer infinities diverge with opposite real signs even though their paired sum is O(k^(-2)).

The exact nonphysical source f(u)=u on integer [1,N] gives I_f(k)=-kappa(N-1)/k and I_f(-k)=kappa(N-1)/k for every k>=1. This is an exact one-sided-divergence witness, not merely an asymptotic example.

## Symmetry Audit

The algebraic identity P_T^0(x)=-P_H^0(-x)+u_(N,x)(x-u_N) is exact, but it is only a polynomial identity.

The reflection x->-x sends u to a^2/u because x=log(u/a). It sends [1,N] to [a^2/N,a^2], which is not [1,N]; for a>=N>=2 the image lies at or above N.

Under u=a^2/v, alpha_Plog u-ru becomes 2alpha_Plog a-alpha_Plog v-r a^2/v. Its v-derivative cannot equal alpha_P/v-s for a constant Fourier mode s on an interval when alpha_P>0.

For real S, conjugate(I_P^(alpha_P)(r))=I_(conjugate P)^(-alpha_P)(-r). Conjugation reverses alpha_P as well as r, so it does not pair two modes inside the fixed physical alpha_P family.

The Hermitian terminal Abel coefficient still vanishes at q=N, but that finite identity neither changes Delta_H nor regularizes either one-sided outer series.

## Coefficient-Blind Guard

For the nonphysical Fourier source f_K(u)=e(Ku) on integer [1,N], I_K(r)=N-1 when r=K and is zero for every other integer r; the starred carrier is also M=N-1.

If K lies in mathcal R_N, then mathscr O_N=0 and M+mathscr O_N=N-1.

If K lies outside mathcal R_N, then mathscr O_N=N-1 and M+mathscr O_N=2(N-1). Negating f_K reverses both values.

Therefore the outer pullback, finite Abel identity, and Fourier algebra alone provide neither universal cancellation nor a sign. Any useful result must exploit the physical ideal-cubic coefficient and its joined terminal structure.

## Route Decision

Retire the requested three-way split into independent negative, low-positive, and high-positive infinite tails. Keep the finite low modes and the absolutely convergent remote pairs as the canonical decomposition of the symmetric outer sum.

The complete ideal targets are mathcal J_H^0=mathscr M_N[P_H^0]+mathscr N_N[P_H^0]+mathscr C_N[P_H^0] and mathcal J_T^0=mathscr M_N[P_T^0]+mathscr N_N[P_T^0]+mathscr C_N[P_T^0]+2iu_Nmathcal T_N^[N][B_(T,N)^0].

Keep K_near,N intact and test a cellwise Euler/Abel composition of the starred carrier with its mean-one finite kernel. The remote paired tail may be attached afterward through its absolute C^2 bound; the transpose terminal survivor remains explicit.

Prove or rigorously obstruct a source-specific signed estimate for the carrier-plus-near package, including endpoint half-weights. Do not assign budgets to the divergent one-sided outer tails.

## Pi Provenance

Every pi in this gate comes from the fixed Fourier character e(x)=exp(2pi i x) and kappa=1/(2pi i). The factor 2pi k is the derivative of that character with respect to u. No circle, polygon, fitted constant, or plotted symmetry is inserted.

## Proof Boundary

No signed complete ideal join, physical derivative-size bound, completed-current theorem, Phi_B theorem, contact exclusion, retained aggregate or Xi theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved. This gate is not a proof of RH.
