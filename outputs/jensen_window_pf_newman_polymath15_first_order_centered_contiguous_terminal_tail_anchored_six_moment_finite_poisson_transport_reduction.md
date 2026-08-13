# Six-Moment Finite-Poisson And Stationary-Transport Reduction

Date: 2026-08-01

Status: exact reduction note with one imported weighted B-process theorem,
`0 explicit uniform remainder constants`, `0 signed flow bounds`, and this is not a proof of RH.

Primary analytic source:

```text
https://arxiv.org/abs/1205.0090
```

## Domain And Moments

Use the physical q=1, L>=50, fixed-N terminal/bulk chart, B=N-M-1, and T_0-epsilon<=xi<=T_0. Hold sigma, t, the normalizer, all terminal jets, and the Section 11.165 coefficient rows fixed along the auxiliary xi segment.

For alpha=xi/(2*pi), define A_j(u)=(log u)^j exp[t(log u)^2/4-sigma log u] and F_(j,xi)(u)=A_j(u)e(alpha log u). Then the six bare sums H_j(xi)=sum_(1<=n<=B)F_(j,xi)(n), 0<=j<=5, differ from Ghat_j only by the fixed normalized source factor.

Writing lambda=log u, partial_lambda A_j=jA_(j-1)+(t*lambda/2-sigma)A_j, with A_(-1)=0.

Since e(alpha log u)=exp(i*xi*log u), partial_xi F_(j,xi)=iF_(j+1,xi) and therefore H_j'(xi)=iH_(j+1)(xi) for 0<=j<=4.

At u=1, A_0(1)=1 and A_j(1)=0 for 1<=j<=5. Thus the lower hard endpoint occurs only in the zeroth moment.

## Full Poisson Formula

For a C2 amplitude on [1,B], full Poisson summation in the symmetric limit gives sum^*_(1<=n<=B)F(n)=lim_(R->infinity)sum_(r=-R)^R integral_1^B F(u)e(-ru)du. The star gives half weight to an integral primal endpoint. This is equation (10) in Vandehey's paper.

The physical bulk sum uses both endpoints with full weight. Therefore H_j=E_j+lim_(R->infinity)sum_(r=-R)^R I_(j,r), where I_(j,r)=integral_1^B A_j(u)e(alpha log u-ru)du and E_j=[F_(j,xi)(1)+F_(j,xi)(B)]/2.

E_j=(1/2){delta_(j,0)+A_j(B)e(alpha log B)}. In particular, E_j'(xi)=iE_(j+1)(xi) for 0<=j<=4.

For every integer r, differentiation under the finite integral gives I_(j,r)'(xi)=iI_(j+1,r)(xi). Applying the sourced Poisson identity separately at orders j and j+1 avoids an unsupported termwise differentiation of a conditionally convergent series.

The mode sum is a symmetric R->infinity limit. It is not silently replaced by an absolutely convergent sum.

## Reciprocal Main

For positive mode r, phi_r(u)=alpha log u-ru has the unique stationary point u_r=alpha/r. Its curvature is negative. Modes r<=0 have no stationary point on [1,B].

The saddle lies in [1,B] exactly when alpha/B<=r<=alpha.

At u_r=alpha/r, phi_r=alpha[log(alpha/r)-1], phi_r''=-r^2/alpha, and 1/sqrt(|phi_r''|)=sqrt(alpha)/r=u_r/sqrt(alpha).

Negative curvature contributes e(-1/8)=exp(-i*pi/4). This restores the stationary signature suppressed when Section 11.165 only located the saddle.

The order-j interior main is M_(j,r)=[sqrt(alpha)/r]A_j(u_r)e(alpha[log(u_r)-1]-1/8).

## Stationary Transport

On a fixed dual chart, M_(j,r)'(xi)=iM_(j+1,r)+C_(j,r), where C_(j,r)=[sqrt(alpha)/r]e(alpha[log(u_r)-1]-1/8)/(2*pi*alpha) times {jA_(j-1)+(1/2-sigma+(t/2)log(u_r))A_j}.

The iM_(j+1,r) term differentiates the phase. C_(j,r) is the saddle-motion and stationary-prefactor transport. Dropping it is not a valid first-frequency derivative.

For R_(j,r)=I_(j,r)-M_(j,r), the exact fixed-chart recurrence is R_(j,r)'=iR_(j+1,r)-C_(j,r), 0<=j<=4.

The first-frequency remainder through order four requires only value remainders R_0,...,R_5 and transport terms involving A_0,...,A_4. No seventh logarithmic moment is introduced.

## Observation Handoff

The four physical observations V,N,A,Q and their four xi derivatives are fixed real linear forms of G_0,...,G_5. The physical x differentiation needed to define N and Q has already been performed exactly in their coefficient rows. Therefore one needs six endpoint-complete value transforms, not eight separate value/physical-derivative transforms. Uniform dependence on the fixed physical chart parameters is still required.

The fixed factor eta/S_a and the external common phase c_xi remain outside every Poisson integral. Endpoint carries c_xi, Hermitian products cancel it, and transpose products carry c_xi^2 exactly as in Section 11.165.

## Cutoff Charts

Define P_j=E_j+sum^*_(alpha/B<=r<=alpha)M_(j,r). The dual star gives half weight when alpha/B or alpha is an integer, matching the standard B-process endpoint convention.

With Delta_j=H_j-P_j, on every open xi interval where neither alpha nor alpha/B crosses an integer, Delta_j'=iDelta_(j+1)-sum^* C_(j,r), 0<=j<=4.

As alpha increases through an integer q, the upper condition r<=alpha admits mode q. The one-sided dual-main jump is +M_(j,q); the starred value at alpha=q is its midpoint, and Delta_j has the opposite jump.

As alpha increases through Bq, the lower condition r>=alpha/B removes mode q. The one-sided dual-main jump is -M_(j,q); the starred value is the midpoint, and Delta_j has the opposite jump.

The sharp starred main is piecewise smooth, not globally C1. A proof across a crossing must either retain these exact jumps with the hard endpoint/adjacent recurrence or replace the sharp main by a sourced uniform endpoint transition.

Since xi=2*pi*alpha, upper crossings occur at xi=2*pi*q and lower crossings at xi=2*pi*Bq. Their possible absence on the short auxiliary segment may not be assumed from interval length.

## Sourced Remainder Architecture

After conjugating the negative-curvature phase, Vandehey Theorem 1.1 gives on a regular block [A,D] the starred stationary main plus an error of size O((V+|g(A)|){M/sqrt(T)+log(f'(A)-f'(D)+2)}), with the implicit constant depending on the derivative-comparability constants. This is an imported qualitative bound, not an explicit finite-height constant.

On [U,2U], |f''| lies between alpha/(4U^2) and alpha/U^2, |f'''|<=2alpha/U^3, and |f''''|<=6alpha/U^4. Thus the sourced weighted theorem applies with M=U and T=alpha, after conjugation, uniformly in the block location.

On [1,B], max|f''|/min|f''|=B^2. Therefore a theorem requiring uniform curvature comparability cannot be applied globally with B-independent implicit constants. A dyadic or smooth localization is compulsory.

Hard dyadic blocks introduce artificial half-endpoints. They must be rejoined algebraically so internal endpoint terms cancel, or replaced by a smooth partition of unity. Bounding each artificial endpoint separately changes the signed problem.

Theorem 1.1 retains a logarithmic derivative-range loss and inexplicit constants. It validates the transform architecture but is not yet the h^2/400 remainder theorem. Vandehey's sharper endpoint analysis is a candidate for removing that loss.

## Pi Provenance

The phase e(alpha log u-ru) uses e(y)=exp(2*pi*i*y) and alpha=xi/(2*pi). Negative stationary curvature contributes e(-1/8)=exp(-i*pi/4). These are the standard Poisson and stationary-signature normalizations recorded by Vandehey and the existing reciprocal-symbol guard; no new pi is inserted.

## Handoff

Choose a smooth dyadic localization or Vandehey's uniform endpoint transition, derive explicit value-remainder constants for H_0,...,H_5 uniformly on the full physical chart, rejoin all internal boundaries, and compose the two real hard endpoints with the terminal recurrence. Insert those six transforms into the eight observation rows, then apply the Section 11.166 Hermitian pairing before any absolute value.

Prove Psi(T_0)-epsilon*integral_0^1 Psi'(T_0-theta*epsilon)dtheta<=h^2/400, preferably h^2/800.

## Proof Boundary

This proves the specialization of full starred Poisson summation to the six smooth bare moments, restoration of both full primal endpoints, modewise frequency closure, the negative-curvature stationary main, five exact stationary-transport and residual recurrences, two active-window jump laws, the reduction from eight observations to six value transforms, and dyadic admissibility for a sourced weighted B-process theorem. It proves no explicit uniform Poisson/B-process remainder constant, smooth endpoint-transition estimate, imaginary Hermitian coefficient bound, signed offset-sum gain, signed flow estimate, Phi_B bound, contact exclusion, retained aggregate or Xi theorem, Q209, cofinal descendant theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion.

## Reproduce

```text
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_finite_poisson_transport_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_finite_poisson_transport_reduction.py
```
