# Hermitian Reciprocal Pairing Reduction

Date: 2026-08-01

Status: exact reduction artifact for the interior Hermitian reciprocal main;
`0 Poisson remainder bounds`, `0 signed flow bounds`, and this is not a proof
of RH.

## Domain

Use the physical q=1, L>=50, fixed-N terminal/bulk chart and the auxiliary segment from Sections 11.164-11.165. The present reduction concerns only the interior Hermitian reciprocal stationary main at xi=T_0; all endpoint and remainder terms remain explicit obligations.

## Swap Pairing

Let Z_(k,l)=A_(k,l)e(alpha*mu_(k,l)), where e(t)=exp(2*pi*i*t), mu_(k,l)=log(l/k), and the interior Hermitian stationary coefficients obey A_(l,k)=conj(A_(k,l)). Then Z_(l,k)=conj(Z_(k,l)) and mu_(l,k)=-mu_(k,l).

Writing A_(k,l)=R_(k,l)+iI_(k,l), the two ordered flow terms satisfy i*mu*Z_(k,l)-i*mu*conj(Z_(k,l))=-2mu{R sin(2*pi*alpha*mu)+I cos(2*pi*alpha*mu)}.

After the existing outer factor 1/2, one unordered pair k<l contributes P_H(k,l)=-mu{R sin(2*pi*alpha*mu)+I cos(2*pi*alpha*mu)}.

The ordered flow pair is exactly (1/(2*pi))*partial_alpha [Z_(k,l)+conj(Z_(k,l))] when the stationary coefficient is held fixed. Equivalently, since xi=2*pi*alpha, it is the xi phase derivative.

At k=l one has mu=0, so the full Hermitian differentiated stationary pair vanishes, including every complex correction.

The real channel is -mu*R*sin(2*pi*alpha*mu), while the imaginary correction is -mu*I*cos(2*pi*alpha*mu). Thus the real leading kernel gains a sine zero at integral phase, but an imaginary correction generally gains only the single mu factor.

## Leading Reciprocal Main

At n_*=a^2/k and m_*=a^2/l, put mu=log(l/k) and tau=log(kl/a^2). Then u_n-u_m=-mu and u_n+u_m=tau.

The Hermitian Hessian is diag(-k^2/a^2,+l^2/a^2), has signature zero, and contributes the positive leading factor a^2/(kl); no Maslov phase is introduced by this mixed signature.

For the ideal q=1 kernel H_0^+=-tau^2/8 and positive real normalized saddle amplitudes W_k,W_l, the interior stationary coefficient is A^0_(k,l)=-C_(k,l), where C_(k,l)=a^2*tau^2*W_kW_l/(8kl)>=0.

The corresponding outer-half unordered Hermitian flow main is P_H^0(k,l)=C_(k,l)*mu*sin(2*pi*a^2*mu). It has the exact diagonal zero and a second zero whenever a^2*mu is integral, but it has no fixed sign between those phases.

## Offset And Lattice Defect

For k<l write l=k+d with integer d>=1. Then mu_(k,d)=log(1+d/k), F_k(d)=alpha*mu_(k,d), F_k'(d)=alpha/(k+d), and F_k''(d)=-alpha/(k+d)^2. The two-dimensional ordered main is therefore reduced to one triangular family of one-dimensional logarithmic phases.

For x=d/k>=0, x-log(1+x)=integral_0^x t/(1+t)dt. Hence 0<=x-log(1+x)<=x^2/2 and 0<=mu<=d/k.

Let rho_k=dist(alpha/k,Z). Choosing an integer q with |alpha/k-q|=rho_k gives dist(alpha*log(1+d/k),Z) <= d*rho_k+alpha*d^2/(2k^2).

The integer dq may be subtracted from the phase. The triangle inequality leaves d|alpha/k-q| plus alpha[d/k-log(1+d/k)], and the displayed logarithmic defect bounds the latter term.

Because |sin(2*pi*y)|<=2*pi*dist(y,Z), the real leading channel satisfies |P_H^0(k,k+d)| <= 2*pi*C_(k,k+d)*(d/k)*[d*rho_k+alpha*d^2/(2k^2)].

For a full stationary coefficient R+iI, the exact swap formula gives |P_H(k,k+d)| <= 2*pi*|R|*(d/k)*[d*rho_k+alpha*d^2/(2k^2)]+|I|*d/k. The first term exposes reciprocal-lattice suppression; the second is the correction obligation that cannot be discarded.

A proof may split by rho_k. Small rho_k is controlled by the displayed sine defect if R and I are bounded sharply. On the complement, the monotone derivative alpha/(k+d) identifies the remaining one-dimensional exponential-sum problem. This is a route decomposition, not a completed estimate.

## Cutoff And Guards

Let K={ceil(a^2/B),...,floor(a^2)}. The interior two-variable stationary main is indexed by KxK, which is invariant under (k,l)<->(l,k). Its |K| diagonal modes vanish and its |K|(|K|-1) off-diagonal modes form exactly |K|(|K|-1)/2 unordered swap pairs.

There is no unpaired mode in the square interior stationary main. This combinatorial fact does not pair the one-variable endpoint family, hard-cutoff half-weights, nonstationary modes, Poisson remainders, or the adjacent terminal recurrence.

For R=1,I=0,mu=1/4 the outer-half pair is -1/4 at alpha=1 and +1/4 at alpha=3, so swap symmetry has no sign. For R=0,I=1,mu=1/2,alpha=2 the phase is integral but the pair is -1/2, so phase resonance does not kill the imaginary defect.

## Pi Provenance

The sine is sin(2*pi*a^2*mu) because Poisson summation uses e(t)=exp(2*pi*i*t) and T_0=2*pi*a^2. The flow multiplier is i*mu because partial_xi e[xi*mu/(2*pi)]=i*mu e[xi*mu/(2*pi)]. Thus the factor 1/(2*pi) in the alpha-derivative identity is forced by the same two inherited normalizations; no new pi is inserted.

## Handoff

Derive the endpoint-complete discrete Poisson/B-process formula with uniform value and first-x derivative remainders. In its Hermitian main, bound the real leading offset sums using the reciprocal-lattice split, prove a source-specific bound for the imaginary coefficient I, and retain the endpoint, transpose, hard-cutoff, and adjacent-recurrence terms before taking absolute values.

Prove Psi(T_0)-epsilon*integral_0^1 Psi'(T_0-theta*epsilon)dtheta<=h^2/400, preferably h^2/800.

## Proof Boundary

This proves the exact Hermitian reciprocal swap pairing, its outer-half sine/cosine formula, the real-leading reciprocal coordinates and stationary factor, the offset phase derivatives, the logarithmic and reciprocal-lattice defect bounds, complete pairing of the square interior dual main, and three algebraic nonpromotion guards. It proves no discrete Poisson or B-process formula, stationary-phase remainder, hard-cutoff boundary estimate, imaginary-correction bound, signed offset-sum gain, signed flow estimate, Phi_B bound, contact exclusion, retained aggregate or Xi theorem, Q209, cofinal descendant theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion.

## Reproduce

```text
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_hermitian_reciprocal_pairing_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_hermitian_reciprocal_pairing_reduction.py
```
