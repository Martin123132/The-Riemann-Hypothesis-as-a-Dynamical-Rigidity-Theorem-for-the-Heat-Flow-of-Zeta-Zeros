# Jensen-Window PF Newman Polymath-15 Critical Scaled Coercivity Target

Date: 2026-07-17

Status: exact corrected finite target for the remaining critical scaled
layer. This is not a proof of its sign, `Lambda <= 0`, or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_critical_scaled_coercivity_target.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_critical_scaled_coercivity_target.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_critical_scaled_coercivity_target.py
```

The finite approximation and endpoint correction are imported from
D. H. J. Polymath, [Effective approximation of heat flow evolution of the Riemann xi function](https://arxiv.org/abs/1904.12438),
especially the A+B-C approximation in Corollary 6.4.

## Critical Region

```text
L=log(x/(4*pi))->infinity, c=t*L in [0,c_*+o(1)], c_*=4911678521/1933561194=2.540223984760008..., N=floor(sqrt(x/(4*pi)+t/16))
```

The oscillatory zeta handoff closes every fixed `c > c_*`, where
`c_* = 4911678521/1933561194 = 2.540223984760008...`. Thus only
`0 <= c <= c_* + o(1)` remains asymptotically open. The earlier
absolute zeta-moment argument stops at `c=4`; published exponent-pair
cancellation is what lowers the boundary from `4` to `c_*`.

## Finite Coercivity

Set

```text
D_(N,t)(x)=sum_(n<=N)exp((t/4)log(n)^2)*n^(-s_*(x)), P_(N,t)=2Re(exp(i*beta_t)D_(N,t))
exp(i*beta)D=X+iY, exp(i*beta)D'=U+iV1, exp(i*beta)D''=W+iV2, B=beta'
```

Direct symbolic differentiation, without dividing by `D`, gives

```text
C[P]/4=B^2(X^2+Y^2)+2B(X*V1-U*Y)+U^2-X*W+beta''*X*Y+V_t*X^2
X*V1-U*Y=Im(conj(D)*D')
```

This remains valid at zeros of the complex finite sum. It exposes the
precise competition between the `beta'^2*|D|^2` phase floor, the
symplectic cross term, and the first two Dirichlet jets.

## Endpoint Correction

The uncorrected A+B theorem treats the leading Riemann-Siegel endpoint
term as error. That is too coarse in the critical layer. Define

```text
T'=x/2+pi*t/8, a=sqrt(T'/(2*pi)), N=floor(a), p=1-2(a-N), U_RS=exp(-i((T'/2)log(T'/(2*pi))-T'/2-pi/8))
C0(p)=(exp(pi*i*(p^2/2+3/8))-i*sqrt(2)*cos(pi*p/2))/(2*cos(pi*p)), with removable values at p=+-1/2
C_t(x)=2*(-1)^N*exp(t*pi^2/64)*Re(M_0(i*T')*C0(p)*U_RS*exp(pi*i/8))
Q_(N,t)=C_t/A_t, J_(N,t)=P_(N,t)-Q_(N,t), A_t=|B_t|
```

Its exact curvature contribution is

```text
C[J]-C[P]=-2P'Q'+Q'^2+P*Q''+Q*P''-Q*Q''+V_t*(-2P*Q+Q^2)
```

## Refined Transfer

After extracting the correction, the paper replaces the `e_C0` factor
`1+tilde_epsilon(s_-)+tilde_epsilon(s_+)` by
`tilde_epsilon(s_-)+tilde_epsilon(s_+)`. Since the latter is `O(1/N)`,
the refined remainder gains one full saddle-index factor:

```text
H_t/A_t=J_(N,t)+r_ref; the published e_C replaces the e_C0 factor 1 by tilde_epsilon(s_-)+tilde_epsilon(s_+)=O(1/N)
On radius rho=1/L collars, the refined normalized C2 remainder is o(L^2) after including one adjacent-N correction
```

A radius `1/L` keeps the normalizer ratio bounded, while Cauchy costs
only the powers of `L` naturally present in the curvature. The refined
remainder is then lower order even under crude absolute main-jet caps.

## Live Target

The remaining theorem obligation is now concrete:

```text
Prove C_t[J_(N,t)] greater than its explicit refined C2 remainder budget uniformly for 0<=t*L<=c_*+o(1)
```

This is an oscillatory arithmetic inequality for a finite sum of
length about `sqrt(x/(4*pi))`, with an explicit endpoint correction
and explicit error. Absolute tails cannot settle it; proving its sign
would be the genuinely new RH-level step.
