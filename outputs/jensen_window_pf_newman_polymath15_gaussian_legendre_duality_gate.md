# Jensen-Window PF Newman Polymath-15 Gaussian/Legendre Duality Gate

Date: 2026-07-17

Status: exact equivalence and route nonpromotion gate. The required
cancellation improvement remains open. This is not a proof of `Lambda <= 0` or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_gaussian_legendre_duality_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_gaussian_legendre_duality_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_gaussian_legendre_duality_gate.py
```

## Exact Gaussian Rewrite

For every finite cutoff, completing the Gaussian square gives

```text
exp(t*y^2/4)=(pi*t)^(-1/2)*integral_R exp(-u^2/t+u*y)du; hence D_(N,t)(s_*)=(pi*t)^(-1/2)*integral_R exp(-u^2/t)S_N(s_*-u)du, S_N(z)=sum_(n<=N)n^(-z)
The same finite interchange gives D_k=(pi*t)^(-1/2)*integral_R exp(-u^2/t)S_(N,k)(s_*-u)du for k=0,1,2
```

This is tempting as an alternative to dyadic exponent pairs, but its
scaled saddle contains exactly the same optimization.

## Legendre Equivalence

```text
For t=c/L and L=2log(N)+o(1), exp(-u^2/t)=N^(-2u^2/c+o(1))
If Phi(r) bounds the phase sum on M=N^r, define E(q)=sup_(0<=r<=1)(Phi(r)-q*r)
For fixed r, sup_u{u*r-2u^2/c}=c*r^2/8 at u=c*r/4
With sigma=1/2+c/4, sup_u{E(sigma-u)-2u^2/c}=sup_r{Phi(r)-r/2-c*r*(2-r)/8}
```

The equality follows by interchanging the two suprema and maximizing
the displayed concave quadratic in `u`. Therefore

```text
A Gaussian partial-sum proof using the same pointwise envelope has strict decay exactly when the weighted dyadic frontier has strict decay; neither formulation is stronger
```

## Exact Contact

At the current threshold, the exposed pair-2/pair-1 corner maps to

```text
r_*=125662/155153
q_*=4800718975/7734244776=0.620709469901577...
u_*=31657/61548
E(q_*)=1989040967/9549356844
2u_*^2/c_*=1989040967/9549356844
net exponent=0/1
```

Thus `c_*` is a true equality point for this proof architecture, not an artifact
of whether the estimates are written blockwise or as a
Gaussian heat average.

At `c=2`, the same corner gives

```text
q_2=92322/155153=0.595038445921123...
E(q_2)=11029137521/48144906818
Gaussian cost=3947734561/24072453409
deficit=3133668399/48144906818=0.065088263870695...
```

This reproduces the cancellation/zero-free wall gate independently.

## Nonpromotion And Handoff

```text
At c=c_* the pair-2/pair-1 corner gives equality, not decay. Rewriting the current eleven-pair hull as a Gaussian heat average therefore cannot lower c_*
Lowering c_* requires a genuinely smaller partial-sum profile E(q) near q=0.6207094699015768 and across the exposed q-range, or an argument exploiting cancellation beyond pointwise partial-sum majorants; reaching c=2 still leaves the fixed c<2 phase wall
```

A useful next theorem must therefore improve the actual partial-sum
profile, exploit cancellation between estimates that the pointwise
profile discards, or attack the corrected Xi phase directly. Merely
rewriting the same hull as a heat-semigroup estimate does not move the
boundary.
