# Jensen-Window PF Uniform First-Summand Heat-Tilt Asymptotic Theorem

Date: 2026-07-13

Status: uniform compact-heat first-summand heat-tilt asymptotic theorem.
This is not a proof of uniform order-four positivity, PF-infinity,
RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_uniform_first_summand_heat_tilt_asymptotic_theorem.json
python work/rh_compute/scripts/jensen_window_pf_uniform_first_summand_heat_tilt_asymptotic_theorem.py
python work/rh_compute/scripts/check_jensen_window_pf_uniform_first_summand_heat_tilt_asymptotic_theorem.py
```

## Published Input

[Cormac O'Sullivan, Zeros of Jensen polynomials and asymptotics for the Riemann xi function](https://arxiv.org/abs/2007.13582), Theorem 5.2, proves
arbitrary-order saddle asymptotics for suitable multipliers. Section 5
explicitly treats the Gaussian-log multiplier used by the heat flow, and
Section 7 supplies the Xi first-summand integral reduction.

```text
I_alpha(f;n)=integral_1^infinity (log t)^n*exp(-alpha*t)*f(t)dt
I_alpha(f;n)=sqrt(2*pi)*u^(n+1)*f(exp(u))*exp(u-n/u)/sqrt((1+u)*n)*(1+sum_(r=1)^(R-1)a_r(f;u)/n^r+O(u^(R*(1+2*lambda_suit))/n^R)), u=W(n/alpha)
```

## Uniform Suitability

For `0<=T<=100`, set

```text
f_T(t)=exp(-T*(log t)^2/16), 0<=T<=100
f_T(t*(1+x))/f_T(t)=exp(-T*((v/8)*log(1+x)+log(1+x)^2/16)), t=exp(v)
```

On the published complex disk, the exponent is uniformly bounded. Cauchy
estimates therefore give

```text
f_T(t*(1+x))/f_T(t)=sum_(r=0)^(K-1) f_(T,r)(v)x^r+O_K(v^K|x|^K), uniformly for 0<=T<=100
sup_(0<=T<=100)|f_(T,r)(v)|=O_r(v^r)
```

Thus the complete compact family is suitable with uniform constants.

## Heat-Tilt Expansion

The exact `t=exp(4u)` change of variables gives

```text
M_k^(1)(T)=C_k*(2*pi^2*I_pi(t^(5/4)*f_T(t);2k)-3*pi*I_pi(t^(1/4)*f_T(t);2k))
```

Applying Theorem 5.2 to both terms and combining their all-order
expansions yields

```text
log M_k^(1)(T)=L_0(k)-T*w^2/16+sum_(r=1)^(R-1) P_(r,T)(w)/k^r+O_R(w^(3R)/k^R), uniformly for 0<=T<=100
log R_T^(1)(k)=-T*w^2/16+sum_(r=1)^(R-1) Q_(r,T)(w)/k^r+O_R(w^(3R)/k^R)
```

## Seven Difference Orders

For `w=W(2k/pi)`,

```text
w'(k)=w/(k*(1+w))
d^m/dk^m w(k)^2=O(w(k)/k^m), m=1,...,7
Delta^m F(k)=integral_[0,1]^m F^(m)(k+s_1+...+s_m) ds
```

Choose the all-order expansion parameter `R>m`. The explicit correction
terms gain `m` inverse powers under a fixed local difference, while
`w^(3R)/k^R=o(w/k^m)`. Consequently

```text
Delta_k^m log R_T^(1)(k)=O(log(k)/k^m) uniformly for 0<=T<=100 and m=2,...,7.
```

The same estimate holds for backward stencils. This closes the sole
first-summand heat-tilt target in the uniform order-four tail reduction.

```text
outputs/jensen_window_pf_uniform_superpolynomial_first_summand_dominance.md
outputs/jensen_window_pf_compound_order4_uniform_heat_eventual_tail_reduction.md
outputs/signed_hankel_jensen_dependency_graph.md
```
