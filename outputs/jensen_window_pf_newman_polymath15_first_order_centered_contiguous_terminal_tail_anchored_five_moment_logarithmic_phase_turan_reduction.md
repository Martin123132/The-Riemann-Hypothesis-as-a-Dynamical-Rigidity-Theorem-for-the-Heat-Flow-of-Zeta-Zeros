# Logarithmic-Phase Jet And Turan Reduction

Date: 2026-08-01

Status: exact structural reduction with `0 signed joint bounds` and `0 Phi_B
bounds`. This is not a proof of RH; the physical endpoint/correction estimate
is open.

## Physical Carrier

Use the fixed physical q=1, L>=50 terminal/bulk chart with B=N-M-1 and the correction-free carrier w_n=(eta/S_a)exp[t*log(n)^2/4-s_*log(n)]. The source gives eta/S_a nonzero; no assumption that S_a is positive or real is made.

```text
Write s_*=sigma+i*tau and polarize the nonzero common factor eta/S_a=omega*|eta/S_a|, |omega|=1. Then w_n=omega*A_n*exp(-i*tau*ell_n), where ell_n=log(n) and A_n=|eta/S_a|*exp[t*ell_n^2/4-sigma*ell_n]>0.

F(z)=sum_(n<=B)A_n*exp(-i*z*ell_n), Z(z)=omega*F(z), and G_r=i^r*Z^(r)(tau) for 0<=r<=4.

For every polynomial P of degree at most four, sum_(n<=B)P(ell_n)w_n=[P(i*partial_z)Z(z)]_(z=tau).
```

Thus the ten real bulk coordinates are the order-four real/imaginary jet of
one positive-amplitude logarithmic-frequency polynomial, not ten freely
chosen parameters.

## Pair Phases

```text
w_n*conj(w_m)=A_nA_m*exp[-i*tau*log(n/m)], while w_nw_m=omega^2*A_nA_m*exp[-i*tau*log(nm)].
```

The ratio family has the exact Bochner factorization

```text
For real x_j and complex v_j, sum_(j,k)v_j*conj(v_k)F(x_j-x_k)=sum_(n<=B)A_n*|sum_j v_j*exp(-i*x_j*ell_n)|^2>=0.
```

Bochner positivity controls the unweighted phase-difference matrix F(x_j-x_k). It does not control the phase-sum family w_nw_m, the polynomially weighted Hermitian kernel, the endpoint-linear sum, or their joint real projection.

## Shifted Trace

```text
Put u_n=log(a/n)=L_a-ell_n and f(y)=Re[e^(-i*y*L_a)Z(tau-y)]=Re sum_(n<=B)w_n*exp(-i*y*u_n); define g(y) as the imaginary part of the same sum.

For S_r=sum_(n<=B)u_n^r*w_n, f(0)=Re(S_0), f'(0)=Im(S_1), f''(0)=-Re(S_2), and g(0)=Im(S_0).

At C=1, D=0, s_*'=-i/2, s_*''=0, b=-1/2, b_x=0, and chi_N=-i*u_N/2, the real bulk observations are V=f, Q=A=f'/2, and N=f''/4+(u_x/2)g at y=0.
```

Consequently,

```text
P_bulk^(0)=V*N-A*Q=[f*f''-(f')^2]/4+(u_x/2)*f*g at y=0.
```

On every interval where f is nonzero, f*f''-(f')^2=f^2*(log|f|)''. The division-free left side is the primary identity.

Most importantly,

```text
If f(0)=0, then P_bulk^(0)=-(f'(0))^2/4<=0; the u_x*f*g term vanishes without division.
```

## Exact Nonpromotion Guard

```text
n=1,2; a=4; ell_1=0, ell_2=log(2); u_1=2log(2), u_2=log(2)
Take tau=pi/log(2), omega=-1, A_1=1/2, A_2=1. Then exp[-i*tau*log(2)]=exp(-i*pi)=-1, so w_1=-1/2 and w_2=1.
f(y)=cos(log(2)*y)-(1/2)cos(2log(2)*y)
f(0)=1/2, f'(0)=0, f''(0)=log(2)^2, g(0)=0
P_bulk^(0)=log(2)^2/8>0 for every real u_x
```

This exact two-atom example has positive amplitudes, integer logarithmic support, and one common logarithmic phase, but it is not an attained large-q=1 Xi chart. It proves that those abstract properties and Bochner positivity alone cannot sign the leading bulk current.

### Pi Provenance

The pi in this guard is chosen solely to make a half-turn on the log(2) frequency: tau*log(2)=pi and exp(-i*pi)=-1. It is not an inserted physical constant and it proves only a nonimplication. The separate physical u_x=h^2/(8*pi) is inherited from the certified growing-tail source.

## Full-Current Decomposition

```text
Let E=Re sum_(n<=B)L_nw_n and let H_0^+,T_0^+ be the ideal kernels H_0^+(n,m)=-(u_n+u_m)^2/8 and T_0^+(n,m)=-(lambda_n-lambda_m)^2/8-i*u_x/2. Then exactly h^2*Phi_B=E+P_bulk^(0)+R_corr, where R_corr=(1/2)Re sum_(n,m<=B){[H^+-H_0^+]w_nconj(w_m)+[T^+-T_0^+]w_nw_m}.
```

Every kernel difference in R_corr vanishes at the ideal specialization, but no quantitative bound on E or R_corr is proved here. They may not be discarded or absorbed into an unnamed error.

## Open Physical Target

```text
Use the actual physical tau, amplitudes A_n, terminal aggregate, and correction coefficients to prove E+P_bulk^(0)+R_corr<=h^2/400, preferably h^2/800. A viable argument must control the phase-sum family jointly with the weighted phase-difference family and endpoint.
```

## Boundary

This proves the positive-amplitude logarithmic Fourier-jet representation, Bochner positivity of its difference kernel, the exact leading Turan-current identity, zero-fibre sign, an integer-log common-phase nonpromotion guard, and an exact correction decomposition. It proves no signed joint bound, upper bound on Phi_B, contact exclusion, retained aggregate sign, Xi residual transfer, Q209, cofinal descendant theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion.

## Reproduce

```text
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_five_moment_logarithmic_phase_turan_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_five_moment_logarithmic_phase_turan_reduction.py
```
