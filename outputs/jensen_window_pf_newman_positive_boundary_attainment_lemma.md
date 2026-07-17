# Jensen-Window PF Newman Positive-Boundary Attainment Lemma

Date: 2026-07-17

Status: exact positive-boundary attainment and arbitrary-multiplicity
energy lemma. This is not a proof of RH or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_newman_positive_boundary_attainment_lemma.json
python work/rh_compute/scripts/jensen_window_pf_newman_positive_boundary_attainment_lemma.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_positive_boundary_attainment_lemma.py
```

Current result:

```text
validated Jensen-window PF Newman positive-boundary attainment lemma: 10 rows, 0 issues, 2 published compactness inputs, 1 finite-boundary attainment theorem, 1 positive-time simplicity equivalence, 1 arbitrary-multiplicity Hermite split, 9 exact Hermite checks, 1 cluster-energy blow-up, 1 open Xi endpoint handoff
```

## Boundary Attainment

The Polymath positive-time theorems have the uniformity needed for
a compactness argument:
Primary source: https://arxiv.org/abs/1904.12438.

```text
Since H_0 has no zero with |Im z|>1, every zero of H_t for 0<t<=1/2 satisfies |Im z|<=sqrt(1-2t)<=1.
There are absolute C,c>0 such that, uniformly for 0<t<=1/2, if x>=exp(C/t) and H_t(x+iy)=0 then y=0; each sufficiently high reference disk contains exactly one zero, which is real.
```

Assume `Lambda>0` and choose `t_n<Lambda` tending to `Lambda` with
`t_n>=Lambda/2`. Each `H_(t_n)` has a nonreal zero. The two
published estimates force every such choice into

```text
|Re z_n|<X_Lambda=exp(2*C/Lambda),  |Im z_n|<=1.
```

After taking a subsequence, local uniform convergence gives a real zero
`z_*` of `H_Lambda`. Both `z_n` and its distinct conjugate converge to
that same point, so Rouche zero counting makes `z_*` multiple. Hence

```text
If Lambda>0, H_Lambda has a finite real multiple zero c with |c|<=X_Lambda=exp(2*C/Lambda). Thus a positive Newman boundary cannot be realized solely by zero collisions escaping to infinity.
```

This closes the high-index escape loophole specifically in the
hypothetical positive-boundary regime. It also yields the exact
reformulation

```text
Using simplicity for every t>Lambda and the published bound Lambda<=1/5, Lambda<=0 if and only if H_t has only simple zeros for every 0<t<=1/5.
```

## Multiplicity Cluster

At a multiplicity-m zero c at time t_*, the nearby zeros for tau=t-t_*>0 are c+sqrt(2*tau)*lambda_a+O(tau), where lambda_a are the simple real roots of the probabilists' Hermite He_m.
The Hermite ODE then gives

```text
He_m''(x)-x*He_m'(x)+m*He_m(x)=0
sum_(b!=a) 1/(lambda_a-lambda_b)=lambda_a/2
sum_(b!=a) 1/(lambda_a-lambda_b)^2=(4*(m-1)-lambda_a^2)/12
sum_a lambda_a^2=m*(m-1)
sum_(a!=b) 1/(lambda_a-lambda_b)^2=m*(m-1)/4
```

Therefore, with `tau=t-Lambda`,

```text
x_b(t)-x_a(t)=sqrt(2*tau)*(lambda_b-lambda_a)+O(tau)
sum_(a!=b) 1/(x_a-x_b)^2=m*(m-1)/(8*tau)+O(tau^(-1/2))
The Rodgers-Tao linear reference corrections are O(1)+O(sqrt(tau)), so the ordered renormalized cluster has the same leading term.
integral_0^epsilon E_cluster(tau) dtau=+infinity for every m>=2
```

The earlier double-zero coefficient is the `m=2` case. No separate
multiplicity assumption is needed for the endpoint energy obstruction.

## Live Handoff

Prove either Xi positive-time simplicity for every 0<t<=1/5, or an Xi-specific finite-truncation energy estimate that remains integrable down to every hypothetical positive boundary. The published Rodgers-Tao estimate starts strictly after a negative boundary and does not supply this endpoint control.

The route is now global enough to cover every hypothetical positive
Newman boundary, but the Xi-specific simplicity or endpoint-integrability
theorem remains the genuinely missing step.
