# Jensen-Window PF Kernel Mellin Upper-Wall Certificate

Date: 2026-07-10

Status: interval theorem certificate. This is not a proof of full cone entry,
Jensen-window PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_kernel_mellin_upper_wall_certificate`.

Proof boundary: this certificate proves the all-`k` upper wall `x_k<=1`.
It does not prove the remaining adjacent-`k` wall `x_(k+1)>=x_k`.

Machine-readable result:

```text
work/rh_compute/results/jensen_window_pf_kernel_mellin_upper_wall_certificate.json
```

Generator and checker:

```text
python work/rh_compute/scripts/jensen_window_pf_kernel_mellin_upper_wall_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_kernel_mellin_upper_wall_certificate.py
```

Current result:

```text
validated Jensen-window PF kernel Mellin upper-wall certificate: 8 rows, 0 issues, 200 positive compact intervals, 1 positive analytic ray, 1 remaining open cone clause, 0 ready-to-apply rows
```

## Exact Mellin Reduction

With `y=u^2` and `g(y)=Phi(sqrt(y))`,

```text
M_k(lambda) = integral_0^infinity y^(k-1/2)*exp(lambda*y)*g(y) dy
A_k(lambda) = sqrt(pi)*4^(-k)*M_k(lambda)/Gamma(k+1/2)
x_k = (A_(k+1)/A_k)/(A_k/A_(k-1))
```

The geometric factor `4^(-k)` cancels from `x_k`. The Berwald-Borell
inequality says that the Mellin transform of a log-concave measure on
`[0,infinity)` divided by `Gamma` is log-concave. The primary references
used for that theorem are:

- Christer Borell, *Complements of Lyapunov's inequality*: https://doi.org/10.1007/BF01362702
- Bo'az Klartag and Joseph Lehec, *Poisson processes and a log-concave Bernstein theorem*: https://www.math.tau.ac.il/~klartagb/papers/log_concave_bernstein.pdf

Thus it remains to prove that `g(y)` is log-concave. Multiplication by
`exp(lambda*y)` preserves log-concavity for every real `lambda`.

## Compact Certificate

Exact full-kernel evenness supplies an analytic series
`g(y)/g(0)=sum r_j*y^j`. On `0<=y<=1/25`, coefficient balls through
`j=20` are combined with a degree-21 Cauchy tail and its first two
`y` derivatives. Arb proves

```text
Q(y)=g'(y)^2-g(y)*g''(y) > 0
```

on `200/200` rational subintervals.

- Minimum outward-rounded `Q` lower bound: `9.78226700758389066271907567855175219791775557224228073147730E-1`.
- Minimum outward-rounded `g` lower bound: `1.89290858572780829556206684120636372818546839098206570604994E-1`.
- Cauchy value-tail radius upper: `1.08985800012647478785720452718545299281878305522247205821276E-5`.
- Cauchy first-derivative radius upper: `5.82614702941173926343971117454215819245366688811936451808941E-3`.
- Cauchy second-derivative radius upper: `2.97248926238969247259718700214806234738405009635974326235041E+0`.

## Analytic Ray

For `u>=1/5`, write `Phi(u)=t_1(u)*(1+R(u))` and
`L=log(Phi)`. Log-concavity in `y=u^2` is equivalent to

```text
D(u)=L'(u)-u*L''(u) > 0.
```

The `n=1` contribution is increasing on the ray. Explicit geometric
bounds for the `n>=2` log-derivative perturbation give:

- `n=1` margin lower at `u=1/5`: `5.61192710776685055199113708934062556596874903381393579023416E+0`.
- Full perturbation upper: `4.66053137595613410944899328277945128651787525109002037816389E-5`.
- Full ray `D` lower: `5.61188050245309099065004259940779777145588385506142489003037E+0`.

The compact interval and ray meet at `y=1/25`, so
`y -> Phi(sqrt(y))` is strictly log-concave on the full half-line.
The kernel's double-exponential tail makes all fixed-real-`lambda`
Mellin moments finite.

## Cone Consequence

Berwald-Borell now gives, for every real `lambda` and `k>=1`,

```text
A_k(lambda)^2 >= A_(k-1)(lambda)*A_(k+1)(lambda)
x_k(lambda) <= 1.
```

The separate Cauchy-Schwarz moment lemma already gives the lower wall
`(2*k-1)/(2*k+1)<=x_k`. The only remaining ratio-cone clause is
`x_(k+1)>=x_k`. That clause is still open and is not implied by ordinary
log-concavity or by Berwald-Borell.

## Integration

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_full_kernel_evenness_cauchy_lemma.md
outputs/jensen_window_pf_heat_flow_boundary_threshold_lemma.md
outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md
outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

A compact Arb/Cauchy certificate and an analytic n=1-dominant ray prove that y->Phi(sqrt(y)) is strictly log-concave. Berwald-Borell then proves the upper ratio-cone wall x_k(lambda)<=1 for every real lambda and k>=1. The adjacent-k monotonicity x_(k+1)>=x_k remains open, so this is not cone entry or Lambda<=0.
