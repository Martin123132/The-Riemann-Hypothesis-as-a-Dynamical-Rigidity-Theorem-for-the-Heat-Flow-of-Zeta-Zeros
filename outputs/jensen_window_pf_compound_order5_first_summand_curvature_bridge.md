# Jensen-Window PF Compound Order-Five First-Summand Curvature Bridge

Date: 2026-07-13

Status: exact order-five first/full curvature transfer with one open
continuous first-summand ceiling. This is not a proof of order-five
entry, PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order5_first_summand_curvature_bridge.json
python work/rh_compute/scripts/jensen_window_pf_compound_order5_first_summand_curvature_bridge.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order5_first_summand_curvature_bridge.py
```

## Second Stable Coordinate

For the continuous first-summand family retain `B,d,ell,g,h=log(g)`
from the completed order-four analysis and put

```text
R(t)=3*B(t)-h(t-1)+2*h(t)-h(t+1)
f(t)=g(t)^2-x(t)^3*g(t-1)*g(t+1)=g(t)^2*(1-exp(-R(t)))
q(t)=log(f(t)/d(t))=2*h(t)-ell(t)+log(1-exp(-R(t)))
```

Then `F_n=f(n+3)` and, with `k=n+4`,

```text
C_n=Delta^2 q(k), k=n+4; C_n^(1)=integral_[-1,1](1-|s|)*q_1''(k+s) ds
```

Differentiation gives the cancellation-preserving formula

```text
q''=2*h''-ell''+phi(R)*R''-chi(R)*(R')^2; phi(z)=1/(exp(z)-1), chi(z)=exp(z)/(exp(z)-1)^2
```

The negative square term is retained; no sign is discarded in the
identity itself.

## Positive Floors

The proved first-summand `J` floor and the moment perturbation give

```text
U_j=2*a_j+L_(j-1)+2*L_j+L_(j+1); |J_j-J_j^(1)|<=U_j<=1/(56*j), j>=319
J_j^(1)>=1/(7*j) and U_j<=1/(56*j) imply min(J_j,J_j^(1))>=1/(8*j)
```

The shifted numerator proving `U_j<=1/(56*j)` has degree
`29` and all
`30` coefficients
positive after `j=319+m`.

At the second stable layer the completed first and full order-four
penalty bounds separately imply

```text
min(R_j,R_j^(1))>=7/(5*j), j>=320
first numerator after j=320+m: m^2+597*m+88622>0
full numerator after j=320+m: 53*m^2+31570*m+4674200>0
```

Thus neither logarithmic Lipschitz step divides by an uncertified
near-zero raw determinant.

## Full-Kernel Transfer

Starting from the proved `0<=delta_j<=2/j^6`, define the explicit
error chain

```text
a_j=2*((j-1)^(-6)+2*j^(-6)+(j+1)^(-6)); |B_j-B_j^(1)|<=a_j
L_j=4*j*a_j; |ell_j-ell_j^(1)|<=L_j
V_j=2*L_j+8*j*U_j; |h_j-h_j^(1)|<=V_j
W_j=3*a_j+V_(j-1)+2*V_j+V_(j+1); |R_j-R_j^(1)|<=W_j
E_j=2*V_j+(5*j/7)*W_j+L_j; |q_j-q_j^(1)|<=E_j
```

The centered difference therefore satisfies the unconditional theorem

```text
|C_n-C_n^(1)|<=E_(k-1)+2*E_k+E_(k+1)<=37/k^2, k=n+4>=321
```

After `k=321+m`, its cleared reserve numerator has degree `52`,
all `53` coefficients are positive, its leading
coefficient is `259`, and its constant
coefficient is `65843676401696298803198228145370180796791894597273289782492818922099539802371335936561917542865953379935593154628824339578880000000`.
At the splice, the scaled error and reserve are

```text
321^2*error_321=36.5675727543293229942652731877
37-321^2*error_321=0.432427245670677005734726812273
```

## Remaining Continuous Target

The exact tent identity shows that it is now sufficient to prove

```text
q_1''(t)<=60/t^2 for every real t>=320
```

Indeed,

```text
q_1''(t)<=60/t^2 => C_n^(1)<=60*[-log(1-1/k^2)]<=60/(k^2-1)<63/k^2, k>=321
q_1''(t)<=60/t^2 on t>=320 => C_n<=100/k^2 on k>=321
```

The constants split the original budget exactly as `37+63=100`.
The displayed continuous theorem is open.

## High-Precision Scout

These finite-upper mpmath saddle quadratures are diagnostics, not
interval enclosures or a uniform theorem.

| t | t J(t) | t R(t) | t^2 q''(t) | margin below 60 |
|---:|---:|---:|---:|---:|
| `321` | `1.004429730197858772672835059` | `1.503847671635953898950880809` | `3.591271183292988374104627833` | `56.40872881670701162589537217` |
| `400` | `1.05780877651205911009808887` | `1.584285208344855018533951904` | `3.884744206666446525850337294` | `56.11525579333355347414966271` |
| `1000` | `1.236184194845727462948537548` | `1.853078071669055362675437732` | `4.792885305725572997992827484` | `55.20711469427442700200717252` |
| `5000` | `1.424423201466112164417804212` | `2.136359935913358137704723217` | `5.497325847320247832090541516` | `54.50267415267975216790945848` |
| `100000` | `1.59329498799636376185408111` | `2.389927924288453642993272314` | `5.823082437117652234154696579` | `54.17691756288234776584530342` |

Observed scaled range: `['3.591271183292988374104627833', '5.823082437117652234154696579']`.
The data are consistent with a limit of `6` from below, leaving
a factor-ten reserve under the proposed analytic constant `60`.

```text
outputs/jensen_window_pf_compound_order5_m100_tail_curvature_reduction.md
outputs/jensen_window_pf_compound_order4_m100_entry_certificate.md
outputs/jensen_window_pf_compound_order4_first_summand_curvature_bridge.md
outputs/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.md
outputs/signed_hankel_jensen_bridge_target.md
outputs/formal_core.md
```
