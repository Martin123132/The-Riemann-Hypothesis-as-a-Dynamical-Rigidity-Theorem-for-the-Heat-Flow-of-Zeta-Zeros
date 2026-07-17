# Jensen-Window PF Newman Polymath-15 Oscillatory Zeta Handoff Theorem

Date: 2026-07-17

Status: quantified asymptotic strict first-Laguerre positivity for the
exact Newman heat flow above an explicit oscillatory threshold. This is
not a proof of `Lambda <= 0` or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_oscillatory_zeta_handoff_theorem.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_oscillatory_zeta_handoff_theorem.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_oscillatory_zeta_handoff_theorem.py
```

The finite heat-flow approximation is imported from D. H. J. Polymath,
[Effective approximation of heat flow evolution of the Riemann xi function](https://arxiv.org/abs/1904.12438).
The cancellation input is the exponent-pair definition and explicit
known-pair hull in Trudgian and Yang,
[Toward optimal exponent pairs](https://arxiv.org/abs/2306.05599); the central
pair is originally due to [Bourgain](https://arxiv.org/abs/1408.5794).

## Scaled Geometry

```text
L=log(x/(4*pi)), c=t*L, N=floor(sqrt(exp(L)+t/16)), L=2*log(N)+o(1), sigma=Re(s_*)=1/2+c/4+o(1), tau=|Im(s_*)|=2*pi*N^2*(1+o(1))
D_k=sum_(n<=N)exp((t/4)log(n)^2)log(n)^k*n^(-s_*), Z_k=sum_(n>=1)log(n)^k*n^(-s_*), k=0,1,2
```

The ordinary zeta tail is already absolutely small for `c>2`. The new
task is to show that the heat weight minus one contributes `o(1)`.

## Exponent-Pair Envelope

On a dyadic block `M=N^r`, the published exponent-pair definition gives

```text
For M=N^r and exponent pair (kappa,lambda), sup_I|sum_(n in I)n^(-i*tau)| <<_eta N^(2*kappa+(lambda-kappa)*r+eta)
B_(kappa,lambda,c)(r)=2*kappa+(lambda-kappa)*r-(1/2+c/4)*r+c*r^2/8
```

The exponent pair may change with `r`. Exact rational optimization of
eleven available pairs gives

```text
c_*=4911678521/1933561194=2.540223984760008...; it occurs at r_*=125662/155153 where exponent-pair rows 2 and 1 cross
```

| pair change | r | required c | maximum |
|---:|---:|---:|:---:|
| 10->9 | 0.461538461538462 | 2.426666666666667 |  |
| 9->8 | 0.502948240762826 | 2.459744358354176 |  |
| 8->7 | 0.565028002489110 | 2.486656146960047 |  |
| 7->6 | 0.584722760759985 | 2.498315468187929 |  |
| 6->5 | 0.597504888587748 | 2.503774097755819 |  |
| 5->4 | 0.622271377649349 | 2.506330852075120 |  |
| 4->3 | 0.755407653910150 | 2.520153596080002 |  |
| 3->2 | 0.788496732026144 | 2.535377425991032 |  |
| 2->1 | 0.809923108157754 | 2.540223984760009 | yes |
| 1->0 | 0.857142857142857 | 2.527777777777778 |  |

For G(r)=a+b*r, d=b-1/2>=0, the derivative numerator of 8*(G(r)-r/2)/(r*(2-r)) is d*r^2+2*a*r-2*a. It is strictly increasing, so any interior critical point is a minimum; the maximum on each active interval is at an endpoint.
Thus the table is a proof certificate for the full continuous scale
interval, not a sampled optimization.

## Zeta Handoff

Fix `epsilon>0`. Choose the low/high split just below the absolute
summability endpoint. Then

```text
For c>=c_*+epsilon choose 2-4/c_*<rho<2-4/(c_*+epsilon). On n<=N^rho, exp((t/4)log(n)^2)n^(-sigma)<=n^(-1-delta_epsilon), so sum (w_n-1)log(n)^k*n^(-s_*)=O_epsilon(1/L).
On N^rho<n<=N, dyadic Abel summation and the 11-pair envelope give sum (w_n-1)log(n)^k*n^(-s_*)=O_epsilon(N^(-delta_epsilon)).
A_k(u)=(exp(q)-1)u^(-sigma)log(u)^k, q=t*log(u)^2/4; |A_k'(u)|<=u^(-sigma-1){exp(q)*(t/2)*log(u)^(k+1)+sigma*(exp(q)-1)*log(u)^k+k*(exp(q)-1)*log(u)^(k-1)}. Since t*log(u)=O(1) on u<=N, one dyadic total variation is O_epsilon(log(N)^(k+1)*max(exp(q)u^(-sigma))).
sum_(n>N)log(n)^k*n^(-s_*)=O_epsilon(N^(1-sigma)log(N)^k)=o_epsilon(1)
```

Partial summation costs only the total variation of the smooth heat
weight on each dyadic block. Its logarithmic derivative is bounded
for bounded `c`, logarithmic moment factors are absorbed in the strict
power margin, and there are only `O(log N)` blocks. Consequently

```text
D_k-Z_k=o_epsilon(1), k=0,1,2, uniformly for c_*+epsilon<=c<=25; hence D and its first two x-jets converge to zeta(s_*) and its corresponding x-jets
sigma>1 and the Euler product give |zeta(s_*)|>=1/zeta(sigma_0)>0 for a fixed sigma_0>1 depending only on epsilon
```

## Curvature Transfer

Writing `D=a*exp(i*psi)` and `theta=beta+psi`, exact differentiation
gives

```text
For D=a*exp(i*psi), theta=beta+psi, u=(log(a))', v=theta': C[P]=4*a^2*(v^2+(V-u')cos(theta)^2+(v'/2)sin(2*theta))
a>=m_epsilon>0, u'=O_epsilon(1), v=-L/4+O_epsilon(1), v'=O_epsilon(1), V=o(1)
C_t[P_(N,t)]>=c_epsilon*L^2 for all sufficiently large L when c_*+epsilon<=t*L<=25
```

The Euler-product floor makes division by `D` legitimate for all
sufficiently large `L`. The `v^2` term is of order `L^2`; all possible
negative logarithmic corrections remain bounded.

On radius-`1/L` cutoff collars the published scalar remainder and one
adjacent cutoff block satisfy

```text
sum_(n<=N)w_n*n^(-sigma)<=N^(1/2-c/8+o(1)); the e_A+e_B bracket is O(N^-2), so e_A+e_B=O(N^(-3/2-c/8+o(1))); e_C0 and one adjacent-N block are O(N^(-1/2-c/8+o(1))).
On radius-1/L cutoff collars, the normalized Polymath-15 A+B remainder plus one adjacent-N block is E=O_epsilon(N^(-1/2-c_*/8+o(1))); its j-th jet is O_epsilon(L^j*E), j=0,1,2, and its curvature cost is o_epsilon(L^2).
```

The exact normalized-curvature perturbation identity therefore
preserves the main sign. Combining with the existing `c>=25` ray gives

```text
For every epsilon>0 there exists L_epsilon such that 0<t<=1/2, L>=L_epsilon, and t*L>=4911678521/1933561194+epsilon imply L_t(x)>0
```

## Remaining Gap

The asymptotically unresolved layer is now

```text
0<t*log(x/(4*pi))<=4911678521/1933561194+o(1)
```

The improvement uses established cancellation estimates; it is not a
finite-height verification and does not prove `Lambda <= 0` or RH.
