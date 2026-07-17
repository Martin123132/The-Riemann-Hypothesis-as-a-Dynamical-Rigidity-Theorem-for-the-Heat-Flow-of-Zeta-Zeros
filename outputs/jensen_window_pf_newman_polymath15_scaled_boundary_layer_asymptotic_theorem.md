# Jensen-Window PF Newman Polymath-15 Scaled Boundary-Layer Asymptotic Theorem

Date: 2026-07-17

Status: historical quantified asymptotic strict first-Laguerre
positivity for the exact Newman heat flow above the scaled threshold
four. The theorem remains valid, but its corpus-wide boundary has
been superseded. This is not a proof of `Lambda <= 0` or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_scaled_boundary_layer_asymptotic_theorem.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_scaled_boundary_layer_asymptotic_theorem.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_scaled_boundary_layer_asymptotic_theorem.py
```

The effective finite approximation and error bounds come from D. H. J.
Polymath, [Effective approximation of heat flow evolution of the Riemann xi function](https://arxiv.org/abs/1904.12438),
Theorem 1.3. The c>=25 branch is supplied by the explicit global-ray
certificate; this note's absolute-value argument treats the historical
bounded-c branch 4+epsilon<=c<=25.

## Scaled Regime

Fix `epsilon>0` and put

```text
L=log(x/(4*pi)), c=t*L, 4+epsilon<=c<=25, t=c/L and L->infinity
D_(N,t)=sum_(n<=N)w_n*n^(-s_*), w_n=exp((t/4)log(n)^2), N=floor(sqrt(exp(L)+t/16))
```

The exact critical-axis formula for the published exponent gives

```text
Re(s_*)=1/2+c/4+o_epsilon(1)>=3/2+epsilon/8 for sufficiently large L
```

At the largest retained index, `log(N)=L/2+o(1)`. Hence

```text
w_n*n^(-Re(s_*))<=n^(-p_epsilon), p_epsilon=1+epsilon/16, 1<=n<=N
```

The number four appears here structurally: the worst exponent tends
to `1/2+c/8`, which is summable exactly when `c>4`.

## Zeta Limit

For `k=0,1,2`, compare the finite weighted moments with the ordinary
Dirichlet moments of zeta. Since `exp(q)-1<=q*exp(q)` for `q>=0`,

```text
D_k=sum_(n<=N)w_n*log(n)^k*n^(-s_*), Z_k=sum_(n>=1)log(n)^k*n^(-s_*); |D_k-Z_k|<=25/(4L)*sum_(n>=1)log(n)^(k+2)n^(-p_epsilon)+sum_(n>N)log(n)^k*n^(-Re(s_*))=O_epsilon(1/L), k=0,1,2
```

The first infinite sum is finite for every fixed epsilon; the second
is exponentially small because `Re(s_*)>1`. Therefore

```text
D_(N,t)^(j)=(zeta(s_*(x)))^(j)+O_epsilon(1/L), j=0,1,2
```

The Euler product supplies a noncircular lower bound independent of
the spectral phase:

```text
|zeta(s_*)|>=1/zeta(3/2+epsilon/8)>0
```

Indeed `|1-p^(-s)|<=1+p^(-Re(s))`, and the resulting product is at
least `1/zeta(3/2+epsilon/8)`.

## Curvature

Write the nonzero finite main sum as `D=a*exp(i*psi)` and set
`theta=beta+psi`. Exact symbolic differentiation gives

```text
For D=a*exp(i*psi), theta=beta+psi, u=(log a)', v=theta': C[P]=4a^2*(v^2+(V-u')cos(theta)^2+(v'/2)sin(2theta))
```

The zeta lower bound and moment convergence imply

```text
a>=m_epsilon>0, u'=O_epsilon(1), v=-L/4+O_epsilon(1), v'=O_epsilon(1), V=o(1)
```

so the `v^2` term is of order `L^2`, while every possible negative
term is bounded. Thus, for a positive epsilon-dependent constant,

```text
C_t[P_(N,t)]>=c_epsilon*L^2 for all sufficiently large L
```

## Exact Heat Flow

On a radius-`1/4` complex collar, the published `eA+eB` term is
exponentially small because its bracket is `O(L^2*exp(-L))`; the
`eC0` term remains `O(exp(-L/4))`; and an adjacent cutoff contributes
one `n^(-p_epsilon)` block at `n` of order `exp(L/2)`. Restoring the
raw normalizer costs only `exp(L/16+O(1))`. Consequently

```text
On fixed radius-1/4 collars, the Polymath-15 raw remainder plus one adjacent-cutoff block is O_epsilon(exp(-kappa_epsilon*L)); its normalized C2 curvature cost is o_epsilon(L^2)
```

Cauchy differentiation and the exact normalized-curvature perturbation
identity preserve the main sign. Combining this bounded-c branch with
the already certified c>=25 ray proves

```text
For every epsilon>0 there exists L_epsilon such that 0<t<=1/2, L>=L_epsilon, and t*L>=4+epsilon imply L_t(x)>0
```

## Historical Method Boundary

This absolute-value argument stops at

```text
0<t*log(x/(4*pi))<=4+o(1)
```

At this threshold the uniform coefficient exponent reaches one, so
absolute Dirichlet control genuinely stops. Any continuation through
that layer must use cancellation rather than another absolute-tail
majorant.

The later oscillatory-zeta handoff does exactly that and lowers the
current corpus-wide boundary to

```text
The stronger oscillatory-zeta handoff lowers the corpus-wide asymptotic boundary to c_*=4911678521/1933561194; the threshold four remains only the boundary of this absolute-value method
```

Thus the strip at and below four is not the live global gap; it is the
historical boundary of this particular method.
