# Jensen-Window PF Newman Theta Modular-Blend Gate

Date: 2026-07-16

Status: exact positive-time-compatible modular theta blend with a
termwise spectral non-promotion guard. This is not a proof of RH or
`Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_newman_theta_modular_blend_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_theta_modular_blend_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_theta_modular_blend_gate.py
```

## Modular Blend

For the standard theta summands, set

```text
phi_n(u)=(2*pi^2*n^4*exp(9u)-3*pi*n^2*exp(5u))*exp(-pi*n^2*exp(4u))
omega(u)=(1+erf(3*sinh(4u)))/2
b_n(u)=omega(u)*phi_n(u)+(1-omega(u))*phi_n(-u)
```

The switch is entire and satisfies `omega(-u)=1-omega(u)`, so every
`b_n` is even and entire. Theta modular evenness gives the exact partition

```text
sum_(n>=1)b_n(u)=Phi(u)
```

## Strict Positivity

For `u>=0`, `phi_n(u)>0`. If the reflected summand is negative,
write `q=pi*n^2`; then `q*exp(-4u)<3/2` and

```text
phi_n(u)/(-phi_n(-u))=exp(10u)*(2q*exp(4u)-3)/(3-2q*exp(-4u))*exp(-q*(exp(4u)-exp(-4u)))
(1-omega(u))/omega(u)<=erfc(3*sinh(4u))<=exp(-9*sinh(4u)^2), u>=0
log(phi_n(u)/(-phi_n(-u)))>10u-(3/2)*(exp(8u)-1)>-9*sinh(4u)^2
10u-(3/2)*(exp(8u)-1)+9*sinh(4u)^2=10u+(3/4)*exp(8u)+(9/4)*exp(-8u)-3>0 in the negative case
```

Thus the positive forward summand dominates the reflected negative
piece by more than the exact blend ratio, proving

```text
b_n(u)>0 for every n>=1 and every real u
```

## Positive Newman Time

The unwanted reflected tail is multiplied by a double-exponential
switch:

```text
For u->+infinity, omega*phi_n(u) is O_n(exp(9u-pi*n^2*exp(4u))) and (1-omega)*phi_n(-u) is O_n(exp(-5u-9*sinh(4u)^2)).
For every finite T, exp(t*u^2)*b_n(u) is Schwartz uniformly for 0<=t<=T.
```

Positivity and Tonelli then give the exact normal-moment identity

```text
sum_n integral_R |u|^k*exp(t*u^2)*b_n(u)du=integral_R |u|^k*exp(t*u^2)*Phi(u)du<infinity
```

Consequently, for every fixed spectral derivative and bounded time
interval,

```text
H_t^(k)(x)=sum_(n>=1)B_(n,t)^(k)(x), k=0,1,2,...
L_t(x)=sum_(m,n>=1)M_(m,n,t)(x)
```

Both series are absolutely and uniformly convergent. This removes the
cusp obstruction of hard-reflected finite summands and the exponential
tail obstruction of the smooth Euler-cell blocks.

## Tail Enclosure

For a finite blended sum, define

```text
S_(N,t)(x)=sum_(1<=n<=N)B_(n,t)(x), R_(N,t)(x)=H_t(x)-S_(N,t)(x)
mu_(N,j)(T)=integral_0^infinity u^j*exp(T*u^2)*sum_(n>N)b_n(u)du, j=0,1,2
```

Positivity of every omitted block gives a frequency-uniform jet bound
and hence the exact a posteriori criterion

```text
|R_(N,t)^(j)(x)|<=mu_(N,j)(T) for every real x and 0<=t<=T
|L[H_t]-L[S_(N,t)]|<=2*|S'|*mu_1+|S|*mu_2+|S''|*mu_0+mu_1^2+mu_0*mu_2
If L[S_(N,t)] exceeds the displayed error, then L_t(x)>0.
```

This turns the decomposition into a usable finite-to-infinite bridge,
but the partial-sum inequality and tail moments still require rigorous
uniform enclosures on any claimed region.

## Frequency Tail

The three-moment bound cannot itself be global. For every fixed `N`,
the omitted positive tail makes all three moments strictly positive,
whereas Riemann-Lebesgue gives

```text
S_(N,t)(x), S_(N,t)'(x), and S_(N,t)''(x) tend to zero as |x| tends to infinity.
The moment-only error majorant tends to the positive constant mu_1^2+mu_0*mu_2 while L[S_(N,t)] tends to zero.
For every fixed N, the three-moment criterion is necessarily a compact-frequency certificate and cannot close the global tail.
```

The exact repair is to retain oscillatory decay. Since the entire
blended remainder is Schwartz after bounded Newman time, define

```text
d_(N,j,m)(T)=(1/2)*sup_(0<=t<=T)||partial_u^m[(i*u)^j*exp(t*u^2)*r_N(u)]||_L1(R)
|R_(N,t)^(j)(x)|<=d_(N,j,m)(T)/|x|^m for x!=0, j=0,1,2
epsilon_(N,j,m)(x,T)=min(mu_(N,j)(T),d_(N,j,m)(T)/|x|^m)
```

Replacing each `mu_j` by `epsilon_j` in the Laguerre error is valid
and removes the positive floor. The remaining hard task is quantitative:
obtain useful rigorous derivative budgets and compare them with a
global lower profile for the finite blended Laguerre expression.

## Spectral Guard

Positive smooth blocks still do not have nonnegative cosine transforms.
Independent coarse/fine quadrature gives

```text
t=0, n=1, x=35, B=-2.68918284116059830e-05
t=0, n=2, x=25, B=-3.41318721640947986e-06
t=0, n=3, x=15, B=-9.21647952068369118e-09
t=0.5, n=1, x=35, B=-2.57791389573067733e-05
t=0.5, n=2, x=25, B=-3.68009874039114352e-06
t=0.5, n=3, x=15, B=-9.35716644247212534e-09
```

The stronger individual Laguerre-Polya shortcut also fails numerically.
For the first block, form `G_t(z)=B_(1,t)(i*sqrt(z))`; its degree-nine
shift-zero Jensen polynomial has the stable nonreal roots

```text
t=0, degree=9, scaled root=-8.20514016059978646e-01 + 2.53798454211509213e-02*i
t=0.5, degree=9, scaled root=-8.20182070086582282e-01 + 1.44435467728611482e-02*i
```

So termwise Fourier positivity and individual block Laguerre-Polya
promotion are closed. The surviving use of this
decomposition is a coupled matrix estimate, a relative coercivity
bound, or a transform-tail dominance theorem that retains the infinite
arithmetic cancellation. None of those signs is proved here.
