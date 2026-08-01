# Newman Theta Adaptive Modular C1 Remainder Contract

Date: 2026-07-24

Status: exact reduction with one open quantitative first-jet target.
This is not a proof of `Lambda <= 0` or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_theta_adaptive_modular_c1_remainder_contract.json
python work/rh_compute/scripts/jensen_window_pf_newman_theta_adaptive_modular_c1_remainder_contract.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_theta_adaptive_modular_c1_remainder_contract.py
```

## Cancellation-Preserving Blocks

Set

```text
omega(u)=(1+erf(3*sinh(4u)))/2
b_n(u)=omega(u)*phi_n(u)+(1-omega(u))*phi_n(-u)
b_n is even, entire, and strictly positive on R
sum_(n>=1)b_n(u)=Phi(u)
```

For every finite T, exp(tu^2)b_n(u) is Schwartz uniformly for 0<=t<=T
Consequently,

```text
B_(n,t)(x)=integral_0^infinity exp(tu^2)b_n(u)cos(xu)du
S_(N,t)(x)=sum_(1<=n<=N)B_(n,t)(x)
R_(N,t)(x)=H_t(x)-S_(N,t)(x)
H_t^(j)(x)=sum_(n>=1)B_(n,t)^(j)(x), j=0,1,2,..., absolutely and uniformly on R_x times [0,T]
```

Unlike a hard half-line theta truncation, every retained block is even, so the modular endpoint cancellation is preserved before transformation.

## Decaying Remainder

For `r_N=sum_(n>N)b_n`, define

```text
mu_(N,j)(T)=integral_0^infinity u^j*exp(Tu^2)*sum_(n>N)b_n(u)du
d_(N,j,m)(T)=(1/2)*sup_(0<=t<=T)||partial_u^m[(iu)^j*exp(tu^2)r_N(u)]||_1
epsilon_(N,j,m)(x,T)=min(mu_(N,j)(T),d_(N,j,m)(T)/|x|^m)
```

Repeated Fourier integration by parts gives

```text
|R_(N,t)^(j)(x)|<=epsilon_(N,j,m)(x,T), j=0,1,2
```

Every d_(N,j,m)(T) is finite because the modular tail is uniformly Schwartz at bounded Newman time.

## Direct C1 Contract

With

```text
J_(N,t)(x)=16*x^4*S_(N,t)(x)
J_(N,t)'(x)=64*x^3*S_(N,t)(x)+16*x^4*S_(N,t)'(x)
```

the full first-jet errors satisfy

```text
|J_t-J_(N,t)|<=16*x^4*epsilon_(N,0,m)
|J_t'-J_(N,t)'|<=64*x^3*epsilon_(N,0,m)+16*x^4*epsilon_(N,1,m)
```

Therefore the exact sufficient condition is

```text
|J_(N,t)|>16*x^4*epsilon_(N,0,m) OR |J_(N,t)'|>64*x^3*epsilon_(N,0,m)+16*x^4*epsilon_(N,1,m)
```

The sufficient disjunction implies (H_t(x),H_t'(x))!=(0,0).

## Adaptive Scale

The component saddle obeys

```text
F_(a,n,t,x)(u)=tu^2+(a+ix)u-pi*n^2*exp(4u), a in {5,9}
4y=atan((x+2ty)/a); n_*^2=sqrt(a^2+(x+2ty)^2)/(4pi)
n_*(a,t,x)=sqrt(x/(4pi))*(1+O(x^-1)) uniformly for bounded t
N_kappa(x,t)=ceil(max(n_*(5,t,x),n_*(9,t,x)))+kappa
```

The diagnostics suggest kappa=2, but no fixed-collar theorem is proved.
The ten stored diagnostics have first scale-level counts

```text
[3, 4, 4, 5, 6, 3, 4, 4, 5, 6]
```

and all lie within the empirical two-index collar. This is supporting evidence, not a collar theorem.

## Stress Gate

The fixed three-block high-frequency scout has

```text
rows=10
minimum x=200 cancellation digits=21.089970237240061894000000000000000000000000000000000000000000000000000000000000
```

Thus the sign is carried by the coupled arithmetic remainder, not by a fixed base margin. The separate exact fixed-block theorem also proves that static absolute bars cannot be cofinal.

## Open Quantitative Target

```text
Find explicit integers m>=5 and kappa>=0, or another controlled adaptive N(x,t), and rigorously prove the direct C1 sufficient disjunction for every 0<t<=1/5 and x>38.
For m>=5, the explicit J error multipliers are x^(4-m)*d_(N,0,m) and 4*x^(3-m)*d_(N,0,m)+x^(4-m)*d_(N,1,m), up to the common factor 16.
```

Derive rigorous N-dependent derivative L1 budgets and a matching lower separation profile for the finite adaptive first jet; intervalize transition cells.

Together with the certified |x|<=38 core and positive-boundary attainment, this target would close the remaining positive-time first-jet obstruction.

The stronger global monotonicity condition is not being smuggled back in:

```text
The stronger global sign -L_t'>0 is false at the Xi Lehmer stress point and is not part of this target.
```

## Proof Boundary

Neither finite diagnostics nor the exact error contract prove the open adaptive separation inequality.
No new `Q_j`, `Lambda<=0`, RH, or Clay-prize result follows.

Current validation:

```text
validated Newman theta adaptive modular C1 remainder contract: 10 rows, 0 issues, 2 exact modular theorems, 3 exact C1 inequalities/compositions, 1 fixed-block guard, 1 adaptive saddle theorem, 10 cancellation diagnostics, 1 open quantitative target
```
