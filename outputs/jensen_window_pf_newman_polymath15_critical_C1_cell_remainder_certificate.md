# Jensen-Window PF Newman Polymath-15 Critical C1 Cell Remainder Certificate

Date: 2026-07-17

Status: explicit corrected `C1` remainder on critical fixed-cutoff
cells. This is not a proof of `Lambda <= 0` or RH; corrected
transversality remains open.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_critical_C1_cell_remainder_certificate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_critical_C1_cell_remainder_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_critical_C1_cell_remainder_certificate.py
```

The scalar approximation comes from the A+B-C theorem in
[Polymath 15](https://arxiv.org/abs/1904.12438). The endpoint lift and its `C0` strip
bounds are supplied by the two preceding local certificates.

## Region

```text
L=log(x/(4*pi))>=50, 0<t<=1/2, 0<t*L<=25, rho=1/L, and the radius-rho disk stays in one prescribed-N cell
On |z-x|<=1/L: X=Re z satisfies L-1<=log(X/(4*pi))<=L+1 and 0<=|Im z|<=1/L
```

## Published Errors

The coefficient envelope and explicit exponential bracket give

```text
sum_(n<=N)(1+|gamma|N^|kappa|n^y)b_n^t*n^(-Re s_*) <=50*exp(L/4)
e_A+e_B<1000*exp(-3L/4)
saved eAB constant = [908.5187784038783415789288628637307876477653783004372588533489246507593 +/- 3.72e-68] < 1000
```

The refined endpoint remainder gives

```text
e_C<100*exp(-3L/4)
saved eC constant = [10.58500008309006664389050936138499443143182591970360549094533621789808 +/- 1.14e-69] < 100
```

Using `e_C` is essential: `e_C0` retains a leading factor one and
would lose the extra saddle decay.

## Analytic Collar

The exact endpoint lift satisfies

```text
|C_hat_(N,t)(z)-C_paper(z)|/A_t(x)<100*exp(-3L/4)
saved lift source constant = [175.3719734790537217943900391351181260828946707735163262435446987656219 +/- 2.82e-68] < 300
sup_collar |B_t(z)|/A_t(x)<2
R_hat=H_t-(A_(t,N)+B_(t,N)-C_hat_(N,t)) is holomorphic and sup_collar |R_hat|/A_t(x)<2500*exp(-3L/4)
```

The saved combined constant is
`2300.000000000000000000000000000000000000000000000000000000000000000000 < 2500`.

## C1 Transfer

Cauchy's estimate on radius `1/L`, followed by differentiation of
the positive real normalizer, yields

```text
Z_t=J_hat_(N,t)+r, where J_hat agrees on the real axis with the corrected main J=P-Q
|r|<2500*exp(-3L/4), |r'|<5000*L*exp(-3L/4)
r^2+(r'/L)^2<32000000*exp(-3L/2)
```

At `L=50` the saved error norm is already
`[8.571638277785849421902128643448128630314615245610054185620196635488349e-26 +/- 1.37e-96] < 10^-24`.
Thus the fixed-cell collision condition is excluded whenever

```text
T_L[J]>32000000*exp(-3L/2) excludes a double zero on a fixed-N cell
```

## Remaining Work

The error side is now explicit on fixed-cutoff cells. Two obligations
remain separate: compare adjacent corrected lifts on cutoff-crossing
disks, and establish the arithmetic lower bound for `T_L[J]`, which
is the remaining RH-level obligation.
Neither is supplied by this certificate.
