# Symmetric finite-Poisson interchange for the Kummer roster

Date: 2026-08-10

Status: exact height-uniform interchange theorem validated; not a proof of the stationary-phase remainder

Write the contiguous odd roster as `alpha=A+2u`, `u=0..L`, and set

```text
f_x(u)=(A+2u) exp(i*pi*x*(A+2u)^2/4),
W_t(x)=exp[i(t/2)log((1-x)/x)]/[x(1-x)]^(1/4).
```

For a symmetric Fourier cutoff `M`, finite Poisson gives the partial
reconstruction

```text
K_(t,M)=int_0^1 W_t(x)
  {[f_x(0)+f_x(L)]/2+sum_(m=-M)^M I_m(x)} dx,
I_m(x)=int_0^L f_x(u)exp(-2*pi*i*m*u)du.
```

Two integrations by parts, followed by pairing `m` and `-m`, give

```text
I_m+I_(-m)=[-2 Delta f_x'+R_m+R_(-m)]/(2*pi*i*m)^2.  (SP1)
```

The apparent `Delta f_x/(2*pi*i*m)` endpoint current cancels exactly in
(SP1).  Consequently

```text
|I_m+I_(-m)|
 <= (|Delta f_x'|+int_0^L |f_x''(u)|du)/(2*pi^2*m^2).
```

The exact derivatives are

```text
f_x'=exp(i*pi*x*alpha^2/4)[2+i*pi*x*alpha^2],
f_x''=exp(i*pi*x*alpha^2/4)
       [6*i*pi*x*alpha-pi^2*x^2*alpha^3].
```

Since `|W_t(x)|=[x(1-x)]^(-1/4)`, beta integration yields, uniformly for
every real `t` and every integer `M>=1`,

```text
|K_t-K_(t,M)| <= C_AB/(2*pi^2*M),                       (SP2)

C_AB=B(3/4,3/4)
 [4+pi(5B^2-A^2)/4+7pi^2(B^4-A^4)/160].                (SP3)
```

Thus no auxiliary endpoint regulator is required: the Kummer endpoint
singularity is absolutely integrable, and the paired Fourier tail supplies
an integrable majorant independent of height.  This closes the interchange
part of the previous obligation.

For the source roster `A=159577`, `B=5122421`, `L=2481422`, interval arithmetic gives

```text
C_AB/(2*pi^2) = [25519454103834020117297802.049102364899541820144413028107520521758596699227478761409410417863347068009183926703 +/- 8.72e-85].
```

This last number is a structural warning.  The raw `C^2` triangle bound would
need a cutoff larger than

```text
[5103890820766804023459560409.8204729799083640288826056215041043517193398454957522818820835726694136018367853405 +/- 1.07e-82]
```

merely to fall below `0.005`.  It proves convergence and interchange, but is
not a viable quantitative Hardy-error estimate.  The next step must exploit
the joint Kummer/Poisson phase, retain the modewise boundary/Fresnel pairing,
and obtain phase-adapted bounds for the interior, turning, and nonstationary
mode ranges before taking absolute values.

Proof boundary: exact finite-roster Poisson reconstruction and a uniform
interchange/tail theorem only.  No useful stationary-phase constant,
source-aligned hybrid remainder, `Lambda<=0`, RH, or prize-level conclusion is
proved.
