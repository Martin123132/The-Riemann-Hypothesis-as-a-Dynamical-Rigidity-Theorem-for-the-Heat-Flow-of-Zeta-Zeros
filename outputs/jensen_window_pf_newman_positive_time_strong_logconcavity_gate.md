# Jensen-Window PF Newman Positive-Time Strong Log-Concavity Gate

Date: 2026-07-17

Status: exact positive-time strong-log-concavity theorem with an
Xi-specific nonpromotion gate. This is not a proof of RH or
`Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_newman_positive_time_strong_logconcavity_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_positive_time_strong_logconcavity_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_positive_time_strong_logconcavity_gate.py
```

Current result:

```text
validated Jensen-window PF Newman positive-time strong-log-concavity gate: 9 rows, 0 issues, 1 published input, 2 exact curvature/admissibility theorems, 1 Arb threshold certificate, 1 Prekopa correlation theorem, 1 square-transform identity, 1 Xi-specific nonpromotion gate, 1 target-window margin, 1 weighted-hierarchy handoff
```

## Published Input

Csordas and Varga proved

```text
q(r)=log(Phi(sqrt(r))) is strictly concave for r>0
```

Primary source: https://doi.org/10.1007/BF02075457.

## Uniform Curvature

Writing `u=sqrt(r)` gives

```text
q''(u^2)=(u*(log Phi)''(u)-(log Phi)'(u))/(4*u^3)<0
With kappa=-(log Phi)''(0)=-Phi''(0)/Phi(0), (log Phi)''(u)<=-kappa for every real u
(log phi_t)''(u)<=(log Phi)''(0)+2*t=-(kappa-2*t)
```

Thus `phi_t` is positive, even, strictly decreasing,
super-exponentially decaying, and uniformly strongly log-concave
whenever `0<=t<kappa/2`.

## Arb Threshold

The origin values use the exact theta-summand formulas

```text
phi_n(0)=a*(2*a-3)*exp(-a)
phi_n''(0)=a*(32*a^3-224*a^2+330*a-75)*exp(-a)
a=pi*n^2
Phi(0)=[0.4466969004671234440869846670547091132204243670948249747 +/- 3.09E-56]
Phi''(0)=[-33.46100154940651373358024870154595919815511266249142184 +/- 5.37E-55]
kappa=[74.90761971801328467366520780869287945907515433156308220 +/- 1.63E-54]
kappa/2=[37.45380985900664233683260390434643972953757716578154110 +/- 8.13E-55]
kappa/2-1/5=[37.25380985900664233683260390434643972953757716578154110 +/- 8.13E-55]
```

The finite sum ends at `n=8`; explicit decreasing-ratio geometric
bounds enclose every omitted `n>=9` term.

## Correlation Propagation

On the positive half-line in `s`, the logarithmic Hessian has
quadratic form

```text
D^2 log integrand on (s,v) has quadratic form a*(p+q)^2+b*(p-q)^2-2*n*p^2/s^2, a=(log phi_t)''(s+v), b=(log phi_t)''(s-v)
D^2 log integrand <=-2*(kappa-2*t)*I for 0<=t<kappa/2
Prekopa marginalization gives (log K_(n,t))''(v)<=-2*(kappa-2*t) for every n>=0
```

This follows by Prekopa marginalization. In particular, throughout
`0<=t<=1/5`, every correlation is much more strongly log-concave
than a merely qualitative admissibility argument records.

## Nonpromotion Gate

The same shape package does not imply Fourier zero-freeness:

```text
Fourier[K_(0,t)](xi)=2*H_t(xi/2)^2
K_(0,0) is positive, even, strongly log-concave, positive definite, and has the same Xi super-Gaussian correlation tail
Hardy's theorem gives infinitely many real zeros of H_0, so the Fourier transform of K_(0,0) has real double zeros
```

Uniform strong log-concavity, admissibility, positive definiteness, and the exact Xi correlation tail do not force Fourier zero-freeness: K_(0,0) has all of them and still has double Fourier zeros. Any proof for K_(1,t) must use the s^2 weighting or a coupling in the correlation hierarchy.

## Live Handoff

Exploit structure that distinguishes K_(1,t) from K_(0,t), such as a relative F_2/F_1 curvature inequality, an s^2-weighted modular square, or a hierarchy coercivity estimate. Shape-only Fourier criteria are closed.
