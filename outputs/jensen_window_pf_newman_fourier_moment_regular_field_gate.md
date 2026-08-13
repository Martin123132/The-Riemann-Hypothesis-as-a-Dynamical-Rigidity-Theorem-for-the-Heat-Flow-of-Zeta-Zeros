# Newman Fourier-Moment Regular-Field Gate

Date: 2026-08-03

Status: exact Xi Fourier/score-moment representation and generic positive-kernel nonclosure; the simultaneous Xi theta-contact theorem and degree-uniform Jensen remainder remain open; not a proof of RH.

## Xi Moment Coordinate

Write

```text
f_t(u)=exp(tu^2)Phi(u),
H_t(x)=integral_0^infinity f_t(u)cos(xu)du,
M_k(c,t)=integral_0^infinity u^k f_t(u)exp(icu)du.
```

The theta tail permits differentiation under the integral, giving

```text
H_t^(k)(c)=Re[i^k M_k(c,t)].
```

At a zero `c>0` of exact multiplicity `m`,

```text
Re[i^j M_j]=0                  (0<=j<m),
Re[i^m M_m]!=0,
B_c=Re[i^(m+1)M_(m+1)]/((m+1)Re[i^m M_m]),
ell=(2cB_c-m)/4.
```

With `C_k=int u^k f_t(u)cos(cu)du` and `S_k=int u^k f_t(u)sin(cu)du`, this becomes

```text
m even: B_c=-S_(m+1)/[(m+1)C_m],
m odd:  B_c= C_(m+1)/[(m+1)S_m].
```

The parity signs are part of the theorem and cannot be replaced by absolute values.

## Positive Score Coordinate

The proved Xi score probability is

```text
dnu_t(u)=-f_t'(u)du/f_t(0),
S_t(x)=E_nu[sin(xU)],
H_t(x)=f_t(0)S_t(x)/x.
```

Thus `S_t` has the same multiplicity at `c`, and

```text
B_c=S_t^(m+1)(c)/[(m+1)S_t^(m)(c)]-1/c,
ell=c*S_t^(m+1)(c)/[2(m+1)S_t^(m)(c)]-(m+2)/4.
```

For even `m`, the score ratio is `Cnu_(m+1)/Snu_m`; for odd `m`, it is `-Snu_(m+1)/Cnu_m`. This is the exact Xi-only quantity a field theorem would have to constrain.

## Positive-Kernel Countermodel

Let `sinc(z)=sin(z)/z`. For `pi<y<2pi`, define

```text
H_(m,y,c)(x)=sinc(pi*x/c)^m sinc(y*x/c).
```

This is the characteristic function of a sum of `m` independent uniform variables on `[-pi/c,pi/c]` and one on `[-y/c,y/c]`. Its frequency kernel is therefore even, nonnegative, compactly supported, and log-concave. At `x=c`,

```text
B_c=(y*cot(y)-m-1)/c,
ell=(2y*cot(y)-3m-2)/4.
```

The map `y*cot(y)` decreases continuously from positive infinity to negative infinity on `(pi,2pi)`. Hence this positive-Fourier-kernel family realizes every real `ell`.

Each sinc factor is Laguerre-Polya. Compact support makes

```text
H_tau(x)=integral p(u)exp(tau*u^2)exp(ixu)du
```

defined for every real `tau`. The factorized operator `exp(-tau D_x^2)` preserves the Laguerre-Polya class for `tau>=0`, while the negative-time local polynomial `exp(D_z^2)z^m` has a nonreal pair. The exact heat threshold is therefore zero.

## Strongly Log-Concave Score Guard

Multiply the transform by `exp(-q*x^2/2)`. In frequency this convolves the compact law with a Gaussian of variance `q`. If its compact summand lies in `[-R,R]`, then

```text
(log p_q)''(u)=Var(Y|X=u)/q^2-1/q
              <=-(q-R^2)/q^2.
```

Choosing

```text
q=((m+2)^2*pi^2+1)/c^2
```

makes `q>R^2` uniformly for every `pi<y<2pi`. The kernel is now smooth, positive, even, strictly decreasing, and uniformly strongly log-concave; its score law has a continuous positive density. Its field is

```text
ell=[2y*cot(y)-3m-2-2q*c^2]/4,
```

which still realizes every real value. Thus Fourier positivity, score-probability positivity, and uniform strong log-concavity, even together with a local first-loss heat split, do not constrain the field.

This smoothed model has Gaussian rather than Xi double-exponential tails. Separately, the promoted weighted countermodel proves that Xi-like strong/root-variable log-concavity and theta-type decay do not imply the needed weighted-correlation sign. These statements are kept separate: no single countermodel here is claimed to combine the exact contact with every Xi tail property.

## What Remains

The surviving target is narrower than before. It must use the simultaneous Xi theta/modular source and the displayed multiplicity-m moment equations, or bypass the field through global first-jet winding. Finite theta blocks and termwise spectral squares are already closed by their endpoint defects. The cofinal Jensen route also still requires a remainder uniform in degree.

Pi provenance is explicit: `pi` appears in the countermodel because `sin(pi)=0`, the first nonzero zero of the uniform characteristic function. No `pi` is inserted into the Xi moment identities themselves.

## Audit

The builder and independent checker cover 15 multiplicities, 15 exact sinc contacts, 5 Gaussian-smoothed contacts, and 15 forward/backward Hermite splits.

```text
python work/rh_compute/scripts/jensen_window_pf_newman_fourier_moment_regular_field_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_fourier_moment_regular_field_gate.py
```

Current result:

```text
validated Fourier-moment regular-field gate: 20 rows, 5 sources, 15 multiplicities, 15 direct moment checks, 15 score moment checks, 15 score-jet checks, 15 sinc models, 5 Gaussian models, 5 strong-curvature checks, 15 Hermite checks, 2 open handoffs, 0 Xi field bounds, 0 degree-uniform bounds
```
