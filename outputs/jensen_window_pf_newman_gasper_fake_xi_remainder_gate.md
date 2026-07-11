# Jensen-Window PF Newman Gasper Fake-Xi Remainder Gate

Date: 2026-07-11

Status: exact fake-Xi comparison with scalar and positive-convolution
obstructions. This is not a proof of RH or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_newman_gasper_fake_xi_remainder_gate`.

```text
work/rh_compute/results/jensen_window_pf_newman_gasper_fake_xi_remainder_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_gasper_fake_xi_remainder_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_gasper_fake_xi_remainder_gate.py
```

Current result:

```text
validated Jensen-window PF Newman Gasper fake-Xi remainder gate: 8 rows, 0 issues, 2 exact transform identities, 1 established real-zero benchmark, 1 scalar algebra theorem, 2 interval scalar witnesses, 2 high-precision cross-checks, 1 exact positive-convolution obstruction, 1 sign-aware handoff
```

## Gasper Benchmark

With

```text
Psi(u)=4*pi^2*exp(-2*pi*cosh(4u))*cosh(9u)
P_t(z)=integral_0^infinity exp(tu^2)*Psi(u)*cos(zu)du
P_0(z)=Xi_star(z/2)/8
P_0(z)=pi^2/2*(K_(9/4+i*z/4)(2*pi)+K_(9/4-i*z/4)(2*pi))
```

Gasper's integral-of-squares theorem proves that the model at `t=0`
has only real zeros. The classical heat universal-factor theorem
preserves this for `t>=0`. This is a valid comparison model, not an
approximate proof for Xi.

## Scalar Obstruction

For `E_alpha=H-alpha*P`,

```text
L[H]=alpha^2*L[P]+B[alpha*P,E_alpha]+L[E_alpha]
B[alpha*P,E_alpha]=alpha*B[P,H]-2*alpha^2*L[P]
L[E_alpha]=L[H]-alpha*B[P,H]+alpha^2*L[P]
```

The absolute-margin strategy would require `R_alpha<1`, where

```text
R_alpha=(abs(B[alpha*P,E_alpha])+abs(L[E_alpha]))/(alpha^2*L[P])
rho=B[P,H]^2/(L[P]*L[H])
If L[P]>0, L[H]>0, B[P,H]>0, and 0<rho<2 at one x, then inf_(alpha>0) R_alpha=3-rho>1.
```

Thus the optimization over the normalization is exact, not a grid
search. Arb midpoint quadrature with interval second-derivative errors
and explicit tails certifies both endpoint times at `x=25`:

```text
Arb t=0: rho=[0.41 +/- 7.05e-3], inf R_alpha=[2.59 +/- 7.05e-3]
Arb t=1/2: rho=[0.43 +/- 3.08e-3], inf R_alpha=[2.57 +/- 3.08e-3]
```

Independent 80-digit quadrature cross-checks give:

```text
t=0: rho=0.414796465135351145115779744580258235307799494729391878976676896781241963434
     alpha_star=3.47657819605454804669500385759263056682910105814050513244566880659488779575
     inf R_alpha=2.58520353486464885488422025541974176469220050527060812102332310321875803657
t=1/2: rho=0.430761090017545245374029004216750108892360442539619846237137865414549035433
     alpha_star=3.3363259884026044810382951892480128214468378718873867757317456041503241612
     inf R_alpha=2.56923890998245475462597099578324989110763955746038015376286213458545096457
```

The corresponding Arb balls lie strictly above `1`; the displayed
digits are independent cross-checks. No scalar multiple of
the fake-Xi block can support the proposed absolute triangle bound.
This does not reject a sign-aware estimate, because the discarded
mixed and residual terms can cancel.

## Convolution Obstruction

The direct positive-convolution transfer would require

```text
M(u)=Phi(u)/Psi(u) for u>=0
M(u)=integral_R cosh(u*s) dmu(s), mu>=0 finite
```

The kernel asymptotics instead give

```text
Phi(0)>pi*(2*pi-3)*exp(-pi)>4*pi^2*exp(-2*pi)=Psi(0), so M(0)>1
lim_(u->infinity) M(u)=1
```

If a positive finite cosh transform has a finite limit as u->infinity, monotone convergence forces mu(R\{0})=0, so the transform is constant.
Therefore the nonconstant ratio `Phi/Psi` cannot be a positive
cosh-mixture multiplier. This rules out the direct Cardon/Polya
one-factor transfer exactly.

## Live Handoff

Retain the Gasper block only inside a sign-aware or multi-block decomposition. The next admissible target is to construct blocks P_j with a proved common real-zero/interlacing structure and control the signed mixed Laguerre forms after the full modular endpoint cancellation. A scalar absolute budget or a single positive-cosh multiplier is now a closed route.

The useful content of the Gasper model is now sharply located:
its positive square can remain as one block, but the Xi correction
must be handled with its sign and modular cancellation intact.

References: https://arxiv.org/abs/0801.2996 and
https://mathdept.byu.edu/cardon/papers/convolutions_and_zeros.pdf.
