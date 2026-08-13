# Exact paired half-domain target residual

Date: 2026-08-13

Status: exact-lemma certificate; not a proof of quantitative pair-block bounds

After the odd-alpha reflection, define

```text
K_m(t)=integral_0^(1/2) W_t(x)I_m(x)dx,
K_H(t)=integral_0^(1/2) W_t(x)H_x dx.                 (PR1)
```

For every symmetric cutoff `M>=39895` the source half-sum is

```text
Q_half,M=K_H+K_0+sum_(m=1)^M (K_m+K_-m),
Q_main,M=2 Re[exp(-i*pi/8)Q_half,M].                  (PR2)
```

Let `tau_m` be the already normalized complex classical carrier, so that

```text
T_target=2 Re sum_(m=622)^39894 tau_m,
tauhat_m=exp(i*pi/8)tau_m.                            (PR3)
```

Pure finite algebra now gives the cancellation-preserving residual

```text
Q_main,M-T_target=2 Re exp(-i*pi/8) {
 K_H+K_0
 +sum_(m=1)^621                 (K_m+K_-m)
 +sum_(m=622)^39694(K_m+K_-m-tauhat_m)
 +sum_(m=39695)^39894   (K_m+K_-m-tauhat_m)
 +sum_(m=39895)^M             (K_m+K_-m) }.   (PR4)
```

The four ranges contain respectively `621`, `39073`,
`200`, and `M-39894` positive modes.  Together with
`K_0` they contain exactly `M+1` nonnegative indices, while every negative
index occurs exactly once in its symmetric pair.  The target is subtracted
exactly once over its `39273` modes.

This form preserves the endpoint-current cancellation:

```text
K_m+K_-m=integral_0^(1/2)W_t(x)(I_m+I_-m)dx,
I_m+I_-m=[-2 Delta f_x'+R_m+R_-m]/(2*pi*i*m)^2.       (PR5)
```

The `1/m` current is absent before any absolute value.  The proved symmetric
Poisson interchange supplies the limit `M->infinity`, so (PR4) is also the
exact limiting residual when its outer block is understood symmetrically.

Two guards are essential.  First, `tauhat_m` in (PR3)--(PR4) is only an
algebraic allocation of the real target projection.  It does not assert that
the full-line Gamma bulk has equal modewise `x<1/2` and `x>1/2` pieces.
Second, (PR5) must retain the canonical quadratic chirp: discrete odd-square
parity at `x=1/2` does not make `I_m+I_-m` vanish there.

The tangent B profile therefore cannot be promoted by itself.  Its grouped
`q`-density must be combined with the low pairs and with the nonpositive and
outer-positive completion in (PR4).  The ordinary block remains separate
from the 200 fold-owned target pairs, matching the certified ownership
ledger.

Pi provenance: `pi/8` comes from exact odd-square Kummer reflection; all
other `pi` factors come from the equation-(9) phase or Fourier-Poisson
characters.  No fitted or geometric constant is inserted.

Proof boundary: exact finite and limiting paired residual normal form only.
No quantitative block estimate, nonlinear B-crossing theorem, A-fold splice,
complete `T_upper`, height-uniform theorem, `Lambda<=0`, PF-infinity, RH, or
prize-level conclusion is proved.
