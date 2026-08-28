# Global R_Dir projector ownership ledger

Date: 2026-08-23

Status: exact saved-height ownership and post-A target certified; joined
remainder still open

For a common finite cutoff `M>=39936` and the inherited Abel weight, the
positive/negative pair identity is

```text
I_m+I_-m-chi_T(m)P_m
 =(1-chi_T(m))P_m+P_-m+A_m+A_-m+B_m+B_-m.          (OL1)
```

The six positive-label classes are `1..621`, `622..39852`,
`39853..39894`, `39895..39936`, and `39937..M`, together with the zero
mode.  Their counts satisfy

```text
621+39231+42+42+(M-39936)=M.                        (OL2)
```

The certified A extraction removes `A_m+A_-m` on all 84 A-window modes and
also removes `P_m` on the 42 outside-target modes `39895..39936`.  It removes
no `P_-m` or B endpoint atom.  Thus

```text
R_KGamma=E_Btr,win+E_outer+A_transition+R_after_A,
R_Dir=A_transition+R_after_A,                        (OL3)

A_transition=[-0.036674124797773646717966091098920512012405425979411807743120038220889293386300000 +/- 3.04e-12].            (OL4)
```

The two certified one-sided continuation targets are

```text
R_after_A < 0.0368147039947
  ==> R_Dir < 0.00014057919999999995,

R_after_A < 0.0368061039947
  ==> R_Dir < 0.00013197919999999995.        (OL5)
```

The first threshold is about
`[261.87874162536145527885406424330743056371932807536 +/- 8.57e-49]` times the original
`R_Dir` working scale.  This is only target arithmetic.  It does not say the
remaining joined current is small.

The half-current, zero mode, negative full-line bulk, remaining endpoint
atoms, and remote-positive modes stay inside one `R_after_A`.  The exact
double-Weber reconstruction shows that the zero, direct, modular, and
endpoint-half sectors are one copy of the finite source, while symmetric
pairing cancels the apparent one-over-`m` current.  Taking separate norms on
those labels would destroy the proved cancellation.

The B trace window is already outside `R_Dir`.  The completed analytic B
outer block is also outside it.  The finite completed B block is not a third
numerical contribution: it is the common Dirichlet/Abel representation from
which the A positive-bulk modes `39895..39936` are now reassigned exactly.

Pi provenance: every `pi` in (OL1)--(OL5) is inherited from the equation-(9)
Fourier kernel, the Gaussian Abel regulator, and the previously certified
Gamma/A endpoint identities.  No geometric or fitted occurrence is added.

Proof boundary: exact mode/component ownership, no-overlap guards, and
post-A sufficient target arithmetic at `t=10^10` only.  No bound for
`R_after_A`, complete `R_Dir`, `Q_K-T`, all-height theorem, `Lambda<=0`,
PF-infinity, RH, or prize-level conclusion is proved.
