# Quartic Outer-Branch Alternate Length-13 Survivor Gate

Date: 2026-07-25

Status: exact alternate length-thirteen survivor and fixed-tail length-fourteen obstruction; not a proof of a uniform finite-tail theorem, `Lambda <= 0`, or RH.

```text
work/rh_compute/results/jensen_window_pf_quartic_outer_branch_alternate_length13_survivor_gate.json
python work/rh_compute/scripts/jensen_window_pf_quartic_outer_branch_alternate_length13_survivor_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_quartic_outer_branch_alternate_length13_survivor_gate.py
```

## Exact Alternate Tail

Keep only the exact `x_2,...,x_5` outer-contact data from the
quartic nonpromotion gate. Define

```text
G_j=d_(j+2)^2-x_(j+2)^2*d_(j+1)*d_(j+3)
C_k=G_(k-4)^2/(x_(k-2)^3*G_(k-5))
G_(k-3)=y_k*C_k.
```

Use the rational coordinates

```text
y_6=37/100
y_7=69/70
y_8=13/15
y_9=7/10
y_10=7/8
y_11=9/20
y_12=1/9
y_13=1/2
```

They give the following exact rationals, displayed only by
certified decimals and hashes because the last fractions contain
thousands of digits:

```text
x_2 ~= 0.9661855600000000000000  sha256=b1f86a9049fd276e84195aee1f7863878317cf5ecfa6595abf38b44d9b2874d2
x_3 ~= 0.9723789559857728810569  sha256=1cb61a8acfc8ea4d18c71e3ab5ad279f73db7ca56e6e7a693077b0cbeb43cbd0
x_4 ~= 0.9772556571312076870058  sha256=4355001a931cdfc4bd92bc348371756757d4c7cc3ab0cfd81aff445576065bcd
x_5 ~= 0.9812000000000000000000  sha256=738abe731814bcedc7f70298453c98b459bb7884c48b28d770e09f4adfca317a
x_6 ~= 0.9840907545684242204979  sha256=525d1ff26bc15db3b4332090aac1f39d6bc0f2bbd232c16f08c7511c6052f349
x_7 ~= 0.9861671591927467742128  sha256=5300aeec28b195d0d50f7797090317a3afc1ee4319f7664f57bafc741df6cc53
x_8 ~= 0.9876510450612590914713  sha256=82820007a1618fcbed54bcb509b9024de6bb2c10262ca9f4e1011f34e217248d
x_9 ~= 0.9887017877414108395887  sha256=d6592d4e6b73989cf9e5c17a44fec250d950f05f0bb22178d0250dcacf05ebf4
x_10 ~= 0.9894261038014610393844  sha256=8c138d5a38b8e902995f3aaaf137e9bc7af8fbd9bcbab44e4ee52656f510ddb1
x_11 ~= 0.9898913826267730819105  sha256=27e1a64af528385634a723a9975a7ba8284cd558d6776b8c726275674de8609a
x_12 ~= 0.9901378098880659919795  sha256=cacdba5f7a722164ea0c87d639b8212245453fe03be73f2e114bbb1f37b7fa91
x_13 ~= 0.9901856014084241800054  sha256=3fcc92ba4f581eaca9bd479297a2dea2cbefa2fd353a5f8a62e5b4ec72645ed5
```

All contractions increase, remain below one, clear the
pointwise walls, have increasing scaled defects, satisfy the
reciprocal-square-root increment bound, and satisfy every
adjacent cubic frontier.

## Finite Signed Layers

Exact direct enumeration over `A_1,...,A_14` gives

```text
order 2: 364 strict signed minors
order 3: 715 strict signed minors
order 4: 792 strict signed minors
```

The checker rebuilds the rationals, hashes every determinant,
and verifies all signs rather than trusting floating point.

## Survival And Death

For the increasing compatibility margin

```text
Delta_k=C_k-[d_(k-1)^2-x_(k-1)^2*d_(k-2)*d_(k-1)],
Delta_13 ~= 4.736242267635457629108e-7 > 0,
Delta_14 ~= -0.000001421531812277118711728 < 0.
```

Thus the same outer quartic contact does have a different tail
surviving through `x_13`; the earlier length-13 obstruction is
genuinely prefix-specific. This selected alternate tail cannot
reach an increasing `x_14` with the next order-four sign.

A bounded 800,000-point scout plus three deterministic
optimizations also found no positive `Delta_14`, but that
observation is numerical evidence only. It is not promoted to a
uniform theorem.

## Corrected Handoff

The live exact target is now a uniform semialgebraic statement
over all admissible earlier corridor choices, or another
explicit tail that survives the length-14 boundary. Neither
finite outcome supplies PF-infinity, an Xi-specific bridge,
`Lambda <= 0`, or RH.
