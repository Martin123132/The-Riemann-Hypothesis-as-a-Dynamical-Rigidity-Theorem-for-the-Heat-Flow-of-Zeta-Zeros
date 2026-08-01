# Q207-to-Q208 Forward Adiabatic Successor Certificate

Date: 2026-07-25

Status: rigorous finite forward-successor certificate from the Q207 old bottom edge; not a proof of `Lambda<=0` or RH, and not an all-`j` theorem.

## Direction Audit

The earlier hybrid collar transported from cells on the new Q208
bottom edge. That proves the complete collar, but it is a
reverse-endpoint calibration. The conditional induction lemma instead
asks whether cells already known on the old Q207 edge absorb the
downward heat displacement.

For `x>=38`, the Q207 source stores origin-free cells for

```text
J=16*x^4*H, J_x.
```

On each old cell this certificate applies the exact
positive-determinant transform

```text
F=(1+x^(-4))*J,
F_x=(1+x^(-4))*J_x-4*x^(-5)*J,                  (1)
```

so that `F=16*(1+x^4)*H`, the regular proxy used by the collar
derivative bounds.

## Forward Transport

The old and new times are

```text
t_207=1/1035, t_208=1/1040,
delta_207=1/215280.
```

The 414 old Q207 parent cells cover
`38<=x<=245`. The stored rigorous heat-jet displacement bounds are
partitioned into 303 half-unit panels
and 222 quarter-unit panels. Every
one of the 525 transport panels satisfies

```text
delta_207 sup||(F_t,partial_t F_x)||_2
    < dist(0,C_207,k).                            (2)
```

The largest rigorous ratio is
`[0.763025371788542828925319376890660876392947077093285710521869220214188543194542700548851194754377387126764 +/- 3.80e-106]` on the panel
`[189,379/2]`.

## Finite Composition

The independent compact transversality certificate already covers
`0<=x<=38` for the complete time interval. Equation (2) transports
the Q207 old edge across the new collar on `38<=x<=245`, and the
independent derivative-positive theorem covers the new right strip

```text
[1/1040,1/5]x[245,246].
```

Together with the contact-free Q207 base, these are exactly the old-edge
collar and right-strip hypotheses of the successor lemma. They give a
genuine forward finite Q207-to-Q208 successor and preserve the zero
first-jet degree. This independently recovers the finite Q208
contact-free theorem.

## Proof Boundary

This certificate proves one finite successor using the Q207 old-edge
cells. It does not provide a Q208-to-Q209 estimate, a uniform all-`j`
collar bound, a uniform right-strip cone, a cofinal boundary theorem,
`Lambda<=0`, RH, PF-infinity, or a Clay-prize proof.

Machine-audited files:

```text
work/rh_compute/results/jensen_window_pf_newman_q207_q208_forward_adiabatic_successor_certificate.json
work/rh_compute/scripts/jensen_window_pf_newman_q207_q208_forward_adiabatic_successor_certificate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_q207_q208_forward_adiabatic_successor_certificate.py
```
