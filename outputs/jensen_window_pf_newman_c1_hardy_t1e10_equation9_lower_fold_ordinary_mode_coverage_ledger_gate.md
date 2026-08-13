# Lower-fold/ordinary mode-and-height coverage ledger

Date: 2026-08-13

Status: exact ownership ledger complete; not a proof of the quantitative ordinary-Morse splice

At the saved height, the source classical target is exactly `622..39894`.
The finite-Poisson geometry partitions it as

```text
622..39852       lower ordinary interior (39231 modes),
39853..39894     lower endpoint transition (42 modes).
```

The 399-event atlas supports `39695..40093`.  Its exact overlap with the
lower ordinary interior is `39695..39852`
(158 modes).  A disjoint target ownership suitable for a
matched proof is therefore

```text
622..39694       ordinary core (39073 modes),
39695..39894     fold-owned target (200 modes).
```

The overlap `39695..39852` is retained as a validation region, not counted
twice.  Modes `39895..40093` belong to exact Poisson/selector bookkeeping but
not to the source `T_upper` target.

The 399 closed event buffers have radius `pi/16` and are pairwise disjoint.
They do **not** cover the selector cell by themselves.  Their complement has
400 open ordinary-height corridors;
the smallest inter-event corridor has exact width `7*pi/8`, enclosed by
`[2.74889357189106908365481296036956502367252322445321759335307651826933935550043287272786082 +/- 1.00e-80]`.  Assigning every event
face to its event buffer makes the 399 buffers plus 400 corridor interiors a
disjoint and exhaustive height ownership table.

The first missing normalized inequality is

```text
sup_n |K_m(tau_n+pi/16)-S_m(tau_n+pi/16)| <= eps_Morse,
```

where `K_m` is the exact grouped finite-Poisson mode and `S_m` is the ordinary
Morse approximation with its classical carrier.  The corresponding physical
unit is `sqrt(2)/pi` times the normalized value.  The already certified fold
error leaves local prototype headroom
`[1.68476298697649586850672926815291725591022432864404247317448332382987150000000000000000000e-5 +/- 3.39e-21]` normalized and
`[7.63109302667009851879293693481828311436503431622493713695322190907971470000000000000000000e-6 +/- 1.70e-21]` physical.  These are local
splice diagnostics, not a complete aggregate budget.

The next aggregate obligation is to bound the ordinary remainder on
`622..39852`, preserve the 158-mode overlap before taking absolute values,
and then prove that the disjoint ordinary-core plus fold-owned target equals
the source `T_upper` in the exact phase normalization.

Proof boundary: exact integer ownership, event ordering, and height coverage
only.  No ordinary stationary-phase remainder, quantitative Airy/Morse
overlap, complete `T_upper`, `Lambda<=0`, PF-infinity, RH, or prize-level
conclusion is proved.
