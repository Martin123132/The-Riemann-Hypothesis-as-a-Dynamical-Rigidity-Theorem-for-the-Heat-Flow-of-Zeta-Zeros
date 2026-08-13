# Hardy block-20 exceptional-mode recombination full gate

Date: 2026-08-07

Status: rigorous all-374 finite-input identity and quadrature certificate; not a height-uniform bound and not a proof of RH

The four-branch exceptional-mode construction was replayed serially at 70 and
110 decimal digits on all 374 saved calls.  Every completed call was fsynced
to an append-only cache.  The independent paper-model identities

```text
W2 = M_b,lin(N)-M_b,quad(N),
W3 = conjugate(M_c,lin(0)-M_c,quad(0)),
W4 = Z_quad(0)+i*exp(-2*pi*i*F(N))/(2*pi*F'(N))
```

and the resulting singularity-free recombinations all enclose equality.  The
quadratic modes use the paper's global `Phi2`, while the exact modes retain
the complete cubic phase.  Thus no `1/delta` or `1/eta` term remains in the
objects that now require uniform bounds.

```text
W3 calls                              = 165
W4 calls                              = 209
maximum paper-model identity gap      <= 1.85675202920745081348041965437585405169225241384569E-106
maximum complete recombination gap    <= 5.23698368491970890712533526212268952804151922464371E-14
maximum exceptional-mode tail         <= 5.86083000123935813095248830503277364858077013704702E-31
maximum four-piece triangle/correction ratio <= 6.48364866601286819371125038802718964175668274915030E+3
```

The worst triangle ratio occurs at chain
`361`.  It
shows that removing the formal poles does not remove the need to preserve
coupling between the upper and lower endpoint families.  The ratio is a
finite-roster diagnostic, not a fitted theorem constant.

This gate does not provide a height-uniform finite-difference majorant,
recursive accumulation, outer Hardy control, `Lambda<=0`, PF-infinity, RH, or
a prize-level conclusion.
