# Phase-cell 12.5--13 anchored Lipschitz gate

Date: 2026-08-28

Status: rigorous two-configuration closed-cell sign certificate; not an
all-height theorem and not a proof of RH

Put `a=10^10+12.75` and split

```text
I_L=[10^10+12.5, 10^10+12.75],
I_R=[10^10+12.75, 10^10+13].                    (PC1)
```

Both intervals lie strictly inside the same fixed-roster event cell.  The
joined packet therefore has no roster jump on either half.

## Reused grouped atlas

For the grouped 752-label collar, a panel atlas is evaluated once at `a`.  If
`h=t-a`, the exact height transport is

```text
F_t(y)=F_a(y) exp(-i h log y).                  (PC2)
```

Factoring the common rotation at log-center `s_0` and one panel rotation at
`s_j` leaves the rigorous errors

```text
|E_j| <= |h| epsilon_j B_j,
|E_j^log| <= log(C)|h| epsilon_j B_j.           (PC3)
```

Here `epsilon_j` is the panel log half-width and
`B_j=752(y_right-y_left)/sqrt(y_left)`.  Production uses 192 forward panels
at 384 bits; the independent path uses 256 reverse panels at 448 bits.  No new
occurrence of pi is introduced in this phase-cell argument.

## Anchor and derivative norms

The production anchor upper endpoint is

```text
Q(a) <= -0.001981574179808376356959342956542968750000. (PC4)
```

The independent anchor upper endpoint is

```text
Q(a) <= -0.001981309902475913986563682556152343750000. (PC5)
```

The transported complete derivative enclosures give

```text
production:  sup_I_L |Q'| <= [0.005257618472683134314660212298670289683352480764044756256 +/- 1e-62]
             sup_I_R |Q'| <= [0.005314597268124998663333858992045042182533146762901665738 +/- 1e-62]
independent: sup_I_L |Q'| <= [0.005249439532462929492406776448526857066164980764044756256 +/- 1e-62]
             sup_I_R |Q'| <= [0.005307681886961290729013873075404611666475726897140852145 +/- 1e-62]. (PC6)
```

For every `t` in either half, `|t-a|<=0.25`, so the real mean-value theorem
gives

```text
Q(t) = Q(a) + integral_a^t Q'(u) du
     <= Q(a) + 0.25 sup_I |Q'|.                 (PC7)
```

The resulting worst-case upper endpoints are

```text
production left:  [-0.0006671695616340400646154893809460407078239891839888109360 +/- 3e-63]
production right: [-0.0006529248627735739774470777076023525830288226842745835655 +/- 3e-63]
independent left: [-0.0006689500193530761861043874421619182407829735589888109360 +/- 3e-63]
independent right:[-0.0006543894307284858769526132854424795907052870257147869637 +/- 5.01e-59]. (PC8)
```

All four are strictly negative.  Since `I_L union I_R` is the full closed
interval in (PC1),

```text
Q_K(t)-T(t) < 0
for every t in [10^10+12.5, 10^10+13].          (PC9)
```

This proof does not assume or prove uniqueness of the stationary point seen
near offset 12.76.  Any additional hidden critical points are absorbed by the
rigorous derivative norm.  The certificate covers only this one half-unit
phase cell and makes no claim about the rest of the upper event cell.
