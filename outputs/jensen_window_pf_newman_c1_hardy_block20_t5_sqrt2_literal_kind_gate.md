# Hardy block-20 t5 `sqrt(2.0)` literal-kind gate

Date: 2026-08-06

Status: rigorous finite exact-argument correction for the default-real square-root factor; not a proof of the contour estimate or RH

## Source issue

Both intrinsic real-erfc arguments in `t5` contain the unsuffixed expression
`sqrt(2.0)`.  Under the admitted gfortran build, `2.0` is default real and the
square root is therefore first rounded to binary32:

```text
bits       = 0x3FB504F3
exact      = 11863283 / 8388608
decimal    = 1.41421353816986083984375000000000000000000000000000E+0
```

It is then promoted into the binary128 expression.  This is distinct from the
principal positive mathematical root `sqrt(2)`.  No circle or polygon supplies
this constant: it is the algebraic root of `x^2=2`.  The source's `sqrt(pi)`
factor remains fixed, so this gate isolates only the literal-kind error.

For each of the 913 admitted retained calls, write the source argument as
`x_src` and reconstruct its pre-`sqrt(2)` factor `B` from exact binary128
inputs.  Define

```text
delta x = B * (sqrt(2) - sqrt2_binary32),
x_corr  = x_src + delta x,
delta q = W * (erfc(x_corr) - erfc(x_src)).             (1)
```

Here `W` is the independently emitted exact binary128 source weight.  The
previous intrinsic-erfc gate corrected function evaluation at `x_src`; (1)
corrects the generation of the argument and does not count that residual twice.

## Result

```text
maximum source-model/emitted argument gap       <= 8.58286045401223318363760908112818905616834268487077E-33
maximum literal-kind argument shift             <= 5.24634741956476672149566508307948573194842907486158E-7
maximum erfc-value shift                        <= 8.28139260697765971686085943033039417379886144581537E-9
maximum one-call q shift                        <= 3.74455843536869250487608333967682430390781023092894E-8
maximum correlated per-chain q shift            <= 3.74455843053224604527289594075761568462959794392346E-8
minimum residual after sqrt(2) correction       >= 4.66664105707053750046113925417788214645456571865052E-3
minimum residual after sqrt(2)+finite-ip stack  >= 4.57823940647513633873846091754273294458559986416668E-3
stacked residuals excluding zero                   374 / 374
```

The stacked residual is `R_q-delta_q_sqrt2-T_ip`, where `T_ip` is the separately
admitted displayed-W1 finite-index completion.  Every operation is repeated at
180 and 260 decimal digits with overlap required.

## Boundary

This closes one source literal-kind normalization on one finite block.  It does
not prove that the displayed W1 approximation equals the exact equation-(69)
integrals, bound higher saddle phase or vertical-leg remainders, reconstruct
the unavailable Maple code, control W2--W4 or outer Hardy errors, establish a
height-uniform recurrence, or prove `Lambda<=0`, PF-infinity, RH, or a
prize-level theorem.
