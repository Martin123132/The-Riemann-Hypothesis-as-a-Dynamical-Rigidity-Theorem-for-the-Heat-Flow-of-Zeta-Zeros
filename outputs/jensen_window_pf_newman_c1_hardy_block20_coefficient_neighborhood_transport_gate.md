# Hardy block-20 coefficient-neighborhood transport gate

Date: 2026-08-09

Status: rigorous_full_cell_endpoint_transport_closes_finite_block20_scale_other_columns_open; finite block-20 theorem only, not a proof of the complete evaluator or RH

## Cellwise contract

Each of the 374 occupied recursive selector cells is used
as an exact box in `(Phi1, xr=2*Phi2, Phi3, fracL)`.  The fixed selector word is
replayed through the cubic Legendre child map.  The transformed child kernel,
the source-anchored recurrence multiplier, the parent phase, the exact
equation-(69) saddle family, and the retained-decay endpoint functional are
recorded in separate columns.

For the saddle family, `|exp(iu)-exp(iv)| <= |u-v|` gives the source-normalized
explicit bound

```text
|Delta sum I_n| <= tpp_src*J*(r1*N^2/2+r2*N^3/3+r3*N^4/4),
```

where `J` is the fixed number of admitted saddle indices and `tpp_src` is the
accepted binary128 value produced by the evaluator's `8*atan(1)`.  Mathematical
`2*pi` and its independently bounded source-constant transport remain distinct.

For the endpoint functional, every positive numerator coefficient is replaced
by its cell maximum, while every gap and every retained quadratic/cubic decay
coefficient is replaced by its cell minimum.  This gives a rigorous upper
envelope without assuming monotonicity of the full composite expression.

## Finite result

```text
maximum center endpoint budget          4.87965760106136006892983729083960062751062901236036E-3
maximum cell endpoint budget            5.04364075180962585766496963351820390820391778527801E-3
maximum cell/center inflation            1.03360546254569138061884060750227216546123793201016E+0
outputs below the 0.005 endpoint scale   15 / 15
maximum parent phase Lipschitz budget    1.31567078180280431776255181133003553399427810050745E+0
maximum saddle-family Lipschitz budget   1.02028243715753055370660743686206641497000105278261E+2
maximum child-kernel variation           1.70998521521687507629394531250000000000000000000000E+0
maximum multiplier relative variation    1.33542570701701948427107236811008282536158363334553E-1
```

The endpoint figure uses the already certified exact output weights.  The
child, multiplier, saddle, and `fracL` arithmetic columns are not silently
added to that endpoint-only budget; they identify the remaining complete
recurrence obligations.  Each is a separate open budget.

## Boundary

Rigorous coefficient-cell transport for the finite block-20 recursive roster. The fixed-weight endpoint column is separate from child-kernel, multiplier, saddle-family, PSI/ERF arithmetic, Legendre-tail, other-block, outer-Hardy, and height-uniform obligations. No complete evaluator theorem, Lambda <= 0, PF-infinity, RH, or prize-level conclusion follows.
