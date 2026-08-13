# Hardy q Selector-Cell Gate

Date: 2026-08-06
Status: rigorous local selector cells; not a proof and physical coverage/analytic q open

## Cell Contract

Each cell is a nonzero-radius box in `(a1, xr=2*a2, a3, fracL)` with fixed discrete labels. Arb interval evaluation preserves the transformed-child integer/parity selectors, q sign and loop branches, all six `PSI` branch paths, both `ERF` branch paths, denominator nonvanishing, and the cubic stationary branch.

All `64` saved chains admit such a cell. The largest shrink count is `11`; the halving histogram is `{"10": 2, "11": 6, "5": 24, "6": 19, "7": 9, "8": 2, "9": 2}`.

## Tightest Certified Margins

- Child linear selector: `8.547542970516879599739464984702844682173638137E-4`
- Child quadratic selector: `1.470126303399505237160160033779479896970093617E-3`
- `PSI` branch boundary: `1.463522510131825601699615809066267978961124023E-3`
- `ERF` branch boundary: `3.022609278559684753417968750000000000000000000E-3`
- Kernel denominator: `9.996752090179856404359680414684543839434763753E-1`
- Stationary discriminant: `9.987014683902089018374681472778320312500000000E-1`

## Critical Scope Boundary

No physical rae/t trajectory is yet enclosed in these derived-input boxes, and no source PSI/ERF approximation error or analytic q value is bounded.

The next proof obligation is to map actual `rae` intervals into these boxes. Only then is it meaningful to replace the source `PSI` and six-term `ERF` approximations with rigorous special-function balls and enclose the correlated sum `q=t3+t5+conjg(t1+t2+t4)`.

## Proof Boundary

This gate proves sixty-four nonzero-radius local boxes in derived q-input coordinates on which the recorded integer, parity, orientation, PSI, ERF, loop, denominator, and cubic stationary selectors remain fixed. It does not prove that a physical coefficient trajectory enters or remains in any box, bound source special-function approximation errors, enclose analytic W1-W5 or qq, control the outer Hardy representation, evaluate a physical carrier, prove a determinant sign, Lambda<=0, PF-infinity, RH, or a prize-level conclusion.
