# Newman Q208 Top Phase-Cell Certificate

Date: 2026-07-25

Status: rigorous complete fixed-time top certificate.
This is not a proof of Q208, `Lambda<=0`, RH, or the Clay prize.

```text
work/rh_compute/results/jensen_window_pf_newman_q208_top_phase_cell_certificate.json
work/rh_compute/results/jensen_window_pf_newman_q208_top_phase_cell_certificate.jsonl
python work/rh_compute/scripts/jensen_window_pf_newman_q208_top_phase_cell_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_q208_top_phase_cell_certificate.py
```

## Contract

The six-term, 352-bit fixed-time engine certifies
`F_t+iF_t'` at `t=1/5` on `0<=x<=246`, using the same
axis-safe proxy and raw omitted-tail theorem as the bottom chain.

## Progress

```text
panels=492/492
certified cells=492
unresolved=0
branches={'value': 319, 'derivative': 173}
minimum closest-point norm lower=[3.09024996991711826010647932123357337545865067077163797156897935865388249628291778385662831894748761048840827080658661707e-29 +/- 1.33e-134]
maximum subdivision depth=0
stop reason=complete
```

## Exact Open Phase Chain

```text
cells=492
witnesses=493
forward positive-ray crossing count=-19
```

The Q208 boundary traverses this top chain in reverse.
Its endpoints must still be matched to the transformed
right edge and positive-real symmetry axis before the
closed winding is evaluated.

## Scope

The bottom chain, right/top corner witnesses, and cyclic exact
winding are separate obligations. No incomplete chain is
promoted to Q208 or a cofinal theorem.
