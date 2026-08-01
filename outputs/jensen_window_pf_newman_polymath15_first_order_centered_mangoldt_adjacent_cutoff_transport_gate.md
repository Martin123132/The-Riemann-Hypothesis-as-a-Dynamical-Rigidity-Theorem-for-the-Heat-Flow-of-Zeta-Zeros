# Centered Mangoldt Adjacent-Cutoff Transport Gate

Date: 2026-07-30

Status: exact componentwise adjacent transport and square/wing transition audit;
not a proof artifact. This is not a signed centered-jet lower bound, Abel
gap, contact exclusion, `Lambda<=0`, PF-infinity, RH, or a Clay-prize result.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_mangoldt_adjacent_cutoff_transport_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_mangoldt_adjacent_cutoff_transport_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_mangoldt_adjacent_cutoff_transport_gate.py
```

## Pi Provenance

The `pi` in the endpoint recurrence and `a^2=x/(4*pi)+t/16` remains the
completed-zeta and Riemann-Siegel normalization. The Mangoldt cutoff,
integer square root, and row transfer introduce no new `pi`.

## Component Jumps

```text
At one fixed (t,x), eta, |f_1|, s_*', and L_a=log(a) are cutoff independent. For n=N+1, Delta e=kappa_NJ_a/|f_1|, Delta g=kappa_N[J_(a,x)+mu_aJ_a]/|f_1|, Delta W_0=z_n+Delta e, and Delta Z_Lambda=log(n)z_n.

For E_N=g_N-s_*'L_a e_N, Delta E_N=Delta g-s_*'L_a Delta e.

The endpoint-value jump cancels exactly: Delta mathfrak B_N=Delta E_N-s_*'Delta Z_Lambda+s_*'L_aDelta W_0=Delta g+s_*'log(a/n)z_n.
```

The two appearances of `Delta e` cancel before an absolute value is taken.

```text
With ell=log(n/a), z_n=f_n/|f_1|, Delta mathfrak B_N={-s_*'ell f_n+kappa_N[J_(a,x)+mu_aJ_a]}/|f_1|=Delta W_A. Thus Re(Delta mathfrak B_N)=Delta A_a/|f_1|.
```

Thus the new centered coordinate has exactly the old physical adjacent
jet; it is not a second transport problem.

## Projection Semantics

```text
The certified adjacent theorem gives |Re(Delta W_0)|<2200exp(-5L/4) and |Re(Delta mathfrak B_N)|<10000exp(-7L/4).

No complex absolute-value bound for Delta W_0 or Delta mathfrak B_N at the displayed real-projection scales is proved. The exact complex cancellation must precede real projection and absolute values.
```

The real trace is the proof-facing quantity. A small complex adjacent
jump is neither needed nor claimed.

## Ordinary Cutoffs

```text
If n=N+1 is not a square, D=floor(sqrt(N))=floor(sqrt(n)); the complete square is unchanged and the wing gains exactly the factor pairs of n. Hence Delta Z_square=0 and Delta Z_wing=log(n)z_n.
```

## Perfect Squares

```text
If n=r^2, D changes from r-1 to r. Put R_r=sum_(d<r)[Lambda(d)+Lambda(r)]z_(dr), D_r=Lambda(r)z_(r^2), and A_n=sum_(d|n,d<r)[Lambda(d)+Lambda(n/d)]z_n. Then Delta Z_square=R_r+D_r and Delta Z_wing=A_n-R_r.

At n=r^2 the transferred row R_r cancels, and D_r+A_n=log(n)z_n by the symmetric Mangoldt divisor identity. Perfect-square changes introduce no residual term.
```

The square/wing boundary moves at a perfect square, but the movement is
an internal transfer. It leaves no extra term in `Delta Z_Lambda`.

## Exact Audit

The independent audit checks `79` successive
cutoffs, comprising `71` ordinary and
`8` perfect-square transitions.
It uses rational-complex carriers and formal prime-log coefficient
vectors, with `0` mismatches.

## Route Decision

```text
The centered jet and balanced hyperbola are now compatible with adjacent cutoffs. Any Type-I/II or Vaughan estimate must retain the ordinary factor-pair insertion, the perfect-square row transfer, the recurrent endpoint slope, and real-projection semantics before absolute values.
```

The next estimate may now be written on the balanced hyperbola without
silently changing the object at a cutoff. It must still be tested at
`W_0=0`, `mathsf_X=0`, `H_a=0`, `q=1`, prime edges, ordinary transitions,
perfect-square transitions, and adjacent real projection.

## Boundary

This gate proves the centered component jumps, exact endpoint-value
cancellation, recurrence match, inherited real-projection bound, ordinary
wing insertion, perfect-square row transfer and cancellation, and finite
coefficient audit. It proves no signed lower bound, Abel-scalar gap,
horizontal successor winding cap, contact exclusion, Q209, cofinal
descendant theorem, `Lambda<=0`, PF-infinity, RH, or prize-level result.
