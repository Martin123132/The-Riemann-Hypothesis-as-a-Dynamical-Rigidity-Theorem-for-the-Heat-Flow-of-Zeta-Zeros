# Source transition and endpoint scheduling gate

Date: 2026-08-12

Status: exact finite-roster decomposition and interval cutoff certificate;
the source transition quadrature remains an approximation and is not promoted
to a theorem-level identity

At the earlier source transition entry,

```text
t_tr=[10000002368.452855115796615822593132249655275276958737919804918428369710000000000000000000000000000000000000000 +/- 1.01e-24].                         (SE1)
```

the pinned implementation changes its ideal odd roster from

```text
R_-={C,C+2,...,B_-},      C=159577, B_-=5122421,
R_+={C+2,C+4,...,B_+},             B_+=5122423.               (SE2)
```

Both rosters contain `2481423` points.  For any assigned
single-alpha summand `F_t(alpha)`, finite cancellation gives the exact identity

```text
sum_(alpha in R_+) F_t(alpha)-sum_(alpha in R_-) F_t(alpha)
  =F_t(B_+)-F_t(C).                                               (SE3)
```

If `T_C(t)` denotes the source transition contribution, the complete idealized
source scheduling change is therefore

```text
Delta_source=[T_C-F_t(C)]+F_t(B_+).                              (SE4)
```

This is a decomposition, not a cancellation theorem.  The first bracket is
the transition approximation defect and the second term is the newly scheduled
upper endpoint.  Interval arithmetic also gives

```text
N_-(B_-)=[621.55615085858927834075625234156252 +/- 6.28e-32],
N_-(B_+)=[621.55590806021966312398073510640727 +/- 6.48e-32],                          (SE5)
```

so the Riemann-Siegel floor is `621` on both sides.  No RS
term absorbs `F_t(B_+)`.

The 64-node source-form quadrature was independently repeated at higher
precision.  As a diagnostic only,

```text
T_C=-0.0002523826968032865590791462135397093267584338551153148935,
F_t(C)=-0.0331488605021633922263382127558995188665993092092979214,
F_t(B_+)=0.0007465905980142806096884382130377093742147968970444658956,
Delta_source=0.0336430684033743862769475047553975189140556722512270724.                        (SE6)
```

The 64- and 96-node values agree beyond 40 decimal places, but (SE6) is not an
interval enclosure of the complete compiled source and is not used as proof.

The primary-source status agrees with this separation.  The 2015 paper calls
the hybrid an approximation, describes direct numerical integration in the
transition region, and notes unproved numerical correlations.  The 2026 paper
states that its hybrid representation is not exact and again sends part of the
transition range to numerical integration.

For theorem-level height transport, fix the analytic endpoint at
`B_fix=5122423`.  The fixed chart `159577,159579,...,5122423` has
`2481424` points and obeys

```text
source-pre  = fixed-chart-F_t(B_+),
source-post = fixed-chart+[T_C-F_t(C)].                           (SE7)
```

Thus the arbitrary implementation threshold and moving endpoint belong in an
explicit source-approximation ledger.  They must not define the exact
`C -> C+2` selector theorem.  The exact selector change remains the one-cell
Poisson identity at `t*=pi C^2/8`, with fixed `B_fix=5122423`.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_source_transition_endpoint_scheduling_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_source_transition_endpoint_scheduling_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_source_transition_endpoint_scheduling_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_source_transition_endpoint_scheduling_gate.py
```

This gate does not bound the complete source approximation error, the
quantitative 399-mode chart, `T_upper`, `Lambda<=0`, RH, or a prize-level
conclusion.
