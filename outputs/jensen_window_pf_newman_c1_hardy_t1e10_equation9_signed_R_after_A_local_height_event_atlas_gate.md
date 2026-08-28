# Signed-closure local height event atlas

Date: 2026-08-28

Status: maximal same-roster cell and all active discrete guards certified;
height-uniform packet bounds and sign transport remain open

Keep the exact odd window fixed at

```text
A=159577, B=5122421, M=2481423.
```

For a fixed odd label `alpha`, a saddle reaches mode `m` at

```text
tau_m(alpha)=pi*m*(alpha-2m).
```

The nearest two active events around `t0=10^10` are therefore

```text
t_-=pi*39852*(159577-2*39852)
   =[9999999793.16411593523938104429435216216653734175894249111793767794172126423867869493726330 +/- 4.40e-81],

t_+=pi*39936*(159577-2*39936)
   =[10000000057.0578988367820130751563963576447796143211720386268266398370670180952568229788192 +/- 1.09e-80].
```

Hence the maximal open interval containing the saved height on which the
same `A` transition roster is valid is

```text
t_- < t < t_+,
t0-t_-=[206.835884064760618955705647837833462658241057508882062322058278735761321305062736695607437 +/- 1.28e-88],
t_+-t0=[57.0578988367820130751563963576447796143211720386268266398370670180952568229788191891474885 +/- 4.41e-89].
```

The walls are not decimal guesses.  Exact complementary-root identities give

```text
t=t_-: y_-(A)=39852,   y_+(A)=39936.5=U,
t=t_+: y_-(A)=39852.5=a, y_+(A)=39936.
```

Inside the open cell the complete active ownership data stay fixed:

```text
ordinary modes:       622..39852 (39231 modes),
transition modes:     39853..39936 (84 modes),
classical split:      39853..39894 | 39895..39936,
B lower/upper floors: 621 and 2560588,
upper complement:     230 stationary labels,
                      2561211.5..2561440.5.
```

The rational boundaries are proof choices, not height events:

```text
L=621.5,
B_core=621+9/16=621.5625,
C_far=621+11/16=621.6875.
```

Their exact Fourier orders remain 2 and 16 independently of `t`.
`240=15*16` and `752=47*16`, so both grouped endpoint cancellations remain
exact.  On the whole closed event cell, interval arithmetic also proves

```text
finite 752-launch margin  > 0: [30815.1704312692524963299920423911992806154103163862600922584533691406250000000000000000000 +/- 36.7],
remote 752-launch margin  > 0: [9197.87508098997582464077110708799978056049440056085586547851562500000000000000000000000000 +/- 35.4],
previous grid margin      < 0: [-182218.903647299554239258190919091562331999512025504373013973236083984375000000000000000000 +/- 19.3],
previous count margin     < 0: [-7260.98895762644301386020545541200021943950559943914413452148437500000000000000000000000000 +/- 34.3].
```

Thus the chosen far launches remain admissible and their immediate tested
predecessors remain inadmissible everywhere in the cell.  The dormant
unequal-truncation donor split also stays `622..860 | 861..39936`, but it is
not part of the final signed route.

What this certifies is the discrete and exact-algebra layer.  With `A,B`
fixed, the finite RSI identity, branch correction, Dirichlet-cell partition,
pole-safe contour join, entire-kernel finite difference, and complementary
half-lattice identities all retain the same form in this cell.

What it does **not** certify is just as important: the complex balls for
`V_L`, `C_U`, `T_upper`, `T_lower`, `O_join`, `T_A^join`, `J_Z`, the signed B
trace, or the outer allowance have not yet been transported in height.  The
saved inequality `Q_K-T<0` therefore remains a theorem only at `t=10^10`.
Crossing either wall requires the already-known exact one-mode ownership
handoff plus matched interval estimates; this atlas does not silently cross
it.

Next obligation: derive interval-valued transport bounds for the joined
complex packets, signed B trace, and outer term on a first certified subcell,
retaining the exact cancellations before every norm.
