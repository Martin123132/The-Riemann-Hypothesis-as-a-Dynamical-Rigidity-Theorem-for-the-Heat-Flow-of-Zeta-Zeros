# Source-aligned hybrid-error ledger at t=10^10

Date: 2026-08-10

Status: exact finite interval identity validated; not a proof of the height-uniform signed theorem

The primary-source route uses the endpoint cutoff from equations (124)--(127)
and the Fortran control flow.  The diagnostic nearest-root midpoint is not used
in this ledger.

Write `Z_src=L_src+Q_src+rho_add`, `Z_exact=L_exact+U_exact`, and let
`T_upper` be the exact classical upper-main sum induced by the 36 source
cutoffs.  For block increments

```text
e_k = (Q_k-Q_(k-1)) - (T_k-T_(k-1)),  k=0,...,35,
```

the independently checked finite identity is

```text
Z_src-Z_exact
  = (L_src-L_exact) + sum_(k=0)^35 e_k
    - (U_exact-T_upper) + rho_add.
```

All 15 output identities and all 540 block
rows are interval certified.  At the central output,

```text
E_total    = [-0.0069830167012009191664618470140087218990724443773612412740312566610214870418568438133304931550212789889083421600 +/- 1.14e-100]
E_lower    = [-7.8608916441819954504123019651585832479714632820450071521859677983236566209728860491955219480001263547591187164e-17 +/- 2.35e-98]
sum e_k    = [-0.0069830180084699455319725803829133704069288214989573701554341205283354305539931196430885535174224228367531020376 +/- 9.46e-96]
R_upper    = [-1.3072691049744271751888591526308660287074286085960356843173854653718142590663242700912616358030642389169740701e-9 +/- 2.41e-97]
rho_add    = [1.0000000000000000000000000000000000000000000000000000000000000000000003153063508356442351959090963753340375797e-41 +/- 7.08e-111]
sum |e_k|  = [0.14466461698436266881072695415053723304287230319610991863848904276790811213302768585631571385112232790152607378 +/- 9.46e-96]
```

The triangle sum rejects `0.005` on all 15
outputs.  Its range is

```text
[0.13303034632772588889577778961552984056331811542139107007600333460949171898645177828910589990288151071169985270 +/- 1.56e-95]
to
[0.15526911168332570866524206560505587017936704396427745771141840173422307665265110948460749888620594547545527539 +/- 1.56e-95].
```

Only [0.046676278888853190883810738547243994065873870564711966046631841291001587350584384053636955042345017900784087929 +/- 1.06e-94] to
[0.049885739432795757895306957463493408724662413234293278453121755392441111230634118209931142229329333425497311454 +/- 1.23e-94] of that absolute mass remains in the
signed block sum.  Equivalently, the certified cancellation fraction is
[0.95011426056720424210469304253650659127533758676570672154687824460755888876936588179006885777067066657450268855 +/- 1.23e-94] to
[0.95332372111114680911618926145275600593412612943528803395336815870899841264941561594636304495765498209921591207 +/- 1.06e-94].  A proof that takes absolute
values block by block therefore destroys the mechanism the data says matters.

The first genuinely unbounded obligation is

```text
E_hyb_main(t)=sum_k e_k(t)=Q_src(t)-T_upper(t).
```

It needs a constant-bearing, height-uniform signed estimate using the published
endpoint convention.  The existing corrected internal model is useful
diagnostic machinery, but it is not this source-aligned outer error and cannot
be substituted for it.

Proof boundary: rigorous finite decomposition for one fifteen-point window
near `t=10^10`.  It proves no height-uniform hybrid estimate, no cofinal
Jensen/PF theorem, no `Lambda<=0`, no RH, and no prize-level conclusion.
