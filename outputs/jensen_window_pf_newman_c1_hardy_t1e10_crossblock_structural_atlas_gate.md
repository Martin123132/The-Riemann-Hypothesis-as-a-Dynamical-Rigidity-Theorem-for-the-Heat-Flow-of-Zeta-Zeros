# Hardy t=1e10 cross-block structural atlas

Date: 2026-08-09
Status: rigorous_t1e10_blocks_20_35_binary_structural_atlas_with_route_transition; finite structural atlas only, not a proof of RH

## Route transition

The one-CPU observer replay leaves the checkpoint state and all fifteen
displayed Hardy values exactly unchanged. It captures 424 cubic calls in
each block from 20 through 35.

```text
block  parent N  recursive  direct  child lengths
   20       104        374      50  {"1": 369, "2": 5}
   21       119        307     117  {"1": 307}
   22       135        240     184  {"1": 240}
   23       154        185     239  {"1": 185}
   24       175        134     290  {"1": 134}
   25       200         90     334  {"1": 90}
   26       228         61     363  {"1": 61}
   27       259         21     403  {"1": 21}
   28       295          2     422  {"1": 2}
   29       336          0     424  {}
   30       382          0     424  {}
   31       435          0     424  {}
   32       495          0     424  {}
   33       563          0     424  {}
   34       641          0     424  {}
   35       729          0     424  {}
```

There are `1414` recursive and `5370` direct calls. The final recursive calls occur in block `28`; blocks `29` through `35` are entirely direct.

## Exact structural margins

Every recursive exact-binary parent preserves its fundamental a1 cell,
dual floor, positive curvature, and real cubic stationary radical.

```text
minimum a1 selector margin       4.80545278534125145753580876110041565510735519096879E-5
minimum dual fractional margin  6.05083558137936839160960298770743230800767291851878E-5
minimum curvature floor         1.87637390985599111106941593730696398655084352742418E-3
maximum child |z|               1.29757399980552990362988882312330031833776619363892E-3
minimum child (1+z)             9.98702426000194470096370111176876699681662233806361E-1
maximum stored-frac rounding    2.36220350894242615583191923888001801441985797652293E-34
```

## Boundary

This identifies the finite scaling domain and its route transition. It
does not promote the block-20 endpoint estimate to blocks 21-28, certify
the later direct-kernel errors, or prove a height-uniform route theorem,
outer Hardy estimate, `Lambda<=0`, PF-infinity, RH, or a prize result.
