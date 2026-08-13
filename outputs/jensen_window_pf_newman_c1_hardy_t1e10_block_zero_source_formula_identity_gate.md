# Block-zero source-formula identity gate

Date: 2026-08-09

Status: finite source localization validated; not a proof of RH

Source-attribution notice: the `0.0161` quantity below is the residual of the
published/source hybrid with its equation-(124)/(125) endpoint convention.
The later midpoint comparison is an introduced diagnostic alternative, not a
cutoff rule stated in the paper.  It therefore cannot be subtracted here and
called a repair of equations (126)--(127) without an independent theorem.

The source starts block zero at `M2=159577` and tests two transition collars.
Arb certifies

```text
g  = [4.0771450318377521111281510822877889442952918335329105911531010857192387762405895732739914062283432528838596665 +/- 1.63e-104]
gm = [88.754631640417825737073375936101142586731690666691962161860756087976318617616306950525190011723084343947006498 +/- 1.63e-104]
```

Both are strictly above `3.2`, so neither branch executes and the transition
contribution is exactly zero in this run.

Put `s=alpha/a`, `b=sqrt(s^2-1)`, and `pc=(s+b)^2`.  Exact polynomial
reduction modulo `b^2=s^2-1`, together with positivity, proves

```text
sqrt(8/a)/sqrt(b) = 2sqrt(2)/(alpha^2-a^2)^(1/4),
b(b-s)+log(s+b)+1 = (log(pc)+1/pc+1)/2.
```

Thus each normalized source `ter` summand is exactly the left side of paper
equation (62) in ideal arithmetic.  The source normalization and phase are not
independent approximations or mismatched conventions.

An independent 70/110-digit Arb reconstruction of all
101 odd-alpha summands at each of the fifteen outputs differs from
the saved binary128 block-zero checkpoint by at most

```text
[2.6810248695185322548867869598850334460176993162941896121550992990274733891489100884338174774987073724697183402e-23 +/- 2.81e-98]
```

which is below `1e-20` on every output.  Yet the reconstructed equation-(62)
sum also exposes the source's unsuffixed default-real `0.01` shift literal.
Replacing that rounded literal by the intended exact hundredths changes the
normalized block by at most

```text
[1.0923997865621384193611837742222089156391867565862224691262334146055244706244195109960493222392047031538459296e-12 +/- 5.61e-98]
```

which is below `2e-12` on every output.  The reconstructed equation-(62)
sum is below its exact equation-(124) classical target by more than `0.005` on
all fifteen, with absolute residual between

```text
[0.016120329380248886474641597661340281075448198691814306353041446207144330032694216229842091826990537370778407959 +/- 3.85e-98]
[0.016159333143851197781225297199688843438084920962450971072587075793684815859389715498352412732688245283294464802 +/- 3.85e-98]
```

The actual source discrepancy is therefore not created by the transition,
the displayed amplitude normalization, phase translation, or observable
binary128 departure from their ideal formula.  It is the source-aligned
equation-(62)/(124)--(127) hybrid residual at these saved heights.  This gate
does not decompose that residual into a proved approximation error and a
proved cutoff correction.  It supplies no height-uniform remainder bound,
does not repair the hybrid formula, and has no RH implication.

The paper's equation (57) supplies only the unsigned global estimate
`E(t)<6.15*t^(-1/12)`.  At `t=10^10` its right side is
`[0.90269654958757276766611802561236967169932023405818200750825515008228403031135997702665091078937250555354130956 +/- 1.35e-110]`,
or `[180.53930991751458109357217820259899551671445294687146088108420372009277343750000000000000000000000000000000000 +/- 3.14e-14]` times
the saved tolerance.  It falls below `0.005` only for
`t>1230^12=11991163848716906297072721000000000000`.
That global estimate is not a signed or cellwise block-zero remainder.
