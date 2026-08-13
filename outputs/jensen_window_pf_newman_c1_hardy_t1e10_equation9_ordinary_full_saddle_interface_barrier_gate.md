# Bare full-saddle interface barrier

Date: 2026-08-13

Status: exact leading-current obstruction certified; not a proof of a full mode-error lower bound

At the reduced Morse saddle `x_m`, retain the complete grouped alpha current

```text
G_m=m*sqrt(2)*x_m^(-3/2) I_0+x_m^(-1) I_1,
I_0=integral_(q_A)^(q_B) exp(i*pi*q^2/2)dq,
I_1=integral_(q_A)^(q_B) q exp(i*pi*q^2/2)dq.
```

The full-line current is
`G_m^full=m*sqrt(2)*x_m^(-3/2)(1+i)`.  Therefore the exact frozen-saddle
ratio is

```text
R_m=I_0/(1+i)+sqrt(x_m) I_1/[m sqrt(2)(1+i)].
```

At the proposed atlas edge `m=39695`,

```text
q_A*=[-0.953603923431587592051739380591228181809442280280507503436494725937069927895274812195661274 +/- 7.75e-86],
|R_m-1|=[0.209940701887517246514837777791161734974948370397338498881498249385791505992322864589536336 +/- 1.38e-79],
|R_m-1|/sqrt(m)=[0.00105372853687510519659312819162210253339777400964453292653697415633008201329808665214798012 +/- 6.91e-82].
```

That normalized leading defect is more than
`[55.4593966776371224465380172765756826812333368792451437911950051784515380859375000000000000 +/- 9.89e-15]` times the local `1.9e-5` target.
At the true lower-interior edge `m=39852`,

```text
q_A*=[-0.000825592621781281439668889857774576378487547421998158096359755897023555569088073283763241127 +/- 7.74e-86],
|R_m-1|=[0.499584469761557615034723230613278321250069134670473874645303894285295991696587384062868091 +/- 9.91e-80],
|R_m-1|/sqrt(m)=[0.00250255636853014738065931472089966134099853708205620682777109341960461775869947317427617325 +/- 4.97e-82],
physical=[0.00112654616534528370979617395741384362056872756635101558084164397324115190070051366314018162 +/- 2.24e-82].
```

The mode is a half-saddle to leading order: the normalized defect exceeds
the local target by `[131.713493080534088902349206729505129084834003094783838605508208274841308593750000000000000 +/- 2.35e-14]` and the
physical defect exceeds its target by `[130.993740156428353952502932333590551089441689214254438411444425582885742187500000000000000 +/- 2.59e-14]`.

This rigorously rules out a **termwise splice whose ordinary overlap object is
only the bare full-line carrier**.  It does not lower-bound the complete
`x`-integrated exact-mode error, because amplitude variation and neighboring
modes may cancel.  The admissible next object must retain the incomplete
Fresnel plus endpoint current through the ordinary Morse step, or aggregate
the signed endpoint defects before taking absolute values.

Pi provenance: every occurrence of `pi` comes from the equation-(9) Kummer
quadratic phase, its Fourier-Poisson dual mode, and the already derived
equation-(9) normalization.  No geometric fit is introduced here.

Proof boundary: exact grouped alpha-current ratios at four saved-height
interface witnesses only.  No full `x`-integral lower bound, uniform ordinary
Morse remainder, complete `T_upper`, `Lambda<=0`, PF-infinity, RH, or
prize-level conclusion is proved.
