# Exact y=64 endpoint-current roster reassembly

Date: 2026-08-11

Status: exact endpoint-current reassembly and upper lattice character
validated; not a proof of the grouped bulk integrals or `T_upper`

The exact scales satisfy

```text
h/sigma=pi,
y_B=(B-C)/sigma,
h*y_B=pi(B-C).                                         (ER1)
```

Both source endpoints are odd.  At the saved endpoints,

```text
B-C=4962844=4*1240711,                                 (ER2)
```

and `1240711` is odd.  Therefore, for every integer mode,

```text
exp(-i d_m y_B)
 =exp(i*pi*C(B-C)/4)=-1,                               (ER3)

D_84(y_B)=-84.                                         (ER4)
```

This is an exact lattice character, not a numerical large-argument phase
evaluation.  At the artificial split endpoint the complete kernel instead
satisfies

```text
|D_84(64)|=[0.611603386768977255057884711066356804326409974292972528833133038880343364332918738988039595 +/- 2.43e-86]<0.612.   (ER5)
```

On every one of the 169 open moving-roster cells, the nonstationary,
Fresnel, and Morse mode blocks reassemble in order to `39853..39936`.
Consequently their endpoint currents sum pointwise to the complete kernel:

```text
sum_chart D_chart(U)=D_84(U),       U=64 or y_B.        (ER6)
```

The `U=64` current cancels exactly against the opposite current from
`K(0,64)` by the primitive telescope of Section 11.320.  The `U=y_B`
current is the genuine source `B` endpoint current and must remain paired
with the grouped Fresnel/bulk contribution; (ER4) evaluates its mode
character but does not bound or discard it.

Proof boundary: exact saved-height scale algebra, endpoint characters, and
169-cell endpoint reassembly only.  No grouped bulk cell integral,
lower-interior join, source `T_upper` identification, complete error theorem,
`Lambda<=0`, RH, or prize-level conclusion is proved.
