# Exact quadratic-domain gap for the local B modes

Date: 2026-08-13

Status: exact-coordinate interval certificate; not a proof of the complete
B-face estimate

Use the global outer Morse coordinate and the exact endpoint Fresnel
coordinate

```text
y=sqrt(t)s/2,
u=q_B(m,x)=sqrt(x/2)(B-2m/x).                         (QD1)
```

The complete two-dimensional phase defect is exactly

```text
-y^2+pi*u^2/2,   grad=(-2y,pi*u).                    (QD2)
```

There is no phase Taylor remainder.  For mode 621 the residual domain is on
the side `u<=q_621(y)` and `q_621(0)<0`; for mode 622 it is on the side
`u>=q_622(y)` and `q_622(0)>0`.  Thus both exact domains exclude the Gaussian
saddle `(0,0)`.

At the outer saddles,

```text
q_621(0)=[-50.45038224483949231296052328957914548227870992779848838385067395149814685990991781494022124314306671 +/- 2.73e-96],
q_622(0)=[40.28709467245251579647477943943994079183154412092771226898320202040835568222202451154780901776790755 +/- 2.72e-96].            (QD3)
```

On the nearest boundary arc the two gradient components vary oppositely, so
the minimax point is the unique balance `2|y|=pi|q_m(y)|`.  Arb brackets the
mode-621 balance by

```text
x in [0.00024238523658, 0.00024238523660],
y=[-28.59244965223986411011826751846067122192751380964729937828486428817916676244193176320533710198901799 +/- 1.30e-92], q=[-18.20251765244168790060854554459839943494074153108967340862475627448088187783655432292728466633309851 +/- 1.49e-96],
max(2|y|,pi|q|)>[57.18489517860224228123608513176879468775159735030455129286760351526122730991374785326555563017240276 +/- 2.34e-92]>57.17. (QD4)
```

For mode 622,

```text
x in [0.00024291643709, 0.00024291643710],
y=[22.82323784031335307266117728952123076428586438553172233267885523561626258219308071540168992338063121 +/- 1.63e-92], q=[14.52972449859799944546643221864543005065614613927855929407578834014820721998978647178471799411522715 +/- 1.48e-96],
max(2|y|,pi|q|)>[45.64647362149232831697694842358886205840350547562184863884365750848658363591046193969512118503030248 +/- 2.94e-92]>45.63. (QD5)
```

These are the exact nonstationary margins for the only two B-local integer
modes.  They are stronger than the common face-coordinate bound because
they use the true quadratic phase scales.  The tangent half-plane model can
be nearly degenerate only after replacing the curved boundary by a remote
line; neither exact local domain contains a joint saddle.

Pi provenance: `pi` comes from the equation-(9) outer and endpoint phases.
No fitted constant is used.

Proof boundary: exact phase coordinates and local-domain gradient margins
only.  No transformed-amplitude derivative bound, integration-by-parts
constant, complete B estimate, A-fold splice, complete paired residual,
`T_upper`, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion is
established.
