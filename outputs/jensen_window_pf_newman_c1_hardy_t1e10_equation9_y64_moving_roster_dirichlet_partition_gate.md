# Moving y=64 roster and contiguous Dirichlet blocks

Date: 2026-08-11

Status: exact moving-roster partition validated; not a proof of the grouped
chart integral bound

The mode dependence of the boundary function is affine:

```text
g_m(z)=h(m-C/4)+G(z),
G(z)=beta*tanh(z/beta)-64x(z)/beta,                     (MR1)

g_(m+1)(z)-g_m(z)=h>0.                                 (MR2)
```

Consequently, at every fixed `z`, the 84 modes split into three contiguous
integer blocks:

```text
m<=mu_-(z):             nonstationary,
mu_-(z)<m<mu_+(z):      endpoint/Fresnel,
m>=mu_+(z):             ordinary Morse,                (MR3)

mu_+-mu_-=2delta/h
 =[2.26146711747207085973182466566297593219491961543696826066809410081306820866740683338407090 +/- 1.62e-89].                (MR4)
```

Since `2<2delta/h<3`, at most three integer modes can lie in the Fresnel
block at one `z`.  Thus the per-mode transition envelope from Section 11.320
is never paid 84 times simultaneously.

The 168 entry and exit events are all distinct.  Their minimum separation is

```text
[0.0141200837569204213851998615434370871677886538849057911923673746488517470700008067048429198 +/- 1.06e-87]>0.014,                    (MR5)
```

so they produce 169 open `z` cells with fixed rosters.  The cell histogram is

```text
Fresnel count 0: 2 cells,
Fresnel count 1: 2 cells,
Fresnel count 2: 83 cells,
Fresnel count 3: 82 cells. (MR6)
```

At `g=-delta` equality belongs to the nonstationary block; at `g=+delta`
it belongs to the Morse block.  The isolated event points have zero `z`
measure, and the half-open assignment makes the full partition exact.

Each nonempty block retains an exact finite Dirichlet kernel.  For consecutive
`a<=m<=b`,

```text
D_[a,b](y)
 =exp[-i(d_a+d_b)y/2]
  sin((b-a+1)hy/2)/sin(hy/2),                          (MR7)
```

with removable denominator zeros interpreted by the original finite sum.
This supplies a cancellation-preserving route for the nonstationary and
Morse blocks rather than summing their per-mode absolute envelopes.

Proof boundary: exact saved-height moving-roster algebra, event atlas, and
contiguous-kernel decomposition only.  No grouped chart integral estimate,
completed z integral, lower-interior join, `T_upper` theorem, `Lambda<=0`,
RH, or prize-level conclusion is proved.
