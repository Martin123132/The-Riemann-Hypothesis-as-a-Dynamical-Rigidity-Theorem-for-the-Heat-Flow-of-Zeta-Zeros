# Newman Theta Switch-Defect Explicit Constant Gate

Date: 2026-07-24

Status: explicit all-`N` modular-tail constants and `m=9`
value/first-jet compiler. This is not a retained-separation theorem
and not a proof of `Lambda<=0`, RH, or a Clay-prize result.

## Direct Tail

For `N>=2`, set `K=N+1`, `y=exp(4u)`, `Q=pi*n^2`, and
`s=erfc(3*sinh(4u))/2`. The proof bounds the direct tail

```text
r_N(u)=sum_(n>=K)[(1-s(u))*phi_n(u)+s(u)*phi_n(-u)].
```

The derivative recurrence and switch Laurent recurrence give four
positive templates for each Leibniz index: forward order zero,
forward positive switch order, reflected order zero, and reflected
positive switch order.

## Integral Split

The forward template uses

```text
u^p<=p!*exp(u),  u^2<=exp(4u)/16,  log(y)<=y-1,
```

and becomes an explicit Gaussian arithmetic tail. The reflected
template is split at `y=sqrt(2)`. Its compact part is another
Gaussian arithmetic tail. On the large part,

```text
Q/y+h >= c0*n^(4/3),
h >= 9*y^2/16,
c_safe=c0/4=3*3^(2/3)*pi^(2/3)/16.
```

One quarter of the phase supplies the discrete stretched
exponential. The remaining three quarters produce a directed upper
incomplete-gamma integral.

## m=9 Witnesses

The following rows sample the symbolic all-`N` formula.

| N | K | log10 d0 upper | log10 d1 upper |
|---:|---:|---:|---:|
| 2 | 3 | `[33.86062629149568464788141732532915624154703589717920610 +/- 2.08e-54]` | `[33.86106718410240630488429883620977944207742894522514746 +/- 2.37e-54]` |
| 3 | 4 | `[33.86062620818788484717412495383867591519106126723342742 +/- 9.55e-55]` | `[33.86106710002014110139629660849706260911013560671280552 +/- 3.64e-54]` |
| 4 | 5 | `[33.86061470708735460169848598593588456523519766339618837 +/- 2.74e-54]` | `[33.86105554244522033300505397382589033628545915280295727 +/- 4.91e-54]` |
| 7 | 8 | `[33.84016969318321122904921195841699659329098195763910340 +/- 3.10e-54]` | `[33.84058788955254255784135346089558171664084429322396969 +/- 4.48e-55]` |
| 10 | 11 | `[33.55592802806234903207680969784594268241754725479111137 +/- 4.56e-55]` | `[33.55624660088155362470088521525144514452260813384927643 +/- 3.46e-54]` |
| 16 | 17 | `[30.99118706319523547854393850584271059829890614035627249 +/- 2.26e-54]` | `[30.99132626756835600364322278790023633950989458266325543 +/- 1.67e-54]` |
| 25 | 26 | `[23.27462150892396851425064946051192768565271291748185619 +/- 4.74e-55]` | `[23.27467964091017543121017907148467094518068308593321462 +/- 4.19e-54]` |
| 64 | 65 | `[-34.25354290581904070514668649483765497349530856098097924 +/- 6.91e-55]` | `[-34.25353396573068080244017421106192314750194676950623556 +/- 4.34e-54]` |

Every stored number is an Arb enclosure at 256-bit precision.
The independent checker regenerates the kernel and switch
recurrences and recompiles every witness.

## Remaining Boundary

The switch-tail constants are no longer an effective-constant
placeholder. The main cofinal obligation is now a rigorous
retained value-or-derivative lower margin across all adaptive
transition cells, strong enough to dominate these upper budgets.
