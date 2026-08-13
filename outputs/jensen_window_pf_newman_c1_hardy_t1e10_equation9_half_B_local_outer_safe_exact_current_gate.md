# Exact outer-safe local B current for modes 621 and 622

Date: 2026-08-13

Status: cancellation-safe saved-height interval certificate; not a complete
B-face estimate

For a positive mode, retain the exact continuous target-step/tail current

```text
R_B,m^+=U_B,m+[1_(m<=621)-1_(q_B<0)]P_bulk,m.        (LC1)
```

Because `B` is odd, quadratic completion gives its exact face coefficient

```text
C_exact,+,m=exp(-i*pi*q_B^2/2)R_B,m^+/x.             (LC2)
```

The jump of `U_B` at `q_B=0` is cancelled by the step in LC1 before any
integral or absolute value is taken.  The negative partner is uniformly far
from its crossing and is represented by its first three direct Fresnel
terms.  On a positive normal-safe arc those same three terms are used for
both signs; on the complementary arc LC2 replaces only the positive terms.

The mode-621 balance point lies at
`xi=[-86.0499957141587921738879935416963255158122899231504821548268 +/- 1.71e-53]`, below the Gaussian window, so
mode 621 uses the exact outer-safe representation throughout `|xi|<=70`.
Mode 622 switches to its normal-safe direct expansion at
`xi=[68.7277355764349373185736515590694064198847585223466476875555 +/- 1.71e-53]`.  Every panel is also cut at the
exact mode crossing:

```text
xi_621=[-63.2531757872986894942948562936101789012866815877696305213028 +/- 1.71e-53],
xi_622=[50.5107901702604273895499366426400915767532393635485956339688 +/- 1.71e-53].               (LC3)
```

An append-only `563`-panel Arb quadrature evaluates the two
mode sum before projection.  It gives

```text
mode 621: [-0.000214200534446594413739745925668072724364927098790966476246701 +/- 4.24e-20],
mode 622: [7.59442119956973905382042540949979525519504283761901964803434e-5 +/- 7.78e-20],

principal sum:
[-0.000138256322450897023201541671573074771812976670414776279766358 +/- 1.21e-19].          (LC4)
```

The existing positive normal-safe remainder allowance is
`[1.94796259922431190234924394619204508559231052593413807257246e-10 +/- 4.30e-70]`.  The negative
partners contribute at most
`[8.42338703285006033858970859018653425557076129774347057108145e-36 +/- 2.07e-90]` beyond their
three displayed terms.  Adding both adversely still proves

```text
E_local,B <= [-0.000138256127654637100770351436648671729221435261123385096643960 +/- 1.21e-19]
           < -1.369e-4.                              (LC5)
```

Thus the local pair is not a positive budget cost: its exact physical
projection is rigorously negative by more than `1.369e-4`.  This sign is
obtained only after the crossing step and Fresnel tail remain joined.

Pi provenance: all `pi` factors come from the equation-(9) Fresnel and
outer phases, odd-endpoint completion, or the paper normalization.  No
fitted constant is used.

Proof boundary: modes 621 and 622 on the saved-height `|xi|<=70` B window
only, with a deliberately overinclusive positive normal-remainder
allowance.  No outside-window grouped tail, complete B estimate, A-fold
splice, complete `T_upper`, height-uniform theorem, `Lambda<=0`,
PF-infinity, RH, or prize-level conclusion is proved.
