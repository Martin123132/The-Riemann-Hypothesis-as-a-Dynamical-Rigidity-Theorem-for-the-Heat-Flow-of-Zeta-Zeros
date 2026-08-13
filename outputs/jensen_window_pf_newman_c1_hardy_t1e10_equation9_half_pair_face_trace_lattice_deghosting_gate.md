# Face-trace lattice geometry and tangent deghosting

Date: 2026-08-13

Status: interval certificate; not a proof of the B face remainder bound

On the triangular face `z=x`, the two mode phases cancel exactly:

```text
Phi_(m,D)(x,x)
 =pi*D^2*x/4+(t/2)log((1-x)/x)=Phi_D(x).              (FT1)
```

The trace is independent of `m`.  Its two stationary points are

```text
x_D^pm=[1 pm sqrt(1-8t/(pi D^2))]/2,                  (FT2)
```

and their lattice centers are

```text
m_D^pm=D*x_D^pm/2
      =[D pm sqrt(D^2-8t/pi)]/4.                      (FT3)
```

These are exactly the two inverse roots of `alpha_m=D`.  At the lower trace
saddle, the normal derivative is

```text
partial_z Phi_(m,D)
 =pi*D^2[1-(m/m_D^-)^2]/4.                            (FT4)
```

Hence a genuine face critical point requires `m=m_D^-`, not merely a null
tangent in a mode-centered approximation.

For `B=5122421`, rigorous intervals give

```text
x_B^-=[0.0002426805620273058955478049164473530951179361333034812162106778015490549729094654724584894389903118538 +/- 5.37e-102],
m_B^-=[621.5560036102371463889412039565834444235567629462757775115581974443558616929385174481439652810961182 +/- 1.38e-95].                             (FT5)
```

Thus only modes `621` and `622` are adjacent to the B face event.  Their
lattice distances exceed `0.55` and `0.44`, and even their exact normal
derivatives have magnitudes above `3.68e10` and `2.94e10`.  Every other
integer mode is at least `1.44` lattice units from the event.

This resolves the apparent tangent singularity near mode `257`.  Although
the old mode-centered tangent defect is close to zero there, the exact B-face
normal derivative is

```text
[17084945357014.21122465308928149375680640293227324524062743569348004304019859504342864206455476380390 +/- 1.62e-85],              (FT6)
```

with magnitude above `1.70e13`; its outer saddle is also more than `0.0002`
away from `x_B^-`.  Mode `257` is not a joint or face critical point.  Its
large tangent correction came from extending a local line to a remote
stationary region.

For the lower A face,

```text
m_A^-=[39852.39138623123896184313678682614784426989001770014282206708149972339563001678357751562596063941840 +/- 4.59e-94],                             (FT7)
```

so modes `39852,39853` are the adjacent lower-face pair.  The reflected root
is `[39936.10861376876103815686321317385215573010998229985717793291850027660436998321642248437403936058160 +/- 4.60e-94]` and is already represented by half-Kummer
conjugation.  The separate `39894,39895` pair belongs to the `x=1/2` outer
saddle corner.

The quantitative theorem can therefore use a two-mode B face chart plus
normal integration by parts for all other B modes.  The old 51-mode tangent
collar remains a diagnostic partition, not the natural exact local roster.

Pi provenance: all `pi` factors come from the exact equation-(9) pair phase
and its boundary trace.  No fitted constant is used.

Proof boundary: exact trace cancellation, stationary roots, and saved-height
lattice/normal margins only.  No normal integration-by-parts constants,
two-mode B face estimate, A-fold splice, complete paired residual,
`T_upper`, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion follows.
