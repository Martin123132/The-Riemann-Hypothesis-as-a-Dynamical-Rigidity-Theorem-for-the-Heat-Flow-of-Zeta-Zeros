# Equation-(62) diagnostic midpoint-partition gate

Date: 2026-08-09

Status: finite diagnostic cutoff comparison validated; not a proof of RH

Paper equation (124) relates each alpha endpoint to a continuous lower root,
and equations (125)--(127) use an endpoint cutoff.  Neither the 2026 paper nor
the source prescribes the nearest-root Voronoi midpoint rule tested here.
This gate introduces that rule as a diagnostic alternative, extends it through
all 101 alpha terms, and certifies every root with Arb.

At `t=10^10`, the first nearest roots are `39852` and `39691`; their midpoint
is `39771.5`.  The introduced diagnostic therefore needs a lattice convention
at this height.  The artifact retains both floor and ceiling completions at
every half-integral boundary.  Its final boundary is not ambiguous: the roots for
alpha `159777` and `159779` are `37946` and `37936`, so their midpoint is the
integer `37941`.

The executable derives its complementary classical tail directly from the
terminal alpha endpoint, in agreement with equations (124)--(127).  It starts
at `37946`.
At the centre output the three relevant values are

```text
alpha prefix - source-endpoint target = [-0.016139907963525290656361230144305829820358816624850213147648360240463455883775822774757447583675703572749896110 +/- 2.80e-98]
alpha prefix - midpoint-cell target   = [0.0014681344930085758486665616493525399471890491221522584292615533140900122642128490317202407635381062158288896025 +/- 2.80e-98]
midpoint target - source target       = [-0.017608042456533866505027791793658369767547865747002471576909913554553468147988671806477688347213809788578785713 +/- 1.26e-98]
```

The identity `source residual = midpoint residual + boundary correction` is
certified at every prefix and output.  Across all fifteen final prefixes, the
source-endpoint residual has magnitude above `0.005`, while the unambiguous
midpoint-cell residual has magnitude below `0.005`.  Its absolute range is

```text
[0.0014665922903781999713192554563354664691760077706987656158526509641960447747584751382173915132979645096353043616 +/- 3.86e-98]
[0.0014696132280717262580692086493457472254471606662875120158175036046911087973257960881473892184137031273117779471 +/- 3.86e-98]
```

Thus replacing the source cutoff by this diagnostic midpoint changes the
comparison by the five terms `37941..37945` and leaves a residual near
`0.00147`.  This is an exact finite comparison between two conventions.  It
does not show that the source omitted those terms or that the midpoint is the
published hybrid's correct boundary.

The smaller residual motivates a possible new cutoff theorem, but that theorem
must be derived independently and must transport the complete hybrid error.
The prefix atlas is diagnostic evidence only and has no RH implication.
