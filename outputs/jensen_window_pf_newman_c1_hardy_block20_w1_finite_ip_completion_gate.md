# Hardy block-20 displayed-W1 finite-ip completion gate

Date: 2026-08-06

Status: rigorous evaluation of every displayed equation-(81) endpoint term omitted by the source index clipping; not a proof of the exact-contour theorem or RH

## Finite index set

The paper is exact through equation (69).  Its `W1` approximation in equation
(81) contains lower- and upper-endpoint correction terms over the finite
secondary index set.  The paper notes that all terms could be summed, but uses
a practical endpoint cutoff.  The accepted source sets `ip=3`; when the
transformed length is 1 or 2 it clips this to `ip1=1`.

For each of the 374 recursive block-20 calls this gate constructs the complete
finite set `jbot <= n <= L`, subtracts the indices actually visited by each
source loop, and evaluates every remaining displayed `W1` term with Arb at 180
and 260 decimal digits.  There is no extrapolation beyond the finite sum.  To
isolate the source's index clipping, the two erfc arguments retain the actual
default-real `sqrt(2.0)` value, binary32 `0x3FB504F3`; its replacement by the
mathematical square root is a separate correction gate.

```text
calls whose two source loops already cover every displayed index = 204
calls with at least one omitted displayed endpoint term          = 170
omitted lower endpoint terms                                    = 5
omitted upper endpoint terms                                    = 170
total omitted displayed terms                                   = 175
```

Each omitted term is split exactly into its algebraic inverse-cube part and
its erfc/exponential part.  The finite inverse-cube sums are the terms that can
equivalently be written as differences of order-two polygamma values; no
unavailable Maple expression is assumed.

## Comparison with the required native correction

Let `T_ip` be the sum of all omitted displayed terms and let `R_q` be the exact
native correction required by the independently summed parent and child.  The
remaining target is

```text
R_after = R_q - T_ip.
```

The finite roster gives

```text
maximum |T_ip|                              <= 4.55725961949695486116266009508599260515896660844525E-4
maximum inverse-cube part                   <= 4.55725961949695486116266009508599260515896660844525E-4
maximum erfc/exponential part               <= 9.29814720943049516762698215552852683656251689692566E-187
minimum |R_after|                           >= 4.57824824240049749388074215994808155173220435825063E-3
maximum |R_after|                           <= 2.05024273257495403610363658542613157352603105485991E-1
calls with R_after excluding zero              374 / 374
rigorously improved / worsened / unresolved    59 / 111 / 204
```

The retained-term replay gap is at most
`6.35886281329559141846144294780269937968237171321547E-26`.  This checks the
formula transcription against the independently admitted source component
probe; it is not used as an analytic contour-error bound.

## Boundary

This gate decides only the finite index truncation inside the already
approximate displayed `W1` formula.  It does not turn equations (72), (75), or
(77) into equalities, bound omitted higher saddle phase, validate an unknown
Maple implementation, control `W2`--`W4`, establish a height-uniform recurrence,
or prove `Lambda<=0`, PF-infinity, RH, or a prize-level theorem.
