# Positive-real-tail K on the first height subcell

Date: 2026-08-28

Status: rigorous positive-real-tail Hardy-operator bound; complete `K_T`
is not yet assembled

For `s=1/2+it`, `R=2*pi*U`, and

```text
F_s(u,t)=u^(-s) exp(i*u^2/(4*pi)) Delta_q(u),
P_tail(t)=E_s(t) integral_R^infinity F_s(u,t)du,
E_s=exp(pi*t/2+i*pi/4)(2*pi)^(s-1),
```

the fixed-roster event atlas keeps `R`, `q_-`, and the finite label count
constant on `I_1=[10^10-10^-4,10^10+10^-4]`.  Therefore differentiation
under the absolutely convergent positive-real tail gives exactly

```text
(d_t+i theta'-H'/H)P_tail
 =E_s integral_R^infinity
  [pi/2-H'/H+i(theta'-log(u/(2*pi)))]F_s(u,t)du.     (PT1)
```

For real `u>=R`, the finite geometric difference is positive and

```text
0<Delta_q(u)<=exp(-q_- u)/(1-exp(-R)).              (PT2)
```

Writing `u=Rv`, the logarithmic weight is bounded without discarding its
endpoint cancellation:

```text
|theta'-log(u/(2*pi))|
 <=|theta'-log U|+log v,
log v<=v-1,                u^(-1/2)<=R^(-1/2).      (PT3)
```

Equations (PT1)--(PT3) yield the elementary closed bound

```text
|L_t[P_tail]|
 <= |E_s| R^(-1/2) exp(-q_-R)/(1-exp(-R))
    *[A/q_-+1/(R q_-^2)],
A=|pi/2-H'/H|+|theta'-log U|.                      (PT4)
```

Arb uses the direct and duplication formulas for `H'/H` as an overlap guard
and stable positive `H` transport.  It obtains

```text
log10 |L_t[P_tail]| <= [-1873216239.351917559937632990376398424766109340948599520 +/- 4.90e-46],
log10 |Hardy_t[L_t[P_tail]]/H|
                      <= [-1873216239.050887564273651795162795246915207837994461268 +/- 3.99e-46].
                                                               (PT5)
```

The checker raises precision, splits `I_1` into two changed half-boxes, uses
the duplication formula for `H'/H`, and reconstructs (PT4).  Its altered
five-label packet independently differentiates `E_s T_R` numerically and
agrees with the weighted integral, catching the `pi/2`, `log(2*pi)`, and
`-log u` signs.

Pi provenance: `pi` comes only from the Mellin--Fresnel kernel, the exact
quarter-disk scale `R=2*pi*U`, the Gamma normalization in `E_s`, and the
Riemann--Siegel phase.  No fitted geometric constant supplies `pi`.

Proof boundary: (PT1)--(PT5) certify only the positive-real-tail contribution
on `I_1`.  The exact tiny target correction and final assembly with the
transition--upper-arc and lower-plus-ordinary packets remain open.  No
complete `K_T`, wider `Q_K-T` sign interval, event-cell theorem, wall handoff,
all-height theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion
is proved.
