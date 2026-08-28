# Alternating boundary-jet transition packet

Date: 2026-08-27

Status: rigorous grouped Fresnel/Gamma/A transition packet at `t=10^10`;
the full `P_W`/unowned-cell join remains open.

On the contiguous transition interval

```text
[a,b]=[39852.5,39936.5],       b-a=84,
F(y)=y^(-s)exp(i*pi*(159577*y-y^2)),
```

the finite source is

```text
sum_(j=0)^(M-1) integral_a^b F(y)exp(2*pi*i*j*y)dy,
M=2481423.
```

Because both endpoints are half integers and their difference is integral,
`exp(2*pi*i*j*a)=exp(2*pi*i*j*b)=(-1)^j`.  Repeated integration by parts gives
the exact finite expansion

```text
I_transition=I_0
 +sum_(r=0)^(R-1) (-1)^r [F^(r)(b)-F^(r)(a)]
    A_(r+1)(M-1)/(2*pi*i)^(r+1)
 +R_R,

A_p(N)=sum_(j=1)^N (-1)^j/j^p.
```

The remainder is bounded without summing the source roster:

```text
|R_R| <= zeta(R)/(2*pi)^R integral_a^b |F^(R)(y)|dy.
```

With `R=12` and 128 derivative slabs,

```text
|R_R| <= [6.449047667892289955345308360003388431686827404967077819e-17 +/- 6.45e-73].
```

The complete transition source integral is

```text
[0.07154278585722209135245427134057774180276338970003299759032738714187419 +/- 6.45e-17]
+ i [0.06062036407180469713273252142459330655801741589770030209519090444070046 +/- 6.45e-17].
```

Subtracting the exact extended Gamma lift and the grouped natural A lift before
projection gives

```text
T_A^join = [0.04434802829267891969187949473383939830430149015295483198684460561537668 +/- 2.13e-12]
          +i [0.04625765106805546494431672804809789345992318548847537361039770496291719 +/- 2.13e-12],

Hardy_t[T_A^join] = [0.1270842386787259385836210648632645330567928594437503538301626196809571 +/- 5.97e-12].
```

Pi provenance: every `pi` comes from the finite Fresnel source, its exact
integer Fourier spacing, the Gamma normalization, or the Riemann-Siegel
phase.  No fitted constant is introduced.

Proof boundary: rigorous transition block at the one saved height only.  This
gate does not join `P_W` to cells `1..621` and `39937..infinity`, bound the
ordinary packet, enclose complete `J_Z` or `D_K`, or prove a non-A, all-height,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion.
