# All-order signed-pair endpoint jets and the turning-ratio route guard

Date: 2026-08-23

Status: exact all-order pair identity and exact turning-ratio obstruction
certified; raw higher-jet route rejected, phase-adapted joined bound open

Let

```text
f_x(u)=(A+2u)exp[i*pi*x*(A+2u)^2/4],
I_m(x)=integral_0^L f_x(u)exp(-2*pi*i*m*u)du,
D_m=2*pi*i*m.
```

Repeated integration by parts at any even order `2r` and pairing `+m` with
`-m` before taking absolute values gives the exact identity

```text
I_m+I_-m
 =-2 sum_(k=0)^(r-1) Delta f_x^(2k+1)/D_m^(2k+2)
  +[R_(m,2r)+R_(-m,2r)]/D_m^(2r),                    (AJ1)

R_(+/-m,2r)
 =integral_0^L f_x^(2r)(u)exp(-/+2*pi*i*m*u)du.       (AJ2)
```

Thus every even endpoint derivative cancels from the signed pair; only odd
endpoint jumps remain.  The elementary absolute remainder is

```text
|remainder|<=2 integral_0^L |f_x^(2r)(u)|du
                    /(2*pi*m)^(2r).                   (AJ3)
```

Put `y=A+2u` and `q=i*pi*x`.  The endpoint derivatives have the exact
polynomial recurrence

```text
p_0(y)=y,
p_(n+1)(y)=2p_n'(y)+q*y*p_n(y),
f_x^(n)(u)=p_n(y)exp(q*y^2/4),                        (AJ4)

p_n(y)=q^n y^(n+1)+lower powers.                     (AJ5)
```

For endpoint `D` in `{A,B}`, the ratio of successive highest odd-jet
monomials is therefore exactly

```text
rho_D(x,m)^2=[x*D/(2m)]^2.                            (AJ6)
```

This is the same normal turning coordinate already resolved by the exact
Fresnel/Morse transition charts.  At a normal endpoint crossing `x=2m/D`,
`rho_D=1`; repeated endpoint extraction is not a descending geometric
hierarchy there.

The saved integers make the obstruction exact.  At `x=1/2` and the top
target mode,

```text
rho_B=32.1001967714443239287903>32,
rho_A=1.00000626660650726584834>1.                                    (AJ7)
```

For the first outer mode `m=39895`,

```text
rho_A=0.999981200651710699922603<1,                                    (AJ8)
```

so the A half-boundary lies precisely between modes 39894 and 39895.  The B
normal crossing for mode 622 is at

```text
x=1244/5122421=0.000242853916146290978650066,    rho_B=1.                 (AJ9)
```

An 80-digit floating scout, used only to choose the route, evaluates the
separate endpoint terms of (AJ1) through derivative order seven:

```text
B_normal_crossing_m622: rho^2=1.0; observed ratios=1.000000000000000000145, 1.000000000000000000524, 1.000000000000000001143
A_half_boundary_m39894: rho^2=1.000012533252284845511; observed ratios=1.000012533252284845547, 1.000012533252284845642, 1.000012533252284845797
A_first_outer_m39895: rho^2=0.9999624016568369774895; observed ratios=0.9999624016568369775258, 0.9999624016568369776208, 0.9999624016568369777758
```

The exact conclusion is limited but decisive.  Raising the raw endpoint-jet
order cannot give a uniform descending absolute hierarchy across the B
crossing or the A half-boundary.  This does not rule out a signed resummation,
a contour argument, or the already-certified transition functions.  It says
that the next `R_after_A` estimate must retain those phase-adapted A/B
extractions and act on the remaining joined projector current; it must not
replace them by a global high-order Bernoulli or endpoint-derivative norm.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_all_order_endpoint_jet_route_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_all_order_endpoint_jet_route_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_all_order_endpoint_jet_route_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_all_order_endpoint_jet_route_gate.py
```

Pi provenance: every `pi` in (AJ1)--(AJ9) is inherited from the Kummer
quadratic phase and integer Fourier character.  No fitted or geometric
occurrence is introduced.

Proof boundary: exact all-order paired integration by parts, derivative
recurrence, highest-monomial ratio, and rational turning geometry only.  The
floating order scout is diagnostic.  No lower bound on the true signed
remainder, no exclusion of cancellation-aware resummation, no quantitative
`R_after_A` or `R_Dir` bound, no complete `Q_K-T`, all-height theorem,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
