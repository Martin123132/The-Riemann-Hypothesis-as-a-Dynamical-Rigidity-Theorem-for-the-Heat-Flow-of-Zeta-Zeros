# Jensen-Window PF Newman Polymath-15 Real-Edge Projective-Current Gate

Date: 2026-07-31

Status: exact division-free edge current and rigorous negative `q=1`
cofinal leading symbol. This is not a proof of a finite-height edge theorem,
an Abel gap, `Lambda <= 0`, or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_real_edge_projective_current_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_real_edge_projective_current_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_real_edge_projective_current_gate.py
```

## Exact Edge Current

```text
V=e+z_N, G=g+s_*'*u_N*z_N, C_edge=Re(V), D_edge=Re(G)-alpha*C_edge, alpha=Re(s_*')*u_N

J_edge=C_edge*D_(edge,x)-D_edge*C_(edge,x)=C_edge*Re(G_x)-Re(G)*C_(edge,x)-alpha_x*C_edge^2
```

The `alpha*C*C_x` terms cancel. The exact derivative retains

```text
G_x=g_x+(s_*''u_N+s_*'u_(N,x)+s_*'u_N*z_(N,x)/z_N)z_N; g_x retains H_(a,xx), mu_(a,x), the moving T_0+i factor, and the real radial kappa rate
```

No ratio is taken at `C_edge=0`, `c_0=0`, or `H_a=0`.

## Correct Terminal Phase

On

```text
theta=(1-p)/2, a=N+theta, q=2tL^2=1, S_a=kappa*T_0, r_S=S_(a,x)/S_a, fixed p in [-1,1], N->infinity
```

the endpoint and terminal limits are

```text
F=C_0(p), H_a->F, a*H_(a,x)->-F'/(4pi), a^2*H_(a,xx)->F''/(16pi^2), a*mu_a->0, a^2*mu_(a,x)->0

Q(p)=exp[-pi*i*(p^2/2-p+3/8)], z_N/S_a->-Q(p), a*u_N->theta, a*(z_(N,x)/z_N-r_S)->-i*theta/2
```

The phase guard is essential:

```text
Q=conjugate(R)*exp[-2pi*i*(1-p)] for R=exp[pi*i*(p^2/2+p+3/8)]. The extra factor is one at p=0 and p=1 but not in the interior; endpoint-only tests cannot determine the terminal phase.
```

At `p=29/30`, the wrong endpoint-interpolated phase gives the Arb current
ball

```text
[0.000106841426805094672587728330430117455888265649129042645 +/- 4.85e-58]
```

whereas the corrected source gives

```text
[-0.000807457738109035704810346742034299023984120305185651123 +/- 6.39e-58]
```

The full finite source was used to reject the wrong branch before promotion.

## Three-Jet Collapse

Define

```text
A=Re(F-Q); B=Re[-F'/(4pi)+i*theta*Q/2]; M=Re[F''/(16pi^2)+Q{i/(16pi)+theta^2/4}]
```

Then

```text
C_edge/S_a->A, a*D_edge/S_a->B, a(C_(edge,x)-r_S*C_edge)/S_a->B, and a^2(D_(edge,x)-r_S*D_edge)/S_a->M

B=-A'/(4pi), M=A''/(16pi^2)

a^2*J_edge/S_a^2->K_edge(p)=[A*A''-(A')^2]/(16pi^2)
```

This is division-free. In particular,

```text
The polynomial limit is primary. At A=0 it equals -(A')^2/(16pi^2), with no division by C_edge or A.
```

## Whole-Cell Certificate

The real source is

```text
A(p)=-cos(pi*(p^2/2-2p+3/8))/(2cos(pi*p))=-Re C_0(p-2), with removable values A(-1/2)=5/4 and A(1/2)=-3/4
```

At the two apparent poles the checker uses

```text
At p=1/2+y, A=(y-3)sinc(pi*y*(y-3)/2)/[4sinc(pi*y)]; at p=-1/2+y, A=(5-y)sinc(pi*y*(y-5)/2)/[4sinc(pi*y)].
```

The complete Arb cover has `3072` rational intervals at
`192` bits. Its weakest saved lower endpoint is

```text
[0.06339593028046580621590837488799260877676282746790246852103520367109893 +/- 3.49e-72]
```

on `4093/4096..1`. Therefore

```text
(A')^2-A*A''>3/50 on -1<=p<=1; since pi^2<10, K_edge(p)<-3/8000.
```

This includes every zero of `A`; no tangent chart is needed.

## Handoff

```text
The complete q=1 cofinal real-edge leading symbol is strictly clockwise on the full saddle cell, including both removable C_0 points and every A=0 fibre.

Derive an explicit uniform finite-a remainder smaller than 3/8000 for a^2*J_edge/(kappa*T_0)^2 on q=1, then extend the signed edge estimate to q>=1 and splice it through adjacent cutoffs before returning to the cumulative contact scalar.
```

The missing finite-`a` error bound is now quantitative and has a fixed
`3/8000` reserve. The Abel gap, winding cap, contact exclusion, and RH
remain open.
