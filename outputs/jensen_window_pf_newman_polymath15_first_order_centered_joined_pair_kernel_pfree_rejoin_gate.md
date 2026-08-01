# Joined Pair Kernel And P-Free Rejoin Gate

Date: 2026-07-30

Status: exact carrier-pair and p-free rejoin audit. This is not a proof
of an Xi cross-current sign, an Abel gap, `Lambda<=0`, PF-infinity, or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_joined_pair_kernel_pfree_rejoin_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_joined_pair_kernel_pfree_rejoin_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_joined_pair_kernel_pfree_rejoin_gate.py
```

## Physical Pair Kernel

On one fixed ray cell put

```text
w_n=eta q_n,
Gamma_n=x[-s_*'log(n)+delta_n]+i Omega_eta
       =a_n+i omega_n.
```

Then `D_j w_n=Gamma_n w_n`. For `n<m`, exact rational
polarization gives the division-free pair current

```text
K_(n,m)
 =Im[(D_jw_n)conjugate(w_m)+(D_jw_m)conjugate(w_n)]

 =(a_n-a_m)Im[q_n conjugate(q_m)]
  +(omega_n+omega_m)Re[q_n conjugate(q_m)].
```

With `s_*'=c+ib` and `delta_n=xi_n+i chi_n`,

```text
a_n-a_m
 =x[c log(m/n)+xi_n-xi_m],

omega_n+omega_m
 =2 Omega_eta+x[-b log(nm)+chi_n+chi_m].
```

The ordered logarithmic rates therefore multiply unrestricted sine and
cosine correlations. Rate ordering alone supplies no sign.

## Moment Collapse

Let

```text
H_0=sum_n q_n,
H_1=sum_n log(n)q_n,
D_0=sum_n delta_nq_n.
```

Summing every diagonal and pair row gives exactly

```text
K_bulk
 =Omega_eta|H_0|^2
  -xc Im[H_1 conjugate(H_0)]
  -xb Re[H_1 conjugate(H_0)]
  +x Im[D_0 conjugate(H_0)].
```

The checker reconstructs both sides over exact rational complex
arithmetic. The recurrent endpoint and its cross terms still have to be
added before an endpoint-complete estimate.

## Heat-Shifted P-Free Telescope

For

```text
A_t(n;s)=exp[t log(n)^2/4-s log(n)],
h=log(p),
M_k=floor(N/p^k),
sigma_k=s-t k h/2,
```

put

```text
T_k=A_t(p^k;s)S_(M_k,t)(sigma_k),
L_k=A_t(p^k;s)O_(M_k,t)^(p)(sigma_k).
```

The exact heat-shift identity

```text
A_t(p^k;s)A_t(p;s-t k h/2)=A_t(p^(k+1);s)
```

gives

```text
L_k=T_k-T_(k+1),
sum_k L_k=T_0=S_(N,t)(s).
```

This is a genuine linear telescope. The phase current is quadratic:

```text
K(T_0)
 =sum_k K(L_k)+sum_(k<l)K_(k,l).
```

Thus retaining every cross term reconstructs the original full-prefix
moment current. P-free reindexing is exact, but it does not create a
new algebraic sign.

## Exact Cancellation Guard

One two-layer exact jet has

```text
L_0=i,       D L_0=-1,
L_1=-i/3,    D L_1=1.
```

It telescopes linearly to `T_0=2i/3`, `D T_0=0`, while

```text
sum_k K(L_k)=4/3,
K_(0,1)=-4/3,
K(T_0)=0/1.
```

The cross current cancels the entire positive diagonal reserve. This is
the same cubic-contact algebra seen by the preceding polarization gate,
now written as an exact telescoping-layer audit. It is not an Xi
counterexample.

## Route Decision

The old prime symmetry remains useful as an exact organization of the
coefficients, but its linear telescope does not reduce the quadratic
joined current below the original `H_0,H_1,D_0` correlation. Any further
gain must come from a new Xi-specific estimate for those oscillatory
correlations together with the recurrent endpoint. Without such a
theorem, the more robust target is the direct endpoint-complete,
division-free Abel-scalar gap.

## Pi Provenance

No `pi` enters the pair-kernel or p-free telescope identities. Any `pi`
in the Xi saddle scale remains the ordinary constant inherited from
completed-zeta normalization.

## Boundary

This gate proves the physical division-free pair kernel, exact bulk moment collapse, correction-free heat-shifted p-free telescope, quadratic rejoin, and strict cross-cancellation guard. It proves no Xi-specific pair-correlation sign, recurrent-endpoint absorption, Abel-scalar gap, horizontal successor winding cap, contact exclusion, Q209, cofinal descendant theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion.
