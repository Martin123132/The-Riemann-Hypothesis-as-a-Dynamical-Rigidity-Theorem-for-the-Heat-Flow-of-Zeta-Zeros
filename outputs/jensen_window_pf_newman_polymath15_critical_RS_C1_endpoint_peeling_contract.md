# Jensen-Window PF Newman Polymath-15 Critical C1 Endpoint Peeling Contract

Date: 2026-07-26

Status: exact first-coefficient extraction and critical-line
cancellation. The uniform second-order heat remainder is open.
This is not a proof of `Lambda <= 0` or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_critical_RS_C1_endpoint_peeling_contract.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_critical_RS_C1_endpoint_peeling_contract.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_critical_RS_C1_endpoint_peeling_contract.py
```

## Source Normalization

[Polymath 15](https://arxiv.org/abs/1904.12438), Proposition 6.2, writes

```text
R_0,N(s)=Dirichlet_N(s)+(-1)^(N-1)U*a^(-sigma)*(sum_(k=0)^K C_k(p,sigma)/a^k+RS_K(s))
```

The current [FLINT Riemann-Siegel documentation]
(https://flintlib.org/doc/acb_dirichlet.html) and implementation expose the coefficient
assembly suppressed in the Polymath notation:

```text
C_k(p,sigma)=pi^(-2k)*sum_(j=0)^floor(3k/2)(pi/(2i))^j*d_j^(k)(sigma)*F^(3k-2j)(p)
d_0^(0)=1, d_1^(0)=0 => d_0^(1)=1/12, d_1^(1)=(1-2sigma)/2
```

Substitution gives the exact first omitted coefficient

```text
C_1(p,sigma)=F'''(p)/(12*pi^2)+i*(2sigma-1)*F'(p)/(4*pi), with F=C_0
```

The identity is source-derived, not fitted from numerical data.

## Critical Cancellation

At the real Newman boundary the zeta argument has `sigma=1/2`, so

```text
C_1(p,1/2)=F'''(p)/(12*pi^2)
C_1(p,1/2+sqrt(t)*v)=F'''(p)/(12*pi^2)+i*sqrt(t)*v*F'(p)/(2*pi)
For dmu(v)=pi^(-1/2)exp(-v^2)dv, integral C_1(p,1/2+sqrt(t)*v)dmu(v)=F'''(p)/(12*pi^2)
```

Thus the first coefficient's fluctuating `F'` term is odd under
the centered heat Gaussian. This cancellation applies to the
frozen-prefactor coefficient layer. The exact gamma/Stirling and
`log M_0` defects must be retained at the next order.

## Explicit Peel

Define the added critical endpoint block by

```text
DeltaC_t(x)=2*(-1)^N*exp(t*pi^2/64)*Re(M_0(i*T_0)*U*exp(pi*i/8)*F'''(p)/(12*pi^2*a))
DeltaQ=DeltaC_t/A_t; Q_[1]=Q_[0]+DeltaQ; J_[1]=P-Q_[1]=J_[0]-DeltaQ
Z=J_[0]+r_[0]=J_[1]+r_[1], hence r_[0]=-DeltaQ+r_[1] and r_[0],x=-DeltaQ_x+r_[1],x
```

The last line is the key contact-normal handoff: relative to the
old corrected main, the explicit signed remainder model is
`r_0=-DeltaQ`, with derivative `r_0,x=-DeltaQ_x`.

## Order Gain

```text
a=sqrt(T_0/(2*pi)); a^(-1)=sqrt(2*pi/T_0), a^(-2)=2*pi/T_0
After retaining C_1, the fixed-sigma Riemann-Siegel tail begins at a^(-2)=2*pi/T_0; the heat-integrated proof must also retain the O(T_0^-1) gamma/Stirling and log-M_0 Taylor terms
```

The old endpoint error starts at `a^-1`. Keeping `C_1` removes
that term. The published fixed-sigma expansion then starts at
`a^-2`, the same order as the first gamma/Stirling and
`log M_0` Taylor defects. This is an asymptotic order audit, not
yet the uniform heat-integrated bound needed by the proof.

## Cutoff Parity

The critical endpoint coefficient has the parity needed at a
Riemann-Siegel cutoff:

```text
F(-p)=F(p) => C_1(-p,1/2)=-C_1(p,1/2)
At a=m, p_left=-1, p_right=1, and (-1)^(m-1)C_1(-1,1/2)=(-1)^m C_1(1,1/2)
For E_N(T)=(-1)^N*a^(-1)*C_1(p_N(T),1/2), p_N'(T)=-a/T and [E_N']_(right-left)=-(-1)^m*F''''(1)/(6*pi^2*T)
```

Thus switching from the left cutoff to the right cutoff creates
no `a^-1` value jump. The first one-sided derivative mismatch is
already proportional to `T^-1`, at the desired second-order
scale. A holomorphic collar still needs an explicit bound for it.

## Remaining Theorem

```text
Prove explicit value and x-derivative bounds for r_[1] on L>=50 and 0<t*L<=25, uniformly across fixed-N cells and adjacent-cutoff splices, then insert r_0=-DeltaQ into the signed contact-normal peeling inequality
```

Three pieces must be proved together: explicit Gaussian-integrated
constants, the `x` derivative of the new remainder, and a
second-order adjacent-cutoff derivative splice. Only then can the smaller
remainder be inserted into the signed contact equation.
