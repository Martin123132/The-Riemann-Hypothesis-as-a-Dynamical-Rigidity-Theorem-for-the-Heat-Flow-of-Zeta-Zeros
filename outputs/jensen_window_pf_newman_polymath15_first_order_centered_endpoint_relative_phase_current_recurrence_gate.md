# Jensen-Window PF Newman Polymath-15 Endpoint-Relative Phase-Current Recurrence Gate

Date: 2026-07-31

Status: exact endpoint-relative phase-current identities, a rigorous
`q=1` terminal sign-reversal guard, and an exact recurrence-defect
reduction. This is not a proof of RH or `Lambda <= 0`.
This is not an Abel gap or contact-exclusion theorem.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_endpoint_relative_phase_current_recurrence_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_endpoint_relative_phase_current_recurrence_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_endpoint_relative_phase_current_recurrence_gate.py
```

## Endpoint Self-Current

For the corrected complex endpoint,

```text
T_(0,x)=1/2, kappa in R, H_a in C, e=kappa(T_0+i)H_a

I_e:=Im(e_x*conj(e))=kappa^2{(T_0^2+1)Im(H_(a,x)*conj(H_a))-|H_a|^2/2}

|e|^2=kappa^2(T_0^2+1)|H_a|^2
```

The `kappa_x` term is radial and cancels. On the ordinary fibre,

```text
On H_a!=0, partial_x arg(e)=Im(H_(a,x)*conj(H_a))/|H_a|^2-1/{2(T_0^2+1)}
```

The equivalent moving-frame form is

```text
For e_x=g+(lambda_a-rho_1)e, lambda_a=u_a+i*v_a, g=kappa(T_0+i)J_a^(der), I_e=kappa^2(T_0^2+1){Im(J_a^(der)*conj(H_a))+v_a|H_a|^2}

Im(mu_a)+v_a=-1/{2(T_0^2+1)}; the radial normalizer cancels from I_e
```

Thus the result retains the full complex `H_a` and `J_a^(der)`.

## Relative Carrier Current

Write

```text
For z_n=r_n*zeta_n and delta_n=Im(d_(n,x)/(1+d_n)), Im(z_(n,x)*conj(z_n))=(v_a+b*u_n+delta_n)|z_n|^2

Omega_(0n):=|e|^2 Im(z_(n,x)*conj(z_n))-|z_n|^2 Im(e_x*conj(e))
```

Direct cancellation of the common frame phase gives

```text
Omega_(0n)=kappa^2(T_0^2+1)|z_n|^2{|H_a|^2(b*u_n+delta_n)-Im(J_a^(der)*conj(H_a))}
```

When both atoms are nonzero this is

```text
When e*z_n!=0, Omega_(0n)/(|e|^2|z_n|^2)=partial_x arg(z_n/e)
```

The division-free polynomial identity remains primary on every zero
fibre.

## Exact Sign-Reversal Guard

Take the physical cofinal family

```text
Fix p in [-1,1], theta=(1-p)/2, a=N+theta, q=2tL^2=1, and n=N-m with fixed m>=0; let N->infinity.
```

The already certified source estimates imply

```text
H_a->C_0(p), a H_(a,x)->-C_0'(p)/(4*pi), a*mu_a->0, b->-1/2, a*u_(N-m)->m+theta, and a*delta_(N-m)->0.

a*Omega_(0,N-m)/(|e|^2|z_(N-m)|^2) -> Phi(p)/(4*pi)-m/2-(1-p)/4
```

For the terminal carrier `m=0`, the two exact witnesses are

```text
C_0'(0)=0, so the m=0 limit at p=0 is -1/4.

At p=1, Phi(1)/pi=1-sqrt(4+2sqrt(2))/4, so the m=0 limit is {1-sqrt(4+2sqrt(2))/4}/4>0.
```

Therefore

```text
Along the physical q=1 family the endpoint-terminal relative current has both signs for all sufficiently large cutoffs. A uniform one-sided terminal orientation is false.
```

This rejects a uniform terminal relative-current sign before any broad
numerical calibration. It does not say that contact conditioning is
irrelevant.

## Tail Composition Before Projection

The one-step recurrence is

```text
Q_(N-1):=z_N+e_N-e_(N-1), hence e_N+z_N=e_(N-1)+Q_(N-1)

Delta W_(A,N-1):=s_*'u_N z_N+g_N-g_(N-1), hence g_N+s_*'u_Nz_N=g_(N-1)+Delta W_(A,N-1)
```

Iteration gives

```text
e_N+sum_(j=0)^m z_(N-j)=e_(N-m-1)+sum_(k=N-m-1)^(N-1)Q_k

g_N+s_*'sum_(j=0)^m u_(N-j)z_(N-j)=g_(N-m-1)+sum_(k=N-m-1)^(N-1)Delta W_(A,k)
```

The value recurrence does not telescope its phase current:

```text
If H_+=J-H, then I(H_+)=I(J)+I(H)-Im(J_x*conj(H)+H_x*conj(J)); the mixed current does not telescope.
```

Moreover,

```text
Existing adjacent-chart certificates control the needed real traces, not a complex norm/current of Q_k. Re Q=Re Q_x=0 alone permits either current sign: E=1, Q=i*y, Q_x=i*sigma gives I(E+Q)=sigma.
```

So a complex recurrence-defect current bound, or a direct real-edge
theorem, is a genuinely new input.

## Contact-Minor Sum

The corrected endpoint minors satisfy

```text
sum_(n=1)^N K_(0n)=c_0*mathcal C_N-d_0*mathsf X

On mathsf X=0, sum_(n=1)^N K_(0n)=c_0*mathcal C_N.
```

On the ordinary endpoint fibre this recovers the complete contact scalar.
On the exceptional fibre,

```text
When c_0=0=mathsf X, the minor sum is zero even if mathcal C_N!=0. The corrected zero-real-projection fibre therefore still needs a division-free direct theorem.
```

## Route Decision

```text
Retire a uniform endpoint-terminal relative-current sign and any claim that the C_0 recurrence alone telescopes phase currents.

The exact tail composition may still be useful if one proves a complex recurrence-defect current bound, or if the terminal-composed real edge block is treated directly.

Derive the division-free projective current of (C_edge,D_edge)=(c_0+c_N,d_0+d_N), retaining the second endpoint jet and d_(N,x). Test whether it is negative; otherwise return to the complete contact scalar and seek a signed cumulative-minor estimate with a separate zero-real-projection theorem.
```

## Proof Boundary

This gate has 16 rows, two exact endpoint-current forms, one
division-free relative current, two opposite asymptotic sign witnesses,
one terminal-tail recurrence, and one contact-minor sum. It proves zero
uniform terminal signs, zero Abel gaps, and zero winding bounds. The
terminal-composed real-edge current and the complete signed contact
scalar remain open.
