# Newman Endpoint Odd-Fibre Correlation Feasibility Guard

Date: 2026-07-28

Correction notice (2026-07-31): the fixed `(T_0,1)` endpoint direction is
superseded by the rotation into the actual complex endpoint vector in the
complex endpoint source-normalization gate and Formal Core Section 11.146.
The finite-current insufficiency guards remain exact.

Status: exact correlation audit and current-input no-go;
not a proof of an actual Xi pivot, `Lambda<=0`, PF-infinity,
RH, or a Clay-prize result.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_endpoint_odd_fibre_correlation_feasibility_guard.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_endpoint_odd_fibre_correlation_feasibility_guard.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_endpoint_odd_fibre_correlation_feasibility_guard.py
```

## Pi Provenance

The pi in T_0, beta, and the endpoint phase cancellation comes from the completed-zeta normalization and the Riemann-Siegel saddle. The phase restoration uses only the unit factor phi=M_t(s)/|M_t(s)|. The odd/even split, finite difference guard, and correlation rotation introduce no new pi.

## Physical Rephasing

```text
Put phi=M_t(s)/|M_t(s)| and f_n=phi*exp[t*log(n)^2/4-s_*log(n)]*(1+d_n). Then f_1=phi*(1+d_1).
For S_odd=sum_(m<=N,m odd)f_m, the odd numerator satisfies S_odd=phi*O_N.
The phase-cancelled recurrent endpoint is g_0=(-1)^N*beta*(T_0+i)*H_a/|M_t(s)|=phi*R_N. Define the real scalar gamma_N=(-1)^N*beta*H_a/|M_t(s)|, so g_0=gamma_N*(T_0+i).
For p_K=2^K, f_(p_K)=phi*exp(i*omega_*K*h)*C_K and |f_(p_K)|=|C_K|=rho_K.
```

## First Pivot

```text
Delta_K={|S_odd+g_0|^2-|f_(p_K)|^2}/|1+d_1|^2. Thus the first Schur pivot is positive iff the physical odd-index subtotal plus endpoint has modulus larger than the terminal power-of-two carrier.
Re(S_odd*conj(g_0))=gamma_N*(T_0*X_odd+Y_odd).
|S_odd+g_0|^2=(X_odd+gamma_N*T_0)^2+(Y_odd+gamma_N)^2. Hence the pivot numerator is (X_odd+gamma_N*T_0)^2+(Y_odd+gamma_N)^2-|f_(p_K)|^2.
With D=sqrt(T_0^2+1), P_odd=(T_0X_odd+Y_odd)/D and Q_odd=(-X_odd+T_0Y_odd)/D, |S_odd+g_0|^2=(P_odd+gamma_ND)^2+Q_odd^2.
```

If H_a=0, then gamma_N=g_0=0 and the target reduces to |S_odd|>|f_(p_K)|. Endpoint phase cancellation supplies no odd-fibre lower bound in this case.

## Total Versus Odd Fibre

```text
With S_even=sum_(n<=N,n even)f_n, E_[1]=S_odd+S_even+g_0 and S_odd+g_0=E_[1]-S_even. Existing control of the total real value, centered scalar, or adjacent total projections does not by itself bound this difference.
```

## Five-Current Null Family

```text
In the correction-free joined coefficient model with N>=32, perturb q_3,q_6,q_12,q_24 by lambda*(1,-3,3,-1). sum_(k=0)^3 w_k(log(3)+kh)^r=0 for r=0,1,2; the r=3 sum is -6h^3. D_0=D_1=0 remains unchanged, the endpoint and terminal power-of-two carrier are fixed, but Only n=3 is odd, so the k=0 odd fibre changes by lambda.
The complete five-current data H_0,H_1,H_2,D_0,D_1, together with a fixed endpoint and fixed terminal carrier, do not determine the first pivot in the unrestricted correction-free joined class. Choosing lambda to cancel the odd-endpoint subtotal makes the pivot numerator negative; choosing |lambda| sufficiently large makes it positive.
```

## Linked Two-Jet Guard

```text
F_lambda(w)=1+lambda*(w-1)^3 satisfies F_lambda(1)=1, F_lambda'(1)=F_lambda''(1)=0 for every real lambda. Yet lambda=1/4 gives Delta_3=1/2. while lambda=1 gives Delta_3=-1. Thus even a linked value and two z derivatives do not determine a global Schur pivot.
```

The moment-null and linked-two-jet families use unrestricted coefficient perturbations. They prove insufficiency of the currently retained observables, not failure of the actual Xi coefficient family. Positivity, heat-shift, correction, cutoff, and endpoint relations may help only through a new proved Xi-specific inequality.

## Feasibility Verdict

The source audit removes the phase mystery: after restoring phi, the endpoint is the real scalar gamma_N times T_0+i, and the missing correlation is exactly gamma_N*(T_0X_odd+Y_odd). No current theorem controls that odd-fibre projection at the terminal-carrier scale. The five-current closure and the retained linked two-jet are algebraically insufficient without additional Xi coefficient structure. Full Schur disk stability is therefore downgraded from the primary proof route to a falsification coordinate.

## Retained Schur Role

Retain Delta_K for a later one-worker physical scout or for testing any source-derived odd-fibre inequality. Do not expand the second pivot unless a uniform actual first-pivot theorem is proved. A selected positive table cannot reverse this route decision.

## Primary Target

```text
Return to the actual linked point. On q>=1, prove |mathsf_X|<=delta_L implies |mathcal_C_N|>A_L+epsilon_term, including W_0=0, and prove that crossings mathsf_X=0 with mathcal_C_N>0 number less than one successor turn after complete boundary composition. This uses the endpoint-complete first jet without demanding a coefficient fibre dominate on every free disk point.
The q=2tL^2<1 multiplicity-compatible parabolic/Hermite chart and finite connectors remain separate. This feasibility audit does not close them.
```

## Boundary

The physical rephasing, odd-fibre and endpoint subtotal, terminal-carrier identity, Cartesian and rotated correlation forms, total/even decomposition, correction-free five-current null family, and linked-two-jet guard are exact. The verdict is a no-go for the current Schur inputs, not a proof that no future Xi-specific first-pivot theorem exists. No uniform actual Xi pivot, physical pivot failure, signed prefix lower bound, crossing count, strict successor flux theorem, q<1 closure, contact exclusion, Lambda<=0, PF-infinity, RH proof, or Clay-prize conclusion is asserted.
