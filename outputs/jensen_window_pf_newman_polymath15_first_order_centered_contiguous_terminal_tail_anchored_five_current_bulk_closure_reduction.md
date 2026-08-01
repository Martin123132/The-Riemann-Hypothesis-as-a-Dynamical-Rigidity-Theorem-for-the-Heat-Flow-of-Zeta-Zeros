# Terminal-Tail-Anchored Five-Current Bulk Closure

Date: 2026-07-31

Status: exact algebraic closure with `0 signed bulk bounds`. This is not a
proof artifact; `Phi_B` and RH remain open.

## Coefficient Current

Use the fixed q=1, fixed-N terminal/bulk partition of the tail-anchored Abel reduction, with B=N-M-1 and no cutoff change during x differentiation.

```text
Let hat(z_n)=z_n/S_a and write eta_x/eta=i*omega_eta. With kappa_0=i*omega_eta-partial_x log(S_a), the exact coefficient law is partial_x hat(z_n)=[kappa_0-s_*'*log(n)+delta_n]hat(z_n).
```

Define the five truncated complex currents

```text
H_j^(B)=sum_(n=1)^B(log n)^j*hat(z_n), j=0,1,2; D_j^(B)=sum_(n=1)^B(log n)^j*delta_n*hat(z_n), j=0,1
```

Their required first derivatives satisfy

```text
H_(0,x)=kappa_0*H_0-s_*'*H_1+D_0; H_(1,x)=kappa_0*H_1-s_*'*H_2+D_1
```

## Value and Abel Moment

```text
hat(F_B)=H_0, partial_x hat(F_B)=kappa_0*H_0-s_*'*H_1+D_0
hat(R_B)=log(N)*H_0-H_1; partial_x hat(R_B)=kappa_0*hat(R_B)-s_*'[log(N)*H_1-H_2]+log(N)*D_0-D_1
```

The `log(N)` coefficient is constant because the calculation stays inside one
fixed-`N` chart. These formulas are equalities of complex currents and retain
all cancellation.

The equivalent terminal-centered basis is better conditioned:

```text
lambda_n=log(N/n), K_j=sum_(n=1)^B lambda_n^j*hat(z_n) for j=0,1,2, and E_j=sum_(n=1)^B lambda_n^j*delta_n*hat(z_n) for j=0,1; K_0=H_0, K_1=log(N)H_0-H_1, K_2=log(N)^2H_0-2log(N)H_1+H_2, E_0=D_0, E_1=log(N)D_0-D_1
With chi_N=kappa_0-s_*'log(N), hat(F_B)=K_0, hat(R_B)=K_1, partial_x K_0=chi_N*K_0+s_*'*K_1+E_0, and partial_x K_1=chi_N*K_1+s_*'*K_2+E_1
```

## Terminal Join

```text
Insert the displayed hat(F_B),partial_x hat(F_B),hat(R_B),partial_x hat(R_B) into A_B=Re(s_*'R_B)-b*u_N*Im(F_B), its exact first derivative, and Phi_B=h^-2[(X_T+X_B)A_(B,x)+X_B*A_(T,x)-(A_T+A_B)X_(B,x)-A_B*X_(T,x)].
```

The full nonterminal input to Phi_B is determined exactly by the five complex currents H_0^(B),H_1^(B),H_2^(B),D_0^(B),D_1^(B), equivalently the centered currents K_0,K_1,K_2,E_0,E_1, the common scalar coefficients, and the four already certified tail coordinates. No individual prefix is needed after those moments have been formed.

This is stronger bookkeeping than the `O(B)` Abel-prefix form: once the five
moments are supplied, evaluating `Phi_B` takes constant additional work.

## Multiplicative Handoff

The joined dyadic odd-prefix identities remain exact with their summation cutoff replaced by B: unique dyadic valuation gives H_0,H_1,H_2 and D_0,D_1 as heat-shifted odd-prefix sums, retaining every external odd phase and correction current.

The remaining theorem can be sought either as a direct signed bound on Phi_B or as a quadratic inequality for this five-current vector after exact dyadic or p-adic rejoining.

The sharp open inequality is

```text
Prove Phi_B<=1/400, preferably Phi_B<=1/800, for the actual physical five-current vector and endpoint-terminal tail.
```

Finite-dimensional closure supplies sufficient data, not a sign. The existing negative interior currents, amplitude contraction, angular ordering, and separate prime-power block windings do not imply the required quadratic inequality without the joined Xi correlations.

## Boundary

The logarithmic moments and dyadic valuation introduce no pi. All pi factors remain inherited from the completed-zeta and Riemann-Siegel normalization already present in the tail jets.

This proves an exact five-current algebraic closure only. It proves no signed bound on Phi_B, bulk aggregate sign closure, complete retained or Xi-level current theorem, Abel gap, winding cap, contact exclusion, Q209, Lambda<=0, PF-infinity, RH, or prize-level conclusion.

## Reproduce

```text
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_five_current_bulk_closure_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_five_current_bulk_closure_reduction.py
```
