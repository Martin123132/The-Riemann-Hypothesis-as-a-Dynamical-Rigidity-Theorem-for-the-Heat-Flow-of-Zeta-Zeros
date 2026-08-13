# Physical P_lin Joined-Coefficient Envelope Gate

Date: 2026-08-02

Status: five physical carrier-rate bounds, three ideal-rate deviation bounds, and joined alpha/beta propagation through all six centered P_lin coefficients. This is not a proof of a grouped Morse estimate, signed flow, or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_joined_coefficient_envelope_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_joined_coefficient_envelope_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_joined_coefficient_envelope_gate.py
```

## Carrier Rates

L>=50, 0<tL<=25, 0<t<=1/2, 0<=theta<=1, h=1/a<1/72000000000, and 3<pi<22/7.

r_0=-c*u_N, r_1=-s_*', t_0=-c_x*u_N+i*b*u_(N,x), t_1=-s_*'', and delta=chi_N-i*b*u_N.

The source identities u_N=-log(1-h theta), u_(N,x)=h^2/(8pi), s_*'=-i/2-it alpha'/4, and s_*''=-t alpha''/8 give

```text
|r_0|<h^3/3000
|r_1-i/2|<h^2/6000 and |r_1|<3001/6000
|t_0|<h^2/16
|t_1|<h^4/16
|delta|<4h^2
```

Here pi is the inherited completed-zeta/Riemann-Siegel constant in the saddle coordinate. Only pi>3 is used to obtain u_(N,x)<h^2/24; no new circle or polygon is introduced.

For G_0=r_0(r_0+delta)+t_0, G_1=r_1(2r_0+delta)+t_1, and G_2=r_1^2, the ideal point is (-i u_(N,x)/2,0,-1/4), and

```text
|G_0+i*u_(N,x)/2|<h^4/15
|G_1|<17h^2/8
|G_2+1/4|<h^2/3000
```

## Row-Sensitive Errors

For a real row f=(f_V,f_N,f_A,f_Q), define the ideal coefficients

```text
bar(u)_(f,0)=f_V-i*u_(N,x)*f_N/2, bar(u)_(f,1)=i(f_Q+f_A)/2, bar(u)_(f,2)=-f_N/4.
bar(v)_(f,0)=f_Q, bar(v)_(f,1)=i*f_N/2.
```

and the nonnegative error envelopes

```text
e_(u,0)(f)=h^4|f_N|/15+(4h^2+h^3/3000)|f_Q|+h^3|f_A|/3000
e_(u,1)(f)=17h^2|f_N|/8+h^2(|f_Q|+|f_A|)/6000
e_(u,2)(f)=h^2|f_N|/3000
e_(v,0)(f)=h^3|f_N|/3000
e_(v,1)(f)=h^2|f_N|/6000
```

Then |u_(f,j)-bar(u)_(f,j)|<e_(u,j)(f) and |v_(f,j)-bar(v)_(f,j)|<e_(v,j)(f). The row is formed before these moduli are taken.

## Joined Rows

alpha=(N_(p,xi),V_(p,xi),-Q_(p,xi),-A_(p,xi)) and beta=(N_p+A_(T,x),X_T+V_p,-Q_p-X_(T,x),-A_T-A_p).

At the ideal rate point their composed coefficients are

```text
bar(gamma)_0=alpha_V-i*u_(N,x)*alpha_N/2
bar(gamma)_1=i(alpha_Q+alpha_A)/2+i*beta_V+u_(N,x)*beta_N/2
bar(gamma)_2=-alpha_N/4-(beta_Q+beta_A)/2
bar(gamma)_3=-i*beta_N/4
bar(eta)_0=alpha_Q
bar(eta)_1=i(alpha_N/2+beta_Q)
bar(eta)_2=-beta_N/2
```

The actual errors compose as

```text
E_gamma=(e_(u,0)(alpha), e_(u,1)(alpha)+e_(u,0)(beta), e_(u,2)(alpha)+e_(u,1)(beta), e_(u,2)(beta)).
E_eta=(e_(v,0)(alpha), e_(v,1)(alpha)+e_(v,0)(beta), e_(v,1)(beta)).
Gamma_j=|bar(gamma)_j|+E_(gamma,j), H_j=|bar(eta)_j|+E_(eta,j). Then |gamma_j|<Gamma_j and |eta_j|<H_j.
```

No independent numerical bound for an alpha or beta component has been assumed.

## Six Coefficients

Combining these joined rows with the certified C/D box from Section 11.174 gives

```text
Pi_0=E_(gamma,0)+4379h^2Gamma_0+16893h^4H_0
Pi_1=E_(gamma,1)+4379h^2Gamma_1+h^2Gamma_0/4+16893h^4H_1+5h^4H_0/16
Pi_2=E_(gamma,2)+4379h^2Gamma_2+h^2Gamma_1/4+h^2Gamma_0/16+16893h^4H_2+5h^4H_1/16+h^4H_0/16
Pi_3=E_(gamma,3)+4379h^2Gamma_3+h^2Gamma_2/4+h^2Gamma_1/16+5h^4H_2/16+h^4H_1/16
Pi_4=h^2Gamma_3/4+h^2Gamma_2/16+h^4H_2/16
p_5=i*rho_2*beta_N*(s_*')^2 and |p_5|<h^2|beta_N|/63
```

|p_j-bar(gamma)_j|<Pi_j for 0<=j<=3, |p_4|<Pi_4, and p_5 has the exact displayed fibre factor.

Thus the actual polynomial is a row-sensitive O(h^2) perturbation of its ideal cubic in coefficient space, while the top channel keeps its exact total-value factor. This does not say the retained rows or the cubic coefficients are small.

## Basis Guard

The p_j above multiply y=lambda-log(a). The K_j outside-tail template in Section 11.172 was stated for lambda powers. It may be used only after exact binomial translation or after deriving centered tail weights; substituting the p_j directly would hide powers of log(a).

## Next Target

Either prove a source-specific joined envelope for bar(gamma), Gamma, and H from the endpoint-composed retained carrier, or leave those rows inside the grouped Morse sum and prove a coefficient-aware cancellation estimate directly. A scalar Phi_B estimate, even once available, cannot by itself be treated as a gradient-row norm.

## Proof Boundary

This gate proves no numerical retained-observation bound, centered-to-lambda tail bound, grouped c-prime sum, reciprocal cancellation estimate, quadratic residual bound, h^2 flow error, Phi_B upper bound, contact exclusion, retained aggregate or Xi theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion.

## Reproduce

```powershell
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_joined_coefficient_envelope_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_joined_coefficient_envelope_gate.py
```
