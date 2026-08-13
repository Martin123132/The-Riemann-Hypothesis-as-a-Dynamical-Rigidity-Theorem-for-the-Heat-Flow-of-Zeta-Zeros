# Physical P_lin Correction-Coefficient Envelope Gate

Date: 2026-08-02

Status: exact quotient identities and certified physical envelopes for the six terminal-centered C/D coefficients. This is not a proof of a physical P_lin norm, signed Morse sum, or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_correction_coefficient_envelope_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_correction_coefficient_envelope_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_correction_coefficient_envelope_gate.py
```

## Source Separation

L>=50, 0<tL<=25, 0<t<=1/2, h=1/a<1/72000000000, x>a^2=h^-2, and 16892h^2<1.

To avoid the collision with the centered coefficient d_1, write the original first Dirichlet correction as frak_d(ell). The critical-ray gate supplies the alpha, chi, and ell=A correction bounds; the all-carrier gates separately supply the ell=0 denominator and derivative bounds.

```text
|alpha'|<h^2, |alpha''|<h^4, |chi|<2. The same normalized correction majorant applies at ell=A because alpha-A=chi, giving |frak_d(A)|<h^2/4 and |frak_d_x(A)|<h^4/32.
|frak_d(0)|<2189/x<2189h^2<1/2 and |frak_d_x(0)|<4223/x^2<4223h^4.
```

## Exact Quotients

Write the original first Dirichlet correction as frak_d(ell)=1/(6s)+alpha'(s){t/4+t^2[alpha(s)-ell]^2/8}. Put A=log(a), chi=alpha-A, a_0=1+frak_d(0), and B=alpha'(s)t^2/8.

```text
rho_2=B/a_0,                 rho_1=-2*alpha*rho_2,
B_x=-i*alpha''*t^2/16,
rho_(2,x)=B_x/a_0-B*frak_d_(1,x)/a_0^2.
```

a_0=1+frak_d(0), so |a_0|>1/2, |a_0|^-1<2, and |a_0|^-2<4.

Therefore

```text
|B|<h^2/32
|B_x|<h^4/64
|rho_2|<h^2/16
|rho_(2,x)|<h^4/32+(4223/8)h^6<h^4/16
```

The last step uses 16892=4*4223 and 16892h^2<1. Thus both candidate quotient bounds are now certified rather than merely suggested.

## Centered Coefficients

Put A=log(a), chi=alpha-A, and shift ell=A+y. The exact rows are

```text
c_0=[1+frak_d(A)]/a_0,
c_1=-2*chi*rho_2,             c_2=rho_2,

d_0=frak_d_x(A)/a_0-[1+frak_d(A)]frak_d_(1,x)/a_0^2,
d_1=-2*alpha_x*rho_2-2*chi*rho_(2,x),
d_2=rho_(2,x),                alpha_x=-i*alpha'/2.
```

Consequently, without separately majorizing log(a) or log(a)^2,

```text
|c_0-1|<4379h^2<1/3; hence 2/3<|c_0|<4/3
|c_1|<h^2/4
|c_2|<h^2/16

|d_0|<16893h^4
|d_1|<5h^4/16
|d_2|<h^4/16
```

The 16893 constant is intentionally retained. It is small after multiplication by h^4, but it must remain visible when the p_j convolution is bounded.

## Next Target

Insert these six envelopes into the exact gamma/eta recurrence of the physical P_lin gate. Bound p_0,...,p_5 only after the physical observation coefficients, R, R_x, and delta are joined; retain the exact top factor p_5=i*rho_2*(X_T+V_p)*(s_*')^2 and the total-value fibre.

## Proof Boundary

This gate proves no P_lin norm, grouped c-prime sum, reciprocal cancellation estimate, quadratic residual bound, h^2 flow error, Phi_B upper bound, contact exclusion, retained aggregate or Xi theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion.

## Reproduce

```powershell
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_correction_coefficient_envelope_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_correction_coefficient_envelope_gate.py
```
