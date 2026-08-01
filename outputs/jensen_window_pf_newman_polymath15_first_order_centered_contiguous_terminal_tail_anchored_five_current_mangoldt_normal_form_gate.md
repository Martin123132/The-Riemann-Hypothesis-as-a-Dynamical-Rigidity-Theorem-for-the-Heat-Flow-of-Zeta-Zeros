# Terminal-Tail Five-Current Mangoldt Normal Form

Date: 2026-07-31

Status: exact correction-free and Mangoldt normal form with `0 signed Type-I/II
bounds`. This is not a proof artifact; `Phi_B` and RH remain open.

## Exact Correction Transform

Use the fixed physical q=1 terminal/bulk partition with B=N-M-1. All sums below stop at B and share the same eta/S_a normalizer as the certified tail.

```text
c_n=(1+d_n)/(1+d_1)=1+rho_1*log(n)+rho_2*log(n)^2
delta_n=partial_x log(c_n), so delta_n*c_n=rho_(1,x)*log(n)+rho_(2,x)*log(n)^2
G_r^(B)=(eta/S_a)sum_(n=1)^B log(n)^r*exp[t log(n)^2/4-s_*log(n)], 0<=r<=4
```

From 1+d_n=a_0+a_1 log(n)+a_2 log(n)^2 and a_0=1+d_1, put rho_1=a_1/a_0 and rho_2=a_2/a_0. Their x derivatives are exact quotient derivatives; no correction remainder is discarded.

The corrected five-current vector is

```text
H_0=G_0+rho_1G_1+rho_2G_2; H_1=G_1+rho_1G_2+rho_2G_3; H_2=G_2+rho_1G_3+rho_2G_4; D_0=rho_(1,x)G_1+rho_(2,x)G_2; D_1=rho_(1,x)G_2+rho_(2,x)G_3
```

The corrected five-current vector needed by Phi_B is an exact linear image of G_0,G_1,G_2,G_3,G_4. Thus the correction does not create an uncontrolled sixth current.

## Mangoldt Form Through Order Four

```text
For 1<=r<=4, G_r^(B)=(eta/S_a)sum_(dm<=B)Lambda(d)log(dm)^(r-1)q_(dm)^(0)=(eta/(2S_a))sum_(dm<=B)[Lambda(d)+Lambda(m)]log(dm)^(r-1)q_(dm)^(0)
q_(dm)^(0)=q_d^(0)q_m^(0)exp[(t/2)log(d)log(m)]
```

The formal divisor audit checked `205` coefficient
rows and `205` symmetric rows through cutoff
`97`.

For D=floor(sqrt(B)), each symmetric Mangoldt moment splits exactly into the square d,m<=D and the two wings with exactly one variable <=D. The endpoint-terminal block is composed with the resulting moments before any absolute values.

## Exact Guard

For every 1<=r<=4 and any prime sqrt(B)<p<=B, the {1,p} principal minor of the real transpose-symmetric Mangoldt kernel is [[0,log(p)^r/2],[log(p)^r/2,0]] and has negative determinant. Transpose symmetry and entrywise positivity therefore do not supply the signed Phi_B bound.

The four determinants are

```text
{
  "1": "-log(p)^2/4<0",
  "2": "-log(p)^4/4<0",
  "3": "-log(p)^6/4<0",
  "4": "-log(p)^8/4<0"
}
```

Thus the old visual symmetry has a precise role as a Type-I/II coordinate,
but it is not a positivity theorem.

## Open Arithmetic Target

```text
Insert the five-current transform into the exact Phi_B quadratic form and prove Phi_B<=1/400, preferably <=1/800, by a signed endpoint-coupled Type-I/II or Vaughan estimate for G_1,...,G_4 together with the value moment G_0.
```

## Boundary

No pi is introduced by the correction transform, Mangoldt identity, or balanced hyperbola. Existing pi factors retain their completed-zeta/Riemann-Siegel origin.

This proves an exact five-moment correction-free normal form, Mangoldt representations through order four, balanced hyperbola organization, and an indefiniteness guard. It proves no Type-I/II cancellation estimate, signed bound on Phi_B, bulk aggregate sign closure, complete retained or Xi-level current theorem, Abel gap, winding cap, contact exclusion, Q209, Lambda<=0, PF-infinity, RH, or prize-level conclusion.

## Reproduce

```text
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_five_current_mangoldt_normal_form_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_five_current_mangoldt_normal_form_gate.py
```
