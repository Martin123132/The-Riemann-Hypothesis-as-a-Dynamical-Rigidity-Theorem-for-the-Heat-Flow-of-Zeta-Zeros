# Physical P_lin Centered-Basis Morse Conjugation Gate

Date: 2026-08-02

Status: exact centered-to-lambda translation, physical Morse conjugation, joined alpha/beta c-prime formulas, and quantitative outside-tail perturbation separation. This is not a proof of the grouped ideal-cubic remainder estimate, signed flow, or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_basis_conjugation_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_basis_conjugation_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_basis_conjugation_gate.py
```

## Exact Basis Map

Put ell=log(a), x=lambda-ell, and

```text
P(lambda)=sum_(k=0)^5 p_k x^k=sum_(j=0)^5 q_j lambda^j.
q_j=sum_(k=j)^5 binom(k,j)(-ell)^(k-j)p_k, where ell=log(a).
p_j=sum_(k=j)^5 binom(k,j)ell^(k-j)q_k.
```

q_j=P^(j)(0)/j! and p_j=P^(j)(ell)/j! for P(lambda)=sum p_k(lambda-ell)^k.

The six forward rows are

```text
q_0=-ell**5*p_5 + ell**4*p_4 - ell**3*p_3 + ell**2*p_2 - ell*p_1 + p_0
q_1=5*ell**4*p_5 - 4*ell**3*p_4 + 3*ell**2*p_3 - 2*ell*p_2 + p_1
q_2=-10*ell**3*p_5 + 6*ell**2*p_4 - 3*ell*p_3 + p_2
q_3=10*ell**2*p_5 - 4*ell*p_4 + p_3
q_4=-5*ell*p_5 + p_4
q_5=p_5
```

## Tail Translation

Khat_k(ell)=sum_(j=0)^k binom(k,j)K_j ell^(k-j), K=(6,14,55,336,2738,27936).

```text
Khat_0(ell)=6
Khat_1(ell)=6*ell + 14
Khat_2(ell)=6*ell**2 + 28*ell + 55
Khat_3(ell)=6*ell**3 + 42*ell**2 + 165*ell + 336
Khat_4(ell)=6*ell**4 + 56*ell**3 + 330*ell**2 + 1344*ell + 2738
Khat_5(ell)=6*ell**5 + 70*ell**4 + 550*ell**3 + 3360*ell**2 + 13690*ell + 27936
```

Therefore the old lambda-basis estimate gives

```text
|mathcal T[P]|<2h^2 sum_(k=0)^5 Khat_k(ell)|p_k|.
```

For a centered monomial P=p_k(lambda-ell)^k, the old coefficientwise majorant is exactly 2h^2 Khat_k(ell)|p_k|. Thus no smaller universal coefficientwise weight follows from that template.

## Physical Morse Conjugation

On the physical q=1 chart, L>=50, a^2=exp(L)+1/(32L^2), ell=log(a), h=1/a, and h<exp(-25)<1/72000000000. Hence 0<ell<L and h^2<exp(-L).

```text
B_a=sigma-(t/2)ell=1/2-(t/2)delta_a, ell=log(a), delta_a=ell-Re(alpha).
S(ell+x)=S(ell)+S_tilde(x), S_tilde(x)=t*x^2/4-B_a*x.
g(ell+x)=t*x/2-B_a=:g_tilde(x).
```

Thus every explicit ell coefficient cancels from the differential operator:

```text
P'(lambda)+g(lambda)P(lambda)=Ptilde'(x)+g_tilde(x)Ptilde(x).
L_B[P]=P''+(2g_tilde+1)P'+(g_tilde^2+g_tilde+t/2+1/6)P.
L_B[P]=exp(-S_tilde)*(D_x^2+D_x+1/6)[exp(S_tilde)P].
```

x_0=lambda_0-ell=log(alpha_P/(a*r)); log(a) enters the sampled interval through x_0, but not as an extra coefficient in L_B.

## Joined Rows

F=(C,mathcal N,A,Q), F_f=f dot F, and P_lin=alpha dot F+i*x*beta dot F.

```text
(D_x+g_tilde)P_lin=alpha dot[(D_x+g_tilde)F]+i beta dot{x(D_x+g_tilde)F+F}.
L_B[P_lin]=alpha dot L_B[F]+i beta dot{x L_B[F]+2F'+(2g_tilde+1)F}.
```

For the full off-saddle numerator:

```text
With A_z=zJ^2/v, B_z=zJJ'-J, E_x=exp(S_tilde(x)), and E_0=exp(S_tilde(x_0)), factor exp(S(ell)) from the c-prime numerator. The remaining centered bracket is alpha dot K_alpha+i beta dot K_beta, where K_alpha=E_x[A_z(D+g_tilde)F+B_zF]+E_0F(x_0) and K_beta=E_x[A_z{x(D+g_tilde)F+F}+B_zxF]+E_0x_0F(x_0). Thus c-prime equals exp(S(ell)) times this bracket divided by r*sqrt(alpha_P)*z^2.
```

This is the direct grouped formulation promised in Section 11.175: no scalar current has been converted into an unproved row norm.

## Quantitative Tail Separation

R=||alpha||_1+||beta||_1, with alpha and beta the actual joined retained-observation rows. No numerical bound on R is assumed.

The detailed Section 11.175 envelopes imply the safe coarse bounds

```text
|bar(gamma)_j|<=R and |bar(eta)_j|<=R; the Section 11.175 errors give Gamma_j<2R and H_j<2R.
|p_0-bar(gamma)_0|<8764h^2R
|p_1-bar(gamma)_1|<8765h^2R
|p_2-bar(gamma)_2|<8762h^2R
|p_3-bar(gamma)_3|<8760h^2R
|p_4|<(2/3)h^2R
|p_5|<h^2R/63
```

For every positive-coefficient polynomial W of degree at most five, L W'(L)<=5W(L), so exp(-L)W(L) decreases for L>=50. Therefore h^2 Khat_k(ell)<Khat_k(50)/(72000000000)^2.

At L=50 the translated weights are

```text
Khat_0=6
Khat_1=314
Khat_2=16455
Khat_3=863586
Khat_4=45394938
Khat_5=2390362436
```

Let Pbar=sum_(j=0)^3bar(gamma)_j(lambda-ell)^j. Then |T[P_lin-Pbar]|<h^2 R/300000000000, while |T[Pbar]|<2h^2 sum_(j=0)^3 Khat_j(ell)|bar(gamma)_j|.

The certified rational coefficient multiplying h^2 R is

```text
8752727719/2916000000000000000000
<1/300000000000.
```

All correction-induced coefficient departures, including the degree-four and degree-five channels, are negligible in the raw outside tail relative to the retained row scale. The unsuppressed ideal cubic is the remaining coefficient-aware tail theorem.

## Next Target

Estimate the ideal cubic with bar(gamma)_0,...,bar(gamma)_3 left inside the joined Morse/endpoint phase sum. The exact alpha/beta vector kernels above are the admissible starting point. A separate numerical norm on either retained row is optional, not assumed.

The four quadratic observation pairings remain a separate obligation, and the fixed source normalizer nu and common phase c_xi must be restored before comparison with the signed flow benchmark.

## Proof Boundary

This gate proves no numerical retained-observation bound, grouped ideal-cubic c-prime sum, reciprocal cancellation theorem, quadratic residual bound, h^2 flow error, Phi_B upper bound, contact exclusion, retained aggregate or Xi theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion.

## Reproduce

```powershell
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_basis_conjugation_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_basis_conjugation_gate.py
```
