# Full-Support Terminal Endpoint Obstruction Gate

Date: 2026-08-02

Status: exact terminal endpoint specialization and unique order-one obstruction proved; signed kernel open; not a proof of RH.

This is not a proof of RH. It identifies which part of the surviving full-support endpoint trace is genuinely unsuppressed.

## Terminal Specialization

At lambda=log N, R=i b u_N, R+delta=chi_N, and R_x=i(b_xu_N+b u_(N,x)).

C_N=C(-u_N), A_N=i b u_NC_N, Q_N=chi_NC_N+D_N, and mathcal N_N=i b u_NQ_N+i(b_xu_N+b u_(N,x))C_N.

F_(alpha_N)(log N)=mathcal N_(p,xi)C_N+V_(p,xi)mathcal N_N-Q_(p,xi)A_N-A_(p,xi)Q_N.

In the correction-free chart F_(alpha_N)^0(log N)=mathcal N_(p,xi)+V_(p,xi)(-u_N^2/4-i u_(N,x)/2)+i u_N[Q_(p,xi)+A_(p,xi)]/2.

## Row Bounds

|C_N-1|<4380h^2 and |C_N|<3/2.

|D_N|<16894h^4.

|A_N|<8h/5, |Q_N|<8h/5, and |mathcal N_N|<17h^2/10.

|F_(alpha_N)(log N)-mathcal N_(p,xi)|<(8/5)h R_alpha, with R_alpha the l1 norm of alpha_N.

## Opposite Endpoint

0<S_1^[N](alpha_P)<5/N<=(15/2)h.

At u=1, a_1>3alpha_P/5 and 0<b_1-a_1<1+alpha_P/(N+1/2); integrating psi'=zeta(2,.) gives 0<S_1^[N](alpha_P)<5/N<=(15/2)h.

## Conditional Quadratic Traces

For ell_mu=log(mu/a), the first conditional endpoint trace obeys B_mu^(1)[dot P]=i ell_mu B_mu^(1)[P].

At mu=N, ell_N=-u_N with |ell_N|<2h. At mu=1, ell_1=-log a while B_1^(1) already contains S_1^[N](alpha_P)=O(h).

Every S_1-by-S_1 quadratic endpoint term contains at least one terminal u_N factor or one lower-endpoint S_1^[N](alpha_P) factor. This is an h-suppression template, not a complete quadratic bound.

## Handoff

The sole order-one terminal conditional coefficient is mathcal N_(p,xi); all other terms in F_(alpha_N)(log N) are bounded by (8/5)h R_alpha.

Derive the source-specific reciprocal-band/Hermitian/transpose cancellation or sign for the mathcal N_(p,xi) S_1^[N](alpha_P/N) endpoint kernel.

The next calculation should start from `mathcal N_(p,xi)S_1^[N](alpha_P/N)`, not from a four-row endpoint norm. Compose it with the returned `q=N` carrier, the exact terminal recurrence, `H_out`, `T_out`, correction perturbation, and the moved-tail defect before taking a modulus.

## Pi Provenance

No new `pi` is introduced. The factor in `u_(N,x)=h^2/(8pi)` is inherited from the completed-zeta/Riemann--Siegel saddle, and all Fourier factors retain the fixed character `e(x)=exp(2pi i x)`.

## Proof Boundary

This gate proves exact terminal carrier specialization, uniform terminal channel bounds, a row-relative collapse of `F_(alpha_N)(log N)` to `mathcal N_(p,xi)`, an `O(h)` lower-endpoint `S_1` bound, and h suppression of every purely conditional quadratic endpoint product. It proves no sign or cancellation for the surviving kernel, complete outer-complement bound, quadratic residual bound, signed current, contact exclusion, retained aggregate or Xi theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion.
