# Newman Endpoint Schur-Cohn First-Jet Guard

Date: 2026-07-27

Status: exact endpoint stability coordinate and first-jet
recursion with countermodels; not a proof of `Lambda<=0`,
PF-infinity, RH, or a Clay-prize result.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_endpoint_schur_cohn_first_jet_guard.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_endpoint_schur_cohn_first_jet_guard.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_endpoint_schur_cohn_first_jet_guard.py
```

## Pi Provenance

The pi in a^2=x/(4*pi)+t/16 comes from the completed-zeta normalization and the Riemann-Siegel saddle; 4*pi*(2N+1) is the cutoff-cell difference and 2*pi is the period of exp(i*theta). Schur-Cohn uses the unit circle |z|=1 and introduces no new pi.

## Endpoint Polynomial

```text
At xi=omega_* and mu=1 put B_(k,r)=sum_(m<=M_k,m odd)(log(2^k m))^r a_(k,m)exp(i omega_*log m)c_(2^k m), E_(k,r)=sum_(m<=M_k,m odd)(log(2^k m))^r delta_(2^k m)a_(k,m)exp(i omega_*log m)c_(2^k m). Then P_r(z)=sum_k B_(k,r)z^k, D_r(z)=sum_k E_(k,r)z^k, and at z_*=exp(i omega_*log 2), P_r(z_*)=H_r and D_r(z_*)=D_r. Define F_N(z)=r_0+P_0(z)=sum_(k=0)^K b_kz^k with b_0=r_0+B_(0,0) and b_k=B_(k,0) for k>=1.
```

## Outside-Disk Schur Step

```text
For degree d, F#(z)=z^d*conj(F(1/conj(z)))=sum_(k=0)^d conj(b_(d-k))*z^k.
S_dF(z)=conj(b_0)F(z)-b_dF#(z)=sum_(k=0)^(d-1)[conj(b_0)b_k-b_dconj(b_(d-k))]z^k.
Delta_d=|b_0|^2-|b_d|^2.
Delta_d F(z)=b_0 S_dF(z)+b_d z (S_dF)#(z), where the latter # has degree d-1.
```

After trimming exact zero top coefficients, F_d has no zero in |z|<=1 iff Delta_d>0 and S_dF_d has no zero in |z|<=1. Indeed, on |z|=1 one has |F_d#|=|F_d|. Rouche applied to S_dF_d=conj(b_0)F_d-b_dF_d# gives the forward zero count; the inverse identity and |z(S_dF_d)#|=|S_dF_d| give the reverse count. Iteration ends at a nonzero constant.

## Reflection Recursion

```text
For p_d(0)=1 and alpha_d=[z^d]p_d, p_(d-1)=(p_d-alpha_d p_d#)/(1-|alpha_d|^2), and p_d=p_(d-1)+alpha_d z p_(d-1)#.
Normalize each effective degree-d polynomial by its nonzero constant and write alpha_d for its leading coefficient. Then F_N is zero-free in |z|<=1 iff |alpha_d|<1 at every effective stage d=K,...,1. These alpha_d are the exact outside-disk reflection coefficients.
```

## Conditional Boundary Margin

```text
The inverse normalized recursion gives, for |z|=1, |p_d(z)|>=(1-|alpha_d|)|p_(d-1)(z)|. Consequently min_(|z|=1)|F_N(z)|>=|b_0|*product_(d=1)^K(1-|alpha_d|)>0 whenever every pivot is strict. This is a quantified sufficient boundary margin, not a claim that the actual Xi reflection coefficients satisfy the inequalities.
```

## Physical Coefficient Current

```text
Inside one fixed-N chart, (B_(k,r))_x=-s_*'B_(k,r+1)+E_(k,r)-i*k*log(2)*omega_*'B_(k,r). Thus the fixed-z tangent is X(z)=F_(N,x)|_z=r_(0,x)-s_*'P_1(z)+D_0(z)-i*omega_*'log(2)zP_0'(z). At z=z_*, adding z_*'=i*omega_*'log(2)z_* gives X(z_*)+z_*'F_N'(z_*)=r_(0,x)-s_*'H_1+D_0=Z_(0,x).
```

## Centered Companion

```text
Let A=log(a) and G(z)=r_A-s_*'[P_1(z)-A P_0(z)]. Then G(z_*)=Z_A. Its fixed-z x tangent is W(z)=r_(A,x)-s_*''[P_1-AP_0]-s_*'{-s_*'[P_2-AP_1]+D_1-AD_0-A_xP_0-i*omega_*'log(2)z[P_1'-AP_0']}. After adding z_*'G'(z_*), the phase terms cancel and one recovers the exact Z_(A,x) formula containing H_2,D_0,D_1,r_(A,x).
```

## Tangent Recursion

```text
For V=sum v_kz^k, (DS_b[V])_k=conj(v_0)b_k+conj(b_0)v_k-v_dconj(b_(d-k))-b_dconj(v_(d-k)).
D2S_b[U,V]_k=conj(u_0)v_k+conj(v_0)u_k-u_dconj(v_(d-k))-v_dconj(u_(d-k)).
For b_(x,nu)=b+xU+nuV+xnuW, partial_(x,nu)S(b_(x,nu))|_0=DS_b[W]+D2S_b[U,V].
Apply DS at every strict Schur stage to U=X and V=G. Apply DS[W]+D2S[X,G] to the mixed x/companion tangent. This propagates Z_0,Z_A,Z_(0,x),Z_(A,x) through the same exact coefficient recursion before any modulus or triangle inequality is taken.
```

## First Physical Pivot

```text
The first physical pivot is Delta_K=|r_0+B_(0,0)|^2-|B_(K,0)|^2, equivalently |alpha_K|=|B_(K,0)/(r_0+B_(0,0))|<1. Phase-frozen dyadic contraction compares neighboring positive B_k before odd phases, d_n corrections, and endpoint addition. It does not prove this endpoint-sensitive inequality.
```

## Exact Guards

### Degree-one orientation

```json
{
  "interpretation": "The constant-dominant orientation corresponds to a root outside the closed unit disk.",
  "pivot": "Delta_1=9-4=5>0",
  "polynomial": "F(z)=3+2z",
  "root": "z=-3/2"
}
```

### First pivot is not sufficient

```json
{
  "first_pivot": "Delta_2=3/4>0",
  "interpretation": "A strict first pivot is necessary but not sufficient; every recursive pivot is needed for full disk exclusion.",
  "polynomial": "F(z)=1-(9/4)z+(1/2)z^2",
  "reduced_polynomial": "S_2F(z)=3/4-(9/8)z",
  "roots": "z=1/2 and z=4",
  "second_pivot": "Delta_1=-45/64<0"
}
```

### Decreasing moduli with physical-style phases

```json
{
  "boundary_zero": "F(1)=b_0+b_1+b_2=0",
  "coefficients": "b_0=1, b_1=(-23+3i*sqrt(55))/40, b_2=(-17-3i*sqrt(55))/40",
  "first_pivot": "Delta_2=51/100>0",
  "interpretation": "Strictly decreasing coefficient moduli do not replace Schur-Cohn once the coefficient phases are free.",
  "moduli": "|b_0|=1>|b_1|=4/5>|b_2|=7/10"
}
```

### Endpoint pivot collapse

```json
{
  "augmented": "F(z)=-2/3+(2/3)z",
  "boundary_zero": "F(1)=0",
  "endpoint": "r_0=-5/3",
  "interpretation": "The endpoint enters the first pivot and can place a zero on the unit circle even when the prefix is disk-zero-free.",
  "pivot": "Delta_1=0",
  "prefix": "P(z)=1+(2/3)z"
}
```

## Cutoff Update

```text
At N->N+1 write n=2^v m, z_*=exp(i omega_*log 2), and u_n=q_n/z_*^v. The exact polynomial update is F_(N+1)(z)-F_N(z)=j_0+u_n z^v, with u_(n,x)/u_n=-s_*'log n+delta_n-i*v*log(2)*omega_*'. At z=z_* this is q_n+j_0=Q_N/f_1. Hence only the constant and valuation-v coefficients change, but every downstream reflection coefficient may change. The certified adjacent bounds are real projections; no complex pivot-continuity bound is inferred from them.
```

## Route Decision

```text
Retain Schur-Cohn as a sufficient endpoint-complete zero-free route and as an exact coordinate for falsification. Reject coefficient-modulus monotonicity, the first pivot alone, and prefix zero-freeness without the endpoint. The next arithmetic task is to attack Delta_K first using the actual recurrent endpoint and terminal dyadic layer. Only if that survives should the reduced coefficients be expanded into nested odd-prefix sums and tested for subsequent strict pivots.
Full unit-disk stability is stronger than final linked-point nonvanishing. If an actual Xi pivot fails, the recursion still identifies the exact failure stage and its inverse formula. A weaker route must then prove nonvanishing at z_* or a signed three-cylinder crossing budget while retaining the endpoint and complete first-jet tangents; it may not promote selected finite pivots to a disk theorem.
```

## Live Theorem

```text
On q=2tL^2>=1, first prove or refute uniformly |r_0+B_(0,0)|>|B_(K,0)| in every fixed-N chart, with the actual d_n correction and endpoint recurrence. If proved, derive the next normalized coefficient array exactly and seek a structural bound |alpha_d|<=1-epsilon_d whose product margin is strong enough to feed |mathsf_X|<=delta_L => |mathcal C_N|>A_L+epsilon_term and the successor trap 0<=kappa_j<1.
```

The separate small-`q` obligation is

```text
The q=2tL^2<1 layer remains a separate multiplicity-compatible parabolic/Hermite first-jet chart. Schur-Cohn disk stability is not imposed there, and no endpoint simplicity, chart join, or contact exclusion follows from this q>=1 algebra.
```

## Boundary

The endpoint polynomial, outside-disk Schur step and inverse, recursive reflection criterion, conditional boundary margin, fixed-z coefficient current, centered companion, Frechet and mixed tangent recursions, cutoff update, and four algebraic audits are exact. No actual Xi pivot inequality, full unit-disk stability, endpoint-complete Xi lower bound, strict successor flux upper bound, q<1 closure, finite connector or chart join, contact exclusion, Lambda<=0, PF-infinity, RH proof, or Clay-prize conclusion is asserted.
