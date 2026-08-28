# Quadratic-theta Poisson current self-duality

Date: 2026-08-24

Status: exact truncated-theta specialization, normalization, contraction, and
pure first-moment saddle transform certified; endpoint-complete fast evaluator
open

Put `K=L-1=2481421`, `c=A+2s`,

```text
a=xc/2,  b=x/2,

F(K,j;a,b)=K^(-j) sum_(n=0)^K n^j exp(2pi i[an+bn^2]). (PS1)
```

Then exactly `S_0=F(K,0;a,b)` and `S_1=K F(K,1;a,b)`.  This puts the
non-A current in the `j=0,1` case of G. A. Hiary's
[truncated-theta algorithm](https://arxiv.org/abs/0711.5002v4).  The source is
used for the normalization identities and endpoint-complete recursive
architecture; none of its unnamed constants is imported into a corpus bound.

On the physical interval `0<x<=1/2`, `0<b<=1/4` already.  Reducing only `a`
modulo one gives

```text
q=floor(a+2bK)=floor(a+xK)<=(K+1)/2.               (PS2)
```

Every main recursion step therefore contracts by at least one half after the
exact unit/half-period and conjugation normalization.  The eight exact
rational traces terminate in at most 4
recorded rows, with at most 3
main Poisson steps.  This is a parameter-map result; the required remainders
have not been implemented.

The special current has a stronger symmetry than a generic `j=0,1` pair.
Apply Poisson frequency `v` directly to its full quadratic phase:

```text
phi(n;v)=x(c+2n)^2/4-2vn,
n_v=v/x-c/2.                                        (PS3)
```

At the stationary point,

```text
c+2n_v=2v/x,
phi(n_v;v)=cv-v^2/x.                                (PS4)
```

Thus the huge constant phase cancels and the completed stationary channel is
the pure absolute-dual first moment

```text
F_saddle(x,s)
 =2 exp(i*pi/4)x^(-3/2)
   sum_(ceil(xc/2)<=v<=floor(xc/2+xK))
      v exp(i*pi[cv-v^2/x]).                        (PS5)
```

There is no independent zeroth-order dual current in (PS5).  If `v=p+r`, it
reappears only in the tied combination `p F_0+N F_1`, where `N=q-p`; splitting
those two terms before estimating would destroy the derivative-current
cancellation.

In the normalized Hiary coordinate the same fact is

```text
F(K,0)=C_0 F(q,0;a*,b*)+R_0,

F(K,1)=C_0[qF(q,1;a*,b*)-aF(q,0;a*,b*)]/(2bK)+R_1,

cS_0+2S_1
 =C_0[(c-a/b)F(q,0)+(q/b)F(q,1)]
  +cR_0+2KR_1.                                     (PS6)
```

The six floating rows compare (PS5) alone with the validated full-roster
reference current.  The stationary roster uses between
2 and 1240710
terms.  The observed omitted endpoint-plus-nonstationary remainder, divided
by absolute source term mass, ranges from
`4.470607e-07` to
`8.189979e-07`.  These values are route
diagnostics only.  In particular, (PS5) is not promoted as an approximation
until `R_0,R_1` are evaluated with explicit error.

The next implementation stage is narrow: specialize the endpoint-complete
`R_0,R_1` formulas to `j<=1`, preserve the tied current `cR_0+2KR_1`, and
cross-check one complete recursion step against the O(`L`) reference oracle.
The small-`b` branch must retain its Euler--Maclaurin integral and correction
terms rather than silently returning the saddle sum.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_poisson_current_self_duality_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_poisson_current_self_duality_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_poisson_current_self_duality_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_poisson_current_self_duality_gate.py
```

Pi provenance: all `pi` factors in (PS1)--(PS6) come from the inherited
Kummer quadratic phase and the ordinary integer Fourier character in Poisson
summation.  No fitted or geometric occurrence is introduced.

Proof boundary: exact parameter mapping, theta normalization identities,
one-half recursion contraction, and completed stationary-current self-duality
only.  The six remainder values are floating diagnostics.  No
endpoint/nonstationary remainder evaluator, complete fast theta algorithm,
uniform or interval error theorem, physical quadrature, non-A bound, joined
`R_after_A`, `R_Dir`, `Q_K-T`, all-height theorem, `Lambda<=0`, PF-infinity,
RH, or prize-level conclusion is proved.
