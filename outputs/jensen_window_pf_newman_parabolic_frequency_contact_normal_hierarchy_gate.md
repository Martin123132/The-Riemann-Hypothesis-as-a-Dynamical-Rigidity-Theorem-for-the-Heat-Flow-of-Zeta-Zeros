# Newman Parabolic-Frequency Contact-Normal Hierarchy Gate

Date: 2026-07-26

Status: exact contact-normal reduction, exact route guards, and one open
Xi `C1` theorem target. This is not a proof of `Lambda<=0` or RH.

## Corrected Contact Normal Form

Write `H=A Z`, where the Polymath-15 normalizer `A` is positive, and use

```text
Z=J+r,
J=2X-Q, J_x=2(U-BY)-Q_x, T_pf[J]=J^2+s_pf^2*J_x^2,
s_pf=sqrt(2t)/sqrt(1+q), q=2tL^2,
beta=L*s_pf=sqrt(q/(1+q)).
```

At a first-jet contact, `Z=Z_x=0`, so `J=-r` and `J_x=-r_x`. Therefore

```text
for L>=50 and 0<tL<=25: epsilon_0=2500*exp(-3L/4), epsilon_1=5000*exp(-3L/4),
if Z=J+r=Z_x=J_x+r_x=0, |r|<=epsilon_0, |r_x|<=L*epsilon_1, then T_pf[J]<=epsilon_0^2+beta^2*epsilon_1^2.             (1)
a contact requires |J|<=epsilon_0 and |J_x|<=L*epsilon_1.                      (1a)
a contact requires |J|<=2500*exp(-3L/4) and |J_x|<=5000*L*exp(-3L/4).             (1b)
```

The sharp noncircular sufficient theorem is

```text
|2X-Q|>epsilon_0 or |2(U-BY)-Q_x|>L*epsilon_1,                         (2)
```

or, equivalently,

```text
on |2X-Q|<=epsilon_0 prove |2(U-BY)-Q_x|>L*epsilon_1.                    (2a)
```

This requires derivative separation only on the thin value band. The
smooth scaled inequality

```text
(2X-Q)^2+s_pf^2*(2(U-BY)-Q_x)^2 >epsilon_0^2+beta^2*epsilon_1^2                         (2b)
```

is sufficient but stronger than the sharp box exclusion. All three use only
the refined value and first-derivative remainder.

For the endpoint-corrected complex component sum `E_c`, write

```text
for E_c with J=2Re(E_c), put X_c=Re(E_c), U_c=Re(E_c,x); a contact requires |X_c|<=epsilon_0/2 and |U_c|<=L*epsilon_1/2,
U_c=M(h_+-h_-)+(X_c/2)(h_++h_-)+D_0,
set h_+=0 when M_+=0 and h_-=0 when M_-=0; the near-crossing identity then remains valid, including the all-zero-real-part case U_c=D_0.
```

The triangle inequality reduces the frequency-layer band theorem to

```text
|M(h_+-h_-)+D_0|>L*epsilon_1/2+epsilon_0*|h_++h_-|/4,                (2c)
|M(h_+-h_-)+D_0|>625*exp(-3L/4)*(4L+|h_++h_-|).       (2d)
```

This is a weighted arithmetic slope gap only where the corrected value is
inside its remainder band. It is weaker than the previous strip-wide cone
and retains the finite remainder terms omitted by an exact-crossing floor.

At an actual contact the signed remainder gives the sharper phase-matched equation

```text
at a contact, M(h_+-h_-)+D_0=-r_x/2+(r/4)(h_++h_-).          (2e)
```

The magnitude gap (2c)-(2d) is obtained by applying the triangle inequality
to (2e). It is a sufficient fallback, but it discards the phase and derivative
correlation of the actual Riemann-Siegel remainder. The primary
contact-conditioned search should retain (2e) and insert a signed asymptotic
model for `r,r_x`.

More precisely, peeling a signed leading remainder gives the exact reduction

```text
if r=r_0+delta_0, r_x=r_(0,x)+delta_1, |delta_0|<=eta_0, and |delta_1|<=L*eta_1, then contact requires |M(h_+-h_-)+D_0+r_(0,x)/2-(r_0/4)(h_++h_-)|<=L*eta_1/2+(eta_0/4)|h_++h_-|.                (2f)
```

Thus the next analytic deliverable is an explicit Riemann-Siegel
`r_0,r_(0,x)` and smaller certified residual `eta_0,eta_1`, followed by

```text
|M(h_+-h_-)+D_0+r_(0,x)/2-(r_0/4)(h_++h_-)|>L*eta_1/2+(eta_0/4)|h_++h_-|.                 (2g)
```

The unpeeled magnitude target is the special case `r_0=r_(0,x)=0`.

There is also an essential endpoint quantifier guard:

```text
if at t=0 a contact has |J(0)|<epsilon_0 and |J_x(0)|<L*epsilon_1, continuity keeps J,J_x inside the same open box for all sufficiently small t>0; therefore a full-wedge pointwise C1-box exclusion rules out endpoint multiplicity and is stronger than Lambda<=0 or RH.        (2h)
```

Thus a pointwise box exclusion uniform down to `t=0` may still prove the
desired result, but it is stronger than the multiplicity-compatible theorem.
The cofinal boundary-degree and delta-localized transport routes remain live.

## Exact Chart Split

```text
q<=1 iff t<=1/(2L^2) iff c=tL<=1/(2L); q>=1 is the frequency chart and contains every fixed 0<c<=c_*+epsilon as L->infinity.                            (3)
```

Thus `q<=1` is the ultra-small layer `t<=1/(2L^2)`. Every asymptotic ray
with fixed positive `c=tL` eventually lies in the frequency chart. The
proof search should not treat the whole critical layer as parabolic.

For the frequency layer the scale parameter can be removed exactly:

```text
put A0=J^2-epsilon_0^2 and A1=(J_x/L)^2-epsilon_1^2; for q>=1, D_pf=A0+beta^2*A1>=min(A0+A1/2,A0+A1), so positivity of both endpoint deficits suffices.              (3a)
```

It is enough to prove the two endpoint deficits at derivative weights
`1/2` and `1`; no continuum of beta values needs a separate theorem.

## Why Scale Alone Cannot Close Theorem

At `Z=Z_x=0`,

```text
for H=A Z and a=(log A)_x, at Z=Z_x=0: H_xx=A Z_xx, H_xxx=A(Z_xxx+3a Z_xx),
partial_t(H,s_pf H_x)=-A(Z_xx,s_pf*(Z_xxx+3a Z_xx)) at a contact.                     (4)
```

The term `s_(pf,t) H_x` vanishes exactly at a contact. Consequently a
time-dependent scale conditions the approach to a collision but does not
regularize the collision itself. A direct proof of (2) uses a `C1`
remainder, whereas a bare relative heat-jet proof exposes `Z_xx,Z_xxx`
and needs `C3` control or an independently closed factorization.

The exact Newman heat flow

```text
H_tau(y)=y^2-2tau
```

has, at `y=0`,

```text
||partial_tau(H_tau,sH_(tau,y))||/
||(H_tau,sH_(tau,y))||=1/|tau|.
```

The pole is nonintegrable through the double contact for every smooth
positive scale. This is an exact generic guard, not an Xi counterexample.

## Correlation-Hierarchy Translation

The first layers satisfy

```text
F_1(2x)=H_x^2-H H_xx; F_2(2x)=3H_xx^2/4-H_x H_xxx+H H_xxxx/4,
F_2-3*partial_xi^2 F_1=H H_xxxx-H_x H_xxx.                        (5)
```

At a nondegenerate double contact,

```text
H=H_x=0, H_xx!=0 implies F_1=F_1'=0, F_1''=H_xx^2/4, F_2=3H_xx^2/4, partial_t F_1=H_xx^2,
at a nondegenerate double contact, F_1*F_3-F_2^2=-9*H_xx^4/16<0.                    (6)
```

A tempting order-Hankel condition `F_1 F_3-F_2^2>=0` would therefore
exclude contact. It is not a viable generic closure: for the real-rooted
quadratic heat flow,

```text
H_tau(y)=y^2-2tau solves H_tau=-H_yy; F_1=2y^2+4tau, F_2=3, F_3=0, F_1F_3-F_2^2=-9, and at y=0 ||partial_tau(H,sH_y)||/||(H,sH_y)||=1/|tau|.                  (7)
```

Hence correlation-order log-convexity is not even necessary for a
Laguerre-Polya polynomial. Any surviving hierarchy inequality must use
additional Xi/theta structure.

Dimitrov and Xu prove that the first Wronskian is positive definite as a
function because its inverse Fourier correlation is nonnegative. Positive
definiteness is not pointwise strict positivity and permits real zeros; their
density characterization becomes an RH-equivalent condition for the Xi
kernel rather than the missing implication.

## Route Decision

The direct contact route is the signed peeling target (2g), divided into:

```text
q<=1: an ultra-small-time Hermite/parabolic contact theorem;
q>=1 and 0<tL<=c_*+epsilon: prove (2g), with the robust band slope gap
as a stronger sufficient fallback;
bounded L: a terminating finite shoulder theorem.
```

Weighted crossing slopes remain an admissible way to prove the second line
at the contact equations themselves. A strip-wide cone is no longer needed.
An endpoint-uniform floor is not logically required, and no proposed
majorant may contain the inverse unknown first-jet norm.

Primary sources:

```text
https://arxiv.org/abs/1606.05011
https://arxiv.org/abs/1801.05914
https://arxiv.org/abs/1904.12438
```

## Proof Boundary

The exact identities and countermodels above do not prove (2a), the
ultra-small-time Xi theorem, the critical asymptotic Xi theorem, the
bounded-L shoulder, `Lambda<=0`, RH, PF-infinity, or the Clay prize.

Machine-audited files:

```text
work/rh_compute/results/jensen_window_pf_newman_parabolic_frequency_contact_normal_hierarchy_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_parabolic_frequency_contact_normal_hierarchy_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_parabolic_frequency_contact_normal_hierarchy_gate.py
```
