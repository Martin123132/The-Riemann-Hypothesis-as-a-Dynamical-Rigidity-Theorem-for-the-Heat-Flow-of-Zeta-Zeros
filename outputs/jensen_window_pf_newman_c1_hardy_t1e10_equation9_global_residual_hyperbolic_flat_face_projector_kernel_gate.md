# Hyperbolic flat-face projector kernel

Date: 2026-08-13

Status: exact canonical profile, B-dictionary, and joint-fold criterion;
not a bound for the curved-face or amplitude remainder

In the exact bi-Morse variables, normalize the flat tangent face as

```text
Phi(P,S)=(P^2-S^2)/2,       P<a+rho S,
c=rho^2-1.                                           (FF1)
```

For `c!=0`, differentiate the Abel-regularized half-plane integral with
respect to `a`.  Completing the square gives

```text
Phi(a+rho S,S)
 =c[S+a rho/c]^2/2-a^2/(2c),                         (FF2)

F_c'(a)=exp(i sgn(c) pi/4) exp[-i a^2/(2c)]
        /sqrt(2 pi |c|).                              (FF3)
```

The normalization is fixed by central symmetry and the full-plane Gaussian:

```text
F_c(a)=1/2+integral_0^a F_c'(r)dr,
F_c(-infinity)=0,       F_c(+infinity)=1.            (FF4)
```

Thus a straight hyperbolic face is exactly a complex Fresnel smoothing of a
projector step.  As `c` tends to zero, (FF3) is an oscillatory approximate
identity and (FF4) tends distributionally to `1_(a>0)`.  At `c=0` itself the
flat face is characteristic, so curvature is the next nonvanishing datum.

The existing B tangent notation is exactly the same kernel.  Put

```text
P=sqrt(pi)q,      S=sqrt(2)y,
a=sqrt(pi)q_B*,   rho=kappa_B sqrt(pi/2),
c=pi kappa_B^2/2-1=-Delta_B.                         (FF5)
```

Then the complementary upper profile `1-F_c(a)` is precisely

```text
U_0=chi T_sign(Delta_B)(q_B*/sqrt(|Delta_B|))/(1+i), (FF6)
```

with `chi=1` for `Delta_B>0` and `chi=i` for `Delta_B<0`.  This proves an
exact dictionary rather than a resemblance.  The affine q-density current
`U_1` remains a separate mandatory amplitude term.

The apparent B slope degeneracy is remote:

```text
m=257: c=[-0.00096321334567788656289760368248151529242077262178337657733264743981896396891153470 +/- 3.33e-84],
a/sqrt(|c|)=[-1889601.1334804816085384098870053148805374372104252929303464387454285191390799098 +/- 4.30e-74].         (FF7)
```

So the face is over `1.8e6` canonical widths from the saddle there.  It is
not a joint fold.  A true fold requires both face incidence and null slope:

```text
a=0 and c=0
 iff alpha_m=D and m=D/4
 iff t=pi D^2/8.                                     (FF8)
```

The A endpoint is the near realization of (FF8).  At its occupancy edge,

```text
m=39852: a/sqrt(|c|)=[-0.044965822283705807492462957730887240674253336755886155083659219249059977136950333 +/- 4.06e-82],
m=39853: a/sqrt(|c|)=[0.069923033069980478389103207569453898740307915685596552791847283378335128596345547 +/- 4.92e-82]. (FF9)
```

At the half-boundary bracket,

```text
m=39894: c=[6.2665739085848621125804364673098811428099067711555368331915967423319074248146463e-6 +/- 2.94e-86],
m=39895: c=[-1.8799722628127119247836956436522493028139012880850372956003044307809712860626366e-5 +/- 1.24e-85]. (FF10)
```

Hence the B edge admits an ordinary Fresnel face correction, while the A
edge reaches the codimension-two null geometry where curved-face terms must
be matched to the certified fold atlas.  This is the canonical reason those
two endpoints cannot share one crude absolute estimate.

Pi provenance: `pi` is inherited from the equation-(9) Fourier/Kummer phase
and the standardized Gaussian coordinates.  Equations (FF3)--(FF6) derive
it from that normalization; no fitted circle constant is inserted.

Proof boundary: exact Abel-regularized flat-face kernel, tangent-profile
dictionary, saved-height scale audit, and joint-fold criterion only.  No
curved-face or transformed-amplitude remainder, signed projector-defect
bound, `R_Dir` estimate, complete `Q_K-T` or `T_upper`, all-height theorem,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
