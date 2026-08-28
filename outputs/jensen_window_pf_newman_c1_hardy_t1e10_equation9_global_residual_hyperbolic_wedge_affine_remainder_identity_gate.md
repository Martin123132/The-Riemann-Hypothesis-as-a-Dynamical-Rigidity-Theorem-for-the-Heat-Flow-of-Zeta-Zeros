# Exact affine-wedge remainder identity

Date: 2026-08-13

Status: exact transformed-amplitude and signed-domain decomposition; not a
bound for either remainder integral

In the exact bi-Morse variables, including both Jacobians, the endpoint-`D`
pair triangle is

```text
eps_D exp(i Phi_D,m) W0_D,m (2pi)^(-1)
 int_(S>h) int_(P<P_D(S)) A_D,m(P,S)
 exp(i(P^2-S^2)/2)dP dS,                              (AR1)

W0_D,m=4(pi Dm-i)r^(3/4)/(D sqrt(pi t)),
r=t/(2pi m^2).                                        (AR2)
```

The normalized transformed amplitude factorizes exactly:

```text
A_D,m(P,S)=E_D,m(u(P)) O(v(S)),
E_D,m(u)=2u(pi Dm u-i)/[(u+1)(pi Dm-i)],
O(v)=v^(3/4)s/(v-1).                                  (AR3)
```

Both factors equal one at the joint saddle.  Their complete affine tangent is

```text
A_aff(P,S)=1+lambda_P P+mu_S S,
lambda_P=(3pi Dm-i)/[2(pi Dm-i)sqrt(pi mD)],
mu_S=(5/12)sqrt(2/t).                                 (AR4)
```

This retains the mandatory endpoint-density current rather than silently
dropping it.  If `C=C_rho(a,h)` is the wedge from Section 11.427, its affine
moment is still exact:

```text
C_aff=C+i mu_S partial_h C
        -i(lambda_P+mu_S rho)partial_a C.             (AR5)
```

The exact-minus-affine-tangent defect now splits before norms as

```text
R_exact-aff=R_amp+R_face,                              (AR6)

R_amp=(2pi)^(-1) int_(S>h) int_(P<P_D(S))
      [A_D,m(P,S)-A_aff(P,S)]exp(i(P^2-S^2)/2)dP dS,

R_face=(2pi)^(-1) int_(S>h) int_(a+rho S)^(P_D(S))
       A_aff(P,S)exp(i(P^2-S^2)/2)dP dS.              (AR7)
```

The inner integral in `R_face` is oriented, so (AR6) remains valid when the
curved face crosses its tangent.  There is no phase remainder.

At the six saved-height transition modes, interval arithmetic gives

```text
|P_D''(0)|<7e-6,   Re(lambda_P)<2e-5,
|Im(lambda_P)|<1e-14,   mu_S<6e-6.                    (AR8)
```

For orientation, the near-corner values are

```text
A,39894: P_D''(0)=[-2.3569885968163113766078221018357698738126210190650855647052064472303310797623521e-6 +/- 4.18e-87],
A,39895: P_D''(0)=[-2.3571362991659320606778609000058031561461270055195888556246461542298802986455538e-6 +/- 3.34e-86].
```

These small tangent coefficients are promising but are not global integral
bounds.  The next certificate must put a finite local box around the A
corner, bound the second-order amplitude and face strip there, and dispatch
its complement by phase-adapted integration by parts.  The B half-boundary
tail can be treated separately because Section 11.427 puts it over 270598
standardized units away.

Pi provenance: every `pi` in (AR1)--(AR8) comes from equation (9), the exact
Morse scalings, and the Fourier Gaussian.  No fitted constant is introduced.

Proof boundary: exact coordinate Jacobians, amplitude factorization, affine
wedge moments, signed remainder identity, and local tangent coefficients
only.  No global `R_amp` or `R_face` bound, `R_Dir` estimate, complete
`Q_K-T` or `T_upper`, all-height theorem, `Lambda<=0`, PF-infinity, RH, or
prize-level conclusion is proved.
