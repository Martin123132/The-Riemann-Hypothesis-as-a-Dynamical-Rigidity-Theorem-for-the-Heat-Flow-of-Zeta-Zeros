# Beta^-6 common-profile leading coefficient

Date: 2026-08-13

Status: exact formal coefficient and rigorous whole-strip coefficient bound;
the exact post-beta-minus-six Taylor remainder remains open

Put `q=beta^-1`, `X=lambda+y`, and factor the exact transformed integrand by
the cubic carrier `exp(i[z^3/3-Xz])`.  Its exact ratio has the expansion

```text
A(q,z,y)exp(iR(q,z,y))
 =1+q^2 C1(z,y)+q^4 C2(z,y)+q^6 C3(z,y)+O(q^8).       (B6.1)
```

The third amplitude and phase coefficients are

```text
a3=z**4*(390*y - 379*z**2)/1920,
r3=-z**5*(189*y**2 - 306*y*z**2 + 124*z**4)/5670.                                         (B6.2)
```

Exact multiplication gives

```text
C3=a3+i r3-r1 r2-i r1^3/6+a1(i r2-r1^2/2)+i a2 r1.  (B6.3)
```

`C3` has degree fifteen in `z`.  Contracting every monomial with

```text
integral_R^Abel z^k e^(i[z^3/3-Xz])dz
 =2pi i^k d_X^k Ai(-X)                               (B6.4)
```

and reducing by `A_XX=-X A` yields exactly

```text
J6(lambda,y)=2pi beta^-6
 [U3(lambda,y)Ai(-lambda-y)
  +V3(lambda,y)d_X Ai(-lambda-y)],                    (B6.5)

U3=(-49856*lambda**6 + 1344*lambda**5*y + 3360*lambda**4*y**2 - 12920*lambda**3*y**3 - 516535*lambda**3 + 7785*lambda**2*y**4 + 5445*lambda**2*y - 1836*lambda*y**5 + 21195*lambda*y**2 + 1674*y**6 - 21285*y**3 - 265860)/9072000,
V3=-(3584*lambda**7 - 1792*lambda**6*y + 1344*lambda**5*y**2 + 2240*lambda**4*y**3 - 43880*lambda**4 - 1960*lambda**3*y**4 + 18580*lambda**3*y + 1764*lambda**2*y**5 + 5445*lambda**2*y**2 - 567*lambda*y**6 + 14130*lambda*y**3 - 127530*lambda + 189*y**7 - 7605*y**4 + 42570*y)/9072000.                                         (B6.6)
```

A `4096`-panel Arb interval atlas encloses the complete top corridor
`0<=lambda<=pi/(16beta)` and `0<=y<=Y`.  It proves

```text
sup |J6(lambda,y)|
 <[6.9498853069067562883133126672886812503536403230029219969649379870047272287496615e-10 +/- 1.59e-88]<6.951e-10,

sup |g6(lambda,s)|=sup |J6|/h
 <[1.2869285917193211541260434073194526400813577442888722008237867161832441315249874e-8 +/- 3.59e-87]<1.288e-8.           (B6.7)
```

Passing this coefficient through the already certified whole weighted
kernel gives

```text
|Delta W6|
 <[4.9421694684043838438164891397630610956837085959570407024596212982435451005057152e-13 +/- 1.40e-91]<4.95e-13.       (B6.8)
```

Equations (B6.7)--(B6.8) rigorously bound the coefficient of the first
omitted asymptotic order; they are not a bound on the entire exact-minus-
beta-minus-four profile until the `O(q^8)` remainder in (B6.1) is enclosed.
No extrapolation from coefficient size is promoted.

Pi provenance: `2pi` in (B6.4) is the inverse Airy Fourier normalization;
`beta^3=pi C^2/8`, `Y=pi C/(2beta)`, and `hY=2pi` are inherited from the
exact Kummer/Fourier selector geometry.  No fitted or geometric constant is
introduced.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta6_common_profile_leading_coefficient_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta6_common_profile_leading_coefficient_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta6_common_profile_leading_coefficient_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta6_common_profile_leading_coefficient_gate.py
```

No exact Taylor-tail estimate, full exact-minus-beta-minus-four profile
bound, weighted finite-`t` remainder, complete source/initial-data splice,
complete `Q_K-T` or `T_upper`, all-corridor or height-uniform theorem,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
