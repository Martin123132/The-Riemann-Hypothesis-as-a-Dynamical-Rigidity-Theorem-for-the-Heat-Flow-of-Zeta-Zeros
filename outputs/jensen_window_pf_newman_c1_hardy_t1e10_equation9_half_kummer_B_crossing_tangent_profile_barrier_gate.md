# Tangent half-plane model for the B crossing

Date: 2026-08-13

Status: countermodel gate; not a proof or bound for the nonlinear half-Kummer residual

On the half-domain from Section 11.352, scale the universal Morse coordinate
by

```text
y=sqrt(t)s/2,
q_B(s)=q_B*+kappa_B y                              (tangent model),
Delta_B=1-pi*kappa_B^2/2.
```

After factoring the classical full-line density, the complete grouped alpha
current has tangent density

```text
1+c_m q,       c_m=sqrt(x_m)/(m sqrt(2)).                (BT1)
```

The second term is the exact current corresponding to the sine/pivot-
derivative channel.  It must not be dropped.

Put `q=q_B*+kappa_B y` and then `rho=q-kappa_B y`.  The joint phase completes
exactly:

```text
-y^2+(pi/2)(rho+kappa_B y)^2
 =-Delta_B[y-pi kappa_B rho/(2Delta_B)]^2
   +pi rho^2/(2Delta_B).                                (BT2)
```

Let `chi=1` for `Delta_B>0`, `chi=i` for `Delta_B<0`, and let

```text
T_eps(u)=integral_u^infinity exp(i eps pi v^2/2)dv.
```

The excluded upper half-plane has normalized tangent profile

```text
U_0=chi T_sign(Delta)(q_B*/sqrt(|Delta|))/(1+i),

U_1=-c_m chi exp[i pi q_B*^2/(2Delta)]
     /[sqrt(|Delta|) i pi(1+i)].                        (BT3)
```

Relative to the sharp classical step at `621|622`, the exact tangent-model
residual is

```text
R_B,m=1_(m<=621)-(U_0+U_1).                             (BT4)
```

No integer mode hits the tangent degeneracy.  The nearest is `m=257`, with

```text
|Delta_B|=[0.00096321334567788656289760368248151529242077262178337657733264743981896396891153469667311302175762836983551716220 +/- 2.29e-110]>0.0009.
```

Rigorous Arb summation of (BT4) over all positive lower-branch modes
`1..39894`, using the exact Riemann--Siegel theta carrier, gives

```text
2 Re sum carrier*R_B,m
 =[2.1781802282775436547425866559010984359154087466055481780792798870163172547375127981253760181453792560485754925e-5 +/- 2.02e-99],

(sqrt(2)/pi) times this
 =[9.8052559952454123160727648576535368252225665922533861892879373794608287111321279836504055114115651092549781154e-6 +/- 9.08e-100].                (BT5)
```

Both exceed their respective targets `0.000019` and `0.0000086`; the
normalized ratio is

```text
[1.1464106464618652227593065170198563681135117686693547511822544038295745849609375000000000000000000000000000000 +/- 2.05e-16].
```

The mechanism is instructive.  If the `c_m q` channel is illegally omitted,
the same signed diagnostic would be

```text
[8.7255373867829377336033659713278013945430416410994702341135265462385194059652836694219401478408508184464270915e-6 +/- 7.56e-100],                       (BT6)
```

which lies below the normalized target.  The retained q-density contributes

```text
[1.3056264895992498813822500587683182964611045824956011546679272323924653141409844311831820033612941742039327834e-5 +/- 1.27e-99].                         (BT7)
```

Thus dropping the sine/current term creates a false pass.  The complete
termwise triangle is `[0.0027205181080351455065682523237039197741822794469118693188607318141099646014725175445770411697194111079817477324 +/- 1.20e-99]`, so the result
also depends on substantial signed cancellation.

This is a barrier only to promoting the tangent crossing by itself.  It is
not a lower bound for the exact half-Kummer error: nonlinear Morse amplitude,
the `m<=0` and `m>=39895` nonstationary completion, the endpoint half-current,
and the characteristic A/half-boundary chart remain outside (BT4).  They must
be assembled before a final sign or magnitude is claimed.

Pi provenance: `pi` in (BT1)--(BT7) comes from the equation-(9) Kummer
quadratic, the Fourier-Poisson character, Gaussian completion, and
Riemann--Siegel normalization.  No fitted constant is used.

Proof boundary: exact tangent half-plane transform and a saved-height finite
aggregate of that model only.  No nonlinear B-crossing remainder, complete
symmetric half-domain completion, A-fold splice, complete `T_upper`,
height-uniform theorem, `Lambda<=0`, PF-infinity, RH, or prize-level
conclusion is proved.
