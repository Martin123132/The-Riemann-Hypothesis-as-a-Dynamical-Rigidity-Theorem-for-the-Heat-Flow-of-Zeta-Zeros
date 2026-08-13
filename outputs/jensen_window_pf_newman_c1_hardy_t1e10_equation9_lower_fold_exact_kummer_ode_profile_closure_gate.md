# Exact Kummer ODE whole-profile closure

Date: 2026-08-13

Status: rigorous full exact-minus-beta-minus-four profile and weighted
finite-height correction on the top selector corridor; this is not a complete
`T_upper` theorem or a proof of RH

Put `rho=1+y/(2beta^2)`, `t=beta^3-beta*lambda`, and remove the common phase
`exp(i*y^2/(4beta))` from the transformed fold.  The exact logistic change

```text
x=(1-tanh(z/beta))/2
```

gives, without asymptotic expansion,

```text
H_ex=beta*sqrt(2)*rho*exp(-i*kappa/2) B(a,b)
      *1F1(a;3/2;i*kappa),
a=3/4-i*t/2, b=3/4+i*t/2,
kappa=2beta^3 rho^2.                                  (KO1)
```

Kummer's differential equation then reduces exactly to

```text
H_ex''+[lambda+y+y^2/(4beta^2)]H_ex=0.                (KO2)
```

The beta-minus-four Airy profile `H4` has no hidden higher residual.  Exact
symbolic reduction in the basis `A=Ai(-X)`, `A_X=d_X Ai(-X)`, `X=lambda+y`,
proves

```text
H4''+[lambda+y+y^2/(4beta^2)]H4
 =R4=2pi*beta^-6[P3 A+Q3 A_X],                        (KO3)
P3=-y**2*(448*lambda**5 + 280*lambda**2*y**3 + 4565*lambda**2 - 105*lambda*y**4 + 30*lambda*y + 63*y**5 - 405*y**2)/201600,
Q3=-y**2*(-40*lambda**3 + 20*lambda**2*y - lambda*y**2 + 9*y**3 - 27)/6720.                                         (KO4)
```

Thus `E=H_ex-H4` obeys the exact forced equation

```text
E''+(lambda+y)E=-beta^-2*y^2*E/4-R4.                  (KO5)
```

No post-beta-minus-six Taylor remainder occurs in (KO5).  An `8192`-
panel Arb atlas on the complete rectangle proves

```text
sup sqrt(Ai(-X)^2+Bi(-X)^2)
 < 0.71573533490300178527832031250000000000000000000000000000000000000000000000000000,
sup |K_Airy(y,s)| < [1.6093658785480690035013337657559609182474846747270979209822851031879221252732877 +/- 2.66e-80],
Integral_0^Y |R4|dy < [1.0005547158740391427161514925107014037237307524437661321076018598275987301151501e-8 +/- 1.54e-86].  (KO6)
```

The exact completed-point and initial-slope gates supply

```text
sup |E(lambda,0)| < [5.2015699622773613719436081916328840497224413599616158520871834091266840839769547e-15 +/- 3.54e-95],
sup |E_y(lambda,0)| < [5.5505800957709680328443899261493410406837270177161196910388693030314087077616985e-22 +/- 3.94e-102]. (KO7)
```

The Volterra feedback factor is
`[0.045506831649232752925283307320983560612624577161871299500757745043840109061692002 +/- 4.93e-82]`, leaving denominator
`>[0.95449316835076724707471669267901643938737542283812870049924225495615989093830800 +/- 4.01e-81]`.  Therefore

```text
sup_(lambda,y)|J_ex-J4|=sup|E|
 < [1.6870306442946823524180217749542396105370988532743790541849983524665594451204251e-8 +/- 4.12e-88] < 1.69e-8,
sup_(lambda,s)|g_ex-g4|
 < [3.1239191373300502730545514946908700555645929029201112700290193087976528076181091e-7 +/- 3.63e-87] < 3.13e-7. (KO8)
```

Finally, applying the already certified whole Fourier kernel only after the
exact coefficient pairing gives

```text
|Delta W_finite-t|
 < [1.1996732283063585006113383021259468974986624767154439756156448962067273848722271e-11 +/- 3.64e-91] < 1.21e-11. (KO9)
```

This closes the full profile obligation left open in Sections 11.397--11.399;
it does not infer the result from the small formal beta-minus-six coefficient.
Pi provenance: `beta^3=pi*C^2/8` comes from the exact transformed Kummer
geometry, `2pi` in (KO3) is the inverse Airy Fourier normalization, and the
Green Wronskian is `1/pi`.  No fitted constant is inserted.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_kummer_ode_profile_closure_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_kummer_ode_profile_closure_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_kummer_ode_profile_closure_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_kummer_ode_profile_closure_gate.py
```

No complete selected-leading plus source-current inequality, all-corridor or
height-uniform theorem, complete `Q_K-T` or `T_upper`, `Lambda<=0`,
PF-infinity, RH, or prize-level conclusion is proved.
