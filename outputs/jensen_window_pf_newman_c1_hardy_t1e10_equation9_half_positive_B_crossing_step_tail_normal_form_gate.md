# Exact B crossing: step mismatch plus sign-adapted tail

Date: 2026-08-13

Status: exact-lemma certificate; not a proof of the quantitative B remainder

For a positive mode put

```text
q_B=sqrt(x/2)(B-2m/x),
a_m=m sqrt(2/x),
T(r)=integral_r^infinity exp(i*pi*u^2/2)du.            (ST1)
```

Define one sign-adapted grouped tail

```text
U_B=E_B/(i*pi)-sgn(q_B)a_m T(|q_B|).                  (ST2)
```

Using oddness of the Fresnel primitive gives the exact identity

```text
P_B=U_B-1_(q_B<0)P_bulk.                              (ST3)
```

Consequently, over the complete positive source target `1..39894`,

```text
P_bulk+P_A+P_B-1_(m>=622)P_bulk
 =P_A+U_B+[1_(m<=621)-1_(q_B<0)]P_bulk.               (ST4)
```

Neither term after `P_A` is continuous separately.  At `q_B=0`, `U_B`
jumps by `-P_bulk`, while the step mismatch jumps by `+P_bulk`; their sum is
exactly continuous.  This is the nonlinear completion missing from the
isolated tangent profile.

Away from `q_B=0`, one Fresnel integration by parts gives the common leading
term

```text
U_B=E_B*B*x/[i*pi(B*x-2m)]+Fresnel remainder.         (ST5)
```

After the completed-square phase and external `1/x` factor, (ST5) becomes

```text
B exp(i*pi*B^2*x/4)/[i*pi(B*x-2m)],                   (ST6)
```

whose phase is independent of `m`.  Pairing the negative mode before the
same normal integration gives the complete first boundary current

```text
-2(2+i*pi*B^2*x)/[pi^2(B^2*x^2-4m^2)].               (ST7)
```

Thus its all-pair leading sum is cotangent-type.  For `a=B*x/2`,

```text
sum_(m=1)^infinity 1/(B^2*x^2-4m^2)
 =pi*cot(pi*a)/(4*B*x)-1/(2*B^2*x^2).                (ST8)
```

Equation (ST8) is used only after the local poles are removed.  The exact
face ledger proves that near the B trace saddle those poles are precisely
modes `621` and `622`; they must remain in the continuous combination
(ST2)--(ST4).  The cotangent remainder and the Fresnel remainder still need
explicit bounds.

Pi provenance: every `pi` comes from the equation-(9) Fresnel phase,
completed-square phase, or the integer cotangent partial fraction.  No fitted
constant is introduced.

Proof boundary: exact step/tail identity, crossing jump cancellation,
first boundary current, and cotangent algebra only.  No Fresnel-remainder or
cotangent-remainder bound, two-mode B enclosure, A-fold splice, complete
paired residual, `T_upper`, `Lambda<=0`, PF-infinity, RH, or prize-level
conclusion follows.
