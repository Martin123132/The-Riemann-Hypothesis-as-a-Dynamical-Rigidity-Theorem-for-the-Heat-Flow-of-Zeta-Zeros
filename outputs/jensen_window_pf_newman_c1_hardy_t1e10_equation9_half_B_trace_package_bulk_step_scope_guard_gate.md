# B trace-package bulk-step scope guard

Date: 2026-08-13

Status: exact ownership audit; not a proof of the complete endpoint-B window
or RH; the negative trace-package bound survives while that promotion is
rejected

With

```text
tau_m=1_(622<=m<=39894),
ell_m=1_(m in {621,622}),
sigma_m=(1-tau_m)-1_(q_B,m<0),
```

the exact finite-mode identity is

```text
I_m+I_-m-tau_m P_bulk,m
 =[U_B,m+P_B,-m+ell_m sigma_m P_bulk,m]
  +[P_A,m+I_-m-P_B,-m
    +(1-ell_m)sigma_m P_bulk,m].                     (SG1)
```

The first bracket is the certified trace package.  The second bracket keeps
every nonlocal bulk step in the grouped remainder.  The window geometry

```text
c_low=[620.94069454105337461110481316249755282003887157159412506007027746541440167484559 +/- 6.05e-78],
c_high=[622.17131267942091816677759475066933602707465432095742996304611742329732171103145 +/- 2.92e-78]
```

places exactly modes 621 and 622 across a crossing.  Hence on the complete
window the omitted nonlocal step block is exactly

```text
sum_(m=1)^620 P_bulk,m - sum_(m=623)^39894 P_bulk,m. (SG2)
```

Modes 621 and 622 retain their variable steps in the local trace package;
modes at least 39895 have `sigma_m=0`.  Therefore the already-certified
bound

```text
E_Btr,win<-1.3198e-4
```

is preserved, but it is not a bound for the complete endpoint-B sector.
No estimate of (SG2) is asserted; it stays in `R_group` under the common
symmetric Abel prescription.

Proof boundary: exact finite-mode ownership, window sign roster, and a
fail-closed non-promotion guard only.  No independent bound for the omitted
bulk-step block, grouped remainder, complete endpoint B sector, complete
`Q_K-T` or `T_upper`, height-uniform theorem, `Lambda<=0`, PF-infinity, RH,
or prize-level conclusion is proved.
