# Frozen endpoint-current aggregate barrier

Date: 2026-08-13

Status: source-wide frozen-current obstruction certified; not a proof of the non-frozen Morse remainder

For every source target mode `622<=m<=39894`, let `R_m` be the exact
incomplete Fresnel plus endpoint-current ratio evaluated at the outer Morse
saddle.  The frozen endpoint-retained model is

```text
M_frozen=2 Re sum_m exp(i(theta_0-t log m)) R_m/sqrt(m),
T_upper =2 Re sum_m exp(i(theta-t log m))/sqrt(m).
```

Here `theta` is the exact Riemann--Siegel theta from `log Gamma`, while
`theta_0=t/2[log(t/(2*pi))-1]-pi/8`.  Rigorous Arb summation over all
39273 modes gives

```text
M_frozen-T_upper=[-0.0962793319571363535562609595823998120610194967715247820196462852429874782908316754510193840 +/- 9.43e-77],
sum_m (frozen-exact)=[-0.0481396659785681767781304797911999060305097483857623910098231426214937391454158377255096920 +/- 4.72e-77]
                     + i [0.0422998980498563418577271981015563070092528024562295571237143389659141371111654353429100672 +/- 4.72e-77],
equation-(9) physical scale=[0.0288477589447288904293890618231644738617507133657955976409670733658970066622609447255567003 +/- 3.00e-77].
```

The carrier-phase difference is only

```text
theta-theta_0=[2.08333333333333333333454861111111111111111495535714285714285717238653554994821314144622617e-12 +/- 9.12e-80],
two-real phase-only aggregate=[-2.92543149278482427731234983414212707093443047612813103102023794654467405752466458784018981e-12 +/- 7.23e-77].
```

It is therefore far too small to explain the discrepancy.  Signed
aggregation does reduce the termwise triangle
`[0.618694996982026575895484815590940902788552532360688008292681491531466779116611867880276221 +/- 6.01e-77]`, but does not cancel the frozen
defect.

This rejects one specific approximation: replacing the exact global Morse
amplitude by its saddle value, even after retaining the incomplete endpoint
current there.  It does **not** reject the endpoint-retained Morse method and
does not lower-bound the exact equation-(9) error.  The omitted integral

```text
integral exp(-i*t*s^2/4) [A_m(s)-A_m(0)] ds
```

can be the same size near the fold and may supply the missing correction.  It
must be derived and bounded before any comparison with `T_upper` is valid.

Pi provenance: every `pi` here comes from the equation-(9) Kummer/Fourier
phase, the Gaussian Morse normalization, or the Riemann--Siegel theta
normalization.  No geometric fit or free polygonal constant is inserted.

Proof boundary: a rigorous finite saved-height sum for the frozen leading
model only.  No bound on the non-frozen Morse-amplitude integral, complete
`T_upper` theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion is
proved.
