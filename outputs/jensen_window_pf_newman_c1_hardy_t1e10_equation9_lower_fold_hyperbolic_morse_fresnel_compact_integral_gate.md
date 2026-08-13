# Signed hyperbolic Morse--Fresnel compact integral

Date: 2026-08-12

Status: cancellation-preserving compact integral certified at three heights;
this is not a proof of the outer tails or complete one-mode join

For mode `39894`, retain the common carrier and paired endpoint/Fresnel
current.  On `|u|<=200` and `0<=Y<=64/sqrt(2 beta)`, compare

```text
Delta I_core(t)=int int [exp(i Delta_c) A_t(u,Y) exp(i Phi_t(u,Y))
                         -exp(i Phi_0(u,Y))] dY du.       (CI1)
```

The `Y` phase is exactly quadratic and `A_t` is exactly affine in `Y`.
Writing

```text
H_0(a,b;Y)=int_0^Y exp(i[a y^2+b y])dy,
H_1(a,b;Y)=-b H_0/(2a)
            -i[exp(i[aY^2+bY])-1]/(2a),                (CI2)
```

evaluates the inner integral in closed Fresnel form.  The remaining outer
integral is performed directly on the signed difference in `delta=v-1`,
with `du/delta=sqrt(t/2)/[(1+delta)sqrt(2k(delta))]` and

```text
k(delta)=[delta-log(1+delta)]/delta^2.                 (CI3)
```

A degree-24 convergent series with an explicit complex geometric
tail resolves the removable value at `delta=0`.  Append-only interval panels
give

```text
lower_face: |Delta I_core|=[0.000481220264885509953424822963362879081205959376119303285993048782329117862843759212410077453 +/- 1.86e-24], physical=[0.000216625228070934553727621249466418872947310628949712015857590006063376764089412063185591251 +/- 8.36e-25]
event: |Delta I_core|=[0.000481223680878171497985630133101200619883315399246971025963863412971965871633983624633401632 +/- 1.86e-24], physical=[0.000216626765807899084481139466324340252180859900815236934424458739318480882118933550373185426 +/- 8.36e-25]
upper_face: |Delta I_core|=[0.000481227065774524764710297185943367833380477806087325838818835295820908193320519785629585385 +/- 1.86e-24], physical=[0.000216628289546606757841079236725987328246950407668809309633445770067183833162971495767123997 +/- 8.36e-25]
```

uniformly at `t=tau-pi/16`, `tau`, and `tau+pi/16`.  In particular,

```text
|Delta I_core| < 0.00049,
(sqrt(2)/pi)|Delta I_core| < 0.000221.                 (CI4)
```

The displayed `pi` comes from the original Kummer quadratic phase and
Fourier--Poisson character; `sqrt(2)/pi` is the exact coordinate/Jacobian
factor already derived in the overlap chart.

## Boundary

This gate certifies only the signed compact rectangle for one prototype mode
at three heights.  It does not bound `|u|>200`, complete the one-mode
Airy/logistic integral join, propagate through 399 events, establish complete
`T_upper`, or prove `Lambda<=0`, RH, or a prize-level conclusion.
