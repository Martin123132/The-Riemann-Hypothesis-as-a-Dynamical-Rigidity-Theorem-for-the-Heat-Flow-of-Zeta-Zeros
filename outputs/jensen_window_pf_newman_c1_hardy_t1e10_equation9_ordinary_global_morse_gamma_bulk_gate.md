# Exact global-Morse Gamma bulk

Date: 2026-08-13

Status: ordinary full-line bulk closed at the saved height; not a proof of the endpoint charts

After the exact logistic substitution, replace only the finite alpha current
by its full-line bulk part

```text
P_bulk=m sqrt(2/x)(1+i).
```

The complete Kummer weight and Jacobian then simplify exactly:

```text
A_bulk(s)=m sqrt(2)(1+i) r^(3/4)
          v(s)^(-1/4) dv/ds.                           (GB1)
```

No stationary-phase truncation is needed.  With Abel damping first and the
limit taken afterward,

```text
I_bulk(t)=exp(i t/2) Gamma(3/4+i t/2)
          /(i t/2)^(3/4+i t/2).                       (GB2)
```

Divide by the leading Gaussian
`I_G=2 sqrt(pi/t) exp(-i*pi/4)` and write `B=I_bulk/I_G`.
At `t=10^10`, rigorous Arb evaluation gives

```text
log B=[-3.125000000000000000019531250000000000000397135416666666666683573404947839478864078268331316193681969e-22 +/- 1.59e-91]
      + i [2.083333333333333333334548611111111111111114955357142857142857172386532738095238872578386377522340835e-12 +/- 3.12e-90],
theta-theta_0=[2.083333333333333333334548611111111111111114955357142857142857172386532738095238258358396366361908533e-12 +/- 5.31e-90],
|B-exp(i(theta-theta_0))|
 =[3.125000000000000000019042968750000000000391082763671875000016764509677809775087377444165302570213851e-22 +/- 2.25e-91].        (GB3)
```

Thus the imaginary part of `log B` supplies the exact saved-height
Riemann--Siegel phase correction at the working precision; only a relative
amplitude defect of about `3.125e-22` remains.  Summing it over every source
target mode `622..39894` before taking the real projection gives

```text
complex aggregate <= [1.092538584153040444695866547991813837208099771680211742572158912289839289161796654670795293041014340e-19 +/- 7.86e-89],
two-real aggregate <= [2.185077168306080889391733095983627674416199543360423485144317824579678578323593309341590586082028679e-19 +/- 1.58e-88],
equation-(9) physical <= [4.918151566720828980806353638226121738513797734778668279706287331790919038545165004173875668351233940e-20 +/- 3.54e-89]. (GB4)
```

This is more than `[173906901555602.3484378057764843106269836425781250000000000000000000000000000000000000000000000000000 +/- 0.0311]` times below
the normalized target.  The ordinary bulk is therefore not the quantitative
wall at the saved height.

The result does not discard finite-endpoint terms.  Their lower
characteristic current and upper nonstationary current are separate exact
pieces and must be reassembled with the already certified fold and endpoint
charts.  In particular, (GB4) does not turn the frozen endpoint-current
model into a valid approximation.

Pi provenance: `pi` comes from the equation-(9) Kummer/Fourier phase, the
Gaussian Gamma integral, and the Riemann--Siegel theta normalization.  No
geometric fit is introduced.

Proof boundary: exact Abel-regularized full-line bulk at `t=10^10` and its
finite source-mode aggregate only.  No finite-endpoint reassembly, ordinary
endpoint-tail bound, complete `T_upper` theorem, height-uniform theorem, or
`Lambda<=0` theorem is proved.  No claim of PF-infinity, RH, or a prize-level
conclusion is made.
