# Local exact-domain affine A carrier

Date: 2026-08-13

Status: rigorous finite-box exact-domain affine carrier on all 84 A modes;
not an exterior, nonlinear-amplitude, complete A, or global residual bound

For `P_A(S)` the exact A face, define the lower Fresnel tail

```text
T_-(P)=integral_(-infinity)^P exp(iu^2/2)du
 =sqrt(pi/2)exp(i*pi/4)
  [1+erf(exp(-i*pi/4)P/sqrt(2))].                    (LE1)
```

The affine amplitude `1+lambda_P P+mu_S S` can then be integrated exactly
in `P` over the exact domain.  On `y_half<=y<=0.0037` this gives

```text
C_A,loc^exact-aff=(2pi)^(-1) integral exp(-iS^2/2)
 {(1+mu_S S)T_-(P_A(S))-i lambda_P exp(iP_A(S)^2/2)}dS. (LE2)
```

This is the localized object required by Section 11.434.  It is equal to
local tangent plus local face, and contains no tangent exterior or artificial
tangent stationary point.

Restoring `epsilon_A=-1`, the common odd-endpoint carrier, paper rotation,
and equation-(9) normalization, deterministic summation before norms gives

```text
target-side local exact affine=[-0.099047632190178364240298873912985701433413404534447512723999842507 +/- 1.34e-12],
outer-side local exact affine =[-0.027817663785323280322128542936183719145832819877690684388449954561 +/- 1.32e-12],
complete local exact affine   =[-0.12686529597550164456242741684916942057924622441213819711244958558 +/- 2.66e-12],
sum of modewise moduli        =[0.28276880019685993831313200063576830575382989042431761278844878715 +/- 2.66e-12],
signed/modulus ratio          =[0.4486537973325900371387520496058165520025795558467507362 +/- 1.36e-11]. (LE3)
```

The production calculation is resumable through a modewise fsynced JSONL
cache.  An independent replay uses higher precision, more logistic-series
terms, and tighter quadrature tolerance.

Equation (LE3) is only the affine-amplitude contribution on the finite exact
domain.  It excludes the exact exterior `y>0.0037`, the exact-minus-affine
transformed amplitude even inside the box, the positive full-line projector
jump, and all other global sectors.  It is therefore not an `R_Dir` estimate.

Pi provenance: every `pi` in (LE1)--(LE3) comes from the exact Fresnel
primitive, equation-(9) bi-Morse phase, odd-endpoint carrier, and paper
normalization.  No fitted constant is introduced.

Proof boundary: exact affine-amplitude integral on the finite exact A domain
`y_half<=y<=0.0037` for modes 39853..39936 only.  No exact exterior affine
current, exact-minus-affine transformed-amplitude bound, complete A endpoint
theorem, `R_Dir` estimate, complete `Q_K-T` or `T_upper`, all-height theorem,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
