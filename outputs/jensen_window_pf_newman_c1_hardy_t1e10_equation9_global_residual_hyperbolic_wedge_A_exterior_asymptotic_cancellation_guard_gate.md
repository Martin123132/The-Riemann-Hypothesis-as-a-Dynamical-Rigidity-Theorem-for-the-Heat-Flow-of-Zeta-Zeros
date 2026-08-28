# A-exterior asymptotic cancellation guard

Date: 2026-08-14

Status: exact endpoint asymptotics and nonconvergence guard; not an exterior
value, complete A carrier, or global residual bound

The affine/exact-minus-affine identity remains exact at every common finite
cutoff `epsilon>0`.  It cannot be passed term by term to the `x=0` exterior
limit.

Let `H_aff` denote the inner `P` integral of
`1+lambda_P P+mu_S S`.  The lower Fresnel expansion and the exact face limits
give

```text
exp(-iP_A^2/2)H_aff(P_A,S) -> i(mu_S-lambda_P),       (EC1)
sqrt(x)P_A -> -m sqrt(2pi),
sqrt(x)S   ->  m sqrt(2pi),
x^(3/2)S_x-> -m sqrt(pi/2).                           (EC2)
```

Across all 84 modes,

```text
[4.708493731800554843088672207461603754008622591596674582731829708006155e-6 +/- 2.25e-76]
 < Re(lambda_P-mu_S) <
[4.719527148625121397422073414727208945658338217579061687124384607573959e-6 +/- 1.34e-76].              (EC3)
```

Thus the affine exterior current has the common-phase asymptotic

```text
J_aff,m(x)=L_m x^(-3/2) exp(i Phi_A,face(x))+O(x^-1),
Phi_A,face(x)=pi A^2 x/4+(t/2)log((1-x)/x).           (EC4)
```

Summing every mode before the limit does not remove it:

```text
sum_m L_m real=[2.228172230992542958424843286017672358364553425176168256666241366651259e-10 +/- 1.30e-80],
sum_m L_m imag=[-8.912633834684513271822826835935198276103352536501490106728630443053637 +/- 2.04e-70],
|sum_m L_m| =[8.912633834684513271825612068439447132071902804953565775712293440159646 +/- 1.73e-70].                    (EC5)
```

Since `exp(i Phi_A,face(x))=x^(-it/2)(1+O(x))`, the cutoff
primitive has a nonzero `epsilon^(-1/2-it/2)` leading term.  Therefore the
standalone affine exterior has no improper limit.  The exterior
exact-minus-affine amplitude has the opposite leading term and has no
standalone limit either.

The full exact endpoint factor supplies the missing cancellation.  With
`K=pi A m`,

```text
E(P)=-2iK/(K-i) P^(-2)+O(P^-4),
exp(-iP^2/2) integral_(-infinity)^P E(q)exp(iq^2/2)dq
 =-2K/(K-i) P^(-3)+O(P^-5).                          (EC6)
```

Together with `x^(1/4)O(v(x))->sqrt(2)r^(-1/4)`, the full exact current is
only `O(x^-1/4)` and is absolutely integrable at `x=0`.  Hence the admissible
exterior object is

```text
C_A,ext^exact=lim_(epsilon->0)
 [C_A,ext^aff(epsilon)+R_A,ext^amp(epsilon)],         (EC7)
```

with the two terms retained under one cutoff, or equivalently the full exact
endpoint factor integrated directly.  The affine/amplitude split remains
valid on the compact local box.

Pi provenance: every `pi` in (EC1)--(EC7) comes from the exact equation-(9)
bi-Morse map, original triangle face, and endpoint normalization.  No fitted
constant is introduced.

Proof boundary: exact endpoint orders, all-mode coefficient separation, and
the common-cutoff nonconvergence/cancellation guard only.  No numerical exact
exterior value, local or exterior nonlinear-amplitude bound, complete A
endpoint theorem, `R_Dir` estimate, complete `Q_K-T` or `T_upper`, all-height
theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
