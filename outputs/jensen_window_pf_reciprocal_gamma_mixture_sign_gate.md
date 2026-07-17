# Jensen-Window PF Reciprocal-Gamma Mixture Sign Gate

Date: 2026-07-12

Status: exact reciprocal-gamma sign regularity with positive-mixture
obstruction. This is not a proof of RH or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_reciprocal_gamma_mixture_sign_gate`.

```text
work/rh_compute/results/jensen_window_pf_reciprocal_gamma_mixture_sign_gate.json
python work/rh_compute/scripts/jensen_window_pf_reciprocal_gamma_mixture_sign_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_gamma_mixture_sign_gate.py
```

Current result:

```text
validated Jensen-window PF reciprocal-gamma mixture sign gate: 10 rows, 0 issues, 3 exact kernel identities, 1 published all-order theorem, 1 fixed-scale theorem, 2 exact mixture countermodels, 1 tilted-variance equivalence, 1 Xi order-two composition, 1 higher-order handoff
```

## Fixed-Scale Theorem

The coefficient normalization is

```text
m_k(lambda)=mu_(2k)(lambda)=integral_0^infinity t^k*rho_lambda(dt)
gamma_k=k!/(2k)!=sqrt(pi)/(4^k*Gamma(k+1/2))
A_k(lambda)=gamma_k*m_k(lambda)
```

Karlin's Lemma 9.2 states:

```text
For alpha>-nu, every (nu+1)-minor of [1/Gamma(alpha+i+j)]_(i,j>=0) has strict sign (-1)^(nu*(nu+1)/2).
For every t>0 and n>=0, K_(n,t)(i,j)=gamma_(n+i+j)*t^(n+i+j) is strictly sign-regular of all orders with signature epsilon_k=(-1)^(k*(k-1)/2).
```

So the fixed-scale reciprocal-gamma kernel already has the complete
signature observed in the Arb staircase, for arbitrary row and column
sets. The contiguous benchmark is explicit:

```text
det[gamma_(i+j)]_(i,j=0..k-1)=(-1)^(k*(k-1)/2)*pi^(k/2)*4^(-k*(k-1))*prod_(r=1)^(k-1)r!/prod_(i=0)^(k-1)Gamma(k-1/2+i)
```

## Mixture Integral

Multilinearity in the determinant rows gives

```text
R_(k,n)(j_1,...,j_k)=integral_((0,infinity)^k) det[gamma_(n+i+j_l)*t_i^(n+i+j_l)]_(i=0..k-1,l=1..k) prod_(i=0)^(k-1)rho_lambda(dt_i)
```

The scale variables are independent. Karlin's theorem signs the diagonal
`t_0=...=t_(k-1)`, not the full product-measure integrand.

## Exact Nonclosure

At order two the symmetrized integrand is

```text
D_sym(t_0,t_1)=(t_0^2-6*t_0*t_1+t_1^2)/24
D_sym(1,1)=-1/6
D_sym(1,10)=41/24
D_sym<=0 iff 3-2*sqrt(2)<=t_1/t_0<=3+2*sqrt(2)
```

A positive measure supported away from zero gives a complete countermodel:

```text
rho=10*delta_(1/100)+delta_1+delta_2
m_0,m_1,m_2=['12', '31/10', '5001/1000']
A_0,A_1,A_2=['12', '31/20', '1667/4000']
A_0*A_2-A_1^2=5197/2000>0
```

Thus the cone of fixed-scale sign-regular matrices is not closed under
the positive scale mixing used by the Xi moment representation.

## Xi Concentration

For the `n`-tilted moment law,

```text
P_n(dt)=t^n*rho_lambda(dt)/m_n(lambda)
CV_n^2=Var_(P_n)(T)/E_(P_n)(T)^2=m_n*m_(n+2)/m_(n+1)^2-1
A_(n+1)^2>=A_n*A_(n+2) iff m_(n+1)^2/(m_n*m_(n+2))>=(2*n+1)/(2*n+3) iff CV_n^2<=2/(2*n+1)
```

The lambda=-100 full ratio-cone entry plus infinite heat-flow invariance proves this tilted relative-variance bound for every n at lambda=0.
This gives a structural interpretation of the completed all-shift
order-two cone: the tilted Xi measure becomes relatively more concentrated
at the exact rate `2/(2n+1)`.

## Live Handoff

Control the sign-changing row-wise compound integrand after Xi mixing for every k and column set, using a hierarchy of tilted concentration, a coupling that keeps row scales near the diagonal, or a new positive compound-kernel factorization.

The subsequent order-three gate proves that the first contiguous compound sign is exactly C_n=q_(n+1)*q_(n+3)-q_(n+2)^2+1>0.
A strict rational sequence satisfies the full ratio, adaptive-defect, and cubic Jensen cones while having C_n<0.
See `outputs/jensen_window_pf_reciprocal_defect_compound_order3_gate.md`.

The next signed-Hankel antecedent is therefore a higher compound
concentration theorem for the actual Xi mixing measure, not another
application of fixed-scale reciprocal-gamma sign regularity.

Primary source: https://doi.org/10.1090/S0002-9947-1964-0168010-2.
