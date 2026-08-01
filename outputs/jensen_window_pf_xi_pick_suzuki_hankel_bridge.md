# Jensen-Window PF Xi Pick/Suzuki Hankel Bridge

Date: 2026-07-23

Status: exact Phi-to-xi directional Pick reduction, published Suzuki
arithmetic-Hankel equivalence, and exact theorem-mismatch guards.
This is not a proof of LP+, PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_xi_pick_suzuki_hankel_bridge`.

```text
work/rh_compute/results/jensen_window_pf_xi_pick_suzuki_hankel_bridge.json
python work/rh_compute/scripts/jensen_window_pf_xi_pick_suzuki_hankel_bridge.py
python work/rh_compute/scripts/check_jensen_window_pf_xi_pick_suzuki_hankel_bridge.py
```

## Exact Xi Coordinate

The workspace normalization gives

```text
H_0(x)=xi((1+i*x)/2)/8; M(w)=integral_R Phi(u)cosh(u*w)du=2*H_0(i*w)=xi((1+w)/2)/4; F(z)=M(sqrt(z))
```

Put `s=sigma+i*tau`, `delta=sigma-1/2>0`, and
`w=2*(delta+i*tau)`. Away from zeros, the Pick numerator becomes

```text
For s=sigma+i*tau, delta=sigma-1/2>0, w=2*(delta+i*tau), and xi(s)!=0: 2*|w|^2*P_Phi(w^2)/|M(w)|^2=tau*Re(xi'(s)/xi(s))-delta*Im(xi'(s)/xi(s))
```

For `X=delta^2-tau^2` and `Y=2*delta*tau`, the same quantity is

```text
For X=delta^2-tau^2 and Y=2*delta*tau: partial_Y log|xi(1/2+delta+i*tau)| at fixed X=[tau*Re(xi'/xi)-delta*Im(xi'/xi)]/[2*(delta^2+tau^2)]
```

Thus the Pick endpoint is strict increase of `|xi|` along every
first-quadrant branch `X=constant` as `Y` increases. This differs
from ordinary horizontal increase in `sigma`.

## Horizontal-Growth Guard

The Sondow-Dumitrescu theorem gives

```text
Sondow-Dumitrescu: Re(xi'(s)/xi(s))>0 in every open right half-plane lying strictly to the right of all xi zeros; in particular sigma>1
```

but the Pick direction also contains `-delta*Im(xi'/xi)`. The
distinction is forced by an exact polynomial.

For a conjugate zero pair,

```text
If F(z)=(z-(A+i*B))*(z-(A-i*B)) and z=x+i*y, y>0, then -Im(F'(z)/F(z))=2*y*((x-A)^2+y^2-B^2)/([((x-A)^2+(y-B)^2)]*[((x-A)^2+(y+B)^2)])
```

Take

```text
F_*(z)=z^2+6*z+25
M_*(w)=F_*(w^2)=w^4+6*w^2+25
zeros(M_*)={+/-1+/-2*i}
w=11/10+(7/5)*i
z=w^2=-3/4+(77/25)*i
(x+3)^2+y^2-16=-14511/10000<0
partial_Re(w) log|M_*(w)|=9246851240/10870189707>0
-Im(F_*'(z)/F_*(z))=-297959200/10870189707<0
```

Every zero of `M_*` has real part at most one, so the horizontal
logarithmic derivative is positive throughout `Re(w)>1`. The Pick
sign still fails at the displayed point. Horizontal xi-modulus
monotonicity is therefore theorem-mismatched to the required
hyperbolic direction.

## Suzuki Arithmetic Hankel Route

For the completed zeta function define

```text
E_(omega,nu)(r)=xi(1/2+omega-i*r)^nu and Theta_(omega,nu)(r)=E#(r)/E(r)=[xi(1/2-omega-i*r)/xi(1/2+omega-i*r)]^nu
K_(omega,nu)=Fourier_inverse(Theta_(omega,nu)); (K_(omega,nu)[t]f)(x)=1_(x<t)*integral_(-infinity)^t K_(omega,nu)(x+y)f(y)dy
```

For zeta, whose Selberg-class degree is one, Suzuki proves
unconditionally that `nu*omega>1` gives an initial interval on
which

```text
For nu*omega>1, Suzuki constructs tau>0 with det(I+/-K[t])!=0 on 0<=t<tau and H(t)=diag(1/gamma(t),gamma(t)), gamma(t)=[det(I+K[t])/det(I-K[t])]^2
```

The global theorem is the exact equivalence

```text
RH iff there exist omega_n decreasing to 0 and integers nu_n with nu_n*omega_n>1 such that det(I+/-K_(omega_n,nu_n)[t])!=0 for every t>=0 and J_(omega_n,nu_n)(t;r,r)->0 as t->infinity for every r in the upper half-plane
```

Suzuki states two explicit obligations for an arithmetic operator
already constructed without RH. The later causal-multiplier audit
sharpens their cofinal logic:

```text
A strictly decreasing cofinal omega_n->0 family with det(I+/-K_(omega_n,nu_n)[t])!=0 for every n and t already implies RH: the compatible finite forms extend to a bounded causal multiplier, and any off-line zero would generate distinct shifted zeros accumulating at itself
```

The real-boundary identity `|Theta(u)|=1` does not supply a
shortcut to the spectral gate:

```text
Theta_*(r)=(r+i)/(r-i) has |Theta_*(u)|=1 for every real u but has a pole at r=i; moving a high-contour inverse Fourier integral to the real axis acquires a nonzero residue proportional to exp(x). Boundary unimodularity alone does not make the high-contour Hankel operator unitary or contractive
```

The Suzuki kernel is defined from a high horizontal contour. Moving
that contour to the real line without residue already requires the
relevant pole-free analyticity, so a boundary-unitary argument would
be circular.

## Signed-Hankel Contact

The classical Fredholm series is

```text
det(I+/-K[t])=sum_(m>=0)(+/-1)^m/m!*integral_((-infinity,t)^m) det(K(x_i+x_j))_(i,j=1)^m dx_1...dx_m
```

so each coefficient is an integral of a continuum Hankel determinant.
This is a genuine contact with the signed-Hankel programme. It does
not make total nonnegativity sufficient:

```text
The constant Hankel kernel K(x+y)=1 on L2(0,1) is totally nonnegative and rank one, but its nonzero eigenvalue is 1, so det(I-K)=0 and det(I+K)=2
```

The surviving first gate is therefore quantitative:

```text
The finite-t spectral target is pointwise-in-t strict contractivity ||K_(omega_n,nu_n)[t]||<1 for every finite t on a cofinal omega_n->0 sequence; no positive gap uniform in t is required or possible in the desired HB case
```

Here `for every finite t` is pointwise in `t`; it does not mean
one epsilon-sized gap valid uniformly as `t` tends to infinity.
The exact truncation-path distinction is developed in
`outputs/jensen_window_pf_suzuki_spectral_frontier.md`.

For the terminal premise, the exact distinction is

```text
For one fixed pair the terminal limit is not a direct consequence of contractivity. Within the cofinal RH criterion, however, determinant nonvanishing first implies RH, after which Suzuki's Theorem 2.3 supplies J(t;r,r)->0
```

The full proof is in
`outputs/jensen_window_pf_suzuki_determinant_only_reduction.md`.
Finite `t` grids, finite-rank discretizations, local Hamiltonian
positivity, or a canonical system built only after assuming the
Hermite-Biehler property still cannot prove the all-time gate.

## Sources

- Masatoshi Suzuki, `Hamiltonians arising from L-functions in the Selberg class`, Theorems 2.1-2.4: https://doi.org/10.1016/j.jfa.2021.109116
- Jonathan Sondow and Cristian Dumitrescu, `A monotonicity property of Riemann's xi function and a reformulation of the Riemann Hypothesis`: https://arxiv.org/abs/1005.1104
- `outputs/jensen_window_pf_phi_pick_kernel_target.md`
- `outputs/formal_core.md`
- `outputs/jensen_window_pf_suzuki_determinant_only_reduction.md`

## Proof Boundary

This artifact proves the exact Phi-to-xi normalization, the directional and hyperbolic forms of the Pick target, and three exact countermodel guards. It records Suzuki's published unconditional local construction and exact global equivalence. It also records the internally audited determinant-only cofinal reduction. It does not prove global Fredholm nonvanishing, the Phi Pick sign, LP+, PF-infinity, Jensen hyperbolicity for zeta, RH, or Lambda<=0.
