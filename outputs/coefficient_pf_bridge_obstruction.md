# Coefficient PF Bridge Obstruction

Date: 2026-07-03

Status: theorem-search note. This records what would be needed to turn the certified finite Toeplitz/PF evidence into a proof. It is not a proof of PF-infinity, RH, or `Lambda <= 0`.

## Normalized Object

The corrected coefficient sequence is:

```text
A_k(lambda) = mu_{2k}(lambda) k! / (2k)!
c_k(lambda) = A_k(lambda) / k! = mu_{2k}(lambda) / (2k)!
```

Thus:

```text
F_lambda(s) = sum_{k >= 0} c_k(lambda) s^k
            = integral_R cosh(u sqrt(s)) exp(lambda u^2) Phi(u) du
```

and, after `s = -z^2`, this is the usual even Fourier/cosine transform form of the de Bruijn-Newman heat-flow object.

The ASW/Edrei route says:

```text
c_k(0) is PF-infinity
  => F_0 is restricted Laguerre-Polya
  => zeros of F_0 lie on the required real/negative axis in s
  => RH in the corresponding z variable.
```

This is the right normalization. It is also almost the whole problem.

## Exact Relation To Jensen Windows

The ordinary generating function of `c` is the exponential generating
function of `A`:

```text
F(z)=sum_k c_k z^k=sum_k A_k z^k/k!.
```

For the diagonal tail operator `T_n[z^j]=A_(n+j)z^j`, one has

```text
T_n[(1+z)^d]
 =sum_(j=0)^d binom(d,j)A_(n+j)z^j
 =P_(d,n)(z).
```

Consequently ASW/Edrei, Polya-Schur, closure of type-I functions under
differentiation, and finite ASW give the exact all-order equivalence

```text
c is PF-infinity
  <=> F is Laguerre-Polya type I
  <=> every shifted Jensen window is PF-infinity.
```

This does not make the route easier: proving all-order PF-infinity of `c` is
already an endpoint theorem. It does correct the programme map. The
coefficient-PF and Jensen-window routes are distinct finite diagnostics but
not distinct all-order targets. The executable gate is
`outputs/jensen_window_pf_coefficient_pf_equivalence_gate.md`.

## Exact Relation To Edrei-Log Hankel Positivity

Normalize `H(z)=F(z)/F(0)` and write

```text
log H(z)=sum_(n>=1)ell_n z^n,
q_n=n ell_n,
p_n=(-1)^(n-1)q_n,
a_r=p_(r+1)=(-1)^r[z^r]H'(z)/H(z).
```

Sokal's logarithmic-derivative criterion gives, for this known entire `H`,

```text
H is Laguerre-Polya type I
  <=> (a_r)_(r>=0) is a Stieltjes moment sequence.
```

Thus the nonlinear Edrei-log Hankel route is another exact language for the
same endpoint. A strict sufficient target is positivity at every order of

```text
det(p_(i+j+1))_(i,j=0)^m,
det(p_(i+j+2))_(i,j=0)^m.
```

This is not ordinary moment positivity of `mu_(2k)` and not the original
signed-Hankel hierarchy of `A_k`. The existing finite Edrei-log determinants
remain diagnostics only. The exact gate is
`outputs/jensen_window_pf_edrei_stieltjes_equivalence_gate.md`.

Krein's Stieltjes/Pick characterization further gives the exact Phi target

```text
P_Phi(z):=-Im(F_0'(z)conj(F_0(z)))>0
for every Im(z)>0.
```

The polarized double-integral formula and positive-resolvent alternative are
in `outputs/jensen_window_pf_phi_pick_kernel_target.md`. The exact radial heat
hierarchy does not provide generic backward invariance; the repeated-zero
countermodel is in
`outputs/jensen_window_pf_edrei_heat_flow_boundary_gate.md`.

The positive-resolvent handoff is no longer merely schematic. Suzuki's
arithmetic-Hankel theorem gives an exact endpoint route through:

```text
det(I+/-K_(omega_n,nu_n)[t])!=0 for every finite t,
J_(omega_n,nu_n)(t;r,r)->0 as t->infinity,
```

on a cofinal sequence `omega_n->0`, `nu_n*omega_n>1`. Its Fredholm series is
an integral expansion in continuum Hankel determinants, so it is genuine
contact with the signed-Hankel programme. It is not yet a bridge: a TN
rank-one kernel can have `det(I-K)=0`, boundary unimodularity cannot justify
crossing possible poles, and Suzuki's local interval does not imply either
global gate. See
`outputs/jensen_window_pf_xi_pick_suzuki_hankel_bridge.md`.

## Why The Finite Evidence Is Not A Bridge

The current finite certificates prove many concrete Toeplitz minor inequalities for `c_k(lambda)`. They do not prove:

```text
all Toeplitz minors of all orders are nonnegative.
```

By ASW/Edrei, proving that infinite statement would imply the required Laguerre-Polya conclusion. Therefore a proof cannot simply say:

```text
many finite minors are positive, therefore c_k is PF-infinity.
```

It must add a new structural reason why every remaining minor is nonnegative.

## Candidate Structural Bridges

### 1. Moment Positivity Is Not Enough

Writing `x = u^2` gives a positive-looking moment shape:

```text
mu_{2k}(lambda) = integral_0^infty x^k rho_lambda(x) dx
c_k(lambda) = mu_{2k}(lambda) / (2k)!
```

If `rho_lambda` is nonnegative, then `mu_{2k}` is a Stieltjes moment sequence and has Hankel positivity. But Toeplitz PF-infinity of:

```text
mu_{2k} / (2k)!
```

does not follow from ordinary moment positivity. Hankel positivity and Toeplitz total positivity are different determinantal regimes.

Executable guard:

```text
python work/rh_compute/scripts/countermodel_gate_examples.py
```

The exact Stieltjes multiplier trap uses the positive measure:

```text
10 delta_0 + delta_1 + delta_2
```

with moments:

```text
mu_0..mu_6 = 12, 3, 5, 9, 17, 33, 65
```

and leading Hankel determinants:

```text
12, 51, 40, 0
```

But for:

```text
c_k = mu_k / (2k)!
```

one gets:

```text
c_0 = 12
c_1 = 3/2
c_2 = 5/24
c_1^2 - c_0 c_2 = -1/4.
```

So even a concrete positive Stieltjes moment source can fail the first
nontrivial coefficient Toeplitz/PF inequality after the `(2k)!`
normalization.

### 2. Factorial Multiplier Route

The sequence:

```text
d_k = 1 / (2k)!
```

has generating function:

```text
sum d_k s^k = cosh(sqrt(s)),
```

whose zeros are real and negative. This makes `d_k` itself compatible with the restricted Laguerre-Polya/PF sequence world.

But this does not prove that:

```text
c_k = mu_{2k} d_k
```

is PF-infinity, because `mu_{2k}` is not known to be a PF-infinity sequence. A multiplier theorem would need hypotheses that the moment sequence actually satisfies.

The exact trap above shows that those hypotheses must be stronger than
ordinary Stieltjes/Hankel positivity.

### 3. Toeplitz Minor Integral Formula

The cleanest missing bridge would be an identity of the form:

```text
det(c_{q_b - p_a})_{a,b=1}^m
  = integral of a manifestly nonnegative determinant/product
```

for every row set `p_a` and column set `q_b`.

This would convert PF-infinity into a Cauchy-Binet or Andreief-type positivity theorem. At present, no such identity is known for `c_k = mu_{2k}/(2k)!`.

### 4. Schur/Jacobi-Trudi Reformulation

Toeplitz minors of a one-sided coefficient sequence are Jacobi-Trudi determinants:

```text
det(c_{lambda_i - i + j})
```

after reindexing. With the normalization:

```text
d_k = c_k / c_0
```

the exact row/column map is:

```text
K = max(0, max_i(r_i - i), max_i(q_i - i))
lambda_i = K + i - r_i
mu_i     = K + i - q_i
det(c_{q_j-r_i}) = c_0^m s_{lambda/mu}(d).
```

The structural nonzero condition `q_i >= r_i` is exactly the condition that `lambda/mu` is a valid skew shape. Therefore the coefficient PF problem can be reformulated as positivity of all skew-Schur evaluations under the specialization:

```text
h_k -> d_k(lambda).
```

This is a promising language because total positivity of Toeplitz matrices is equivalent to positivity of all such Schur-type determinants. The hard missing step is to prove this specialization is Schur-positive from the zeta/heat-kernel structure.

The exact reindexing workbench and note are:

```text
work/rh_compute/scripts/toeplitz_jacobi_trudi_map.py
outputs/toeplitz_jacobi_trudi_bridge_note.md
```

## Direct-Kernel Warning

The direct de Bruijn-Newman kernel positivity route is separate from the coefficient route. A 2026 arXiv preprint claims a certified PF5 failure for the direct kernel `K(u) = Phi(|u|)`.

This does not contradict the clean finite evidence for `c_k = A_k/k!`, because:

```text
translation minors of K(u) != coefficient Toeplitz minors of c_k.
```

But it strongly warns against any manuscript path that tries to prove the original kernel is PF-infinity.

## Current Verdict

The coefficient PF route is alive but not yet a proof strategy. It becomes a proof strategy only if we can supply one of:

```text
1. an all-order determinant positivity identity;
2. a known total-positivity theorem that applies to the heat-kernel coefficient specialization;
3. a stable polynomial approximation with real negative zeros and coefficient convergence;
4. a noncircular heat-flow preservation theorem for coefficient PF-infinity;
5. Suzuki's cofinal all-t Fredholm spectral gate; the determinant-only
   causal-multiplier reduction makes the published terminal premise redundant.
```

Until then, finite Toeplitz/PF certificates are falsification pressure and route-selection evidence, not a bridge theorem.

The explicit theorem targets are recorded in:

```text
outputs/missing_coefficient_pf_theorem.md
```

The route triage and failure gates are recorded in:

```text
outputs/coefficient_pf_decision_tree.md
```

## Sources

- Aissen, Schoenberg, Whitney, "On the generating functions of totally positive sequences I": https://doi.org/10.1007/BF02786970
- Edrei, "On the generating functions of totally positive sequences II": https://doi.org/10.1007/BF02786971
- Groechenig, "Schoenberg's Theory of Totally Positive Functions and the Riemann Zeta Function": https://arxiv.org/abs/2007.12889
- Michałowski, "On the Polya Frequency Order of the de Bruijn-Newman Kernel": https://arxiv.org/abs/2602.20313
