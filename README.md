 Riemann Zeta Zero Dynamics under Heat Flow  
 A Dynamical Rigidity Perspective on the Riemann Hypothesis

This repository contains research notes, draft manuscripts, and a large
collection of standalone Python experiments investigating the dynamics of
Riemann zeta zeros under Gaussian heat-flow smoothing.

The central theme is the emergence of a **low-dimensional, rigid
renormalisation-group (RG) structure** governing zero motion under *forward*
heat flow, and the implications of this structure for the Riemann Hypothesis.

> **Scope note**  
> This repository does **not** claim a proof of the Riemann Hypothesis.
> It presents exact analytic identities together with extensive
> high-precision numerical evidence for a *dynamical rigidity mechanism*
> that appears to obstruct complex zero formation under forward heat flow.

---

## 1. Core Idea 

- Apply **Gaussian smoothing (heat flow)** to the Riemann–Siegel Z-function.
- Track the induced motion of its real zeros as a function of the heat
  parameter σ.
- Observe that zero motion:
  - is **exactly σ-additive**,
  - admits a well-defined generator `dt/dσ`,
  - is governed (numerically) by a **rank-2 invariant manifold**.
- This manifold is spanned by:
  1. an antisymmetric **Hadamard-type repulsive interaction**, and
  2. a smooth **density-balancing mode**.
- Perturbations transverse to this manifold produce **dramatic RG breakdown**
  and exhibit a **non-analytic cusp** in the associated energy.
- Creation of complex zeros would require a **rank increase** in the
  generator, which is empirically obstructed under forward heat flow.

This reframes the Riemann Hypothesis as a statement of **dynamical rigidity
and stability**, rather than a static property of zero locations.

---

## 2. Manuscripts / Formal Write-Ups

### Primary manuscript
- **`The Riemann Hypothesis as a Dynamical Rigidity Theorem for the Heat Flow of Zeta Zeros`**
  - Full formal exposition.
  - Exact zero-tracking lemmas.
  - Hadamard log-derivative structure.
  - Rank-2 RG closure.
  - Transverse instability and cusp behaviour.
  - Relation to the de Bruijn–Newman constant.

### Consolidated empirical summary
- **`part 1`**
  - Narrative summary of all numerical findings.
  - Operator dependence, σ-additivity, kernel dominance, and open questions.

---

## 3. Formal Analytic Backbone (What Is Exact)

The following components are *analytic* and exact:

- Zero tracking via the implicit function theorem.
- Exact generator formula under heat flow:
  

dt_i/dσ = −c · (F_tt / F_t)(t_i)

- Hadamard log-derivative expansion yielding a nonlocal repulsive interaction:

(F_tt / F_t)(t_i) = 2 Σ_{j≠i} 1/(t_i − t_j) + A(σ)

- Repulsive barrier preventing finite-time zero collisions.
- Collision necessity for complex zero creation.
- Square-root bifurcation at the de Bruijn–Newman boundary (Weierstrass normal form).

These results are summarised and unified in the **Formal Core** section of the
main manuscript.

---

## 4. Numerical Discoveries (What Is Empirical)

All RG structure claims are supported by extensive numerical experiments:

- Exact σ-additivity of zero displacements.
- Dominance of a **two-dimensional velocity subspace**.
- Rank-2 RG equation:

d a / d log T = B a

with numerical rank ≈ 2 and relative closure error ~10⁻¹⁴.
- Violent amplification of RG error (>10¹²) under transverse perturbations.
- A non-analytic cusp in transverse “energy”:

E(δ) ~ |δ|^α,   α ≈ 0.6

These findings are **numerical evidence**, not proofs, and are presented as
such.

---

## 5. Code Structure and Experimental Map

Each `part N` file is a **self-contained Python experiment** with inline
commentary and, in many cases, embedded sample output.

### Orientation guide

| File | Purpose |
|----|----|
| `part 1` | Consolidated narrative of all experimental findings |
| `part 2` | Operator tests, convolution drift, scaling laws |
| `part 3` | Global balance tests, phase-only Fourier masks |
| `part 5` | Continuous σ-flow tracking and instantaneous velocity |
| `part 6` | Generator extraction and odd-kernel correlation |
| `part 7` | Multi-σ renormalisation (solve for α) |
| `part 8` | Operator closure failure & cross-σ SVD |
| `part 9` | Flow extrapolation tests |
| `part 10` | Acceleration vs dv/dσ, β-function RG |
| `part 11` | High-T RG stress tests |
| `part 12` | Controlled trajectory corruption (“t-kink”) |
| `part 13` | Dense near-zero cusp and curvature tests |
| `Part 14` | Energy curve collapse across windows |
| `Part 15` | Bracketed splitting rate near double root |
| `Part 16` | Winding number diagnostics under σ-flow |
| `Part 17` | σ-flow universality checks |
| `Part 18` | Local no-escape certification & Hadamard alignment |
| `Part 19` | RG flow of coefficients across height |

> There is no build system; scripts are intended to be read, modified,
> and run individually.

---

## 6. Relation to the de Bruijn–Newman Constant

- The classical Newman deformation corresponds to **backward heat flow**.
- Backward flow can destabilise the real-zero manifold if Λ > 0.
- This work focuses on **forward heat flow**, where:
- zero ordering is preserved,
- collisions are obstructed,
- the generator is empirically confined to a rigid rank-2 manifold.
- From this perspective, the conjectural bound Λ ≤ 0 may be interpreted as
a **dynamical stability condition**.

---

## 7. Limitations and Status

- The RG manifold and transverse instability are established **numerically**.
- A fully rigorous proof would require:
- excluding all finite-σ escape scenarios analytically,
- proving global attraction to the RG manifold.
- The present results should be read as **strong structural evidence**,
not a completed proof.

---

## 8. Intended Audience

This repository is aimed at readers with background in:

- analytic number theory,
- spectral theory,
- dynamical systems,
- or numerical experimentation with special functions.

Sceptical reading is encouraged.

---

