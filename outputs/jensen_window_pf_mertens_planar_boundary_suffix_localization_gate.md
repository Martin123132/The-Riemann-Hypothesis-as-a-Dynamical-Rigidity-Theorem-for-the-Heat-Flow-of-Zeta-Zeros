# Jensen-Window PF Mertens Planar Boundary Suffix-Localization Gate

Date: 2026-07-24

Status: exact bulk-localization, terminal-collar recovery, cofinal
RH-equivalence, and strict finite-dimensional norm-separation guard for
the positive boundary suffix energy. This is not an estimate for that
energy and is not a proof of RH, PF-infinity, `Lambda<=0`, or a
Clay-prize result.

```text
work/rh_compute/results/
jensen_window_pf_mertens_planar_boundary_suffix_localization_gate.json
python work/rh_compute/scripts/
jensen_window_pf_mertens_planar_boundary_suffix_localization_gate.py
python work/rh_compute/scripts/
check_jensen_window_pf_mertens_planar_boundary_suffix_localization_gate.py
```

## Odd-Sine Tail And The Interior Plateau

Retain `0<alpha<1`, dyadic `K>=4`, and the notation of Corollary
11.22Z.9. Put

```text
N:=2K+1,             R:=ceil(K^(1-alpha/2)),
q_r:=2r-1,           a:=pi/N,

S_R(u):=sum_(r=1)^R sin(q_r u)/q_r.                 (PBSL.1)
```

For finite `1<=h<=ell`,

```text
sum_(r=h)^ell sin(q_r u)
 =sin((ell-h+1)u)sin((h+ell-1)u)/sin(u).            (PBSL.2)
```

Thus every consecutive odd-sine partial sum has modulus at most
`1/sin(u)`. Abel summation of the infinite tail of the square-wave
series gives the explicit estimate

```text
|S_R(u)-pi/4|
 <=1/[(2R+1)sin(u)],                   0<u<pi.       (PBSL.3)
```

Define

```text
eta_R:=arcsin(8/[pi(2R+1)]),
L:=ceil(N eta_R/pi),
J:=min(K,K+2-L).                                    (PBSL.4)
```

Then

```text
S_R(u)>=pi/8,                         eta_R<=u<=pi-eta_R,

L<=1+4N/[pi(2R+1)]
 <=1+2N/(pi R).                                     (PBSL.5)
```

The minimum in the definition of `J` matters when `L=1` or `L=2`.

## Exact Slivers And Bulk Coercivity

Write

```text
Delta_i:=w_i-w_(i-1),        1<=i<=K,
```

where `w_K=c`. Product-to-sum followed by one integration gives, for
`1<=i<K`,

```text
Delta_i
 =(2/N)[
   integral_((K-1-i)a)^((K-i)a) S_R(u)du
  +integral_((K-2+i)a)^((K-1+i)a) S_R(u)du
 ].                                                     (PBSL.6)
```

At the terminal anchor,

```text
Delta_K
 =(2/N)integral_(pi-3a)^(pi-a)S_R(u)du
 =(2/N)integral_a^(3a)S_R(u)du.                         (PBSL.7)
```

If `L>=3`, the right sliver in (PBSL.6) lies inside the plateau for
every `i<=J`. If `L<=2`, this remains true for all current indices and
the endpoint interval in (PBSL.7) has a plateau subinterval of length
at least `a`. Therefore

```text
Delta_i>=pi^2/(4N^2),                    1<=i<=J.       (PBSL.8)
```

Set

```text
B_i:=A_i+beta,          1<=i<K,
B_K:=beta,

E_B:=sum_(i=1)^K Delta_i B_i^2,
H_K:=sum_(i=1)^K B_i^2,
Q_J:=sum_(i=1)^J B_i^2.
```

The bulk lower bound gives

```text
Q_J<=4N^2 E_B/pi^2.                                  (PBSL.9)
```

## Terminal-Collar Recovery

The arithmetic coefficients enter only through the deterministic
increment law

```text
B_i-B_(i+1)=u_i,             |u_i|<=1.               (PBSL.10)
```

Put `m=K-J`. For `1<=t<=m`,

```text
|B_(J+t)|<=|B_J|+t.
```

Consequently,

```text
H_K
 <=(1+2m)Q_J+m(m+1)(2m+1)/3

 <=[4N^2(2L+1)/pi^2]E_B+2L^3.                       (PBSL.11)
```

For `R=ceil(K^(1-alpha/2))`, equation (PBSL.5) gives
`L=O(K^(alpha/2))`. Hence an all-epsilon pointwise estimate

```text
E_B=O_delta(K^delta) for every delta>0                (PBSL.12)
```

implies, after choosing `delta<alpha/2`,

```text
K^(-2-alpha)H_K
 =O(K^(-alpha/2+delta))+O(K^(-2+alpha/2)).            (PBSL.13)
```

Both terms are summable over dyadic `K`. By the exact affine-tent
criterion (MATBH.39), (PBSL.12) at one fixed positive `alpha` implies
`R_alpha<infinity`.

The quantifier boundary is essential. On every member of one fixed
cofinal sequence `alpha_j->0+`, (PBSL.12) implies RH through the prior
cofinal reciprocal-tail criterion. Conversely, RH gives

```text
M(x)=O_delta(x^(1/2+delta)).
```

The weight defining `beta_K` is monotone and has uniformly bounded
variation. Partial summation therefore gives

```text
beta_K=O_delta(K^(1/2+delta)),
B_i=O_delta(K^(1/2+delta)),
H_K=O_delta(K^(2+delta)).                            (PBSL.14)
```

The already proved upper comparison

```text
E_B<=pi(1/2+2pi)K^(-2)H_K
```

then gives (PBSL.12). Thus, for any fixed cofinal positive-alpha
sequence,

```text
RH
 iff
E_(B,alpha_j,K)=O_delta(K^delta) for every delta>0,
every j, and every dyadic K.                         (PBSL.15)
```

This is a theorem-scale equivalence, not an estimate of `E_B`.

## Strict Norm-Separation Countermodel

The reverse comparison is not a uniform norm inequality. For each
`K`, choose bounded synthetic coefficients

```text
u_(K-1)=-1,       u_i=0 otherwise,
v_(2K)=1,         v_m=0 otherwise.                    (PBSL.16)
```

Because the future weight at `2K` is one,

```text
beta=1,       B_i=0 for i<K,       B_K=1,
H_K=1,        E_B=Delta_K.                            (PBSL.17)
```

Direct spectral differencing gives

```text
0<Delta_K
 =(4/N)sum_(r=1)^R
   sin((2K-1)theta_r)sin(theta_r)/q_r^2
 <=8pi^2 R/N^3.                                      (PBSL.18)
```

Therefore

```text
K^2 E_B/H_K
 <=8pi^2 R K^2/N^3
 =O(R/K)
 =O(K^(-alpha/2))
 ->0.                                                (PBSL.19)
```

No constant `c>0` can make `E_B>=c K^(-2)H_K` hold uniformly for
bounded coefficient vectors. The two facts are compatible: `E_B`
misses a short terminal collar as a norm, while (PBSL.10) recovers that
collar at the all-epsilon cofinal theorem scale.

## Finite Diagnostic

At `alpha=1/2`, direct computation gives:

| K | R | L | J | min bulk N^2 Delta | K^2 Delta_K | Mobius K^2 E_B/H_K |
|---:|---:|---:|---:|---:|---:|---:|
| 16 | 8 | 2 | 16 | 6.967305490 | 2.551726426 | 2.352630050 |
| 32 | 14 | 2 | 32 | 7.068824342 | 2.596386886 | 2.376051298 |
| 64 | 23 | 3 | 63 | 7.247416957 | 2.497087829 | 2.427715794 |
| 128 | 39 | 3 | 127 | 7.227781985 | 2.333757268 | 2.443250924 |
| 256 | 64 | 4 | 254 | 9.125891467 | 2.078006159 | 2.452812821 |
| 512 | 108 | 4 | 510 | 8.721113843 | 1.842244233 | 2.462827953 |
| 1024 | 182 | 5 | 1021 | 9.268712610 | 1.608344198 | 2.462331675 |


The rigorous bulk threshold is `pi^2/4=2.467401100...`; every displayed
row is above it. The terminal column is the strict synthetic witness.
The Mobius column is diagnostic only and proves no asymptotic estimate.

## Route Decision

The positive boundary/future complement is not an easier independent
gate. Its all-epsilon estimate on a cofinal positive-alpha sequence is
already RH-equivalent. The direct planar route can still be useful for
structural decomposition, but it cannot claim progress from positivity
alone: it must prove an RH-strength anchored cancellation theorem.

The inherited current-interior projection estimate for `Y_(K,r)`
remains separately open. No `E_B` estimate, `Y` estimate, RH,
PF-infinity, `Lambda<=0`, or Clay-prize result follows here.

## Claim Ledger

| ID | Role | Status | Statement | Boundary |
|---|---|---|---|---|
| `pbslg_01_inherited_energy` | `exact_definition` | `available_exact` | Retain N=2K+1, R=ceil(K^(1-alpha/2)), Delta_i=w_i-w_(i-1), B_i=A_i+beta for i<K, B_K=beta, E_B=sum_(i=1)^K Delta_i B_i^2, and H_K=sum_(i=1)^K B_i^2. | Here 0<alpha<1, K is dyadic, and K>=4. |
| `pbslg_02_odd_sine_partial_sum` | `exact_identity` | `available_exact` | For q_r=2r-1, sum_(r=h)^ell sin(q_r u)=sin((ell-h+1)u)sin((h+ell-1)u)/sin(u). | The formula holds for 0<u<pi and finite 1<=h<=ell. |
| `pbslg_03_square_wave_tail` | `exact_bound` | `available_exact` | For S_R(u)=sum_(r<=R)sin(q_r u)/q_r, |S_R(u)-pi/4|<=1/[(2R+1)sin(u)] on 0<u<pi. | This is Abel summation of the infinite odd-sine tail. |
| `pbslg_04_interior_parameters` | `exact_definition` | `available_exact` | Set eta_R=arcsin(8/[pi(2R+1)]), L=ceil(N eta_R/pi), and J=min(K,K+2-L). | The minimum at K is required when L is one or two. |
| `pbslg_05_increment_slivers` | `exact_identity` | `available_exact` | Writing a=pi/N, Delta_i for i<K is 2/N times the S_R integrals over [(K-1-i)a,(K-i)a] and [(K-2+i)a,(K-1+i)a], while Delta_K=(2/N)integral_(a)^(3a)S_R(u)du. | The endpoint formula uses S_R(pi-u)=S_R(u). |
| `pbslg_06_interior_plateau` | `exact_bound` | `available_exact` | S_R(u)>=pi/8 throughout [eta_R,pi-eta_R]. | No arithmetic coefficient enters this spectral bound. |
| `pbslg_07_bulk_increment_lower_bound` | `exact_bound` | `available_exact` | Every 1<=i<=J has an interior sliver of length at least pi/N, hence Delta_i>=pi^2/(4N^2). | For L<=2 the endpoint sliver is included; otherwise J<K. |
| `pbslg_08_localized_mass` | `exact_bound` | `available_exact` | Q_J:=sum_(i=1)^J B_i^2 satisfies Q_J<=4N^2 E_B/pi^2. | This uses only positivity and the bulk increment lower bound. |
| `pbslg_09_suffix_lipschitz` | `exact_bound` | `available_exact` | Because B_i-B_(i+1)=u_i and |u_i|<=1, |B_(J+t)|<=|B_J|+t for 1<=t<=K-J. | The future anchor beta need not be bounded separately. |
| `pbslg_10_exact_collar_recovery` | `exact_bound` | `available_exact` | For m=K-J, H_K<=(1+2m)Q_J+m(m+1)(2m+1)/3. | This is a deterministic bounded-increment inequality. |
| `pbslg_11_coarse_collar_recovery` | `exact_bound` | `available_exact` | Since m<=L, H_K<=[4N^2(2L+1)/pi^2]E_B+2L^3. | The terminal collar is the only loss in the reverse comparison. |
| `pbslg_12_collar_length` | `exact_bound` | `available_exact` | The elementary arcsine chord bound gives L<=1+4N/[pi(2R+1)]<=1+2N/(pi R). | For R=ceil(K^(1-alpha/2)), L=O(K^(alpha/2)). |
| `pbslg_13_fixed_alpha_summability` | `exact_reduction` | `available_exact` | If E_B=O_delta(K^delta) for every delta>0 at a fixed alpha, then sum_(K dyadic)K^(-2-alpha)H_K is finite. | Choose delta<alpha/2; both collar terms are dyadically summable. |
| `pbslg_14_affine_tent_composition` | `theorem_composition` | `available_exact` | By MATBH.39, fixed-alpha all-epsilon subpower control of E_B implies R_alpha<infinity. | One fixed alpha gives only its shifted reciprocal-tail criterion. |
| `pbslg_15_cofinal_forward_implication` | `theorem_composition` | `available_exact` | If the E_B subpower condition holds for every member of one fixed cofinal sequence alpha_j->0+, then RH follows. | No alpha=0 or uniform-in-alpha passage is used. |
| `pbslg_16_rh_anchor_bound` | `theorem_composition` | `available_exact` | Under RH, M(x)=O_delta(x^(1/2+delta)); partial summation gives beta_K=O_delta(K^(1/2+delta)) and B_i=O_delta(K^(1/2+delta)) uniformly in i. | The beta weight is monotone with uniformly bounded variation. |
| `pbslg_17_rh_energy_bound` | `theorem_composition` | `available_exact` | Under RH, H_K=O_delta(K^(2+delta)); PBSE.17 then gives E_B=O_delta(K^delta) for every delta>0. | The epsilon parameter is renamed after the square. |
| `pbslg_18_cofinal_equivalence` | `exact_equivalence` | `available_exact` | Along any fixed cofinal positive-alpha sequence, the all-epsilon E_B subpower condition is equivalent to RH. | This is a composition of the collar theorem and prior RH criteria. |
| `pbslg_19_synthetic_terminal_vector` | `countermodel_definition` | `available_exact` | Set u_(K-1)=-1, all other current u_i=0, put one future coefficient v_(2K)=1, and set all other future coefficients to zero. | All synthetic coefficients have modulus at most one. |
| `pbslg_20_synthetic_energies` | `exact_identity` | `available_exact` | For the terminal vector beta=1, B_i=0 for i<K, H_K=1, and E_B=Delta_K=c-w_(K-1)>0. | The witness isolates the endpoint anchor direction. |
| `pbslg_21_endpoint_increment_bound` | `exact_bound` | `available_exact` | Spectral differencing gives 0<Delta_K<=8pi^2 R/N^3. | The bound uses |sin(2theta)sin(theta)|<=2theta^2 modewise. |
| `pbslg_22_no_uniform_lower_comparison` | `countermodel_guard` | `guard_active` | For the terminal vector, K^2 E_B/H_K<=8pi^2 R K^2/N^3=O(R/K)->0. | Therefore no uniform E_B>=c K^(-2)H_K norm bound exists. |
| `pbslg_23_finite_diagnostic` | `finite_validation` | `validated_finite` | At alpha=1/2 and K=16,...,1024, every bulk increment exceeds the proved threshold and direct Mobius values give K^2 E_B/H_K near 2.35--2.47. | Finite values prove no asymptotic Mobius estimate. |
| `pbslg_24_route_guard` | `route_guard` | `guard_active` | E_B is strictly weaker as a finite-dimensional norm, but bounded coefficient increments recover its missing O(K^(alpha/2)) terminal collar at the cofinal theorem scale. | The boundary route is an RH-equivalent gate, not an easier shortcut. |
| `pbslg_25_open_energy_gate` | `theorem_target` | `open` | Prove E_B=O_delta(K^delta) for every delta>0 on every dyadic scale and every member of one fixed cofinal positive-alpha sequence, without assuming RH. | This target is now explicitly calibrated as RH-equivalent. |
| `pbslg_26_proof_boundary` | `proof_guard` | `guard_active` | The localization, collar recovery, cofinal composition, and synthetic norm-separation guard are exact. | No E_B estimate, Y projection estimate, RH, PF-infinity, Lambda<=0, or Clay-prize result is proved. |

Summary:

- rows: 26
- exact/theorem-composed results: 22
- finite validations: 1
- active route/proof guards: 3
- open RH-equivalent energy gates: 1
