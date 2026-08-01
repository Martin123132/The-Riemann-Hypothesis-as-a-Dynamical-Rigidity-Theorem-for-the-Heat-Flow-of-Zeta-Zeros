# Newman Prime-Power Heat-Block Composition Guard

Date: 2026-07-27

Status: Exact reduction and countermodel-guard artifact; not a proof of
`Lambda<=0`, RH, PF-infinity, or a Clay-prize result.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_prime_power_heat_block_composition_guard.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_prime_power_heat_block_composition_guard.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_prime_power_heat_block_composition_guard.py
```

## Pi Provenance

The pi in a^2=x/(4*pi)+t/16 comes from the completed zeta normalization xi(s)=s*(s-1)*pi^(-s/2)*Gamma(s/2)*zeta(s)/2 and its Riemann-Siegel saddle. The cutoff-cell difference is 4*pi*(2N+1). The 2*pi in a winding is the period of exp(i*theta). No circle, prime-power block, or prefix polygon defines or approximates pi.

## Absolute Carrier Rate

Write the absolute carrier as `z_n=eta q_n=r_n zeta_n`. The exact
phase-rate identity is

```text
nu_n=partial_x arg(z_n)
    =v_a+b u_n+Im[d_(n,x)/(1+d_n)],
u_n=log(a/n).                                      (1)
```

The certified bounds give

```text
nu_n<=8449/x^2-u_n/2.                              (2)
```

Thus every `n<=N-1` carrier rotates strictly in the same absolute
direction. Only `n=N` can be near stationary, and that can happen
only in

```text
u_N<=16900/x^2,
0<=x-x_N<50701/x.                                  (3)
```

This is the terminal collar that must be combined with the already
proved adjacent endpoint recurrence.

## Prime-Power Blocks

Fix a prime `p`, a `p`-free integer `m`, and
`n_k=m p^k<=N`. Put `h=log p` and `u_k=log(a/n_k)`. Before the
`d_n` correction, consecutive amplitudes satisfy

```text
log(a_(k+1)/a_k)
 =-B_a h-(t/2)u_k h+(t/4)h^2.                     (4)
```

Whenever the next power remains below the cutoff, `u_k>=h`. Since
`B_a>49/100`, restoring the absolute `1+d_n` amplitudes still gives

```text
a_(k+1)/a_k<p^(-12/25)<1.                          (5)
```

The residual phase defect relative to one affine block phase is
less than `8756/x`. It is tiny, but a quantitative boundary-modulus
margin is needed before a homotopy may discard it.

## One-Block Theorem

For positive strictly decreasing coefficients, let

```text
Q_R(z)=sum_(k=0)^R a_k z^k,
P_R(z)=z Q_R(z).                                   (6)
```

Enestrom-Kakeya puts every zero of `Q_R` outside the closed unit
disk because `min a_k/a_(k+1)>1`. Thus `Q_R` is nonzero with
internal winding zero, while the shifted pure-power block `P_R`
has winding one. For a general p-free base, its external phase was
factored out and still moves with `x`; this is not an actual
blockwise `x`-winding theorem for the joined Xi sum.

## Gaussian Heat Identity

For a standard normal `Z` and `sigma=sqrt(t/2)`,

```text
exp(t log(n)^2/4)=E exp(sigma Z log n).             (7)
```

Hence a correction-free block is exactly an expectation of finite
geometric blocks. The tempting pointwise argument nevertheless
fails: after exponential tilting by the block base, the conditional
expanding-ratio mass is

```text
barPhi((1/2+(t/2)(u_m-delta_a))/sigma).             (8)
```

For a terminal block at `t=1/2`, this approaches
`barPhi(1)=0.158655...`. The Gaussian identity therefore does not
make every conditional block contract, and expectation does not
preserve zero-free winding.

## Two-Block Guard

Let `r=2^(-1/2)` and

```text
A(z)=z*(1+r*z+r^2*z^2),
B(z)=-(4/5)*z*(1+r*z).
```

Each block is nonzero on the unit circle and has winding one. Their
sum is

```text
A+B=(z/5)(1+r z+5r^2 z^2).                         (9)
```

The inner quadratic has discriminant `-19/2`; its two roots have
modulus `sqrt(2/5)<1`. Consequently the sum has winding three.
More generally this obstruction persists for every opposite block
scale `1/2<c<1`. This is an exact generic composition guard, not an Xi counterexample.

## Joined Structure

The blocks must be rejoined before winding is estimated. For
`s_Z=s_*-sigma Z`, define the `p`-free prefix

```text
O_M^(p)(s_Z)=sum_(m<=M, p does not divide m)m^(-s_Z).
```

Then

```text
sum_(n<=N)n^(-s_Z)
 =sum_k p^(-k s_Z) O_floor(N/p^k)^(p)(s_Z).        (10)
```

This identity retains the cross-block phases that the rejected
blockwise count loses.

## Live Theorem

The next `q>=1` theorem must combine the joined `p`-free prefixes,
the terminal recurrence, the exact endpoint, and the full first-jet
argument. It must establish

```text
|mathsf_X|<=delta_L => |mathcal_C_N|>A_L+epsilon_term,
0<=kappa_j<1,                                      (11)
```

with all connector, chart-join, and shoulder terms retained. The
`q<1` multiplicity-compatible parabolic/Hermite chart remains a
separate obligation.

## Boundary

This artifact proves the absolute saddle-frame carrier rate, the nonterminal rotation and terminal-collar theorem, exact prime-power heat ratios, corrected amplitude contraction, the correction-free normalized-block zero-free and shifted-winding theorem, the Gaussian mixture and joined p-free identities, and an exact generic two-block winding-three countermodel. It does not prove a uniform corrected-block homotopy, the joined Xi prefix lower bound, the strict composed successor flux bound, the q<1 multiplicity-compatible theorem, finite shoulder closure, contact exclusion, Lambda<=0, PF-infinity, RH, or a Clay-prize conclusion.
