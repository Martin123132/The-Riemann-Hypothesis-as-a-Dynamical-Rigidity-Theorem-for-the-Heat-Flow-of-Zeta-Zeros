# Natural A endpoint-integral carrier lift

Date: 2026-08-27

Status: exact transition-cell carrier alignment; no transition coefficient bound.

Let `s=1/2+it`, `c_t=(pi/(32t))^(1/4)`, and let
`C_A,m^exact` denote the convergent exact endpoint-wedge integral from
Formal Core (11.428.1).  For the odd endpoint `A`, endpoint orientation
`epsilon_A=-1`, phase `Phi_A,m`, and prefactor `W0_A,m`, define

```text
X_A,m=c_t exp(-i*pi/8) epsilon_A exp(i Phi_A,m)
      W0_A,m C_A,m^exact,
a_m=2 Re X_A,m.
```

Odd-endpoint parity gives

```text
exp(-i*pi/8)exp(i Phi_A,m)=exp(i theta_0)m^(-it),
theta_0=t[log(t/(2pi))-1]/2-pi/8.
```

Consequently

```text
X_A,m=m^(-s) eta_A,m,
eta_A,m=epsilon_A c_t sqrt(m) exp(i theta_0)
        W0_A,m C_A,m^exact.
```

The natural Hardy lift and its same-cell coefficient are

```text
mathcal A_A,m^nat=H(t)exp(-i theta(t))X_A,m
                 =m^(-s)C_A,m,
C_A,m=H(t)exp(-i theta(t))eta_A,m,
Hardy_t[mathcal A_A,m^nat]=H(t)a_m.
```

It differs from the earlier real canonical lift by a Hardy-null complex
component.  Keeping that component is essential because transition cells now
reduce exactly to the single pre-norm coefficient

```text
m^(-s)[beta_m-C_G-C_A,m], 39853<=m<=39936.
```

The production replay reconstructs saved complex A rows at modes
`39853, 39894, 39936`.  A separate altered witness uses height
`23.75`, odd endpoint
`19`, and modes
`[2, 4, 7]`.

Pi provenance: every `pi` above descends from the equation-(9) endpoint phase,
odd-endpoint Fourier parity, the exact bi-Morse normalization, the paper
rotation, or the Riemann-Siegel phase.  No fitted constant is introduced.

Proof boundary: exact natural complex lift and common `m^(-s)` coefficient
alignment only.  This gate does not provide 84 modewise full-endpoint complex
enclosures, bound `beta_m-C_G-C_A,m`, join `P_W` to unowned cells, enclose
`J_Z` or `D_K`, or prove a non-A, all-height, `Lambda<=0`, PF-infinity, RH, or
prize-level conclusion.
