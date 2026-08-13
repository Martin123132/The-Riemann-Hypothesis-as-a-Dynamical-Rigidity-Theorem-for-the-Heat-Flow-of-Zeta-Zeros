# Newman C1 Carrier-Near Odd-Harmonic Anchor

Date: 2026-08-04

Status: exact-lemma gate extracting two raw endpoint anchors and their phase/cutoff guards; this is not a proof of RH.

```text
work/rh_compute/results/jensen_window_pf_newman_c1_carrier_near_odd_harmonic_anchor_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_c1_carrier_near_odd_harmonic_anchor_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_c1_carrier_near_odd_harmonic_anchor_gate.py
```

## Odd-Harmonic Coordinate

Put

```text
D_N=sum_(m_N<=r<=n_N, r odd)1/r.
```

For the contiguous near kernel, direct integration of each Fourier
character gives

```text
integral_0^(1/2)K_(near,N)(v)dv=1/2+iD_N/pi
integral_(-1/2)^0K_(near,N)(v)dv=1/2-iD_N/pi
```

Thus the starred endpoint halves and the near kernel compose exactly as

```text
V_(near,N)[P]=integral_0^(1/2){f_P(1+v)-f_P(1)}K_(near,N)(v)dv+sum_(q=2)^(N-1)integral_(-1/2)^(1/2){f_P(q+v)-f_P(q)}K_(near,N)(v)dv+integral_(-1/2)^0{f_P(N+v)-f_P(N)}K_(near,N)(v)dv
mathscr N_N[P]=mathscr M_N[P]+iD_N{f_P(1)-f_P(N)}/pi+V_(near,N)[P]
V_(near,N)[P]+mathscr C_N[P]=-V_N[P]
mathcal J_P=2mathscr M_N[P]+iD_N{f_P(1)-f_P(N)}/pi-V_N[P]
```

No remote one-sided series has been separated.

## Physical Size

On the fixed physical q=2tL^2=1 chart with L>=50, the floor
inequalities and an integral comparison over odd integers prove

```text
(n_N+1)/(m_N+1)>a
D_N>=(1/2)log((n_N+1)/(m_N+1))
D_N>(1/2)log(a)>L/4
```

## Hermitian Anchor

The lower ideal endpoint is purely imaginary and the terminal endpoint
vanishes. Therefore

```text
mathcal J_H^0=-A_H+R_H
A_H=D_N*log(N)*log(a^2/N)^2/(4*pi)
R_H=2mathscr M_N[P_H^0]-V_N[P_H^0]
A_H>L^3*(L-2)/(128*pi)>14900 on q=1, L>=50
```

## Transpose Anchor

Keeping the terminal conditional survivor gives

```text
E_N=e(alpha_P*log(N))*exp(S(log(N)))
mathcal J_T^0=-A_T+R_T+Q_T
A_T=D_N*log(a^2/N)*log(N)^2/(4*pi)
R_T=2mathscr M_N[P_T^0]-V_N[P_T^0]
Q_T=-iD_N*u_(N,x)*log(a^2/N)/pi+iE_N*u_N*u_(N,x)*(2D_N-S_(1,N))/pi
A_T>L^2*(L-2)^2/(128*pi)>14300 on q=1, L>=50
```

## What This Buys

The anchors are useful without choosing a phase if the stationary
remainders can be placed in smaller disks:

```text
|R_H|<=A_H-delta implies |mathcal J_H^0|>=delta
|R_T|+|Q_T|<=A_T-delta implies |mathcal J_T^0|>=delta
```

They are not standalone sign terms. The terminal multiplier rotates the
negative real centers, and D_N changes when an odd reciprocal mode crosses
a roster boundary. The complete Hermitian/transpose projection and exact
mode transfer must remain attached.

## Next Stage

Use the exact negative anchors only as centers of phase-independent modulus disks. Derive a reciprocal-stationary decomposition of R_H=2M_H-V_H and R_T=2M_T-V_T that preserves diagonal/adjacent cancellation and the contact-band condition. First test on the fixed q=1 chart whether |R_H|<A_H or |R_T|+|Q_T|<A_T is even scale-compatible. Do not infer a physical-current sign from the unrotated anchor and do not freeze D_N across an odd roster transfer.

## Pi Provenance

The factor pi comes only from e(y)=exp(2pi i y): integrating an odd
Fourier mode across a half-cell gives i/(pi r). The rational estimate
pi<22/7 is used only to display conservative decimal-free anchor floors.
No circle, polygon, or fitted plotting constant is introduced.

## Proof Boundary

This gate proves exact endpoint-anchor decompositions and quantitative raw anchor lower bounds on the fixed q=1 chart. It proves no extension of those floors to other parabolic frequencies, stationary carrier-variation bound, complete physical-current sign, contact-box exclusion, ordinary retained-main separation, Q209 shell, cofinal successor, Lambda<=0, RH, or prize-level conclusion.

validated Newman C1 carrier-near odd-harmonic anchor gate: 18 rows, 0 issues, 96 half-cell identities, q=1 D_N>L/4, 2 explicit raw anchors, 1 phase-rotation guard, 1 live stationary residual
