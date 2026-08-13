# Newman C1 Actual-Source Saddle-Phase Lock and Centre Symmetry Gate

Date: 2026-08-04

Status: exact actual-source normalization, reverse full-turn theorem, and bare-centre symmetry proved; joined current sign open; not a proof of RH.

## Source phase

The common source normalizer is nu=eta/S_a=omega_c/(B_0T_0), where omega_c=(-1)^N eta and eta=phi(1+d_1)/|1+d_1|. At xi=Omega one has c_xi=1, |nu|=(B_0T_0)^(-1), and exp(2i gamma_phys)=eta^2; cutoff parity cancels from the transpose phase.

For Q(v)=exp(i[v log(v/(2pi))-v-pi/4]) and gamma_t=conjugate(M_t)/M_t=exp(-2i beta_t), Section 11.119 supplies the continuous psi=arg(Q(Omega)/gamma_t), with |psi|<1/x. Hence arg Q(Omega)+2beta_t=psi exactly on the chosen lift.

## Reverse turn

```text
2chi_N=F_N(Omega)+pi/4+psi+2delta_1,
F_N(v)=v[1+log(2pi N^2/v)],
Delta chi_N<-2pi.
```

Put u=1/N. The alternating lower bound log(1+u)>=u-u^2/2+u^3/3-u^4/4 yields G_N<=-2-(7/16)u. Thus the ideal increment obeys Delta(2chi_N)<=-4pi-7pi/(8N).

Replacing the ideal endpoint frequencies by Omega=2pi a^2-epsilon contributes less than 1/N^2 to Delta F_N. The lock phases contribute |Delta psi|<2/(9N^2). Since |arg(1+d_1)|<2|d_1|<4378/x, the term 2Delta delta_1 contributes less than 17512/(9N^2). The total perturbation of Delta(2chi_N) is therefore less than 1947/N^2.

The continuous absolute terminal carrier eta exp(i phi_N)=exp(i chi_N) covers the full unit circle on every complete cell, now in the reverse direction. This is source-inclusive and is distinct from the much faster positive winding of exp(i phi_N).

The previous Delta phi_N>2pi and the new Delta chi_N<-2pi imply Delta theta_eta=Delta chi_N-Delta phi_N<-4pi. The transpose phase phi_N+2theta_eta=2chi_N-phi_N has increment below -6pi. Thus the terminal, source, absolute-terminal, and transpose phases all cover the circle; their correlations, not a fixed sector, are the remaining content.

## Centre symmetry

At xi=Omega, write tau_0=rho exp(i[phi_N-pi/2]), rho>0. The actual-source raw-centre current is R_center=-(rho|nu|^2/2)[A_H sin(phi_N)+A_T sin(phi_N+2theta_eta)].

```text
R_center=-rho|nu|^2 A_T[
  cos(theta_eta)sin(phi_N+theta_eta)
  +(u_N/log N)sin(phi_N)].
```

Equivalently, R_center=-rho|nu|^2 A_T[Re(eta)Im(eta exp(i phi_N))+(u_N/b)Im(exp(i phi_N))]. Choosing omega_c=(-1)^Neta instead gives the same product, so this form is parity invariant and branch free.

Across a complete cell, 0<=u_N<=log(1+1/N)<1/N and b=log N>1/2, hence 0<=u_N/b<2/N. The raw centres are therefore an extremely close equal-anchor pair, but this small coefficient does not control the current near zeros of Re(eta) or Im(eta exp(i phi_N)).

D_N changes at internal odd-roster ties. The ratio and factorization hold on both adjacent raw-centre charts, but the centre is not itself a tie-invariant physical observable. Any sign theorem must restore the reciprocal residual, terminal survivor, correction rows, affine row, and moved-tail defect before crossing a tie.

## Gate rows

| id | role | state | claim |
|---|---|---|---|
| asl_01_normalizer | source normalizer | proved | The physical source phase and modulus are fixed exactly. |
| asl_02_lift | phase lift | proved | The source and terminal phases admit continuous unwrapped lifts. |
| asl_03_lock | stationary lock | proved | The reciprocal stationary phase is conjugate-locked to the normalizer. |
| asl_04_identity | absolute phase identity | proved | The source-inclusive terminal phase has an exact scalar normal form. |
| asl_05_endpoint | ideal endpoint | proved | The ideal complete-cell increment is one elementary logarithmic expression. |
| asl_06_decrement | ideal decrement | proved | The ideal absolute phase decreases by more than one full turn. |
| asl_07_scale | physical scale | proved | The actual q=1 cells make every endpoint correction tiny relative to the turn margin. |
| asl_08_error | physical perturbation | proved | All physical phase-lock and first-anchor corrections fit one O(N^-2) budget. |
| asl_09_reverse | reverse full turn | proved | The actual source-inclusive terminal phase loses more than 2pi per complete cell. |
| asl_10_circle | absolute circle | proved | The actual absolute terminal carrier covers every unit direction. |
| asl_11_inventory | phase inventory | proved | The source and transpose phases also make complete turns. |
| asl_12_ratio | anchor ratio | proved | The Hermitian and transpose raw anchors differ by exactly 2u_N/log N. |
| asl_13_current | actual centre | proved | The physical source can be inserted into the bare two-channel centre exactly. |
| asl_14_symmetry | centre symmetry | proved | The bare centre factors into the source and absolute-terminal phases plus one mismatch. |
| asl_15_branch | branch-free form | proved | The symmetry has a parity-invariant branch-free expression. |
| asl_16_mismatch | mismatch bound | proved | The unequal-anchor term is below 2/N relative to A_T. |
| asl_17_ties | roster guard | proved | The raw-centre factorization cannot be promoted through ties by itself. |
| asl_18_joint | joint phase target | open | A signed Xi-specific joint-phase product estimate remains open. |
| asl_19_complete | complete current | open | The joint phase theorem must be recomposed with every retained residual and correction row. |

## Next action

Derive a local scale-separation theorem for the joint path Re(eta)*Im(eta exp(i phi_N)) on roster interiors. Use the exact reverse-turn coordinate chi_N and the much faster terminal phase phi_N to determine whether the dominant product must take both signs, while retaining the O(u_N/log N) mismatch. If that centre-level sign oscillates, promote it as a route obstruction and move directly to the complete reciprocal residual; if a one-sided Xi correlation survives, insert it into the tie-invariant complete (h,t) current before estimating corrections.

## Pi provenance

The pi in Q(v)=exp(i[v log(v/(2pi))-v-pi/4]) is the completed-zeta reciprocal stationary phase, with pi/4 its negative-curvature signature. The factor 2pi in K_N=2pi N^2 and T_0=2pi a^2 is the same completed-zeta/Fourier normalization. The threshold 2pi is one period of exp(i chi_N). Only elementary pi>3 is used in the coarse endpoint-error audit.

## Proof boundary

This gate proves the exact physical source normalizer, a source-inclusive reverse full turn on every complete q=1 cutoff cell, the induced source and transpose winding inventory, and the exact parity-invariant bare-centre symmetry with its O(u_N/log N) anchor mismatch. It proves no sign for the joint phase product, raw-centre sign theorem, signed reciprocal residual estimate, complete actual-source current inequality, correction-row bound, tie-spliced current theorem, all-q transport, contact-box exclusion, Q209 shell, cofinal successor, Lambda<=0, PF-infinity, RH, or prize-level conclusion.

built Newman C1 actual-source saddle-phase lock/centre-symmetry gate: 19 rows, 0 issues, 1 reverse full-turn theorem, 1 source-inclusive circle theorem, 2 exact centre factorizations, 1 live joint-phase target
