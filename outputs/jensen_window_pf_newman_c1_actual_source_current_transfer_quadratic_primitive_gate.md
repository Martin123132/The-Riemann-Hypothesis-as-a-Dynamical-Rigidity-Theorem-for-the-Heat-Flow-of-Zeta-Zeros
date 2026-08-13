# Newman C1 Current/Transfer Ownership and Quadratic Primitive Gate

Date: 2026-08-05

Status: exact ownership correction and determinant primitive; not a proof of the signed current or RH.

## Ownership corrigendum

D_move(xi)=tilde(Psi)'(xi)-Psi_B'(xi) is a pointwise current defect. R_move^tr=-varepsilon integral_0^1D_move(T_0-theta varepsilon)dtheta is its integrated transfer. The two quantities must have distinct names and budgets.

The 1/1000000 theorem remains valid for E_move^tr=2R_move^tr/|nu|^2. It does not bound the pointwise E_move^cur used in Sigma_full.

At the physical point,

```text
At the physical point xi=Omega, Delta omega_C(Omega)=0, so R_move^cur=D_move(Omega)=2[omega_E+omega_B(Omega)+omega_C(Omega)]^TJdot(omega_C)(Omega).
|E_move^cur|<8P_*^2S_0T_1+(24/5)(B_0T_0)P_*T_1<(128/7)hR_*^4K^2+(48/5)(B_0T_0)h^(3/2)R_*^2K^2.
|E_move^cur|<rho A_T/100.
```

The corrected six-package current budget is

```text
|E_R+E_term+E_aa+E_Delta+E_aff^C+E_move^cur|<(289339/858000)rho A_T.
```

## Quadratic primitive

```text
E_quad^cur=(2/|nu|^2){e_Ve_(mathcalN,xi)+e_mathcalNe_(V,xi)-e_Ae_(Q,xi)-e_Qe_(A,xi)}, with e containing only the nonterminal remainder after epsilon=a+rho_obs.
Q_quad=(2/|nu|^2)(e_Ve_mathcalN-e_Ae_Q), and E_quad^cur=dQ_quad/dxi exactly on every fixed chart.
For xi_theta=T_0-theta varepsilon and Omega=T_0-varepsilon, -varepsilon integral_0^1E_quad^cur(xi_theta)dtheta=Q_quad(Omega)-Q_quad(T_0).
```

For normalized complex nonterminal amplitudes U_X, Q_quad=Re{U_V conjugate(U_mathcalN)-U_A conjugate(U_Q)+U_VU_mathcalN-U_AU_Q}; differentiation gives the Hermitian phase-difference and transpose phase-sum currents without discarding either phase.

E_quad contains no terminal a factor. E_aa is the a-a block; E_term and the finite/remote/correction channels contain terminal-nonterminal terms; E_aff is linear; E_move^cur is the terminal-motion current defect.

## Revised frontier

Define E_main^cur=E_FAN+E_quad^cur. With the corrected absolute budget, adverse-witness negativity necessarily requires E_main^cur<-(551501/858000)rho A_T. The stronger E_main^cur<-(1130179/858000)rho A_T is sufficient independently of the absolute package signs.

Choose one coherent route: prove the signed pointwise inequality for E_FAN+E_quad^cur, or use the transfer primitive and bound the endpoint determinant difference jointly with tilde(Psi)(T_0). No separate |E_quad| budget is assumed.

This is not a proof of a complete-current sign, Lambda<=0, PF-infinity, RH, or a prize-level conclusion.

built Newman C1 current/transfer quadratic-primitive gate: 17 rows, 0 issues, pointwise |E_move^cur|<rho*A_T/100, 1 determinant primitive, 1 signed-main reclassification, 1 open theorem row
