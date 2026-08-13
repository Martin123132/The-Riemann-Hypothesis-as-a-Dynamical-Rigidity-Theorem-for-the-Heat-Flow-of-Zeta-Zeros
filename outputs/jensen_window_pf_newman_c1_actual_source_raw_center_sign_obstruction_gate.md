# Newman C1 Actual-Source Raw-Centre Sign Obstruction Gate

Date: 2026-08-04

Status: actual physical source phase proved to force both raw-centre signs on every complete q=1 cell; complete residual current open; not a proof of RH.

## Scale separation

Write chi_N=chi_0+r with chi_0=(1/2)F_N(Omega)+pi/8, F_N(v)=v[1+log(2pi N^2/v)], and r=psi/2+arg(1+d_1).

The cell crosses K_N=2pi N^2 once. From that crossing to the right endpoint, F_N(Omega) decreases by more than 4pi, so chi_0 decreases by more than 2pi.

In any real interval of phase length greater than 2pi, the periodic set sin(theta)>=1/2 has a connected component of intersection of width at least pi/3. Choose a subinterval J on which chi_0 decreases by pi/3 and sin(chi_0)>=1/2.

With phi_N=Omega log N and Omega>=K_N on J, |d chi_0/d phi_N|=log(Omega/K_N)/(2log N)<log(1+1/N)/log N<2/N.

The pi/3 change of chi_0 therefore requires Delta phi_N>pi N/6. Thus theta_0=chi_0-phi_N decreases by more than pi N/6+pi/3>2pi on J, and assumes representatives with cos(theta_0)=+1 and cos(theta_0)=-1.

## Strict sign witnesses

The phase-lock errors give |r|<1/(2x)+4378/x<487/N^2<1/1000 throughout the cell. The exact anchor mismatch e_N=u_N/log N obeys 0<=e_N<2/N<1/1000.

At a point of J with theta_0=0 modulo 2pi, theta_eta=theta_0+r and chi_N=chi_0+r. Hence cos(theta_eta)>999/1000 and sin(chi_N)>499/1000. The normalized bracket B=cos(theta_eta)sin(chi_N)+e_N sin(phi_N) is greater than 497501/10^6>49/100.

At a point of J with theta_0=pi modulo 2pi, cos(theta_eta)<-999/1000 while sin(chi_N)>499/1000. The same normalized bracket satisfies B<-497501/10^6<-49/100.

```text
R_center=-rho|nu|^2 A_T B,
B>49/100 at one point,
B<-49/100 at another point.
```

Since R_center=-rho|nu|^2A_T B with rho|nu|^2A_T>0, every complete physical q=1 cell contains raw-centre points with R_center<-(49/100)rho|nu|^2A_T and R_center>(49/100)rho|nu|^2A_T. The actual physical source phase therefore does not rescue a cell-uniform raw-centre sign.

The selected points may be perturbed off the finite roster-tie set without losing the strict 49/100 reserve. The raw-centre sign obstruction is valid on roster interiors, but only the complete residual-composed current is tie invariant.

With W_H=R_H and W_T=R_T+Q_T, the complete ideal actual-source scalar is Sigma_id=-2rho A_TB+E_W, where E_W=Re(conjugate(tau_0)W_H)+Re(eta^2tau_0W_T), and R_id=(|nu|^2/2)Sigma_id.

At the witness B<-49/100, negativity of the complete ideal current requires E_W<-(98/100)rho A_T. Since |E_W|<=rho(|W_H|+|W_T|), a necessary condition is |W_H|+|W_T|>(98/100)A_T. The residual must therefore be macroscopic and correctly signed on the adverse arc; this may be a uniform negative bias or a phase-correlated response.

## Gate rows

| id | role | state | claim |
|---|---|---|---|
| rco_01_main | locked main phase | proved | The actual absolute phase splits into an explicit main and a uniformly tiny remainder. |
| rco_02_monotone | frequency monotonicity | proved | The physical frequency is strictly increasing across each complete cell. |
| rco_03_turn | main reverse turn | proved | The explicit locked main makes more than one reverse turn after the saddle crossing. |
| rco_04_arc | safe sine arc | proved | One connected locked-phase arc has sine at least one half and width pi/3. |
| rco_05_speed | scale separation | proved | The locked phase moves at less than 2/N of the terminal phase speed. |
| rco_06_sweep | source half-turn pair | proved | The source phase makes more than a full turn inside the safe locked-phase arc. |
| rco_07_remainder | physical remainder | proved | The actual phase and unequal-anchor errors are each below 1/1000. |
| rco_08_positive | positive bracket witness | proved | One roster-interior point has normalized centre bracket above 49/100. |
| rco_09_negative | negative bracket witness | proved | A second roster-interior point has normalized centre bracket below -49/100. |
| rco_10_obstruction | actual-source centre obstruction | proved | The raw-centre current takes both signs at the actual physical source phase. |
| rco_11_ties | tie guard | proved | The strict witnesses can be chosen away from every roster tie. |
| rco_12_compensation | residual compensation | proved | Any negative complete ideal current needs a macroscopic correctly signed residual on the positive-centre arc. |
| rco_13_target | complete current target | open | Only a signed complete residual current can now restore a one-sided theorem. |

## Next action

Stop trying to obtain the physical sign from either raw anchor centre. Compose the full tie-invariant actual-source scalar h+Re(eta^2 t) using R_H, R_T+Q_T, the physical correction polynomials, fixed affine row, and moved-tail defect. Determine whether it has a cell-uniform negative bias or a phase-correlated O(A_T) response on the strict adverse arcs. A lower bound at one adverse witness, or a common-unit residual bound below the sharp 98/100 threshold, would instead prove that the complete ideal block changes sign and retire this entire sign route.

## Pi provenance

The only pi-dependent phase is inherited from the completed-zeta saddle lock in the source artifact. The pi/3 safe arc and 2pi turn are fixed fractions of the period of exp(i theta); they introduce no geometric model constant. The 49/100 reserve is a rational consequence of the physical N floor and the exact phase/mismatch error bounds.

## Proof boundary

This gate proves a two-scale phase theorem and two strict opposite-sign witnesses for the bare Hermitian/transpose centre at the actual physical source phase on every complete q=1 cell. It retires a cell-uniform raw-centre sign, not the complete current. It proves no signed reciprocal residual estimate, correction/affine/moved-tail compensation, tie-spliced complete current sign, all-q transport, contact-box exclusion, Q209 shell, cofinal successor, Lambda<=0, PF-infinity, RH, or prize-level conclusion.

built Newman C1 actual-source raw-centre sign-obstruction gate: 13 rows, 0 issues, 1 safe phase arc, 1 scale-separation theorem, 2 strict sign witnesses, 1 actual-source centre obstruction, 1 macroscopic compensation condition, 1 live complete-current target
