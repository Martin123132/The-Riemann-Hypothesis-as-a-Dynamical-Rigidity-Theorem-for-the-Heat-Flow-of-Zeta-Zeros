# Q208 Base and Q209 Successor-Shell Coverage Gate

Date: 2026-07-28

Status: exact finite-base correction and Q209 shell coverage
classification. This is not a proof of Q209, the uniform
successor theorem, Lambda<=0, or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_q208_base_q209_shell_coverage_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_q208_base_q209_shell_coverage_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_q208_base_q209_shell_coverage_gate.py
```

Current result:

```text
validated Q208-base/Q209 shell coverage gate: 15 rows, 5 coverage regions, 2 open shell regions, 7 oriented boundary arcs, 416 outer old half-cells, 2 new strip half-columns, Q209 not certified
```

## Corrected Finite Base

The cofinal target's Q207/j>=208 handoff is retained as a historical artifact. The later rigorous Q208 closed-boundary certificate supersedes only that finite base: Q208 is proved and the first unresolved linear stage is Q209.

The current theorem-level finite base is

```text
[1/1040,1/4]x[0,246]
Q208 is contact-free and has zero first-jet winding.
```

The first unresolved rectangle in the linear schedule is

```text
[1/1045,1/4]x[0,247]
```

## Exact Q209 Geometry

```text
t_j=1/(5*j), R_j=j+38
t_208=1/1040, t_209=1/1045
delta_208=1/217360
R_208=246, R_209=247
P_209=P_208 union C_208 union S_208
P_209=P_208 union C_208_compact union C_208_outer union S_208
area(P209)=2704/55 and the four rational piece areas sum exactly to it
```

The standard `(x,t)` boundary orientation is

```text
gamma_B(x)=(x,1/1045), 0<=x<=247
gamma_R(t)=(247,t), 1/1045<=t<=1/4
gamma_T(x)=(x,1/4), 247>=x>=0
gamma_L(t)=(0,t), 1/4>=t>=1/1045
twice signed area=13533/110>0
```

## Coverage Table

```text
q209_cov_01_inherited_p208: covered_exact; [1/1040,1/5]x[0,246]; none
q209_cov_02_compact_collar: covered_exact; [1/1045,1/1040]x[0,38]; none
q209_cov_03_outer_collar: open_xi_antecedent; [1/1045,1/1040]x[38,246]; For every one of the 416 old outer half-unit cells I_(208,k), prove delta_208*M_(208,k)<d_(208,k) using enclosures valid on [1/1045,1/1040]xI_(208,k).
q209_cov_04_new_right_strip: open_xi_antecedent; [1/1045,1/5]x[246,247]; Construct one fixed q_208!=0 and eta_208>0 with q_208 dot V>=eta_208 throughout [1/1045,1/5]x[246,247].
q209_cov_05_high_time: covered_exact; (1/5,1/4]x[0,247]; none for t>1/5
```

Exactly two low-time shell regions remain open: the outer
old descendant collar and the new right strip.

## Oriented Boundary Audit

```text
q209_edge_01_bottom_compact: bottom; 0<=x<=38; covered_exact; compact transversality theorem
q209_edge_02_bottom_outer: bottom; 38<=x<=246; open_xi_antecedent; outer-collar transport is unproved
q209_edge_03_bottom_new: bottom; 246<=x<=247; open_xi_antecedent; new right-strip half-plane is unproved
q209_edge_04_right_low: right; 1/1045<=t<=1/5; open_xi_antecedent; new right-strip half-plane is unproved
q209_edge_05_right_high: right; 1/5<t<=1/4; covered_exact; published Lambda<=1/5 composition
q209_edge_06_top: top; 247>=x>=0; covered_exact; t=1/4 is strictly above 1/5
q209_edge_07_left: left; 1/4>=t>=1/1045; covered_exact; high-time theorem plus compact origin theorem
```

The three open boundary rows are traces of the same two open
two-dimensional shell regions; they are not three independent
analytic obligations.

## Conditional Q209 Theorem

```text
For every one of the 416 old outer half-unit cells I_(208,k), prove delta_208*M_(208,k)<d_(208,k) using enclosures valid on [1/1045,1/1040]xI_(208,k).
Construct one fixed q_208!=0 and eta_208>0 with q_208 dot V>=eta_208 throughout [1/1045,1/5]x[246,247].
```

Those two antecedents, the compact-core theorem, and the proved Q208 base imply that P209 is contact-free with zero first-jet degree. The published Lambda<=1/5 theorem then covers 1/5<t<=1/4; equality on the new x interval must come from the new-strip antecedent. Only then is Q209 certified.

Neither antecedent is currently proved.

## Noninheritance Guards

The Q207-to-Q208 derivative enclosures were proved only on [1/1040,1/1035]. The smaller scalar delta_208 does not extend them to the disjoint lower-time interval [1/1045,1/1040]. New interval enclosures are required.

The proved right strip is [245,246]. Oscillatory Xi phase has no established monotone x-translation theorem, so its half-plane certificate cannot be shifted to [246,247].

A finite Q209 computation may calibrate or falsify a candidate estimate, but it is not the parameter-uniform successor theorem.

## Proof-Level Route

The linear Q209 instance is useful only as a bounded calibration.
The preferred cofinal geometry is

```text
c=25, a=100, t_j=25/(100+j), L_j=101+j, R_j=4*pi*exp(L_j)
t_(j+1)*L(x)>=t_(j+1)*L_j=25 on every new strip
closed by the existing dominant-saddle tL>=25 theorem
C_j^out=[t_(j+1),t_j]x[38,R_j]
s_pf=(L^2+(2t)^(-1))^(-1/2)=sqrt(2t)/sqrt(1+2tL^2)
```

Construct an a priori source-level kappa_j(t) such that ||partial_t(H,s_pf H_x)||<=kappa_j(t)||(H,s_pf H_x)|| on every outer descendant collar, with finite one-step integral. Defining kappa as the quotient of these two unknown norms is forbidden because it assumes nonvanishing.

Pi provenance:

Here pi is the usual circle constant already present in the established Xi frequency coordinate L(x)=log(x/(4*pi)); R_j=4*pi*exp(L_j) is exactly the inverse coordinate change. No polygon, curvature image, or arbitrarily chosen circle introduces this pi.

## Boundary

This artifact corrects the finite handoff, proves the exact Q209 geometry, orientation, area decomposition, and coverage classification, and composes only already certified source domains. It does not certify either open Q209 region, compute a Q209 winding, prove a uniform successor theorem, prove Lambda<=0, prove PF-infinity, prove RH, or establish any Clay-prize conclusion.
