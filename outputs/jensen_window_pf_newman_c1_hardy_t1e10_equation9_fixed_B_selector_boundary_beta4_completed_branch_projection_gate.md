# Beta^-4 completed selector strip and branch projection

Date: 2026-08-13

Status: rigorous beta-minus-four all-mode completion and selector-height error
cell; not an ordinary-Morse splice or a proof of `T_upper` or RH

For `A=159577`, `B=5122423`, `beta^3=pi*A^2/8`, and
`lambda=(t*-t)/beta`, the beta-minus-four Airy reduction has endpoint weights

```text
U(lambda,0)=1-13*lambda/(60*beta^2)
 -[448*lambda^5+4565*lambda^2]/(50400*beta^4),
V(lambda,0)=2*lambda^2/(15*beta^2)
 +[40*lambda^3+27]/(1680*beta^4).                  (BP1)
```

The shared-boundary coefficients are the Fourier coefficients of the smooth
profile `g4_lambda`.  Dirichlet-Jordan and the matching endpoint half-current
therefore give the exact completed identity

```text
H4 + sum_(m=39696)^40094 G4_m + C4 = g4_lambda(0),   (BP2)
```

where `C4` is the zero, negative, and outer-positive complement retained as
one object.  No outer coefficient is bounded or discarded separately.

At `lambda=0`, (BP1) gives

```text
J4-J0=-(9/560)*2*pi*Ai'(0)/beta^4
 =[1.213101703019356253530112027628027692795993392788486563216232186640908870652503602640682472794949393e-15 +/- 2.19e-104],
J_exact-J4=[-6.537196550387106204006970496750689646743880172587764643728299595202431491963949412328675057069208644e-22 +/- 2.94e-29].       (BP3)
```

Thus the beta-minus-four term cancels the previously certified exact-minus-
leading discrepancy by a factor greater than
`1855689.951459832838736474514007568359375000000000000000000000000000000000000000000000000000000000000`.  The residual is not
set to zero; it remains an explicitly enclosed finite-height error.

The same subtraction is rigorous on the complete selector cell.  The exact
and beta-minus-four center derivatives differ by
`[-2.479645130705650201894239978980217731353149294528850276077826542470615189927514503257221662043936813e-22 +/- 1.71e-28]`.  Using the old exact-minus-
leading second-derivative majorant together with an interval bound for the
beta-minus-four correction proves

```text
t*-pi/16 <= t <= t*+pi/16,
|J_exact(lambda)-J4(lambda)|
 < [5.201569942677361309405315923722830475260056596675501465715788551493407241442412657681320304231516536e-15 +/- 1.96e-23],
equation-(9) normalized error
 < [4.822320847949362466807646352612142566409230056103210914252969710407027797151012741597178367417706632e-16 +/- 1.81e-24] < 4.9e-16. (BP4)
```

On the ordinary-side half `t*-pi/16<=t<=t*`, the positive-`X` branch basis is
valid away from its regular endpoint limit.  In the old-strip coordinate
`x=beta*y/A`, the 399 block splits exactly into selected and opposite pieces,

```text
B_sel=e^(-i*799*x)C_+
 +e^(-i*x) sin(398*x)/sin(2*x)
   [e^(i*398*x)C_-+e^(-i*398*x)C_+],                 (BP5)

sum_block G4 = Integral exp(i*y^2/(4*beta))
                    (B_sel+B_opp)dy.                 (BP6)
```

The admissible completed remainder is consequently

```text
R_comp=H4+C4+Integral exp(i*y^2/(4*beta))B_opp dy
      =g4_lambda(0)-Integral exp(i*y^2/(4*beta))B_sel dy. (BP7)
```

Equation (BP7), rather than `B_opp` alone, is the cancellation-preserving
object.  It inserts the reflection pairing into the completed Poisson theorem
without reallocating the complement or endpoint half-current.  The adjacent
`A+2` event chart still has different detunings and is not identified here.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_beta4_completed_branch_projection_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_beta4_completed_branch_projection_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_beta4_completed_branch_projection_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_beta4_completed_branch_projection_gate.py
```

No selected-branch/logistic amplitude match, all-corridor estimate, complete
`Q_K-T` or `T_upper`, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion
is proved.
