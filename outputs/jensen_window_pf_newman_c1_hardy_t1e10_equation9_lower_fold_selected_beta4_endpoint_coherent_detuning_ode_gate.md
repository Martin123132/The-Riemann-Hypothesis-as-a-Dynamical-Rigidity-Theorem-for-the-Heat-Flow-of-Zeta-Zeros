# Endpoint-coherent detuning ODE for the selected beta-minus-four branch

Date: 2026-08-13

Status: exact finite-transform ODE, endpoint-source compression, and grouped
forcing structure proved; no quantitative finite-integral splice theorem

Let `X=lambda+y`, `epsilon=beta^-2`, and

```text
C_sigma=(U W_sigma+V W_sigma,X)/2,
W_sigma,XX=-X W_sigma.
```

For the finite selected-branch transform

```text
G_sigma(d)=Integral_0^Y exp(i[y^2/(4beta)-d y])
                         C_sigma(lambda,y)dy,          (ED1)
```

two integrations by parts give the exact detuning equation

```text
G_sigma''/(4beta^2)+i(1+d/beta)G_sigma'
 +(lambda-d^2+i/(2beta))G_sigma
 =R_0,sigma+R_Y,sigma+Q_sigma,                        (ED2)

R_0,sigma=C_sigma,y(lambda,0)+i d C_sigma(lambda,0), (ED3)

R_Y,sigma=e^(i[Y^2/(4beta)-dY])
 [i(Y/(2beta)-d)C_sigma(lambda,Y)-C_sigma,y(lambda,Y)]. (ED4)
```

The beta-minus-four correction does not leave an opaque bulk error.  Exact
Airy reduction gives

```text
Q_sigma=Integral_0^Y e^(i[y^2/(4beta)-dy])
                    [E_0 W_sigma+E_1 W_sigma,X]/2 dy,

E_0=-epsilon*y^2/4
    +epsilon^2(13lambda*y^2/240+y^3/80),
E_1=-epsilon^2*y^2(8lambda^2-4lambda*y+3y^2)/240.    (ED5)
```

Both forcing channels vanish quadratically at the lower endpoint.  On the
whole selector top corridor and `0<=y<=Y`, interval arithmetic proves

```text
|E_0| <= [0.000182274586145251347307784557253423893603914354877482548889983925586740222977950958319805529 +/- 5.47e-4] <7.31e-4,
|E_1| <= [6.64480633892670228382832872955478594264778785475819927630994417508088765688591400024838378e-9 +/- 9.97e-8] <1.08e-7. (ED6)
```

Now remove event zero and put `h=4beta/C`, `d_0=-h/4`.  The remaining roster
is `d_0-hj`, `1<=j<=198`, and `d_0+hj`, `1<=j<=200`.  If
`C_-=a+ib` and `C_-,y=a'+ib'`, direct finite summation gives

```text
sum R_0=398a'+(79601/2)h b
        +i[-2b'+(599/2)h a],                          (ED7)

sum e^(-i common phase)R_Y
 =-398a'-(79599/2)h b
   +i[2b'-(201/2)h a].                               (ED8)
```

Here (ED7) uses the values at `y=0` and (ED8) the values at `y=Y`.  The
upper phase is genuinely common: `hY=2pi`, `C=1 mod 4`, and therefore
`exp(-i d_m Y)=i` for every selected mode.  The interval source aggregates
are

```text
lower: [712.344072818756103515625000000000000000000000000000000000000000000000000000000000000000000 +/- 0.0458],
upper: [111.564610242843627929687500000000000000000000000000000000000000000000000000000000000000000 +/- 0.267]. (ED9)
```

They are coherent boundary channels, not 398 unrelated errors.  The interior
forcing also retains the exact event-ordered projection (11.385.4), with
`C_sigma` replaced by `(E_0 W_sigma+E_1 W_sigma,X)/2`.  Thus no termwise
absolute value is required anywhere in (ED2)--(ED8).

This explains structurally why a bare full-saddle comparison can fail even
when individual interior modes look accurate: it omits coherent finite
endpoint currents.  Equations (ED2)--(ED8) do not yet prove that those
currents cancel the completed Poisson complement, nor do they bound the
variation-of-constants kernel.

Pi provenance: `Y=pi C/(2beta)` and `beta^3=pi C^2/8` come from the exact
Kummer/Fourier selector geometry.  They force `hY=2pi` and
`Y/(2beta)=h/2`; no independent circle or fitted period is inserted.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selected_beta4_endpoint_coherent_detuning_ode_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selected_beta4_endpoint_coherent_detuning_ode_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selected_beta4_endpoint_coherent_detuning_ode_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selected_beta4_endpoint_coherent_detuning_ode_gate.py
```

No completed endpoint cancellation, uniform finite-integral kernel bound,
grouped 398-mode splice, all-corridor continuation, complete `Q_K-T` or
`T_upper`, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
