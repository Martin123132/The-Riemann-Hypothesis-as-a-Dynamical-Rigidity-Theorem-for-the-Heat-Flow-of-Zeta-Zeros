# Global extraction of the certified B trace window

Date: 2026-08-13

Status: exact residual decomposition and quantitative sufficient target; not
a proof of the grouped-remainder bound or RH; the grouped remainder is open

At a finite symmetric cutoff, the bulk-extracted fixed-state integrand is

```text
G_M(x)=H_x+integral_0^L f_x(u)C_M(u)du
       +sum_(m=622)^39894(P_A,m+P_B,m).               (GE1)
```

Put

```text
tau_m=1_(622<=m<=39894),
sigma_m=(1-tau_m)-1_(q_B,m<0),
ell_m=1_(m in {621,622}).
```

Since `P_B,m=U_B,m-1_(q_B,m<0)P_bulk,m`, exact modewise
algebra fixes the allocation

```text
I_m+I_-m-tau_m P_bulk,m
 =[U_B,m+P_B,-m+ell_m sigma_m P_bulk,m]
  +[P_A,m+I_-m-P_B,-m
    +(1-ell_m)sigma_m P_bulk,m].                      (GE2)
```

The first bracket is `Btr_M` after summation.  It is exactly the trace-current
package certified on the window: the local 621/622 steps remain attached to
their crossing tails, while every nonlocal constant step remains explicitly
in the second bracket.  Put

```text
chi_W=1_[x_low,x_high],
xi=sqrt(H_B)(x-x_B),      W={|xi|<=70}.               (GE3)
```

Linearity at every finite cutoff gives the exact extraction

```text
G_M=chi_W Btr_M+[G_M-chi_W Btr_M].                   (GE4)
```

The established symmetric Abel limit and the complete B trace-window theorem
therefore give

```text
R_end=E_Btr,win+R_group,
E_Btr,win<-1.3198e-4.                                (GE5)
```

Here `R_group` is not merely two outside B tails.  It retains every nonlocal
`sigma_m P_bulk,m` step omitted from the trace package, the A endpoint
current, both B outside-window trace pieces, the zero/negative/outer-positive
completion, and the endpoint half-current in the exact global allocation.
In particular, the left B tail cannot be normed separately: the Abel-theta
corner majorant uses the endpoint difference

```text
Delta_g(z)=g_B(z)-g_A(z)=O(z^(1/2)),                 (GE6)
```

whereas separating the endpoints loses that cancellation.

Against the working physical target, the next sufficient theorem is now

```text
R_group<1.4058e-4  =>  R_end<8.6e-6.                (GE7)
```

The stronger target `R_group<1.3198e-4` would prove `R_end<0`.  GE7 is a
target specification, not an estimate of `R_group`; no part of the grouped
remainder has been assigned an independent triangle budget here.

Pi provenance: `pi` comes from the equation-(9) phase, Fourier-Poisson
character, exact half-Kummer projection, and Jacobi transform.  No fitted
or geometric surrogate is introduced.

Proof boundary: exact extraction of the certified saved-height B trace package
from the global fixed-state residual and the resulting sufficient numerical
target only.  No grouped-remainder estimate, complete `Q_K-T` or `T_upper`,
height-uniform theorem, `Lambda<=0`, PF-infinity, RH, or prize-level
conclusion is proved.
