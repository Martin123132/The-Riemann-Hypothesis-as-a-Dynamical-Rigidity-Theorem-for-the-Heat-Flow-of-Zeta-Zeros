# Uniform Gamma-target insertion in the global residual

Date: 2026-08-13

Status: exact common-carrier identification and corrected target arithmetic;
not a proof of the compressed endpoint-defect bound

For the finite Abel weight `w_m=exp(-pi*epsilon*m^2)`, write

```text
S_(M,epsilon)=H_x+integral_0^L f_x(u)D_(M,epsilon)(u)du,
I_(T,epsilon)=integral_0^L f_x(u)G_(T,epsilon)(u)du
             =sum_(m=622)^39894 w_m I_m.             (GI1)
```

The exact current decomposition is

```text
I_m=P_bulk,m+Q_m,       Q_m=P_A,m+P_B,m.             (GI2)
```

Substituting (GI2) before any norm gives

```text
H_x+integral f_x(D_(M,epsilon)-G_(T,epsilon))
   +sum_(m=622)^39894 w_m Q_m
 =S_(M,epsilon)-sum_(m=622)^39894 w_m P_bulk,m.      (GI3)
```

Thus the target kernel in Section 11.423 already inserts one exact full-line
Gamma bulk carrier for every target mode.  The ordinary/fold partition

```text
622..39694      39073 modes,
39695..39894      200 modes                           (GI4)
```

does not create a second bulk carrier.  Its distinction begins only in the
finite endpoint defect `Q_m`.  In particular the exceptional A-face term is
retained in `Q_m`; it is not an exception to (GI3).

After the common Abel limit and physical projection, refine the notation of
Section 11.423 by writing

```text
R_KGamma=Q_K-G=E_Btr,win+E_outer+R_Dir,              (GI5)
Q_K-T=R_KGamma+(G-T).                                (GI6)
```

Here the exact common mode factors are

```text
g_m=B(t) exp(i(theta_0-t log m))/sqrt(m),
tau_m=exp(i(theta-t log m))/sqrt(m),
g_m-tau_m=[B(t)-exp(i(theta-theta_0))]
             exp(i(theta_0-t log m))/sqrt(m).        (GI7)
```

The saved Arb certificate over all 39273 target modes gives

```text
|G-T|_physical <= [4.918151566720828980806353638226121738513797734778668279706287331790919038545165004173875668351233940e-20 +/- 3.54e-89]
                 < 5E-20.                (GI8)
```

Therefore the fully explicit sufficient targets are

```text
R_Dir < 0.00014057919999999995
        ==> Q_K-T < 8.6e-6,

R_Dir < 0.00013197919999999995
        ==> Q_K-T < 0,                               (GI9)
```

using the certified negative B window and the `8e-10` outer allowance.  The
`5e-20` correction is immaterial numerically but necessary for an exact
object identity.

The remaining theorem is now unambiguous: bound the single finite endpoint
defect `R_Dir`, with its A face, B complement, zero, negative, outer, and
half-current sectors still grouped.  No separate ordinary or fold bulk
identification remains.

Pi provenance: `pi` in (GI1)--(GI7) comes from the equation-(9) Gaussian
Abel weight, Fourier character, Kummer phase, and exact Riemann--Siegel phase.
No fitted or geometric constant is introduced.

Proof boundary: exact carrier insertion, object identification, and a
saved-height Gamma-to-classical correction bound only.  No quantitative
bound for `R_Dir`, complete endpoint residual theorem, all-corridor or
height-uniform theorem, `Lambda<=0`, PF-infinity, RH, or prize-level
conclusion is proved.
