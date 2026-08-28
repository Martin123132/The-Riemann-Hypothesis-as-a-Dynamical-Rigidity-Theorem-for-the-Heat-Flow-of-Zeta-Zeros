# Affine-weighted Mordell recursion and full second steps

Date: 2026-08-24

Status: exact-recursion and pointwise interval-certificate note; the general
affine-weighted Mordell step and two full-roster second-step rows are
certified, but recursive parameter boxes and the small-tau neighbourhood
branch remain open

## Affine-weighted current

Define

```text
W_n(p;a,tau)=sum_(k=0)^n (p+k)exp(2pi i[ak+tau k^2]). (AW1)
```

For `a=z+r`, `-1/2<=z<1/2`, and `tau>0`, put

```text
m=floor(2n tau),  w=z/(2tau),  sigma=-1/(4tau),
C=exp(i*pi/4-i*pi*z^2/(2tau))/sqrt(2tau).
```

Differentiating Kuznetsov's exact theta identity with
`D_a=(2pi i)^(-1)partial_a` gives

```text
W_n(p;a,tau)
 =C/(2tau) W_m(2tau p-z;w,sigma)+E_n(p;a,tau),       (AW2)
```

where the endpoint-complete affine current is

```text
E_n=-i/4{
 E_-[(2p-1)h(u_-,-2tau)+h_z(u_-,-2tau)/(pi i)]
 +(-1)^m E_+[(2p+2n+1)h(u_+,-2tau)
                    +h_z(u_+,-2tau)/(pi i)]}.       (AW3)
```

Thus the child affine weight is exactly `p'=2tau p-z`.  The special physical
identity in Section 11.463 is the case `p=a/(2tau)`, for which `p'=r`.

## Exact normalization

For integer `j,l`, the exact symmetries used before each recursive step are

```text
W(p;a+j,tau+l)=W(p;a,tau),
W(p;a+1/2,tau-1/2)=W(p;a,tau),
W(p;a,-tau)=conjugate(W(p;-a,tau)).                  (AW4)
```

The second line follows from `(k-k^2)/2` being integral.  These operations
normalize every point row to `-1/2<=a<1/2`, `0<=tau<=1/4` without changing
`p`.

## Certified rows

Four short rational rows enclose their independent exact-period currents.
For the full `x=2/5,s=1/3` row, the first child

```text
W_992568(31916;-7/6,-5/4)
```

normalizes by conjugation to `W_992568(31916;1/6,1/4)`.  Equation (AW2)
then contracts the active index from `992568` to
`496284` and produces

```text
p'=95747/6,  w'=1/3,  sigma'=-1.                     (AW5)
```

The child is period three.  Its interval second step overlaps the original
period-six transformed current, and after restoring the first-step multiplier
and joined endpoints the complete source difference contains zero.

At `x=1/2,s=0`, the first child normalizes exactly to

```text
W_1240710(39894;-1/2,0),                             (AW6)
```

an alternating affine sum with period two.  The recursion therefore
terminates exactly at the centre; it does not divide by `tau=0`.  A nonzero
box around this centre contains very small positive `tau`, so a rigorous
small-tau neighbourhood branch is still required.

The largest full-row complete/source difference-ball radius is
`5.792175785049672968251976117642415829888e-47`.  Both differences
contain zero.  An altered 384-bit, `Y=11` replay independently reproduces all
short and full rows.

## Decision

The affine weight is now closed under the exact Mordell step, and the first
interior recursion contracts by one half before reaching a period-three
child.  The corner centre has an exact terminal branch.  What remains is to
lift (AW2)--(AW4) from points to branch-partitioned parameter boxes and add a
uniform small-tau enclosure around the terminal face.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_mordell_affine_weighted_second_step_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_mordell_affine_weighted_second_step_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_mordell_affine_weighted_second_step_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_mordell_affine_weighted_second_step_gate.py
```

Primary-source boundary: Kuznetsov's exact theta/Mordell identity supplies
the unweighted transform and endpoint functions.  Equations (AW2)--(AW4),
the affine endpoint coefficients, and their recursion are derived here.
Practical Gauss--Laguerre quadrature is not used.

Pi provenance: every `pi` comes from the inherited quadratic Fourier phase,
the differential normalization `D_a`, or the exact Mordell representation.
No fitted circle or polygon constant is introduced.

This gate proves the affine-weighted exact recurrence, four short interval
rows, one contracting full second step, and one exact tau-zero terminal full
row.  It does not prove recursive parameter boxes, a small-tau neighbourhood
branch, physical quadrature, the non-A bound, joined `R_after_A`, `R_Dir`,
`Q_K-T`, an all-height theorem, `Lambda<=0`, PF-infinity, RH, or a prize-level
conclusion.
