# Affine A-corner carrier and projector orientation

Date: 2026-08-13

Status: exact phase/projector dictionary and rigorous two-mode tangent-wedge
projection; not an exact A-corner or global residual bound

The stationary phase of the endpoint-`D` triangle is

```text
Phi_D,m=pi*m*D-pi*m^2-t/2+(t/2)log[t/(2pi*m^2)].       (PO1)
```

For an odd endpoint `D=2d+1`,

```text
m(D-m)=2md+m(1-m) is even.                              (PO2)
```

Consequently both odd endpoints have the same exact carrier,

```text
exp(i Phi_D,m)
 =exp(i{t[log(t/(2pi))-1]/2-t log m}),

exp(-i*pi/8)exp(i Phi_D,m)
 =exp(i[theta_0-t log m]),
theta_0=t[log(t/(2pi))-1]/2-pi/8.                      (PO3)
```

Thus the affine A endpoint term and its equation-(9) projection are

```text
K_A,m^aff=-exp(i Phi_A,m)W0_A,m C_aff,A,m,
q_A,m^aff=2(pi/(32t))^(1/4)
 Re[exp(-i*pi/8)K_A,m^aff].                            (PO4)
```

The minus sign in (PO4) is `epsilon_A=-1`; it is not supplied by the target
projector.  Indeed, with `chi_T(m)` the positive target indicator,

```text
I_m+I_-m-chi_T(m)P_m
 =(1-chi_T)P_m+P_-m+Q_m+Q_-m.                         (PO5)
```

The paired endpoint term `Q_m+Q_-m` therefore has coefficient `+1` on both
sides of `39894|39895`.  Only the positive full-line bulk coefficient jumps
from zero at mode 39894 to one at mode 39895.

Restoring the common carrier and the paper normalizer gives

```text
mode 39894: q_A^aff=[-0.0058473814010481522425368130416436205985365837218438414497241135099759924844233798 +/- 4.95e-34],
mode 39895: q_A^aff=[-0.0018925895733072011429497060880330715371198896341834248137864234348109381020191090 +/- 8.70e-35],

two-mode affine sum=[-0.0077399709743553533854865191296766921356564733560272662635105369447869305864424888 +/- 5.82e-34],
two-mode scalar sum=[-0.0077336609886276787660889777139776826796924642372691463368894712536065087117744125 +/- 5.82e-34],
affine-minus-scalar=[-6.3099857276746193975414156990094559640091187581199266210656911804218746680762774e-6 +/- 1.17e-33]. (PO6)
```

The two canonical affine A terms reinforce rather than cancel, and the
physical affine correction alone is more negative than `6.3e-6`.  That is
already comparable with the final `8.6e-6` target scale, so it cannot be
dropped.  The roughly `-0.00774` two-mode value is not an error estimate: it
must remain grouped with the rest of the A roster, the exact curved-face and
amplitude remainders, the full-line projector jump, and the other endpoint
sectors.

Pi provenance: `pi` in (PO1)--(PO6) comes from the exact equation-(9)
endpoint phase, Fourier parity, odd-square reflection, Gaussian
normalization, and paper prefactor.  No fitted constant is introduced.

Proof boundary: exact phase reduction, paired-projector orientation, and
rigorous normalized tangent-wedge values at modes 39894 and 39895 only.  No
exact A-corner remainder, complete A roster, `R_Dir` estimate, complete
`Q_K-T` or `T_upper`, all-height theorem, `Lambda<=0`, PF-infinity, RH, or
prize-level conclusion is proved.
