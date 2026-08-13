# Exact completed point versus weighted profile

Date: 2026-08-13

Status: exact logical separation and conditional whole-kernel transfer;
the physical exact-profile identity and uniform profile error remain open

Let `delta J_lambda(s)` denote the reduced exact-minus-beta-minus-four fold
integral before the all-mode evaluation, normalized so that the saved scalar
theorem is `delta J_lambda(0)`.  Since `hY=2pi`, the corresponding canonical
Fourier profile is

```text
delta g_lambda(s)=[Y/(2pi)]delta J_lambda(s)
                  =delta J_lambda(s)/h.                (EP0)
```

Use the selector-centred Fourier convention

```text
Delta G_(39895+n)(lambda)
 =Integral_0^1 delta g_lambda(s)e^(-2pi i n s)ds.       (EP1)
```

Here (EP1) is the common-profile representation that a finite-`t` remainder
theorem would have to establish; this gate does **not** assume that the still
missing physical amplitude/Fresnel remainder already has this form.  The
completed-strip theorems evaluate the all-mode object at `s=0`.  In the raw
reduced-fold and canonical Fourier normalizations they prove only

```text
|delta J_lambda(0)|
 <9.49612500096965347063853751189785591899931291950937710064324691677659606545530133694636546E-12<9.497e-12,

|delta g_lambda(0)|
 <1.75842251125259831742031899764141519403286918141964047104483729464309891537402844719769180E-10<1.759e-10. (EP2)
```

Equation (EP2) does not control the weighted Fourier correction.  For an
arbitrary complex `z`, define the endpoint-preserving counterprofile

```text
phi_z(s)=z[e^(2pi i 198s)-e^(2pi i 199s)].             (EP3)
```

Both `phi_z(0)` and `phi_z(1)` vanish exactly.  Thus adding (EP3) leaves the
Fourier midpoint, endpoint half-current, their completed value, and the
event-zero coefficient unchanged.  Its only nonzero Fourier coefficients
in the selector roster are

```text
Delta G_40093=z,             Delta G_40094=-z.          (EP4)
```

The certified weights satisfy

```text
|w_40093-w_40094|
 =[1.917549273811047783055983018130064010620117187500000000000000000000000000000000000000000000000000000e-7 +/- 4.06e-11]
 >1.91714327381104778305598301813006401062011718750000000000000000000000000000000000000000000E-7>1.917e-7.              (EP5)
```

Consequently (EP3) changes the weighted correction by
`z(w_40093-w_40094)` while preserving all data used by the completed-point
estimate.  This is a logical insufficiency result for the existing scalar
theorem, not a claim that the physical exact profile can be varied freely.

There is a cancellation-preserving sufficient replacement.  If the missing
physical correction is first proved to satisfy (EP1), define the whole
weighted kernel

```text
P_w(s)=sum_(n=-199)^199 w_(39895+n)e^(-2pi i n s),
w_39894=0.                                               (EP6)
```

The two independently enclosed pieces of this same kernel obey

```text
Integral_0^1|P_w(s)|ds
 <=0.000038402825923711157469357553976338890036124815208824976467044011286419031499973287470<3.840283e-5.          (EP7)
```

Therefore one uniform common-profile theorem

```text
sup_(lambda,s)|delta g_lambda(s)|<=epsilon_profile       (EP8)
```

would imply directly, with no coefficientwise estimate,

```text
|sum_m w_m Delta G_m|
 <=Integral_0^1|delta g_lambda(s)P_w(s)|ds
 <3.840283e-5 epsilon_profile.                           (EP9)
```

Equivalently, a raw-fold profile estimate
`sup|delta J_lambda|<=epsilon_J` would give the same weighted correction
below

```text
0.000711115255886724357667503555655620273299661419993347461371735065747527115808614214504886573 epsilon_J
 <7.112e-4 epsilon_J.                                   (EP10)
```

The next obligation is not 398 separate mode estimates.  It is to derive
the exact finite-`t` common profile before the all-mode collapse, prove that
the missing weighted physical remainder is its Fourier pairing with (EP6),
and bound that profile uniformly on the ordinary top corridor.  A completed
cutoff-family estimate is an equivalent admissible route.

Pi provenance: the `2pi` in (EP1), (EP3), and (EP6) is forced by the exact
selector relation `hY=2pi` and the integer Fourier--Poisson character.  This
gate introduces no fitted period, geometric circle, or polygon constant.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_completed_point_weighted_profile_obligation_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_completed_point_weighted_profile_obligation_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_completed_point_weighted_profile_obligation_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_completed_point_weighted_profile_obligation_gate.py
```

No exact physical common-profile identity, uniform exact-minus-beta-minus-four
profile bound, weighted finite-`t` amplitude/Fresnel remainder, completed
source/initial-data splice, complete `Q_K-T` or `T_upper`, all-corridor or
height-uniform theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion
is proved.
