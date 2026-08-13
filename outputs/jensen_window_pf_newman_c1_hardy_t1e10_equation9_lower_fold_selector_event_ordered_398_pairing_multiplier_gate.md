# Event-ordered pairing after removing event zero

Date: 2026-08-13

Status: exact 398-term selector regrouping and leading-multiplier diagnostic;
not a uniform finite-integral splice

At the completed selector boundary, remove only the already transported
selected summand `m=39894`, whose detuning label is `-1`.  The remaining
`398` selected terms split exactly into
`198` consecutive event pairs

```text
(m_-,m_+)=(q-j,q+j),       q=39894, 1<=j<=198,          (EP1)
```

and the selector-exchange edge `(40093,40094)`.  For each event pair,

```text
C_- exp(i(4j+1)x)+C_+ exp(-i(4j-1)x)
 =2 exp(ix) Re[C_- exp(i4jx)].                        (EP2)
```

The two edge terms obey

```text
C_+[exp(-i795x)+exp(-i799x)]
 =2 C_+ exp(-i797x) cos(2x).                          (EP3)
```

Therefore the complete remainder of the selected projection is

```text
B_rem=2*exp(i*x)*sin(396x)/sin(2x)*Re(C_-*exp(i*398x))+2*C_+*exp(-i*797x)*cos(2x).                  (EP4)
```

The sine quotient in (EP4) always means its original finite 198-term sum.
This derivation removes no complement mode and does not split the endpoint
half-current from the completed remainder.

For the exact leading stationary multiplier, define

```text
a_m=[R_A(s_m,mu) exp(i beta^3 F(s_m,mu))-1]/sqrt(m).
```

Across the top corridor, the measured interval bounds are

```text
max |a_m| < 7.760623205399497237522155046463012695312500000000000000000000000000000000000000000000000000000000000e-6,
max |a_(q+j)-conj(a_(q-j))|
 < 3.197604518945240670291241258382797241210937500000000000000000000000000000000000000000000000000000000e-7,
sum pair leakage < 5.141749161404267829918787041654226754872070159763097763061523437500000000000000000000000000000000000e-6,
coherent Abel variation < 1.508412461872843252628541288462571401396417059004306793212890625000000000000000000000000000000000000e-5,
two selector-edge defects < 1.532949141846984275616705417633056640625000000000000000000000000000000000000000000000000000000000000e-5.
                                                               (EP5)
```

Factoring out the common classical carrier leaves

```text
sup_(t*-pi/16<=t<=t*)
 |sum_(m!=q) a_m(t) exp[-it log(m/q)]|
 < 1.625705479568750888574868440628051757812500000000000000000000000000000000000000000000000000000000000e-5 <1.64e-5. (EP6)
```

This is a rigorous 256-panel height enclosure.  The factor `pi/16` is the
exact selector-cell width inherited from the Kummer/Fourier event geometry.

These are leading-carrier coefficients only.  Their exact coherent/leakage
decomposition identifies the finite differences needed by Abel summation,
but (EP5) is not yet multiplied by a rigorously controlled finite-integral
kernel.  It therefore cannot be added to the one-mode transport bound as a
complete splice estimate.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selector_event_ordered_398_pairing_multiplier_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selector_event_ordered_398_pairing_multiplier_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selector_event_ordered_398_pairing_multiplier_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selector_event_ordered_398_pairing_multiplier_gate.py
```

No uniform finite-integral selected-branch theorem, grouped 398-mode
remainder bound, all-corridor continuation, complete `Q_K-T` or `T_upper`,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
