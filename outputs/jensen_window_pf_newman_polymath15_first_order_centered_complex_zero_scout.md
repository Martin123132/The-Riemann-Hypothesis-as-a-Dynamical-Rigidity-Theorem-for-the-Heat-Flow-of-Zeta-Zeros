# Newman First-Order Centered Complex-Zero Scout

Date: 2026-07-26

Status: selected high-precision route diagnostic; this is not an interval certificate
and is not a proof of RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_complex_zero_scout.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_complex_zero_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_complex_zero_scout.py
```

## Selected Point

| quantity | value |
|---|---:|
| x | 519.3859242675413194059971571412331940721769526032205179156613690290221 |
| t | 0.4018521828453945587120659361547194513935678848113752316604652083241877 |
| L | 3.721622951906665326301222553998987712660275250709871084007571745557969 |
| tL | 1.495542306951214314837719794638763945645461375641256609355554496641399 |
| q=2tL^2 | 11.13168915019416477029800738881244263972340424499351800797385267849907 |
| N | 6 |
| p | 0.1381907382835385429460253910721235855959766683609610143926652168171849 |
| Re(E_[1]) | -1.198812800858597535815615189309415059426601115545462756995952138837094e-85 |
| Im(E_[1]) | -4.755323085809729657479879313026051931526110347583474719701934783779065e-85 |
| |E_[1]| | 4.904105400777740719377975962574449307862526489542584846732299853413914e-85 |
| U=Re(E_[1],x) | 0.5561369634808909530620777743972843641869680841933197849058953865185877 |
| V=Im(E_[1],x) | -0.08143499934642560982074511513290129550450580496562479235547682563394576 |
| det D_(x,t)(Re E_[1],Im E_[1]) | 0.3954233157000086639951462017528374382208742150389412209998142939375605 |
| W_[1] | 2.742236260967209843483934071737898452867312320339222853105881406904536e-85 |
| S_a | 0.5561369634808909530620777743972843641869680841933197849058953865185877 |
| A_a | 0.5561367537278263285896083739449502955574373214320552068951007530889460 |
| Re(D_(1,x)) | 0.0000002097530646244724694004523340686295307627612645780107946334296417240085 |

The corrected complex main is numerically zero while the real
crossing slope is about `0.556`. Consequently the Wronskian is
zero at the limiting point, but the centered scalar still
classifies the crossing as upward.

## Precision Replay

```text
|Delta x|<1e-50
|Delta t|<1e-50
|Delta U|<1e-50
```

The low- and high-precision runs agree beyond fifty decimal
places. This is strong numerical stability, not an interval
existence proof.

## Route Consequence

This selected point makes the exceptional class concrete:
`E_[1]=0` does not imply a multiple real-part crossing. A proof
based only on the sign of `W_[1]Im(E_[1])` would delete it, while
`S_a=U-u_aX` and the centered core retain it. The optional next
finite task is a tiny interval Newton box around this point.

This scout supplies a reproducible selected moderate-height diagnostic and a cross-precision route guard. It does not prove existence of the complex zero by interval arithmetic, any L>=50 complex-main zero, the q>=1 centered scalar theorem, contact exclusion, Lambda<=0, RH, PF-infinity, or a Clay-prize conclusion.
