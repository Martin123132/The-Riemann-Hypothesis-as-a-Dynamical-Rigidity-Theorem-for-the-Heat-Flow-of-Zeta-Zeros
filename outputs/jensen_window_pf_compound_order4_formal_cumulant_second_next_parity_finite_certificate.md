# Jensen-Window PF Order-Four Second-Next Finite Certificate

Date: 2026-07-13

Status: rigorous finite interval theorem for the second-next coefficient
functions. This is not a proof of the exact cumulant corridors,
order-four entry, PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_finite_certificate`.

```text
work/rh_compute/results/jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_finite_certificate.json
work/rh_compute/results/jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_finite_chunks.jsonl
python work/rh_compute/scripts/jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_finite_certificate.py
```

## Certified Bounds

A deterministic cover of `3600` adjacent
blocks of width `1/200` proves:

```text
0 < D_2(u) < 1/20
-1/5 < D_3(u) < 0
-3/5 < D_4(u) < 0
-1 < D_5(u) < 0
-2 < D_6(u) < 0
0 < D_7(u) < 9/2
0 < D_8(u) < 25/4
```

Here `D_r` multiplies `q^-4` for orders 2-4, `q^-3` for orders
5-6, and `q^-2` for orders 7-8 in the corridor-scaled
epsilon-nine/ten correction.

| order | certified lower | certified upper | lower margin | upper margin |
|---:|---:|---:|---:|---:|
| 2 | `3.10692696881464805701457372124565217385187169154301723800689E-3` | `2.98952434841370971932054401797338810990917840578813629979004E-2` | `3.10692696881464805701457372124565217385187169154301723800689E-3` | `2.01047565156355291313513277607873591352832159421186370020917E-2` |
| 3 | `-1.67454045550531767050475584301823585595023213457590978212372E-1` | `-1.40688879177776106150884495495281647330514243697303613098380E-1` | `3.25459544494611275221668146963177031623010052924090217875966E-2` | `1.40688879177776106150884495495281647330514243697303613098380E-1` |
| 4 | `-5.02388158281384408452018879272730730007609586573322142737867E-1` | `-4.11807654391921851035963951039726188608243589770064016099887E-1` | `9.76118417185019047102595046975298901095779134266778572620385E-2` | `4.11807654391921851035963951039726188608243589770064016099887E-1` |
| 5 | `-7.35669677437071819624286303460232610005131663930266436742903E-1` | `-6.33016430725545258378882183127855501385259381992735703343463E-2` | `2.64330322562928180375713696539767389994868336069733563257097E-1` | `6.33016430725545258378882183127855501385259381992735703343463E-2` |
| 6 | `-1.34273308005895631140779922344038974548259049272919338786251E+0` | `-1.11209349213355326349997742910554926518235267504107334108596E-1` | `6.57266919941043688592200776559610254517409507270806612137500E-1` | `1.11209349213355326349997742910554926518235267504107334108596E-1` |
| 7 | `3.60733401313211621996599522555305621094368712474878988519416E+0` | `4.47216979012616466138400085027805549968019499212905831883167E+0` | `3.60733401313211621996599522555305621094368712474878988519416E+0` | `2.78302098738353386159991497219445003198050078709416811683388E-2` |
| 8 | `4.83215345133096423098844336353492607344434814520915262414985E+0` | `6.17624150792470266594014955264867600862887657123098305867318E+0` | `4.83215345133096423098844336353492607344434814520915262414985E+0` | `7.37584920752973340598504473513239913711234287690169413268286E-2` |

The same Taylor cells also prove the normalized potential-jet caps

```text
0 < L_3(u) < 6/5
0 < L_4(u) < 3/2
0 < L_5(u) < 2
0 < L_6(u) < 3
0 < L_7(u) < 9/2
0 < L_8(u) < 7
0 < L_9(u) < 12
0 < L_10(u) < 21
0 < L_11(u) < 38
0 < L_12(u) < 71
```


## Taylor Gate

Each block uses a centered sixth-order Arb Taylor model with exact
rational-center derivatives and a full-block seventh-derivative
remainder. The normalized potential geometry is evaluated through `V^(12)`,
so all 30- and 42-monomial cancellations occur before range
enclosure. This is not a point sample or floating-point fit.

The weakest outward-rounded margin is `0.00310692696881464805701457372124565217385187169154301723800689`
for order `2` on the
`minimum_lower_margin` side.

## Remaining Boundary

This closes the explicit epsilon-nine/ten coefficient layer only on
`2<=u<=20`. Its asymptotic companion and the exact central and two-tail
density theorem remain open.

```text
outputs/jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_certificate.md
outputs/jensen_window_pf_compound_order4_exact_cumulant_remainder_budget.md
outputs/formal_core.md
```
