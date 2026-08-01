# Newman Q208 Selected-Boundary Pilot

Date: 2026-07-25

Status: rigorous complete right-strip theorem and selected-bottom
diagnostic, not a proof of Q208, `Lambda<=0`, or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_q208_selected_boundary_pilot.json
python work/rh_compute/scripts/jensen_window_pf_newman_q208_selected_boundary_pilot.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_q208_selected_boundary_pilot.py
```

Current result:

```text
validated Q208 selected-boundary pilot: 8 panels, 46 certified leaves, 0 unresolved, 0 global promotions
```

## Domain

The prior theorem closes

```text
Q_207=[1/1035,1/4]x[0,245].
```

The next linear-exhaustion rectangle is

```text
Q_208=[1/1040,1/4]x[0,246].
```

This pilot samples six half-unit panels across the new bottom
microstrip and both half-unit panels across the new right strip:

```text
bottom time strip: [1/1040,1/1035]
right time strip: [1/1040,1/5]
```

It uses the same 352-bit six-term retained transform, directed
Taylor enclosure, and raw omitted-tail bars as the Q207 theorem.
The raw tail bars remain valid at x=246; no gamma-scale comparison
outside its recorded x<=245 theorem domain is promoted.

## Results

| region | x panel | leaves | value branch | derivative branch | minimum ratio |
|:---|:---|---:|---:|---:|:---|
| low_time_left_strip | [38,77/2] | 1 | 1 | 0 | [414262127671158275821122758051282934226411808540862986017007.461228110532757 +/- 8.44e-17] |
| low_time_left_strip | [137/2,69] | 1 | 1 | 0 | [9767725848992189759352612460270108076119185140404703141.16953146923147950300 +/- 4.12e-21] |
| low_time_left_strip | [199/2,100] | 1 | 0 | 1 | [25318936010231177921008518033650957413518277093606.7352744601831553516098602 +/- 1.78e-26] |
| low_time_left_strip | [299/2,150] | 1 | 0 | 1 | [438171792602138617842916945509958841388939.104706479036433189249311560762718 +/- 4.21e-34] |
| low_time_left_strip | [399/2,200] | 1 | 1 | 0 | [3242412588292633904336114969489200.85225656135869464445046125441700232943627 +/- 4.05e-42] |
| low_time_left_strip | [489/2,245] | 1 | 0 | 1 | [24806486883006562068195320.7897359078477768495015264174607557008043674472940 +/- 2.69e-51] |
| full_time_right_strip | [245,491/2] | 20 | 0 | 20 | [24315642501122114739555076.4042810360224373772504245579408469289728661881210 +/- 1.96e-50] |
| full_time_right_strip | [491/2,246] | 20 | 0 | 20 | [16775977756171684820231611.2572779119257685637018233703176116638156341196057 +/- 4.76e-50] |

Summary:

```text
branch counts={'value': 3, 'derivative': 43}
full-value sign counts={'negative': 24, 'positive': 1, 'not_sign_separated': 21}
minimum certified ratio=[16775977756171684820231611.2572779119257685637018233703176116638156341196057 +/- 4.75e-50]
minimum location={'t_low': '1/1040', 't_high': '57/5200', 'x_low': '491/2', 'x_high': '246', 'branch': 'derivative'}
right-strip theorem={'domain': '[1/1040,1/5]x[245,246]', 'certified_cells': 40, 'full_derivative_positive': True, 'minimum_full_derivative_lower': '[2.74842995422409717150882800350618905119382156826663560755859516271491876751e-29 +/- 9.46e-105]', 'conclusion': "J_t'(x)>0 throughout [1/1040,1/5]x[245,246]; hence the full first jet is nonzero there and its image lies in the open upper half-plane."}
baseline CPU samples=[23.1, 21.4, 14.3, 21.0, 13.7, 21.9]
runtime CPU samples=[8.3, 17.7, 23.2, 20.7, 24.9, 19.8, 27.6, 32.6]
```

## Interpretation

All selected panels are contact-separated with a large interval
margin. The two right panels form a complete cover, not a sample:

```text
J_t'(x)>0 throughout [1/1040,1/5]x[245,246]; hence the full first jet is nonzero there and its image lies in the open upper half-plane.
```

The right-edge first-jet phase is therefore confined to the open
upper half-plane and contributes no full turn. Across the selected
bottom samples both value and derivative branches are needed, so
a one-sign bottom-edge proof is not the expected continuation.
This supports a continuous first-jet phase or winding certificate.

The pilot does not cover the unsampled bottom panels, does not
compute a closed-boundary winding, and does not certify Q208.
The next proof-facing step is a rigorous boundary phase-unwrapping
algorithm, followed by a complete Q208 boundary run if that
algorithm closes.
