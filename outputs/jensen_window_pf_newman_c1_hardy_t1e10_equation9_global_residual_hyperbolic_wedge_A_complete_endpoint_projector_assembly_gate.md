# Complete exact A endpoint and projector transition

Date: 2026-08-14

Status: complete saved-height A endpoint carrier and adjacent positive
full-line projector jump certified; not a complete `R_Dir` bound

The endpoint coefficient does not jump at the target edge.  The exact finite
paired-projector identity is

```text
I_m+I_-m-chi_T(m)P_m
 =(1-chi_T(m))P_m+P_-m+Q_m+Q_-m.                    (AP1)
```

Therefore `Q_m+Q_-m` has coefficient `+1` throughout modes `39853..39936`,
while the positive full-line Gamma bulk has coefficient zero through mode
`39894` and one on modes `39895..39936`.

The complete exact endpoint carrier is assembled before norms:

```text
local exact affine       = [-0.126865295975501644562427416849169420579246224412138197112449585580000000000 +/- 2.67e-12],
local exact-minus-affine = [1.10622962615098585944923918500230782093161208480147221594949779180000000000e-8 +/- 2.20e-14],
full exact exterior      = [-0.00240153401848142095784295609580361814159723755519748799500000000000000000000 +/- 3.41e-13],
complete exact A endpoint= [-0.129266818931686804010411778452581188697765252651214837092727426085022082000 +/- 3.03e-12]. (AP2)
```

For the occupancy jump, the exact global-Morse Gamma factor `B(t)` gives

```text
G_m=2 Re[B(t) exp(i(theta_0-t log m))]/sqrt(m),
sum_(39895<=m<=39936) G_m
   =[0.0925926941339131572924456873536606766853598266718030293496073878641327886137 +/- 2.45e-77],
Gamma-minus-classical sum=[-2.89352169168478616540746435313925180778989081864917336171846434540295366433e-23 +/- 5.23e-80]. (AP3)
```

Combining (AP2) and (AP3) yields the projector-completed transition

```text
A_transition=[-0.0366741247977736467179660910989205120124054259794118077431200382208892933863 +/- 3.03e-12].        (AP4)
```

The endpoint and full-line terms remain separately recorded because they
have different projector ownership.  Their signed combination in (AP4) is
the admissible A-transition object for later `R_Dir` assembly.

Pi provenance: every `pi` in (AP1)--(AP4) comes from equation (9), the exact
bi-Morse/Fresnel maps, odd-endpoint parity, the Gamma/Gaussian bulk identity,
and the paper normalization.  No fitted constant is introduced.

Proof boundary: the exact A endpoint carrier and its adjacent 42-mode
positive full-line projector jump at `t=10^10` only.  The B endpoint, zero
mode, negative and remote-positive paired currents, complete `R_Dir`,
`Q_K-T`, `T_upper`, any all-height theorem, `Lambda<=0`, PF-infinity, RH,
and any prize-level conclusion remain unproved.
