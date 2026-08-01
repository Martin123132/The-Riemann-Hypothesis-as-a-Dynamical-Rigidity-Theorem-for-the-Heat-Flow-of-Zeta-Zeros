# Order-Eleven Lambda-Zero Completion Certificate

Date: 2026-07-22

Status: **rigorous all-shift signed Hankel order eleven and fixed-order arbitrary-column completion at lambda zero**. This is not a proof of RH or `Lambda <= 0`.

## Contiguous Layer

`Q_(11,n)(-100)>0 for every integer n>=0`
`Q_(11,n)(lambda)>0 for every n>=4 and -100<=lambda<=0`
`Q_(11,n)(0)>0 for every integer 0<=n<=3`

Together these prove `Q_(11,n)(0)>0 for every integer n>=0` and
`Q_(m,n)(0)>0 for every integer 1<=m<=11 and n>=0`.

## Arbitrary Columns

`epsilon_k*R_(k,n)(j_1,...,j_k)(0)>0 for 1<=k<=11, n>=0, and 0<=j_1<...<j_k`.

## Boundary

This proves contiguous signed Hankel positivity at lambda zero through order eleven and the corresponding consecutive-row arbitrary-column signs through order eleven. It proves no order above eleven and is not PF-infinity, Lambda<=0, or RH.

This is not a proof of RH.

## Reproduce

```powershell
python work/rh_compute/scripts/jensen_window_pf_compound_order11_lambda0_completion_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order11_lambda0_completion_certificate.py
```
