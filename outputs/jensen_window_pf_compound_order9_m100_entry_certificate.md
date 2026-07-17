# Order-nine lambda=-100 entry certificate

Date: 2026-07-14

Status: rigorous all-shift signed order-nine entry at `lambda=-100`; This is not a proof of PF-infinity, RH, or `Lambda <= 0`; it closes
one fixed compound order at one heat endpoint.

## Analytic tail

`w_1''(t)<=4200/t^2 for every real t>=1250`

feeds the tent and full-kernel bounds

`w_1''(t)<=4200/t^2 => Y_k^(1)<=4200*[-log(1-1/k^2)]<4201/k^2, k>=1251`

`|Y_k-Y_k^(1)|<=F_(k-1)+2*F_k+F_(k+1)<550/k^2, k>=1251`

`w_1''(t)<=4200/t^2 on t>=1250 => Y_k<4201/k^2+550/k^2=4751/k^2<4900/k^2, k>=1251`

and therefore proves

`Q_(9,n)(-100)>0 for every n>=1243`.

## Exact splice

The finite interval certificate proves `Q_(9,n)(-100)>0 for every 0<=n<=1242`.
The two ranges meet exactly between `n=1242` and `n=1243`, giving

`Q_(9,n)(-100)=H_(9,n)(-100)>0 for every integer n>=0`.

The existing cooperative-flow theorem can now propagate this endpoint
sign through `-100<=lambda<=0`; that composition remains a separate
fixed-order-nine artifact.
