# Order-Ten Lambda=-100 Delayed Entry Certificate

Date: 2026-07-16

Status: rigorous delayed signed order-ten entry and complete endpoint
sign chart at `lambda=-100`. This is not all-shift positivity,
and this certificate is not a proof of PF-infinity, RH, or `Lambda <= 0`.

```text
z_1''(t)<=4200/t^2 for every real t>=1251
z_1''(t)<=4200/t^2 => Z_k^(1)<=4200*[-log(1-1/k^2)]<4201/k^2, k>=1252
|Z_k-Z_k^(1)|<10/k^2 for every integer k>=1252
z_1''(t)<=4200/t^2 on t>=1251 => Z_k<4201/k^2+10/k^2=4211/k^2<5500/k^2, k>=1252
Q_(10,n)(-100)>0 for every n>=1243
```

The analytic tail joins the finite positive block exactly after
`n=1242`, proving

```text
Q_(10,n)(-100)>0 for every integer n>=4
Q_(10,n)(-100)<0 for 0<=n<=3 and Q_(10,n)(-100)>0 for every n>=4
```

The rows `n=0,1,2,3` remain negative and are part of the theorem, not
a numerical exception to be discarded.
