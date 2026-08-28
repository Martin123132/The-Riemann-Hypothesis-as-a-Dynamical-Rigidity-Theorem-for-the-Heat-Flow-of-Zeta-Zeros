# Outer differentiated Fresnel remainder for the completed B trace

Date: 2026-08-13

Status: exact remainder-transport and uniform outer derivative-tail
certificate; not a proof of the complete completed-B estimate

Write the sign-adapted Fresnel tail as

```text
z(q)=exp(-i*pi*q^2/2)T(q),
P_3(q)=i/(pi q)+1/(pi^2 q^3)-3i/(pi^3 q^5),
r(q)=z(q)-P_3(q).                                      (OR1)
```

The already certified contour remainder gives
`|r|<=30/(pi^4|q|^7)`.  Differentiating the exact tail equation before
taking norms gives

```text
r' =-i*pi*q*r-15i/(pi^3 q^6),
r''=-(i*pi+pi^2 q^2)r-15/(pi^2 q^5)
     +90i/(pi^3 q^7),                                  (OR2)

|r'| <=45/(pi^3|q|^6),
|r''|<=45/(pi^2|q|^5)+120/(pi^3|q|^7).                (OR3)
```

For either completed outer sign, `R_s=f r(q_s)` with
`f=m sqrt(2)x^(-3/2)`.  On `m>=B=5122421` and
`1e-4<=x<=1/2`,

```text
|q_s|  >=(3/4)sqrt(2/x)m,
|q_s'| <=(5/4)sqrt(2/x)m/(2x),
|q_s''|<=(13/4)sqrt(2/x)m/(4x^2).                     (OR4)
```

Applying the product/chain rules and summing both signs by the decreasing
series integral test gives

```text
sum_s sum_(m>=B)|R_s|    <= [8.177621601255279253174440019086446364896919897059543984434365095403272983840093444249992989005817510e-36 +/- 1.05e-135]
sum_s sum_(m>=B)|R_s'|   <= [6.319730068315199634077323306801988973447805701160105620177883354824859679200973119821076858770145190e-21 +/- 7.19e-121]
sum_s sum_(m>=B)|R_s''|  <= [146.5181078942200066680557500533319227804485186798367466632727819160310424962248423992510269719693360 +/- 5.50e-98].                (OR5)
```

After the Kummer weight and normalized coordinate are included,

```text
|w R|          <= [8.177826054573802780018255559568164861249116756756913378719194779069516460472361039900530729664953229e-35 +/- 1.51e-134]
|d_xi(w R)|    <= [2.168999188133511393447220559901740520365935698946404074621566372066930864893293091184846176067238795e-28 +/- 4.85e-126]
|d_xixi(w R)|  <= [1.725846495155548946454921996039811633372624444049255567562930306831702427118996464651306648287028573e-14 +/- 7.70e-112].                         (OR6)
```

Thus the asymptotic remainder is not being dropped: its complete outer tail
is uniformly controlled through the two derivatives needed for tangential
integration by parts.  The companion rational-current gate and this gate do
not yet close the analytic block: the exact Fresnel-to-boundary dictionary
term `D_m` from Section 11.411 must also be differentiated and summed.

The dictionary term, finite block `m<B`, local phase-coupled tangential
integration, window and endpoint boundary terms, and the common-regulator
`R_join` combination remain open.  This gate establishes none of a complete
B estimate, `T_upper`, a height-uniform theorem, `Lambda<=0`, PF-infinity,
RH, or a prize-level conclusion.

Pi provenance: `pi` is the same equation-(9) Fresnel phase constant used in
the exact tail differential equation; no fitted or decorative value enters.
