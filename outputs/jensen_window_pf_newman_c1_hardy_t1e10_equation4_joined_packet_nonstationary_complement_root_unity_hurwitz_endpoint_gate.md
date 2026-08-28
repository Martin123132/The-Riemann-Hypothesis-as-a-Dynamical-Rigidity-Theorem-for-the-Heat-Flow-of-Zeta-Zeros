# Root-of-unity Hurwitz endpoint algebra for the nonstationary complements

Date: 2026-08-27

Status: every phase-weighted endpoint denominator sum through power
`k=16` is certified for `R_B`, `R_L`, and the lower complementary
tail; no integrated remainder or complete tail enclosure is claimed.

Put

```text
p=t/(2*pi),              Q(y)=y+p/y,
D_q(y)=q-Q(y),           phi_q'(y)=2*pi*D_q(y).     (EH1)
```

For a nontrivial root of unity `omega^h=1` and `M=hC`, define

```text
F_fin(k;delta,omega,M)=sum_(n=0)^(M-1) omega^n/(n+delta)^k.
```

Grouping by residue class gives, for `k>1`,

```text
F_fin=h^(-k) sum_(r=0)^(h-1) omega^r
 [zeta(k,(delta+r)/h)-zeta(k,(delta+r)/h+C)].       (EH2)
```

At `k=1` the bracket in (EH2) is replaced by
`psi((delta+r)/h+C)-psi((delta+r)/h)`.  Direct Arb summation of all 240
terms independently overlaps (EH2) at both `B` and `a` for every
`1<=k<=16`.

For the infinite root-of-unity sum,

```text
F_inf=h^(-k)sum_r omega^r zeta(k,(delta+r)/h),       k>1,
F_inf=-h^(-1)sum_r omega^r psi((delta+r)/h),         k=1.       (EH3)
```

The apparent `k=1` Hurwitz poles cancel exactly because

```text
sum_(r=0)^15 exp(2*pi*i*9*r/16)=0,
1+exp(pi*i)=0.                                      (EH4)
```

Thus `pi` here comes from the inherited Fresnel/Fourier phase and
`p=t/(2*pi)`.  The order 16 is forced by the exact fractional endpoint
`B=621+9/16`; it is not a fitted circle constant.

The production shifts remain strictly positive.  Representative certified
values are

```text
delta_RB(B)    = [27.7560315984190023728366412973428342339549046725188912242283133747992482 +/- 4.08e-71],
delta_RL(L)    = [10.3203234861501386361484712089777238995053828407281147651800660209552015 +/- 3.87e-72],
delta_lower(a) = [0.999772133576505585286587419700654779410511374806786943690874825117027596 +/- 1.07e-73].
```

At the half-integer endpoints the infinite sums are alternating and satisfy

```text
F_inf(k,delta,-1)+F_inf(k,delta+1,-1)=delta^(-k).   (EH5)
```

Arb interval checks of (EH5), and independent agreement with the explicit
alternating eta/digamma form, pass for every power through 16.  The
lower tail additionally carries the exact orientation factor `(-1)^k` because
`D=-(delta+n)`.

Starting from `h_0=y^(-1/2)` and

```text
h_(n+1)=d/dy[h_n/D],                                (EH6)
```

the certified table supplies all endpoint denominator sums needed through
the first eight integration rounds.  It does not bound the integral containing
`h_8`; that is the next quantitative obligation.

Proof boundary: exact root-of-unity/Hurwitz endpoint identities, explicit
`k=1` pole cancellation, production positive-gap checks, and rigorous endpoint
sum balls through power 16 only.  No integrated remainder, `R_B`, `R_L`, lower
complementary tail, complete ordinary or joined packet, `J_Z`, `D_K`, non-A,
all-height, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
