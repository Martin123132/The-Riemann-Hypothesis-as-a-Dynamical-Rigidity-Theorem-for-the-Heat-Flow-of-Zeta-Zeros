# Small-`tau` alternating Taylor corner box

Date: 2026-08-24

Status: certified local interval lemma; recursive interior boxes remain open

## Status

Certified interval lemma at `t=10^10` for the one-sided physical box

```text
x in [1/2-10^-22, 1/2],   s in [0,10^-22].
```

This closes the local small-`tau` neighbourhood branch at the corner. It is
not a physical quadrature, a non-A bound, or an RH-level result.

## Exact normalization

Write

```text
e(u)=exp(2*pi*i*u),
W_n(p;a,t)=sum_(k=0)^n (p+k) exp(2*pi*i*(a*k+t*k^2)).
```

For this corner, `p=39894` and
`n=1240710`. Integer curvature periodicity and
conjugation give the exact identity

```text
(ST1)  W_n(p;w,sigma)=conjugate(W_n(p;-w,-sigma-1)).
```

The normalized parameters are centred at `a_0=-1/2`, `t_0=0`, and the raw
`a` coordinate is retained continuously rather than reduced modulo one.

## Cancellation-preserving Taylor enclosure

Put `alpha=a+1/2`. Since
`exp(2*pi*i*(-k/2))=(-1)^k`, define the exact signed moments

```text
M_j=sum_(k=0)^n (-1)^k (p+k) k^j.
```

Then

```text
(ST2)  W_n(p;a,t)
       = M_0 + 2*pi*i*(alpha*M_1+t*M_2) + R_2,

M_0 = 660249,
M_1 = 794429714775,
M_2 = 985657301007258645.
```

Here `pi` is not an inserted geometric constant: the factor `2*pi` comes
directly from differentiating the defining phase
`exp(2*pi*i*(a*k+t*k^2))`.

For real `d`, the integral Taylor remainder gives
`|exp(i*d)-1-i*d|<=d^2/2`. Therefore, with
`S_j=sum_(k=0)^n (p+k)k^j`,

```text
(ST3)  |R_2| <= 2*pi^2 [delta_a^2*S_2
                         +2*delta_a*t_max*S_3
                         +t_max^2*S_4].
```

All `M_j` and `S_j` are evaluated by exact integer power sums split into the
even and odd residue classes. Only the quadratic Taylor remainder is placed
under an absolute-value bound.

## Certified source reassembly

The normalized child is conjugated back, multiplied by the first exact
Mordell prefactor, and joined to the already certified physical endpoint box:

```text
(ST4)  SourceBox = P_1 * conjugate(TaylorBox) + EndpointBox,
P_1 = 2*x^(-3/2) exp(pi*i*(1/4+r*(A+2s)-r^2/x)).
```

The complete source radius is at most
`0.01018014730652794415377560000024459441192` and the box contains the exact
periodic source current at `(x,s)=(1/2,0)`. The exact normalized terminal has
period `2`.

## Scale audit

| half-width | max |a+1/2| | max t | linear radius | quadratic remainder | propagated child variation |
|---:|---:|---:|---:|---:|---:|
| `1e-18` | `79788/499999999999999999` | `1/499999999999999999` | `13.18266638421671288483594253193587064743` | `57.87134706800636507750823511742055416107` | `401.9421979406969285264494828879833221436` |
| `1e-20` | `79788/49999999999999999999` | `1/49999999999999999999` | `0.1318266638421671244074673268187325447798` | `0.005787134706800636449464114718921337043867` | `0.7784612011105166429203450206841807812452` |
| `1e-21` | `79788/499999999999999999999` | `1/499999999999999999999` | `0.01318266638421671313463612307259609224275` | `0.00005787134706800636232623680221820450242376` | `0.07489979213078218345245318232628051191568` |
| `1e-22` | `79788/4999999999999999999999` | `1/4999999999999999999999` | `0.001318266638421671226727438508419254503679` | `0.0000005787134706800637079656627476120789310698` | `0.007460515933275521129053320379398428485729` |
| `1e-24` | `79788/499999999999999999999999` | `1/499999999999999999999999` | `0.00001318266638421671172517329884144032803306` | `5.787134706800637463261636547587373532314e-11` | `0.00007457274972497223927431903023332893098996` |
| `1e-28` | `79788/4999999999999999999999999999` | `1/4999999999999999999999999999` | `0.00000000131826663842167131677762871339205930088` | `5.78713470680063733535770152731557193237e-19` | `0.000000007457242238793364175131740093256660362186` |

The certified half-width `10^-22` is one million times the earlier direct
first-difference width `10^-28`. This is a local Taylor gain, not a global
small-`tau` theorem.

## Validation

- Production: `320` bits,
  Mordell cutoff `10`.
- The independent checker uses exact direct witnesses for every power-sum and
  parity formula, verifies the three full moments, and replays the endpoint
  and source box at 384 bits with cutoff 11.
- Builder and dependency hashes are pinned in the JSON artifact.

## Proof boundary

A one-sided width-1e-22 small-tau affine-current Taylor box at the physical x=1/2, s=0 corner, its exact period-two terminal join, and one complete source-current box only. No uniform small-tau branch away from this corner, recursive interior parameter boxes, physical quadrature, non-A bound, joined R_after_A, R_Dir, Q_K-T, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.
