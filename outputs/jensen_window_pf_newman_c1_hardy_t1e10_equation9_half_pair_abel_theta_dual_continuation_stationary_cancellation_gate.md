# Modular-dual continuation stationary cancellation

Date: 2026-08-13

Status: exact leading-stationary cancellation and complete saddle-roster
partition; not a proof of a full modular-continuation bound

For a positive dual index, write

```text
alpha=D-2n,  z_*=x alpha/D,  r_*=(D-alpha)/(x alpha). (MC1)
```

At the `z` saddle the dual phase and Hessian are exactly

```text
psi_alpha(x)=pi alpha^2 x/4+(t/2)log((1-x)/x),
Psi_zz(z_*)=-pi D^3/[2(D-alpha)x].                    (MC2)
```

Include the modular factor `r^(-1/2)`, the endpoint driver
`z^(-1/2)(2+i*pi*D^2*z)`, and the leading stationary scale
`sqrt(2*pi/|Psi_zz|)`.  Their product simplifies without an asymptotic
approximation to

```text
M_D(alpha,x)=4 sqrt(x)/D+2 i*pi*alpha*x^(3/2).        (MC3)
```

For each matched continuation pair

```text
alpha_(B,L+k)=alpha_(A,k)=A-2k,                      (MC4)
```

the common stationary phase, Hessian signature, dual `+/-n` multiplicity,
and universal second term in (MC3) agree.  With the compulsory endpoint
signs, the leading difference is therefore exactly

```text
M_B-M_A=4 sqrt(x)(1/B-1/A).                          (MC5)
```

This cancels the larger `2 i*pi*alpha*x^(3/2)` channel.  It does not cancel
the scalar endpoint channel, the stationary-expansion remainder, or any
face and corner boundary term.

The arithmetic threshold is rigorously separated by Arb:

```text
A-2 < sqrt(8t/pi) < A,
sqrt(8t/pi)=[159576.9121605730711759784239737527473903434524659738630663703318682631703597207354005009335629227746 +/- 3.28e-95],
A-sqrt(8t/pi)=[0.08783942692882402157602624725260965654753402613693362966813173682964027926459949906643707722542970855 +/- 3.01e-96],
sqrt(8t/pi)-(A-2)=[1.912160573071175978423973752747390343452465973863066370331868263170359720735400500933562922774570291 +/- 3.01e-96]. (MC6)
```

Consequently the dual continuation has the complete partition

```text
k=0: exceptional A-face/B-interior fold chart;
1<=k<=79788: matched interior-z saddles, but no real x saddle;
k>=79789: alpha<=-1, hence no interior z saddle.             (MC7)
```

For the middle block,

```text
|psi_alpha'(x)| >= 2t-pi(A-2)^2/4
                      > [479304.7041125798458878933024084241487050751473858115512743255266943577473971531088360046745008905107 +/- 3.37e-91]
```

throughout `0<x<=1/2`.  Thus its reduced leading phase is uniformly
nonstationary.  The surviving scalar coefficient in (MC5) is
`4*(159577-5122421)/(159577*5122421)` with Arb ball
`[-2.428538818921985054081296773582498652697662687871837964467097719657559868827644325303470803504936578e-5 +/- 3.46e-105]`.

Pi provenance: `pi` is inherited from the equation-(9) Fourier/Kummer phase
and the canonical Jacobi transformation.  The checker independently bounds
it from Machin's identity using exact rational alternating sums.

Proof boundary: the identity (MC5), the threshold separation, and the
continuation saddle partition are proved only for the leading interior
`z`-stationary contribution.  No uniform stationary-expansion remainder,
summed scalar-channel estimate, exceptional `k=0` fold splice, complete
modular continuation bound, source-minus-target estimate, `T_upper`,
`Lambda<=0`, PF-infinity, RH, and any prize-level conclusion remain open.
