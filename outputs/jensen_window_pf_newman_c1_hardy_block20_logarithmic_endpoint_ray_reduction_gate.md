# Hardy block-20 logarithmic endpoint-ray reduction gate

Date: 2026-08-07

Status: exact analytic reduction with non-rigorous branch diagnostics; not a proof of W2--W4 error control or RH

## Summed ray kernel

For `Im(u)>0` and integer `m>=1`, define

```text
K_m(u)=sum_(n=m)^infinity exp(2*pi*i*n*u)/n
      =-Log(1-exp(2*pi*i*u))-sum_(n=1)^(m-1) exp(2*pi*i*n*u)/n.  (1)
```

The principal logarithm is analytic here because
`|exp(2*pi*i*u)|<1`.  Near the ray origin, `K_m(u)=-Log(-2*pi*i*u)+O(1)`,
so its logarithmic singularity is integrable.  This is the exact summed
version of the absolute-interchange statement in the sector-contour gate.

Put `E_-(z)=exp(-2*pi*i*F(z))`, `E_+(z)=exp(2*pi*i*F(z))`, and

```text
J^-_m(x;theta)=integral_0^infinity F'(x+r*e^(i*theta))
               E_-(x+r*e^(i*theta)) K_m(r*e^(i*theta)) e^(i*theta) dr,

J^+_m(x;theta)=the same expression with E_+ in place of E_-.
```

With `L=floor(xi)`, `m=L+1`, the exact source-oriented nonsaddle pieces are

```text
B^-=J^-_m(0;theta_b)-J^-_m(N;theta_b),
C^-=conj(J^+_1(N;theta_c)-J^+_1(0;theta_c)),           (2)
```

where `theta_b` and `theta_c` are the sign-aware angles from the preceding
gate.  Integer `N` is essential: it makes `K_m(N+u)=K_m(u)`.

## Exact recurrence remainder identity

Let

```text
P^-=sum_(k=0)^N E_-(k),
I69^-=sum_(n=j0)^L integral_0^N exp(2*pi*i*(n*y-F(y)))dy,
j0=ceil(Phi1) in {0,1},
H_L=sum_(n=1)^L 1/n.
```

Then direct integration by parts in `(67a)`, followed by (2), gives

```text
P^-=I69^-+Q_exact^-,

Q_exact^-=(1+E_-(N))/2 + i*(E_-(N)-1)*H_L/(2*pi)
          +B^-+C^-+1_(Phi1>0)*integral_0^N E_-(y)dy.  (3)
```

Thus the infinite nonsaddle families are reduced to four endpoint-ray
integrals, one optional finite zero-mode integral, and explicit endpoint
terms.  No truncation in `n` appears in (3).

## Roster and diagnostics

The exact selector audit covers 374 calls, with
`L=1 or 2`, `m=2 or 3`, and all four `(sign(Phi1),sign(Phi3))` combinations.
The following high-precision mpmath integrations test orientation only:

```text
chain   1: Phi1 positive, Phi3 positive, |gap|=2.0403116723908029322e-20
chain   2: Phi1 negative, Phi3 negative, |gap|=5.0689777221650352658e-49
chain   3: Phi1 negative, Phi3 positive, |gap|=2.7498248862893647909e-34
chain  36: Phi1 positive, Phi3 negative, |gap|=7.0057575059703867783e-32
```

These four values are non-rigorous diagnostics and are not used to certify
the exact theorem.  Their purpose is to catch conjugation, endpoint, and
zero-mode mistakes before implementing Arb quadrature.

## Pi provenance and boundary

Pi in (1)--(3) is inherited from the paper's Fourier exponential and fixes
both the logarithmic kernel and phase orientation.  It is not inserted as a
fitted normalization.

The reduction is exact under the already certified contour and interchange
hypotheses.  It does not yet rigorously enclose the logarithmic endpoint
integrals, compare them with W2--W4 on all 374 calls, establish a
height-uniform recurrence, control the outer Hardy remainder, or prove
`Lambda<=0`, PF-infinity, RH, or a prize-level conclusion.
