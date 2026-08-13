# Hardy block-20 retained-decay endpoint majorant gate

Date: 2026-08-07

Status: rigorous finite-roster application of a selector-gap-finite retained-decay theorem; not a proof of the complete evaluator or RH

## Summed generic modes

For one endpoint family choose the same sign-aware ray as in Sections
11.270--11.272.  Let `alpha=2*pi*sin(theta)`,
`beta=2*pi*B*|sin(2 theta)|`, `gamma=2*pi*C*|sin(3 theta)|`, and

```text
phi(x)=(1-exp(-x))/x,  phi(0)=1.
```

Integrating the linear-to-cubic homotopy parameter before applying the
triangle inequality and then summing the modes geometrically gives

```text
S(B,C,q)=2*pi integral_0^infinity
 (B*r^2+C*r^3) exp(-alpha*q*r)
 phi(beta*r^2+gamma*r^3)/(1-exp(-alpha*r)) dr.       (1)
```

This is no larger than the Hurwitz bound in Section 11.272, but it retains
all common-ray quadratic and cubic decay.  Its apparent origin singularity is
removable.  The gate integrates from `10^-12` to `64`, bounds the omitted
origin explicitly, and bounds infinity by exact exponential moments.

For W2 and W3 the reciprocal boundary has already cancelled.  Retaining the
nonnegative selector-gap decay gives

```text
E_b=2*pi*C integral (3*N*r^2+r^3) exp(-alpha*delta*r-a_b*r^2) dr,
E_c=2*pi*C integral r^3 exp(-alpha*eta*r-a_c*r^2) dr. (2)
```

Both remain finite at `delta=0` or `eta=0`; no reciprocal selector gap has
been reintroduced.

## Finite block result

All 374 calls overlap across the 60/90-digit Arb ladder and dominate their
independently enclosed upper and lower endpoint corrections.  Transport
through the exact source-observed block weights gives

```text
maximum prior transported bound       8.58013843685509986097703809709438792569406786379705E-3
maximum retained-decay bound           2.18921492147663827569784546460747795269403233708512E-3
requested input scale                  5.00000000000000000000000000000000000000000000000000E-3
maximum retained bound/scale ratio     4.37842984295327655139569092921495590538806467417025E-1
outputs below requested scale          15 / 15
minimum improvement factor             3.91927642767378283482263672383620163685420119216970E+0
```

This closes the **local endpoint-model portion** of the finite block-20
nominal budget without using observed signed cancellation.  It does not close
coefficient transport, recurrence arithmetic, Legendre truncation, other
blocks, the outer Hardy representation, or a height-uniform theorem.
