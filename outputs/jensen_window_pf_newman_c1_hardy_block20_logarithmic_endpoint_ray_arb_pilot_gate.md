# Hardy block-20 logarithmic endpoint-ray Arb pilot gate

Date: 2026-08-07

Status: rigorous four-branch finite-input Arb pilot; not a proof of all-call or height-uniform W2--W4 control or RH

## Certified decomposition

This pilot evaluates the exact logarithmic-ray identity from Formal Core
Section 11.264 on one witness from each `(sign(Phi1),sign(Phi3))` branch.
For each endpoint ray it integrates on

```text
epsilon <= r <= R,   epsilon=1.00000000000000000000000000000000000000000000000000E-16, R=6.40000000000000000000000000000000000000000000000000E+1.
```

Near the origin the exact absolute majorant is

```text
|F'(x+r*e^(i theta))| * exp(a*epsilon)
* [epsilon*(1-log(a*epsilon))+a*epsilon^2/2],
a=2*pi*sin(theta).                                      (1)
```

It follows from
`sum_(k>=0) exp(-a*k*r)/(m+k) <= exp(a*r)*[-log(1-exp(-a*r))]`
and `1-exp(-a*r)>=a*r*exp(-a*r)`.  Thus the logarithmic endpoint is
enclosed rather than sampled.

For `r>=R`, the certified phase bound
`Im phase >= ell*r+q*r^2+c*r^3` is replaced by the linear lower bound
`(ell+q*R+c*R^2)r`.  The resulting polynomial-times-exponential integral
is evaluated explicitly, including the geometric tail of the summed
`1/n` kernel.

On the compact pieces the kernel is evaluated as a principal logarithm near
the origin and as

```text
K_m(w)=w^m*2F1(1,m;m+1;w)/m                            (2)
```

for small `|w|`, avoiding cancellation of the first `m-1` terms.

## Four-branch result

The 70- and 110-digit constructions overlap for every ray and complete
nonsaddle value.  Both independent identities contain zero:

```text
Q_exact^- - (P^- - I69^-) = 0,
(Q_exact^- - Q_paper^-) - R_paper = 0.
```

```text
chain   1  Phi1 positive Phi3 positive |formula-target| <= 2.81987901316120380531757283470994934759801253676414E-14
chain   2  Phi1 negative Phi3 negative |formula-target| <= 2.47523777642029717315289016887902562302770093083382E-14
chain   3  Phi1 negative Phi3 positive |formula-target| <= 2.39057651570487463390038418431515765405492857098579E-14
chain  36  Phi1 positive Phi3 negative |formula-target| <= 3.58883332360543979427935923354198166634887456893921E-14
```

Aggregate bounds:

```text
maximum omitted-origin radius       <= 7.33395552186407205313176671400045494913423498494520E-15
maximum far-tail radius             <= 7.48577985778680634169871981009110712649984275922310E-60
maximum |formula-target gap|        <= 3.58883332360543979427935923354198166634887456893921E-14
maximum |correction-residual gap|    <= 3.58883332889939571461873635271899729559663683176041E-14
```

This independently confirms the source orientation, endpoint harmonic term,
positive-`Phi1` zero mode, and all four ray signs with rigorous balls.

## Pi provenance and boundary

Pi is the mathematical constant inherited from `exp(2*pi*i*phase)`.  It
enters the exact kernel, phase decay, and endpoint harmonic coefficient; no
fitted normalization is used.

This is a rigorous pilot on four fixed low-height calls.  It does not yet
enclose the other 370 calls, derive a symbolic height-uniform W2--W4 error,
control recursive accumulation or the outer Hardy remainder, or prove
`Lambda<=0`, PF-infinity, RH, or a prize-level conclusion.
