# Hardy block-20 endpoint-linearization decomposition pilot

Date: 2026-08-07

Status: rigorous four-branch finite-input decomposition; not an all-height proof and not a proof of RH

## Exact decomposition

At endpoint `x`, write `A_x=F'(x)` and

```text
F(x+u)=F(x)+A_x*u+B_x*u^2+Phi3*u^3,
B_x=Phi2+3*Phi3*x.
```

The exact ray is compared with its endpoint-linear model, not with a uniform
quadratic model.  Summing the linear rays analytically gives

```text
J_b,lin(x)=i*E_-(x)/(2*pi) [psi(m)-psi(m-A_x)],
J_c,lin(x)=i*E_+(x)/(2*pi) [psi(1+A_x)-psi(1)].       (1)
```

The four linear rays plus the exact harmonic endpoint term equal the paper's
endpoint half-sum plus W5.  Therefore

```text
Q_exact-Q_paper
 = R_b+R_c+zero_mode-W2-(W3 or W4),                  (2)
```

where `R_b,R_c` are the full-ray minus linear-ray endpoint differences with
the source orientations retained.  W2 and W3 are the paper's exceptional
quadratic-mode replacements; W4 is its positive-Phi1 zero-mode replacement.

## Four-branch certificate

The 70- and 110-digit Arb constructions overlap throughout.  Both (1) and
(2) enclose equality on chains 1, 2, 3, and 36:

```text
chain   1: |generic-W5 gap| <= 1.70024966440998848199119869662965845588927394826072E-108, |correction gap| <= 3.53657251837246049564617011995437678706366568803787E-14
chain   2: |generic-W5 gap| <= 9.32915047949097693369688566540098189229422553891695E-109, |correction gap| <= 3.29901692435632227194675092363240764825604856014252E-14
chain   3: |generic-W5 gap| <= 6.87860748350171474954842254306579245356776838535732E-109, |correction gap| <= 3.29681104939134710096748648311404394917190074920654E-14
chain  36: |generic-W5 gap| <= 1.82388016188312862703883057071963467016256349900634E-106, |correction gap| <= 4.99380188630967080699873950067058103741146624088287E-14
```

Aggregate bounds:

```text
maximum |generic linear aggregate-(half+W5)| <= 1.82388016188312862703883057071963467016256349900634E-106
maximum |decomposition prediction-correction| <= 4.99380188630967080699873950067058103741146624088287E-14
```

This identifies the correct analytic remainder architecture and rejects the
naive replacement of every endpoint ray by one quadratic logarithmic ray.
The next step is an all-374 replay and explicit majorants for `R_b`, `R_c`,
and the exceptional W2--W4 differences.  This pilot proves no height-uniform
bound, recursive accumulation, outer Hardy control, `Lambda<=0`, PF-infinity,
RH, or prize-level conclusion.
