# Newman Theta Adaptive Gamma-Scale Tail Gate

Date: 2026-07-25

Status: explicit cofinal omitted-tail theorem. This is not a
retained-separation theorem and not a proof of `Lambda<=0` or RH.

## Adaptive Cells

Set

```text
N(x)=ceil((1+x)^(3/4)).
C_N=((N-1)^(4/3)-1,N^(4/3)-1].
```

On `C_N`, the all-`N` `d0/d1` budgets are constant. For `x>=245`,
`exp(pi*x/8)/x^5` and `exp(pi*x/8)/x^6` increase, so each
gamma-normalized error is largest at the right endpoint.

## Cofinal Monotonicity

The forward and compact-reflected endpoint bounds have decreasing
geometric ratios, and their Gaussian exponents dominate the
increment of `exp(pi*x_N/8)`. For the large reflected term, write

```text
f_A(y)=y^A*exp(-c_safe*y^(4/3)),
c_safe=3*3^(2/3)*pi^(2/3)/16.
```

The normalized supremum-plus-integral tail decreases whenever

```text
(4/3)*(c_safe-pi/8)*y^(1/3)-A/y>0.
```

The worst case `A=22,y=63` is strictly positive by directed Arb.
All Hermite and Leibniz coefficients are positive, so the complete
`d0/d1` budgets inherit this decrease.

## Boundary Witnesses

| N | right x | J error / gamma | J' error / gamma |
|---:|---:|---:|---:|
| 62 | `[244.3892798001851396940412056992044185827414744060725117092215733554167 +/- 1.78e-68]` | `[1.750167937020677750070835672706186293011151979098644507045101834452345 +/- 2.90e-70]` | `[1.778851899991011897996464262506904294168463610470105220417600863502479 +/- 4.53e-70]` |
| 63 | `[249.6806040974726871588652632281722314445091373693887429032732467596498 +/- 3.82e-68]` | `[0.2131862947908342683091949461015912325087617764755618912931282051638912 +/- 3.40e-71]` | `[0.2166061676824924744964355490728303308241113737132384352775927075800318 +/- 5.55e-71]` |
| 64 | `[255.0000000000000000000000000000000000000000000000000000000000000000000 +/- 3e-72]` | `[0.02554891452578618993699670818864681976237035655250647524825125568094196 +/- 5.26e-72]` | `[0.02595020775031113278741846684752064572228748302940064037767098959645861 +/- 6.83e-72]` |

`N=62` is an explicit non-promotion guard. At `N=63`, both
ratios are below `1/4`, and the monotonicity theorem propagates
that bound through every later adaptive cell.

## Theorem

For every `0<=t<=1/5` and every real `x>=245`,

```text
N=N(x),
|J_t-J_(N,t)|<exp(-pi*x/8)/4,
|J_t'-J_(N,t)'|<exp(-pi*x/8)/4.
```

The remaining cofinal obligation is a retained first-jet lower
margin strong enough to dominate this quarter-gamma square.
The finite bridge `69<x<245` is also not certified here.
