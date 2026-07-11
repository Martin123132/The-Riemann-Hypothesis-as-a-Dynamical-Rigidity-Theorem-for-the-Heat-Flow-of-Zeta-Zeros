# Jensen-Window PF Newman Classical Three-Block Residual Gate

Date: 2026-07-11

Status: exact and interval-certified classical three-block residual
obstruction. This is not a proof or disproof of RH or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_newman_classical_three_block_residual_gate`.

```text
work/rh_compute/results/jensen_window_pf_newman_classical_three_block_residual_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_classical_three_block_residual_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_classical_three_block_residual_gate.py
```

Current result:

```text
validated Jensen-window PF Newman classical three-block residual gate: 10 rows, 0 issues, 2 established classical real-zero benchmarks, 2 positive-kernel residual theorems, 1 compact parameter reduction, 1 exact bivariate Laguerre identity, 3 Acb spectral certificates, 64908 parameter boxes covered, 2 classical residual obstructions, 1 Gasper square-scope guard, 1 coupled handoff
```

## Classical Blocks

In the corpus normalization,

```text
P_a(x)=pi^2/2*(K_(a+i*x/4)(2*pi)+K_(a-i*x/4)(2*pi))
B_P2=P_(9/4)-a*P_(5/4)
B_dB=P_(9/4)+b*P_(5/4)+P_(1/4)
```

Polya's P2 transform and de Bruijn's transform are established
real-zero benchmarks in the cited sources. The present gate checks
their normalization and studies what remains after subtracting them.

## Positive Residuals

For `q=exp(-4u)` and `a=3/(2pi)`, the P2 first-summand
residual is bounded below by

```text
-r**2*(2*pi*r**7 - 3*r**5 + 3*pi*r**2 - 2*pi**2)/(2*pi)
Phi(u)-Psi_P2(u)>0 for every finite u>=0
```

For de Bruijn's block, the fourth Taylor lower bound for `exp(pi*q)`
reduces the same claim to a degree-six polynomial in `sqrt(q)`.
All seven Bernstein coefficients are rigorously positive:

```text
b_0=[2.434802200544679309417245499938075567656849703620395313206674688110022411209602621500886701859276115912012956887011572038886174061 +/- 1.75e-130]
b_1=[2.268135533878012642750578833271408900990183036953728646540008021443355744542935954834220035192609449245346290220344905372219507394 +/- 5.34e-130]
b_2=[2.288903419868520982777183498186338147916161106605918996537311215203591534422688905320226083176646775768795796578817863621811499804 +/- 8.53e-130]
b_3=[2.363899467300498967939259890024765318530094389469279539461987141384173480179685003291020174959220204976664998045780025160841395838 +/- 8.79e-130]
b_4=[2.46600468670120314091076781159928187969908611029512533472462690290257468815116990333458215861161095514110468828535793489817158547 +/- 5.09e-129]
b_5=[2.50752082367418528193075980652650208839536557502710665821634316693363100501414471762808973446188268664544626587578843761168093742 +/- 3.05e-129]
b_6=[0.23561139971388698116944480086655912489054096876672310475806986675497466877709272192333297389053444825315244504448097238314803141 +/- 3.26e-129]
Phi(u)-Psi_dB(u)>0 for every finite u>=0
```

Thus both classical real-rooted bases leave pointwise positive
kernel residuals. Those residual transforms nevertheless fail below.

## Full Positive Family

```text
Psi_(beta,gamma)(u)=4*pi^2*exp(-2*pi*cosh(4u))*(cosh(9u)+beta*cosh(5u)+gamma*cosh(u))
beta<=b=pi-3/(2*pi)<821/308<27/10, using pi<22/7
Phi(0)/(4*pi^2*exp(-2*pi))<61/10, hence beta+gamma<51/10
Every beta,gamma>=0 with a globally nonnegative residual lies in 0<=beta<27/10, gamma>=0, beta+gamma<51/10.
```

The origin bound uses the first eight exact theta summands and an
explicit geometric majorant for `n>=9`. These are necessary
conditions for a globally nonnegative residual and place every
candidate in one rational triangle.

## Interval Theorem

The exact bivariate expression is

```text
L[R_(beta,gamma)]=L[E]-beta*B[E,P5]-gamma*B[E,P1]+beta^2*L[P5]+beta*gamma*B[P5,P1]+gamma^2*L[P1]
```

For each real entire component, central differences with step h enclose f' and f''. An Acb square bounds abs(f) by M on every Cauchy circle of radius r around the difference segment; the added errors are h^2*M/r^3 and 2*h^2*M/r^4.
The rational triangle is covered by closed `1/80` boxes:

```text
I_i=[i/80,(i+1)/80], J_j=[j/80,(j+1)/80], 0<=i<216, 0<=j<408-i
total boxes=64908
x=48 assignments=12929
x=52 assignments=5833
x=86 assignments=46146
assignment sha256=48e36d40cd04261b379508ff7afb58ed7c95c85df9dee3d2cef29d0bb50731e8
```

Every box has a strictly negative Arb enclosure at one assigned
spectral point. In particular,

```text
L[R_P2](86)=[-2.34e-27 +/- 5.35e-30]
L[R_dB](52)=[-4.8e-17 +/- 1.15e-19]
```

So neither classical positive residual is an independent
Laguerre-Polya component, and no nonnegative-coefficient 9/5/1
block with a globally nonnegative residual can repair this route.

## Square Scope

For F_(a,c)(z)=K_(i(z-i*c))(a)+K_(i(z+i*c))(a), a zero gives equality of two Bessel moduli. Gasper converts their difference into one integral of [f_t(y+c)-f_t(y-c)] times K_(i*x)(a/sqrt(t))^2, and convexity forces y=0.
A linear combination with several shifts does not yield that two-modulus equality. Its squared modulus contains mixed products between distinct shifts, which the single-shift identity does not sign.

This leaves a precise open handoff: derive a genuinely coupled
square or matrix identity for signed higher blocks after modular
cancellation. The published single-shift square alone does not do it.

References: https://arxiv.org/abs/0801.2996,
https://arxiv.org/abs/1502.06844,
https://doi.org/10.1215/S0012-7094-50-01720-0, and
https://doi.org/10.1007/BF02565336.
