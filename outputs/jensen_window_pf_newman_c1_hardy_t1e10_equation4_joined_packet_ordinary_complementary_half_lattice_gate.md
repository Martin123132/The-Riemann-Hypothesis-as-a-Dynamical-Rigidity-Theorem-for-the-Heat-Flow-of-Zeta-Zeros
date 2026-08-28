# Complementary half-lattice reduction of the ordinary packet

Date: 2026-08-27

Status: exact ordinary-packet complement and stationary-label census certified;
quantitative complementary-tail enclosure remains open.

For `622<=m<=39852`, put

```text
w_m(u)=(1+u/m)^(-s)exp(-i*pi*u^2),
beta_m=integral_(-1/2)^(1/2) w_m(u)
       sum_(q=q_-)^(q_+-1) exp[2*pi*i(q-m)u]du,       (CL1)

q_-=79788.5,       q_+=2561211.5.
```

The complete half-integer Fourier lattice has the Abel/Poisson identity

```text
sum_(q in Z+1/2) exp[2*pi*i(q-m)u]
 =sum_(k in Z)(-1)^k delta(u-k).                    (CL2)
```

Only `k=0` lies in the cell, and `w_m(0)=1`.  Consequently the complete
lattice contributes exactly one, so the finite coefficient defect is exactly
the negative of the two missing half-lattice tails.  Reassembling the finitely
many cells before taking the Abel limit gives

```text
O_join=sum_(m=622)^39852 m^(-s)[beta_m-C_G]
      =-T_lower-T_upper
       +(1-C_G)sum_(m=622)^39852 m^(-s),             (CL3)

T_lower=sum_(q in Z+1/2, q<=79787.5)
        integral_L^a y^(-s)exp(-i*pi*y^2+2*pi*i*q*y)dy,

T_upper=sum_(q in Z+1/2, q>=2561211.5)
        integral_L^a y^(-s)exp(-i*pi*y^2+2*pi*i*q*y)dy,

L=621.5,       a=39852.5.                            (CL4)
```

No source label is duplicated: the finite roster, lower complement, and upper
complement partition `Z+1/2` exactly.

The Gamma factor also simplifies exactly.  The phase of `B=I_bulk/I_G` is
`theta-theta_0`; combining its modulus with the exact `H(t)` normalization
gives

```text
C_G(t)=[1+exp(-2*pi*t)]^(-1/2),
0<1-C_G<exp(-2*pi*t)/2.                              (CL5)
```

At `t=10^10`, the entire final term in (CL3) has

```text
log10 absolute upper bound
 < [-27287527074.594509719139191054181568140886570852598338579157932597698166119909387 +/- 1.85e-70].
```

It is retained exactly; it is not silently replaced by zero.

For one complementary label the phase is

```text
Phi_q(y)=-t log(y)-pi*y^2+2*pi*q*y,
Phi_q'(y)=2*pi[q-y-p/y],       p=t/(2*pi).           (CL6)
```

Since `y+p/y` decreases on `[L,a]`, the lower complement has no stationary
label.  Its nearest gap is

```text
[0.99977213357650558528658741970065477941051137480678694369087482511702759589363029 +/- 1.46e-82].
```

The upper complement has exactly 230 stationary labels,

```text
q=2561211.5, 2561212.5, ..., 2561440.5,             (CL7)
```

because

```text
y+p/y at L minus q_+
 = [229.67967651384986136385152879102227610049461715927188523481993397904479849613982 +/- 4.85e-78].
```

Their lower saddle roots all lie in the first ordinary cell:

```text
[621.50016499467458464125532308222573458314325504102629961178993300508841771420086 +/- 5.92e-79]
 <= y_q <=
[621.55576081192507870131278288344926918563290127646748069839952373328749026286757 +/- 3.61e-78].
```

The next upper label is nonstationary with endpoint gap
`[0.32032348615013863614847120897772389950538284072811476518006602095520150386017516 +/- 5.22e-83]`.  Thus the old 39,231-cell
ordinary wall is reduced to one 230-label endpoint-saddle packet plus two
nonstationary complementary tails and the explicit negligible Gamma defect.

Pi provenance: each `pi` comes from the inherited Fresnel phase, the exact
half-integer Fourier spacing, or `p=t/(2*pi)`.  No fitted geometric constant is
introduced.

Proof boundary: exact Abel/Poisson half-lattice completion, exact
complementary-tail identity, exact Gamma coefficient reduction and bound, and
production stationary-label census only.  No quantitative enclosure of the
230-label packet or nonstationary tails, complete ordinary or joined packet,
`J_Z`, or `D_K` is proved, nor is any non-A, all-height, `Lambda<=0`,
PF-infinity, RH, or prize-level conclusion.
