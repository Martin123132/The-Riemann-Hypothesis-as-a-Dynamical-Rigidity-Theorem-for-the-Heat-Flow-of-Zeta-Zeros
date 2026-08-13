# Weber source-roster and same-index theta-tail reassembly

Date: 2026-08-13

Status: exact fixed-height reassembly and source normalization; not a proof
of a bound for the surviving theta tail or zero-mode sector

Write `B_n` and `A_n` for the positive modular-dual endpoint kernels and
`Full_(D,n)`, `Tail_(D,n)` for their Weber completion and `p>=sqrt(x)`
tail.  The licensed cutoff shift gives

```text
sum_(n>=1)(B_n-A_n)
 =sum_(n=1)^L B_n+sum_(k>=1)(B_(L+k)-A_k).           (WR1)
```

Every fixed continuation completion cancels, while
`B_n=Full_(B,n)-Tail_(B,n)`.  Therefore (WR1) becomes exactly

```text
sum_(n>=1)(B_n-A_n)
 =sum_(n=1)^L Full_(B,n)
  +sum_(n>=1)[Tail_(A,n)-Tail_(B,n)].                (WR2)
```

For `alpha=B-2n`, the completed kernel is

```text
Full_(B,n)=2*i*pi*alpha*x^(3/2)
              exp(i*pi*alpha^2*x/4).                 (WR3)
```

The Abel-triangle prefactor is `1/(4i*pi)`, positive dual modes have
multiplicity two, and `a_t(x)=x^(-7/4)(1-x)^(-1/4)` before the common
logistic phase.  Hence

```text
[2/(4i*pi)] a_t(x) Full_(B,n)
 =alpha[x(1-x)]^(-1/4)exp(i*pi*alpha^2*x/4),         (WR4)
```

which is exactly the equation-(9) source integrand.  Thus the finite block
in (WR2) recovers the `L=2481422` labels

```text
B-2,B-4,...,A.                                       (WR5)
```

There is no normalization or phase discrepancy.  The missing upper label
`B`, the dual `n=0` endpoint difference, and the direct `-1` in
`Theta(tau)-1` remain in one exceptional sector; this gate does not split
them.

The remainder in (WR2) can be summed at equal dual index.  Put

```text
s=p^2/x>=1,
tau_x(s)=x(1-1/s),
g_D(x,s)=exp[i*pi*D^2*x(1-s)/4]
           {2+i*pi*D^2*x(1-s)}.                   (WR6)
```

Then, in the inherited cutoff/Abel sense,

```text
2sum_(n>=1)[Tail_(A,n)-Tail_(B,n)]
 =e^(i*pi/4)x int_1^infinity s^(-1/2)
   [g_A(x,s)-g_B(x,s)]
   [Theta(tau_x(s))-1] ds.                           (WR7)
```

This restores the endpoint difference before any absolute value.  At the
new face,

```text
g_A(x,1)-g_B(x,1)=0,
partial_s(g_A-g_B)(x,1)
  =(3i*pi/2)(B^2-A^2)x.                              (WR8)
```

So the apparent termwise tail boundary is absent from the summed object.
The remaining challenge is quantitative: control (WR7) jointly in `s` and
`x` and close the exceptional zero/direct sector against the existing
face/fold atlas.

Pi provenance: the `pi` factors are inherited from equation (9), the Jacobi
modular factor, and the Weber boundary value.  Equation (WR4) explicitly
shows their cancellation into the original source normalization.

Proof boundary: exact source recovery for labels (WR5), same-index tail
reassembly, and face cancellation (WR8) only.  No estimate for (WR7), no
identification of the exceptional sector with the missing `B` source term,
no complete source-minus-target estimate, `T_upper`, height-uniform theorem,
`Lambda<=0`, PF-infinity, RH, and any prize-level conclusion remain open.
