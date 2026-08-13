# Double-Weber exact reconstruction of the finite source

Date: 2026-08-13

Status: exact fixed-height modular/Abel reconstruction; not a proof of a
quantitative source-minus-target estimate

The previous gate reduced all nonzero modular modes to the completed finite
B block plus the same-index theta tail

```text
R_tail=e^(i*pi/4)x int_1^infinity s^(-1/2)
 [g_A(x,s)-g_B(x,s)][Theta(x(1-1/s))-1]ds.           (DW1)
```

Apply Jacobi once more to the theta factor in (DW1).  For a transformed
mode `ell>=1`, put

```text
p=sqrt(s-1),
a=pi*D^2*x/4,
b=pi*ell^2/x.                                        (DW2)
```

Apart from the common modular measure, the endpoint-D integral is

```text
int_0^infinity [2-i*pi*D^2*x*p^2]
 exp(-i*a*p^2-i*b/p^2) dp.                           (DW3)
```

In the convergent half-plane use

```text
I0(a,b)=sqrt(pi)exp(-2sqrt(a*b))/(2sqrt(a)),
I2=-partial_a I0=I0[1/(2a)+sqrt(b/a)].               (DW4)
```

Substitution into (DW3) gives

```text
-2*i*pi*ell*e^(-i*pi/4)e^(-i*pi*D*ell)/sqrt(x).      (DW5)
```

It is not zero endpoint by endpoint.  It is independent of the chosen
endpoint because both `A` and `B` are odd, so
`e^(-i*pi*D*ell)=(-1)^ell`.  Their endpoint difference therefore vanishes
for every `ell>=1`.  The `ell=0` mode vanishes separately, without splitting
divergent pieces, from the Abel primitive

```text
d/dp [2p exp(-i*pi*D^2*x*p^2/4)].                   (DW6)
```

Therefore the common theta tail collapses exactly to its subtracted direct
term:

```text
2sum_(n>=1)[Tail_(A,n)-Tail_(B,n)]
 =Tail_(B,0)-Tail_(A,0).                             (DW7)
```

Adding the modular dual zero mode completes it:

```text
(B_0-A_0)+(Tail_(B,0)-Tail_(A,0))
 =Full_(B,0)-Full_(A,0).                             (DW8)
```

After the Abel prefactor, (DW8) supplies one half of the B endpoint minus
one half of the A endpoint.  The direct `-1` in `Theta-1` is exactly minus
the original Poisson zero mode `K_0`; these cancel before estimation.  The
retained half-Poisson endpoint current is

```text
K_H=1/2[f_x(0)+f_x(L)].                              (DW9)
```

Consequently

```text
K_H+1/2[f_x(L)-f_x(0)]=f_x(L),                       (DW10)
```

which is the missing upper source label `B`.  Together with the already
recovered labels `B-2,...,A`, this reconstructs exactly

```text
B,B-2,...,A                                           (DW11)
```

with the equation-(9) normalization and phase.

This is a strong structural closure: the endpoint-grouped Abel-theta
triangle, modular dual roster, Weber completions, theta tail, zero mode,
direct subtraction, and endpoint half-current form one exact involution.
It explains why attempting to budget those sectors separately kept hitting
large barriers.

It is not yet a small-error theorem.  Equation (DW10) reconstructs the
source itself.  The next step is to splice the exact target-carrier
subtraction into this representation: recover the ordinary and fold target
terms as local evaluations of the same Weber kernels and bound only the
difference between exact finite tails and those target carriers.

Pi provenance: all `pi` factors come from equation (9), the two Jacobi
transformations, and the Weber boundary formula.  The endpoint normalization
is replayed symbolically; no fitted or geometric constant is introduced.

Proof boundary: exact fixed-height source reconstruction and cancellation of
the modular continuation/zero/direct bookkeeping sectors only.  No small
source-minus-target bound, ordinary-carrier remainder, A-fold splice,
complete `T_upper`, height-uniform theorem, `Lambda<=0`, PF-infinity, RH, or
prize-level conclusion is established.
