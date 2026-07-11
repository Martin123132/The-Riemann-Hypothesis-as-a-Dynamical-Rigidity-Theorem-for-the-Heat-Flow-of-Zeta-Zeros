# Jensen-Window PF Newman Root External-Field Lemma

Date: 2026-07-11

Status: exact Newman root external-field reduction with sign
countermodels. This is not a proof of RH or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_newman_root_external_field_lemma.json
python work/rh_compute/scripts/jensen_window_pf_newman_root_external_field_lemma.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_root_external_field_lemma.py
```

Current result:

```text
validated Jensen-window PF Newman root external-field lemma: 10 rows, 0 issues, 5 exact canonical-product identities, 1 pair-flow reduction, 1 gap-stiffness expansion, 2 sign countermodels, 1 cosine equilibrium benchmark, 1 open Xi-balance handoff
```

## Squared-Zero Field

At a real-zero Newman time, write the positive zeros as `x_j` and
suppose `x` is a double zero. In `s=-z^2`,

```text
At a real-zero time, F_lambda(s)=F_lambda(0)*product_(j>=1)(1+s/x_j(lambda)^2), with sum_j x_j^(-2)<infinity.
F_*(s)=(s-rho)^2*U(s), rho=-x^2<0
E_x=U'(rho)/U(rho)=sum_(j!=*) m_j/(x_j^2-x^2)
K_s=-E_x'(rho)=sum_(j!=*) m_j/(x_j^2-x^2)^2>0
```

The second-order Jensen correction from the boundary-layer theorem
therefore depends on a concrete principal root field rather than an
unspecified higher-degree quantity.

## Newman Variables

Returning to the signed zeros of `H`,

```text
If H_*(z)=(z-x)^2*V(z), then B_x=V'(x)/V(x)=1/x-2*x*E_x and H_*'''(x)/H_*''(x)=3*B_x.
K_x=sum_(signed y outside pair)1/(x-y)^2=1/(2*x^2)+2*E_x+4*x^2*K_s>0.
```

For a pair with center `c`, gap `g`, and `q=g^2`, the zero ODE gives

```text
c'=sum_y[1/(x_+-y)+1/(x_--y)] -> 2*B_x
g'=4/g-2*g*sum_y 1/((x_+-y)*(x_--y))
q'=8-4*q*S(c,q), S=sum_y 1/((c-y)^2-q/4)
g(t)^2=8*(t-t_*)-16*K_x*(t-t_*)^2+O((t-t_*)^3)
```

Thus the universal `8*(t-t_*)` square-root birth receives a negative
second-order correction from the strictly positive stiffness `K_x`.

## Sign Guard

Two finite positive-coefficient Laguerre-Polya products at `rho=-1` are

```text
(1+s)^2*(1+s/2): E=1, collision D^-2 coefficient=-3/64
(1+s)^2*(1+2*s): E=-2, collision D^-2 coefficient=9/64
```

So neither LP membership nor coefficient positivity fixes the field
sign. Any usable estimate must exploit Xi-specific zero balance.

## Equilibrium Benchmark

For `F(s)=cosh(sqrt(s))^2`, every squared root is double and

```text
E_x=1/(2*x_k^2), B_x=1/x_k-2*x_k*E_x=0.
```

This is the exact arithmetic-lattice equilibrium. The live target is a
uniform comparison between the Xi field and this paired cancellation,
including the far zero tail and possible escape to infinite height.

Primary zero-dynamics anchor: https://arxiv.org/abs/1801.05914
