# Jensen-Window PF Compound Order-Four Condensation Gate

Date: 2026-07-12

Status: exact contiguous order-four condensation coordinate with 317-row
lambda=-100 prefix and strict lower-order countermodel. This is not a
proof of PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_compound_order4_condensation_gate`.

```text
work/rh_compute/results/jensen_window_pf_compound_order4_condensation_gate.json
python work/rh_compute/scripts/jensen_window_pf_compound_order4_condensation_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order4_condensation_gate.py
```

Current result:

```text
validated Jensen-window PF compound order-four condensation gate: 8 rows, 0 issues, 2 exact identities, 1 exact sign equivalence, 317 positive lambda=-100 prefix margins, 1 strict lower-order countermodel, 1 forbidden promotion, 2 open handoffs
```

## Exact Condensation

Write `H_(m,n)=det[A_(n+i+j)]_(i,j=0..m-1)`. Desnanot-Jacobi gives

```text
H_(4,n)*H_(2,n+2)=H_(3,n)*H_(3,n+2)-H_(3,n+1)^2
```

The completed lower layers have `H_(2,n)<0` and `H_(3,n)<0`. Put
`T_n=-H_(3,n)>0`. Then

```text
H_(4,n)>0 iff T_(n+1)^2>T_n*T_(n+2).
```

In the order-three defect coordinate

```text
G_n=d_(n+2)^2-x_(n+2)^2*d_(n+1)*d_(n+3)>0
T_n=A_n^3*r_n^6*x_(n+1)^3*G_n
T_(n+1)^2/(T_n*T_(n+2))=G_(n+1)^2/(x_(n+3)^3*G_n*G_(n+2))
```

so the exact new target is

```text
G_(n+1)^2>x_(n+3)^3*G_n*G_(n+2) for every n>=0
```

## Repaired Prefix

At `lambda=-100`, 1024-bit Arb arithmetic certifies all `317` available
margins on `0<=n<=316`. The smallest lower bound is

```text
4.71789538689370302040913771438131967598144836682129515711914E-3 at n=316.
```

The largest checked margin is approximately `0.0263397201892780` at
`n=0`; the margins decrease to approximately `0.00471789538689370` at
the current collar. This is rigorous finite evidence, not a tail theorem.

The same balls prove the stronger finite curvature cap

```text
P_n=log(G_n*G_(n+2)/G_(n+1)^2)<2/(n+3)^2, 0<=n<=316.
```

For the analytic tail put `k=n+3>=320`. The completed scaled-defect
theorem and `s_319>251/500` give

```text
-3*log(x_k)>3*d_k>=753/(250*(2*k+1)).
```

Therefore the single sufficient tail estimate

```text
P_n<=4/(n+3)^2, n>=317,
```

would prove the order-four inequality. The rational comparison is exact:
after clearing denominators its numerator is `76466200>0` at `k=320`
and increases thereafter.

## Strict Countermodel

Take scaled defects

```text
s_1,...,s_5=['3/10', '2/5', '41/100', '49/100', '11/20']
d_1,...,d_5=['1/5', '4/25', '41/350', '49/450', '1/10']
x_1,...,x_5=['4/5', '21/25', '309/350', '401/450', '9/10']
cubic frontiers=['-319/15625', '-89697/10937500', '-15257191/2756250000', '-9791/2250000']
G_0,G_1,G_2=['1417/156250', '10943/76562500', '603553/236250000']
```

The contractions satisfy every strict pointwise ratio wall and increase;
the defects decrease; the scaled defects increase; all four cubic
frontiers are strictly negative; and every available order-two,
contiguous order-three, and arbitrary-column order-three determinant has
the required sign. Nevertheless

```text
G_1^2-x_3^3*G_0*G_2=-933340447356927/58618164062500000000<0,
H_(4,0)=-6608596712914764288/582076609134674072265625<0.
```

Thus complete lower compound layers do not imply order four. The next
analytic task is the all-index lambda=-100 G-gap inequality, followed by
a cooperative flow law and arbitrary-column transfer.

```text
outputs/jensen_window_pf_order3_noncontiguous_secant_transfer_lemma.md
outputs/jensen_window_pf_negative_lambda_m100_compound_order3_entry_certificate.md
outputs/jensen_window_pf_compound_order4_first_summand_curvature_bridge.md
outputs/signed_hankel_jensen_bridge_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```
