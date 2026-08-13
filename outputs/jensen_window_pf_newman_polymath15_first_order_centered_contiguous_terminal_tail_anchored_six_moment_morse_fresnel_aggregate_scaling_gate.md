# Six-Moment Morse-Fresnel Aggregate-Scaling Gate

Date: 2026-08-02

Status: exact scaling gate with one absolute-transition route obstruction,
six explicit outside-roster integral-tail bounds, `0 enumerated transition
modes`, `0 signed flow bounds`, and this is not a proof of RH.

## Physical Domain

Use q=1, L>=50, t=1/(2L^2), a^2=exp(L)+1/(32L^2), a=N+theta with 0<=theta<=1, B=N-M-1<a, h=1/a, and alpha=xi/(2pi) on T_0-epsilon<=xi<=T_0. Then alpha>a^2/2>22.

## Lower-Saddle Collar Geometry

Write `s=-log(v)`. The exact geometry is

```text
P(s)=2*((s - 1)*exp(s) + 1)/(1 - exp(s))**2,

Q(s)=-(-2*((s - 1)*exp(s) + 1)*(exp(s) - 1)**2*exp(s) + (exp(s) - 1)**4)/(exp(s) - 1)**5.
```

The identity controlling `P<=1` is

```text
2*(-s + sinh(s))*exp(s)=2 exp(s)[sinh(s)-s]>=0.
```

For `0<=s<=1/10`, exact Taylor-tail arithmetic gives

```text
Q(s)>=16/25,

Q-floor reserve=8364091/29232000,

tail allowance=2187/81200000.
```

## Absolute Transition Aggregation Guard

R_alpha={r integer: ceil(alpha*exp(-1/10))<=r<=floor(alpha)}; then u_0=alpha/r and 0<=log(u_0)<=1/10

partial_y b_(0,r)=(A_0(u)/r){Q(v)+P(v)[t log(u)/2-sigma]}, P=J^2/v, Q=J partial_v J

P<=1, Q>=16/25, sigma<=201/400, and A_0(u)>=1-sigma log(u)>=3799/4000 on the collar

Therefore

```text
partial_y b_(0,r)>=13/(100r),

c_(0,r)(y_1)>=13/(100r).
```

exp(-1/10)<10/11 and alpha>=22 imply #R_alpha>=alpha/22, hence sum_(r in R_alpha)1/r>=1/22

After the exact factor `1/(2pi)` in the transition majorant and `pi<22/7`,

```text
sum_(r in R_alpha) transition_majorant_(0,r)
 >91/96800.
```

At L>=50, a^2>2^50, so the sum of the j=0 termwise absolute transition majorants exceeds 10^12*h^2. This lower-bounds the chosen majorant, not the true signed residual.

This is a lower bound on the nonnegative *majorant being proposed*, not on
the signed residual `sum_r Q_(0,r)`. It retires termwise absolute aggregation
at bare-value `h^2` scale and positively points to grouped cancellation.

## Outside-Roster Integral Tail

For k>=2 and q>=1, zeta(k,q)<=q^(-k)+q^(1-k)/(k-1)<=2q^(1-k).

For x=alpha/u, a_x=x-m+1>=alpha/(2u) and b_x=n+1-x>=alpha_+>=alpha.

Hence

```text
Z_2<=6u/alpha, Z_3<=10u^2/alpha^2, Z_4<=18u^3/alpha^3.
```

and the exact Section 11.168 tail envelope reduces to

```text
The outside-roster integral remainder is at most [4pi^2 alpha]^(-1) integral_0^(log B){6|(D^2-D)A_j|+30|DA_j|+74|A_j|}d lambda.
```

On q=1, L>=50 and 0<=lambda<=log B<log a<L, g'=t lambda/2-sigma<=-49/100, so |A_j|<=lambda^j exp(-49lambda/100).

With |g'|<=51/100 and g''=t/2<=1/10000, |DA_j|<=[j lambda^(j-1)+(51/100)lambda^j]e^(-49lambda/100) and |(D^2-D)A_j|<=[j(j-1)lambda^(j-2)+(101j/50)lambda^(j-1)+(3851/5000)lambda^j]e^(-49lambda/100).

integral_0^infinity lambda^k exp(-49lambda/100)d lambda=k!(100/49)^(k+1).

alpha=a^2-epsilon/(2pi)>a^2/2, so 1/alpha<2h^2. Also pi>3, hence 1/(4pi^2)<1/36.

The resulting explicit constants are

```text
j=0: C_j=234803/1225, K_j=6, tail <6/alpha<=12*h^2
j=1: C_j=1145600/2401, K_j=14, tail <14/alpha<=28*h^2
j=2: C_j=232001200/117649, K_j=55, tail <55/alpha<=110*h^2
j=3: C_j=69600360000/5764801, K_j=336, tail <336/alpha<=672*h^2
j=4: C_j=27840144000000/282475249, K_j=2738, tail <2738/alpha<=5476*h^2
j=5: C_j=13920072000000000/13841287201, K_j=27936, tail <27936/alpha<=55872*h^2
```

These bounds cover only the absolutely convergent integral remainder after calB_j(B)-calB_j(1) is extracted. The two endpoint functionals remain exact and unbounded here.

## Route Decision

Keep the exact Morse-Fresnel chart, but retire summing its per-mode transition majorants by absolute values as a bare h^2-scale theorem. Group roster modes through Hermitian pairing, Poisson/Abel summation, or a smooth aggregate transform before absolute values. The outside-roster integral tail may be bounded absolutely because its six explicit bounds have h^2 scaling.

## Next Target

Construct one grouped roster theorem that keeps the incomplete-Fresnel transition exact while summing its residuals with phase cancellation. Compose calB_j(B)-calB_j(1), the physical endpoint E_j, c_xi, transpose terms, and the terminal recurrence before inserting the six values into the eight observations.

## Proof Boundary

This proves a fixed positive floor for the sum of one particular termwise absolute j=0 transition majorant, and explicit O(h^2) bounds for the six outside-roster integral tails. It does not lower-bound the true signed transition residual, bound either endpoint functional, include the physical observation/current coefficients, prove a grouped roster estimate, signed flow bound, Phi_B bound, contact exclusion, retained aggregate or Xi theorem, Q209, Lambda<=0, PF-infinity, RH, or a prize-level conclusion.

## Reproduce

```text
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_aggregate_scaling_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_aggregate_scaling_gate.py
```
