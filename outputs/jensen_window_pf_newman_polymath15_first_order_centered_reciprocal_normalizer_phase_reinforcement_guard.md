# Newman Reciprocal Normalizer-Phase Reinforcement Guard

Date: 2026-07-28

Status: exact phase-lock and route-rejection gate. This is not a proof
of improved cancellation, `Lambda<=0`, PF-infinity, RH, or a Clay-prize
result.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_reciprocal_normalizer_phase_reinforcement_guard.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_reciprocal_normalizer_phase_reinforcement_guard.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_reciprocal_normalizer_phase_reinforcement_guard.py
```

Primary source: https://arxiv.org/abs/1904.12438

## Exact Normalizer

```text
Use the exact Polymath-15 normalizer M_0(s)=sqrt(2*pi)*s*(s-1)*pi^(-s/2)*exp[(s/2-1/2)Log(s/2)-s/2]/16 and M_t(s)=exp[t*alpha(s)^2/4]M_0(s), with the standard Log branch. On the real axis write phi_t=M_t(s)/|M_t(s)|=exp(i*beta_t) and gamma_t=conj(M_t(s))/M_t(s)=exp(-2i*beta_t).
```

## Pi Provenance

The pi in M_0 and Q(v)=exp(i[v*log(v/(2*pi))-v-pi/4]) is the same completed-zeta and stationary-phase normalization. The pi/4 is the negative-curvature stationary signature. No new geometric construction of pi is introduced.

## Zero-Time Phase

```text
For s=1/2-iT, T>0, direct evaluation of the exact M_0 phase gives gamma_0/Q(T)=exp(i*delta_0(T)), where delta_0(T)=(T/2)log(1+1/(4T^2))+(1/2)atan(1/(2T)).
Since log(1+u)<u and atan(u)<u for u>0, 0<delta_0(T)<3/(8T)=3/(4x). Thus the reciprocal stationary phase Q(T) is already conjugate-locked to the normalizer coefficient gamma_0, with O(1/x) angular error.
```

## Heat Transport

```text
Let omega=-Im(s_*)=T-t*Im(alpha)/2, h=omega-T>0, T_0=T+pi*t/8, a_0^2=T_0/(2*pi), and delta_a=log(a_0)-Re(alpha). For the continuous argument psi_t=arg(Q(omega)/gamma_t) chosen from t=0, psi_t=-delta_0(T)+integral_T^omega log(v/T_0)dv+2h*delta_a. This identity is exact.
On L>=50, 0<=tL<=25, one has 0<h<pi*t/8, 0<=delta_a<1/(4x), and log(T_0/T)<pi*t/(8T). Hence |psi_t+delta_0(T)|<pi^2*t^2/(32x)+pi*t/(16x)<3t/(8x), using pi<22/7 and t<=1/2. Therefore |psi_t|<15/(16x)<1/x.
```

## Conjugate Lock

```text
Combine the phase bound with the previous reciprocal amplitude logarithm eta, |eta|<13/(2x). The normalized leading reciprocal coefficient divided by the conjugate physical carrier is R=exp(eta+i*psi_t), and for x>120 one has |R-1|<8/x.
```

## Reinforcement Guard

```text
The lock is conjugating, not destructively phased. In the allowed real-carrier test z=1, the normalized primal plus reciprocal leading pair is z+R*conj(z)=1+R and has modulus greater than 2-8/x. Thus reciprocal self-duality cannot by itself yield a uniform pairwise cancellation gain; for one legitimate phase it asymptotically doubles the contribution.
For a general carrier z, the locked pair is z+R*conj(z)=2Re(z)+O(|z|/x). It can be small when z is nearly imaginary, but that is the original real-projection oscillation, not a new reciprocal power saving. Any useful gain must control the internal block phase or a genuinely subtractive completed expression.
The c=2 wall needs an exponent gain strictly larger than d_2=3133668399/48144906818. The relative lock R=1+O(1/x) changes a paired leading term only by a constant/conjugate projection and exponentially smaller corrections. It cannot supply an N^(-d_2-eta) gain by self-duality alone.
```

## Endpoint Scope

```text
The C_0+C_1/a endpoint and adjacent recurrence remain essential for a correct discrete transform, but their known role is to repair the cutoff transition and leading last-saddle mismatch. No existing endpoint identity reverses the conjugate lock uniformly across the active r_* block.
```

## Route Decision

```text
Reject raw primal-plus-reciprocal destructive interference as the promised source of the missing power. Retain the reciprocal coordinate as an exact organization of the functional-equation partner and endpoint, but do not expect a power gain from pairing equal leading amplitudes.
```

## Replacement Target

```text
The next cancellation calculation must transform the actual signed difference used in the zeta handoff, sum_(n<=N)(exp[t log(n)^2/4]-1)n^(-s_*)-sum_(n>N)n^(-s_*), before absolute values. Derive its endpoint-complete reciprocal kernel with first x derivative and determine whether the weighted-minus-unweighted tail has a sign-changing or vanishing leading symbol at r_*. If that symbol is also reinforcing, retire the reciprocal route and return to bilinear/additive-energy or direct Xi Abel-phase methods.
```

## Boundary

This artifact proves the exact zero-time normalizer/stationary phase ratio, its heat-transport identity and O(1/x) bound, the combined conjugate lock, and a reinforcement nonpromotion guard. It rejects only a power gain from raw reciprocal pairing. It does not evaluate the signed weighted-minus-unweighted tail kernel, prove cancellation, improve c_*, prove the Abel-scalar gap, close the inner degree, exclude contact, prove Lambda<=0, PF-infinity, RH, or a Clay-prize conclusion.
