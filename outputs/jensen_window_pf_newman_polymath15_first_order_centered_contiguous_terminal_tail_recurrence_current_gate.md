# Contiguous Terminal-Tail Recurrence and Current Gate

Date: 2026-07-31

Status: exact recurrence and Arb-uniform fixed-M leading current theorem. This
is not a proof of a growing-tail estimate, complete aggregate bound,
`Lambda<=0`, or RH.

## Exact Collapse

Fix p in [-1,1], theta=(1-p)/2, q=2tL^2=1, a=N+theta, and M>=0; take N->infinity with M fixed.

C_0(p)=sum_(m=0)^M(-1)^m R(p-2m-2)+(-1)^(M+1)C_0(p-2M-2)

Re[Q(p)r(p)^m]=(-1)^m Re[R(p-2m-2)], Q=exp[-pi*i*(p^2/2-p+3/8)], r=-exp(2*pi*i*p)

Therefore

```text
C_M(p)=Re[C_0(p)-Q*sum_(m=0)^M r^m]=(-1)^(M+1)Re C_0(p-2M-2)
D_M=-C_M'/(4*pi), N_M=C_M''/(16*pi^2)
P_M=[C_M*C_M''-(C_M')^2]/(16*pi^2)
```

This explains why the isolated `N-3` carrier can reverse the edge current
while the canonical edge plus `N-1,N-2,N-3` prefix remains clockwise.

## Uniform Curvature

f(t)=Re C_0(t)=cos(pi*(t^2/2+3/8))/(2*cos(pi*t)), with analytic continuation at half-integers

For M=0, t in [1,3], the imported 3072-box edge theorem gives D(t)>3/50.

For t>=3 and |cos(pi*t)|>6/(5t), direct phase minimization gives 4*cos(pi*t)^4*D >= pi^2*t^2*cos(pi*t)^2-pi^2/2-(pi/2)*sqrt(pi^2+cos(pi*t)^4). Using pi>3, pi<22/7, and t^2*cos(pi*t)^2>36/25 proves D>4/5.

Otherwise choose the nearest alpha in Z+1/2, write y=t-alpha, e=1/alpha, z=alpha*y. Then 0<e<=2/5 and |z|<1/2. With q_e(z)=(1+e^2*z/2)*sinc(pi*z*(1+e^2*z/2))/sinc(pi*e*z), f(alpha+y)=+/-q_e(z)/(2e) and D=(q_z^2-q*q_zz)/(4e^4).

arcsin(x)<=x/sqrt(1-x^2), pi>3, and t>=3 give |z|<=2/sqrt(21)+4/189<1/2.

A 192-bit Arb cover of 1024 rational
boxes proves

```text
q_z^2-q*q_zz>7/5
minimum lower enclosure: [1.423743769665932752336601416713784486545133562574235671638083294451386 +/- 2.00e-70]
```

Since `e<=2/5`, this gives `D>875/64` in every near-removal box. The
direct phase minimum gives `D>4/5` away from those boxes, while the imported
base cell gives `D>3/50`. Hence the global fixed-prefix margin is
`P_M<-3/8000`.

The 30-term sinc derivative series includes an explicit tail ball smaller
than `1e-50`; the argument bound is below `2` on the whole compact box.

## Theorem and Boundary

For every integer M>=0 and every p in [-1,1], D(2M+2-p)>3/50 and the fixed-M contiguous terminal-tail leading symbol satisfies P_M(p)<-3/8000.

The exact endpoint check is `P_M(1)=P_0(1)-M*(M+1)/16.`.

Every pi is inherited from the completed-zeta/Riemann-Siegel normalization and C_0 recurrence. No circle, polygon, fitted period, or plotted symmetry introduces a new constant.

This is not yet a complete aggregate theorem. The near-terminal limits are
proved with `M` fixed before `N -> infinity`; they do not justify a prefix
whose length grows with height. Prove a finite-height estimate uniform for a growing terminal prefix and join it to the nonterminal Abel/cross-current bulk, or prove the endpoint-complete division-free Abel gap directly.

No growing-prefix finite-height bound, complete cross-current estimate, Abel
gap, winding cap, contact exclusion, `Q209`, cofinal descendant theorem,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is claimed.

## Reproduce

```text
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_recurrence_current_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_recurrence_current_gate.py
```
