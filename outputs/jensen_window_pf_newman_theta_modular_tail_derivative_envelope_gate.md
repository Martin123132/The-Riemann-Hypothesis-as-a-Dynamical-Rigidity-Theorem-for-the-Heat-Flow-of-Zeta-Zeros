# Newman Theta Modular-Tail Derivative Envelope Gate

Date: 2026-07-24

Status: exact derivative envelope and effective tail-scale theorem.
This is not a proof of `Lambda<=0`, RH, or a Clay-prize result.

```text
work/rh_compute/results/jensen_window_pf_newman_theta_modular_tail_derivative_envelope_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_theta_modular_tail_derivative_envelope_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_theta_modular_tail_derivative_envelope_gate.py
```

## Explicit Heat-Derivative Envelope

For `v>=0`, define

```text
P_b(T,v)=b!*sum_(ell=0)^floor(b/2) (2v)^(b-2ell)*T^(b-ell)/(ell!*(b-2ell)!)
|partial_u^b exp(tu^2)|<=exp(Tu^2)P_b(T,|u|), 0<=t<=T
```

With `r_N=sum_(n>N)b_n`, put

```text
E_(N,k)(u;T)=exp(Tu^2)*sum_(b=0)^k binom(k,b)P_b(T,u)|r_N^(k-b)(u)|, u>=0
d_(N,0,m)(T)<=integral_0^infinity E_(N,m)(u;T)du
d_(N,1,m)(T)<=integral_0^infinity [u*E_(N,m)(u;T)+m*E_(N,m-1)(u;T)]du
```

The second formula is the exact identity
`D^m[u*f]=u*D^m[f]+m*D^(m-1)[f]`, followed by the
triangle inequality. Evenness turns half of the real-line norm
into the displayed positive-half-line integral.

## Reflected-Blend Phase

Set `y=exp(4u)` and `q=pi*n^2`. For `y>=sqrt(2)`,

```text
9*sinh(4u)^2>=9y^2/16 for y>=sqrt(2)
q/y+9y^2/16
y_*=(8q/9)^(1/3)
(3/2)*(9/8)^(1/3)*q^(2/3)
```

On `1<=y<=sqrt(2)`, the forward `q/y` term is stronger.
Consequently, for every `u>=0` and `n>=2`,

```text
pi*n^2/y+9*sinh(4u)^2>=c_0*n^(4/3), u>=0, n>=2
c_0=3.34638071167606735286846132651
```

Differentiating either blended summand only introduces fixed
polynomial factors. Splitting the exponential phase twice absorbs
those factors, the bounded Newman heat weight, integration, and
the arithmetic tail sum. Thus

```text
For fixed T,p,k there are effective C_(T,p,k)>0 and A_(p,k)>=0 such that integral_R |u|^p*exp(Tu^2)*|r_N^(k)(u)|du<=C_(T,p,k)*N^A*exp(-c_*N^(4/3)) for N>=2.
For every fixed j,m,T, d_(N,j,m)(T)<=C_(j,m,T)*N^A*exp(-c_*N^(4/3)).
c_*=0.836595177919016838217115331627
```

## Correct Absolute-Error Scale

The saddle transition still identifies which arithmetic blocks
matter spectrally, but the proved absolute derivative envelope has
a different scale:

```text
N_sad(x,t)=ceil(max(n_*(5,t,x),n_*(9,t,x)))=Theta(sqrt(x))
The proved absolute derivative estimate at N_sad+O(1) is only polynomial(x)*exp(-c*x^(2/3)).
N_K(x)=ceil(K*(1+x)^(3/4))
d_(N_K,j,m)(T)/x^m<=poly(x)*exp(-c_*K^(4/3)*x)
K>K_gamma=(pi/(8*c_*))^(3/4)
K_gamma=0.567098366759640997160953095281
N(x,t)=max(N_sad(x,t),ceil(K*(1+x)^(3/4)))
```

The saddle-scale statement is a non-promotion guard: a slower
available upper bound does not prove that the exact oscillatory
remainder is large. It does show that a fixed additive collar
does not itself supply the required absolute-error theorem.

## Open Separation Target

Instantiate effective constants for T=1/5 and m=5 or a higher optimized order, then prove the retained direct-C1 value-or-derivative separation on all transition cells.

The derivative envelope does not give a lower separation bound for (S_(N,t),S_(N,t)') and does not intervalize the adaptive transition cells.

No new `Q_j`, strict Laguerre theorem, `Lambda<=0`, RH, or
Clay-prize conclusion is claimed.
