# From four profile points to four complete compact integrals

The remaining source-dependent inverse payment is now controlled over
the entire tested profile interval, not only at x=1/2. This is a finite
continuum certificate plus a uniform reduction theorem, NOT the final
unbounded favourable inequality or an RH proof.

## Four complete finite cases

For the SAME frozen numerator B_N and original nonlinear response Xi_N,
the calculation certifies Re(B_N(x) Xi_N(x))<0 at every

    1/200 <= x <= 2.

It also gives these outward-rounded upper bounds for

    integral_[1/200,2] Re(B_N(x) Xi_N(x)) dx/(4a(x)):

| N | original window centre | complete profile tiles | upper bound |
|---|---|---:|---:|
| 1024 | 1/4 | 660 | < -7.71e-8 |
| 1024 | 3/4 | 338 | < -4.27e-7 |
| 2048 | 1/4 | 3108 | < -1.00e-8 |
| 2048 | 3/4 | 338 | < -6.37e-7 |

The original height window is covered in every case. There are no
uncovered profile gaps. The certificates pay for the inverse correction,
the complete analytic source-expansion tail, profile variation, and
whole-cell quadrature error. No interpolation between sampled signs is
used. These are upper bounds, not measurements of the exact integrals.

They concern the retained FROZEN middle band before the external beta
multiplier. They do not certify the pre-freezing band after its phase
payment, the outside contribution J_out, profile endpoints, or RH counts.

## What made the continuum test possible

The coefficient family has an exact even-power expansion in its fixed
logarithmic coordinate y. Instead of repeating the full source integration
for each x, use the shared functions

    S_k(r)=sum_core y_n^(2k) n^(-1/2-it(r)).

Nine of their positive-weight squared norms, plus an explicitly paid
analytic tail, bound the source's squared norm for EVERY x in this
compact. Minkowski and Jensen are used under the original shared-height
measure. No arithmetic independence is assumed.

The finite tests use nine moments; that number is not extrapolated to
unbounded N. On any FIXED compact I inside (0,infinity), retaining

    L_N=ceil(log N/log log N)

basis terms makes the source error at most N^(-3/2+o(1)). After passage
through the ORIGINAL inverse kernel and multiplication by the inherited
divisor majorant, the full frozen-band substitution costs

    O_I(d^(3/4) N^(-1+o(1))) = o_I(1).

This is an analytic growing-basis reduction, not extrapolation from the
four cases. Each S_k still contains the original arithmetic n sum; this
does NOT claim that its arithmetic evaluation is free or polylogarithmic.
No useful explicit onset for the asymptotic rate is claimed.

## A sharper optional correction

Only the phase sector where the discarded inverse correction can raise
the response needs an adverse payment. With

    q=(-Re(B_N e))_+, m_q=Av(hq), E_q=Av(hq |C|^2),

the upper bound becomes

    Re(B_N Xi_N) <= L0 [Re(B_N Gamma_N)+m_q R(E_q/m_q)].

Zero mass means zero payment. The perspective m R(E/m) is increasing
in both arguments, proving that this is never worse than the previous
whole-window payment. It keeps phase/source dependence, rather than
replacing it by an independent average. The four certificates above use
the original, looser payment; this directional refinement is an analytic
option, not an extra claimed numerical certificate.

## What remains

We still need an arithmetic estimate on the signed numerator and these
gated source moments that stays within the logarithmic budget, together
with J_out and every inherited payment. The order remains: fix I, take
N large, then handle profile endpoints. Neither a compact finite pass
nor a negligible basis-truncation error proves that estimate.

The next proof-facing target is joint control of the signed divisor
sum and the shared positive-weight moments. The directional payment is
available where the whole-window one is too costly. Another isolated
profile scan would not settle this missing uniform arithmetic bound.

## Verification

- 4444 complete profile tiles at EACH of 128 and 192 bits; identical meshes.
- 147456 complete height cells per precision, each certifying all nine moments.
- All 52 compared source-value enclosures overlap.
- 5162 series, tail, norm and shared-measure checks; 6000 exact rational
  directional inequalities; 10 tamper-rejection guards.
- Two failures of the first profile-wrapper path (including its diagnostic)
  and one overstrict verifier attempt are retained. No failed attempt is
  credited as a passing certificate.
- Both source runs use FLINT/Arb, not independent special-function libraries.
  The analytic arguments are documented, not proof-assistant formalized.

See the delivery receipts for the fresh extraction replay, package seal,
historical integrity and process cleanup. No providers, Forge experiments,
GitHub updates or external publication were performed.
