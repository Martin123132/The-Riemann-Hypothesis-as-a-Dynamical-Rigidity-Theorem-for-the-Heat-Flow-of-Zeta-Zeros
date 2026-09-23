# Full-source lattice interaction and coherent cubic heat reduction

22 September 2026. TWO HANDS NETWORK LTD. Private additive RH research.

## Result

Two deductions remove specific obstacles, without proving the favorable
signed bound or RH. All historical evidence is preserved.

1. A growing arithmetic lattice can be removed from the COMPLETE nonlinear
   response, with its interaction with all other terms explicitly paid.
2. Both angular modes of the complete source can be replaced by one coherent
   Gaussian-damped cubic expression at O_I,S((log N)^(-4)) cost, under the
   inherited complete-source second-moment bound. No amplitude is excluded.

The full original profile interval, atomic norm, two height windows, fixed
prime set, global Haar subtraction, and inherited exterior/core payments
remain in the analytic statements. The finite controls below use one prime
and one profile; they do not evaluate the original target.

## 1. The lattice's cross terms are now paid

Let U be the full source and V its exact terms at n=mc. With d=ceil(log(N+1)),
the original positive measure dnu and an arbitrary predetermined retained
mask, an elementary logarithmic-spacing estimate gives

    E_V = integral |V|^2 dnu
        = O_I,S(1/m + N^(1/3)/m^2).

The finite version, including all constants and the exact gate normalizer,
is equation (4) of DERIVATION.md. There is no random-phase substitution.

Writing D for the original complete nonlinear interaction,

    |D(U)-D(U-V)|
    <= (2 alpha M_1^* + 3 beta M_3^*)
       [sqrt(E_U E_V) + E_V/2].

This bounds the effect of deleting V from U, not merely D(V). In the
predecessor's beat N=k^5, the ENTIRE lattice n=c k^2, 1<=c<=k^3, is removable
at O_I,S(k^(-1)) cost when k is coprime to S. More generally, spacing
m/N^(1/6) tending to infinity suffices. The sparse-lattice bounds may not
be summed over arbitrarily many groups without paying their total cost.

## 2. One factor for both angular modes

For sigma>0 and Q_l(z)=|z|^2 exp(i l arg z), l=1,3, set

    T_1(z)=|z|^2 z/sqrt(sigma+|z|^2),
    T_3(z)=z^3/sqrt(sigma+|z|^2).

Both equal Q_l(z) times the SAME factor |z|/sqrt(sigma+|z|^2). The global
error is at most sigma/2, including z=0 and arbitrarily large amplitudes.
The combined integrand's pointwise sign is preserved. Its integral's sign
need not be: attenuation affects favorable and adverse contributions alike.

Choose sigma=a(x)^2/d^3. The Laplace integral for the reciprocal square root,
restricted to the explicitly paid band

    rho_S/(4 N d^8) <= tau <= 8 log(d)/sigma,

produces the cubic numerator

    alpha M_1 |U|^2 U + beta M_3 U^3

with the common weight exp[-tau(sigma+|U|^2)]/sqrt(pi tau).
If G_N^cub is its integral and C=alpha M_1^*+beta M_3^*, then

    |D(U)-G_N^cub|
    <= C { A_I/(2 d^4)
         + [2/(sqrt(pi)d^4) + 1/(8 sqrt(pi)d^8 log d)] E_U },

where A_I=integral_I a(x)/4 <9. Thus the comparison error is O_I,S(d^-4)
when the inherited E_U=O_I,S(1) is used. A higher moment or rare-amplitude
uniform-integrability assumption is not required.

The earlier valid asymmetric construction is retained in sections4-7 and
the SOURCE records. Section8 and the separate CUBIC records give the
stronger common-factor construction; nothing was rescored or overwritten.

## 3. Complete finite checks

Frozen controls: N=512, p=257, x=1, spacing16, BOTH original windows.
There are511 full source terms,32 lattice terms and479 complementary terms.
All gate, kernel, phase and nonlinear cross terms are retained. The
whole exponential prime multiplier is used with its global Haar mean.

Outward enclosures from the160-bit,131072-node coherent replay:

| Quantity | Window0 | Window1 |
| --- | --- | --- |
| Original full response | [0.00041395553,0.00041699554] | [0.00032871593,0.00033173594] |
| Coherent heat response | [0.00041196389,0.00041836390] | [0.00032680182,0.00033318183] |
| Heat response after its finite error payment | [0.00040976389,0.00042056390] | [0.00032460182,0.00033538183] |
| Mixed full-minus-parts response | [-0.00003401058,-0.00002825057] | [-0.00000141157,0.00000432844] |

The first mixed interval excludes zero; the second does not. The theorem
pays both without needing their sign. The finite lattice-energy upper bound
is intentionally coarse (about0.3856 and0.3834, versus measured enclosures
near0.0002422 and0.0002724). This is not a claim of sharp finite lattice cost.

None of these one-prime, one-profile values may be substituted for the
original huge fixed S or full profile integral, or compared with the
original0.00092027989 core threshold as if they were its evaluation.

## 4. Verification and scope

-669 named exact/scalar checks,42 independent checks and65 coherent checks.
-1612 integration/audit checks, including11 tamper-rejection cases.
-308 paired interval comparisons; all768 append-only source-cache blocks
  cover the complete windows at128/65536 and160/131072 precision/grid pairs.
-Eight new full-window builds: four preserved asymmetric and four coherent.
-The complete literal Q response agrees with authenticated predecessor
  enclosures. The two constructions retain byte-identical energy outputs.
-Cartesian derivative and independent Gaussian-quadrature controls use
  mpmath; source enclosures use Arb with explicit whole-interval midpoint
  errors. These are numerical implementation checks, not formal proofs.

The manifest and external SEAL_RECEIPT authenticate the package. A separate
external CLEAN_REPLAY record will report the actual clean-extraction replay;
this report does not anticipate its success. INPUTS binds8243 prior pins and
all copied dependencies. Preseal and final integrity receipts check them.

There were no source-calculation failures. A read-only bookmark parser issue
and the additive mathematical refinement are retained in ITERATIONS.md.
No provider calls, Forge experiments, GitHub updates or external publication.

## 5. Exact remaining obstacle

The signed Gaussian-damped cubic correlation STILL needs an arithmetic bound:

    limsup_N max_(both original windows) G_N^cub <0.00092027989,

with the complete0.0011250999 budget and all inherited payments unchanged.
Its Gaussian depends on the WHOLE source, so this is not an ordinary cubic
moment theorem. Factoring that weight, dropping mixed terms or expanding it
without a paid remainder would recreate the original gap.

The next research target is a bound for this joint weighted correlation,
uniform across the paid heat band. The lattice theorem can remove selected
arithmetic subpopulations first, but cannot dispose of the complete spectrum.
The positive shift also depends on x: any future transport differentiation
must retain its commutator terms. No such differentiation is used here.

Status: FULL_SOURCE_LATTICE_HEAT_REDUCTION_PROVED_SIGNED_BOUND_OPEN.
RH remains open. No new zeta count, crossing certificate or efficacy claim.
