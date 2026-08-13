# Proof Claim Ledger

Date: 2026-08-13

Status: claim-classification ledger. This is not a proof of RH or `Lambda <= 0`; it records which claims are exact, finite, diagnostic, open, rejected, or hygiene gates.

## Purpose

The proof programme now has many artifacts. This ledger keeps the main claims
separated by type so finite evidence and diagnostics cannot be silently
promoted into an all-order theorem.

Machine-readable ledger:

```text
work/rh_compute/results/proof_claim_ledger.json
```

Validator:

```text
python work/rh_compute/scripts/check_proof_claim_ledger.py
```

Current result:

```text
validated proof-claim ledger: 776 claims, 0 issues, 10 open theorem targets
```

The source-normalization corrigendum remains the baseline: the Polymath-15
coefficient `C_0(p)`, and therefore the retained endpoint `H_a`, is generally
complex.  Formal Core Section 11.146 replaces the historical real-`H_a`
endpoint rows.  Section 11.147 now derives the corrected division-free
endpoint and endpoint-carrier phase currents.  Its exact `q=1` witnesses
show that the terminal relative current is eventually negative at `p=0`
and positive at cutoff equality `p=1`; a uniform terminal orientation is
therefore false.  Section 11.148 composes that endpoint with the terminal
carrier before real projection.  The resulting `q=1` cofinal leading symbol
is `[A A''-(A')^2]/(16pi^2)`, and a 3072-box Arb cover proves it is below
`-3/8000` on the full saddle cell.  Section 11.149 now proves the missing
finite-height remainder inside the retained first-order `q=1` model:
`|a^2J_edge/(kappa T_0)^2-K_edge|<1/10000000` for `L>=50`, hence the
normalized model current is below `-3749/10000000`.  The omitted higher-order
Xi approximation error remains open.  Section 11.150 redoes the heat
majorants on the full box `0<tL<=25` and proves the same finite-height sign
throughout the critical range, including every `q>=1` point.  The
adjacent-cutoff projective splice is completed in Section 11.151.  Section
11.152 keeps the global approximation residual on the whole first jet and
shows that the boundary transfer requires only the established `C1`
homotopy.  Section 11.153 proves that an edge plus one isolated near-terminal
carrier can have either current sign, so termwise absorption is unavailable.
Section 11.154 restores the canonical consecutive order: for each fixed
finite `M`, the edge plus `N-1,...,N-M` collapses to one shifted real
`C_0` trace and has leading current below `-3/8000`.  Section 11.155 makes
this effective at finite height on `q=1`.  For
`1<=M<=floor(a^(1/18))-1`, exact carrier-ratio and four-coordinate bounds
give a current defect below `1/400`, hence the retained growing terminal
block has normalized current below `-1/400`.  The remaining nonterminal bulk,
its cross-current with this block, the omitted Xi remainder, and every RH-level
conclusion remain open.  Section 11.156 now performs the exact join without
expanding all carrier pairs.  Preserving the `u_N` shear, truncated Abel
summation reduces the bulk and its first jet to `F_B,R_B` and their first
derivatives, while four-coordinate polarization isolates one signed remainder
`Phi_B`.  The sharp retained-model obligation is `Phi_B<=1/400`; this bound is
not yet proved.  Section 11.157 closes the inputs to `Phi_B` on the five
truncated complex currents `H_0,H_1,H_2,D_0,D_1` and the certified tail jets.
The existing joined dyadic/odd-prefix formulas apply exactly at the bulk
cutoff `B`, but no signed five-current inequality is yet available.  Section
11.158 then removes the quadratic first Dirichlet correction exactly: the
five corrected currents are a linear image of correction-free moments
`G_0,...,G_4`.  Moments `G_1,...,G_4` have symmetric Mangoldt forms on the
balanced hyperbola `dm<=B`.  Their prime-edge minors are negative through
order four, so the live use of that symmetry is a signed endpoint-coupled
Type-I/II estimate, not kernel positivity.
Section 11.159 makes that target fully explicit.  The ten real components of
`G_0,...,G_4` enter `Phi_B` through four real observations, and
`h^2 Phi_B=g^T M g+ell_T^Tg` with `rank(M)<=4`.  A generic exact witness has
inertia `(2,2,6)`, so no coefficient-blind semidefinite promotion is
available.  On the simultaneous retained value/centered-scalar contact,
`Phi_B=-P_T>1/400`; the open arithmetic theorem must therefore exclude that
fibre source-specifically while proving the explicit rank-four inequality.
Section 11.160 expands the rank-four form into its original endpoint-linear,
Hermitian phase-difference, and transpose phase-sum carrier kernels.  The
transpose kernel loses `chi_N` exactly and the Hermitian kernel loses its real
part.  The physical terminal source further gives
`|chi_N+i h theta/2|<3h^2`, with both the terminal and shared-anchor
correction derivatives retained.  The ideal leading kernels are explicit,
but a two-carrier opposite-phase guard reverses the isolated-carrier sign.
The live theorem is therefore one joint endpoint-coupled Type-I/II or Vaughan
bound on both kernel families; no such signed bound is yet proved.
Section 11.161 identifies the exact source restriction behind those kernels.
After polarizing the common factor, `G_0,...,G_4` are the order-four jet of
one positive-amplitude Fourier polynomial on the frequencies `log n`.  Its
phase-difference function is Bochner-positive, while its phase-sum family is
not covered by that positivity.  The ideal bulk current reduces exactly to
`[f f''-(f')^2]/4+[u_(N,x)/2]fg`; hence it is nonpositive on the physical
zero fibre `f(0)=0`.  An exact positive-amplitude `n=1,2` half-turn guard has
current `log(2)^2/8>0`, so the actual phase, amplitudes, endpoint, and
correction remainder are indispensable.  The full signed bound remains open.
Section 11.162 now restores the exact physical `q=1` laws.  It proves
`Omega=-Im(s_*)=T_0-epsilon` with
`3/(8L^2x)<epsilon<7/(8L^2x)`, the strictly ordered profile
`A_n=A_a exp[B_a u_n+u_n^2/(8L^2)]` with `49/100<B_a<1/2`, and the
parity-fixed common phase `omega_c=(-1)^N eta`.  The five unnormalized
moments transfer from `Omega` to `T_0` with an explicit positive-mass error.
Ordered amplitudes sign every two-carrier pure Turan numerator, so the old
reverse-order witness is physically excluded.  On the ordinary fibre the
leading current is equivalently the branch-free complex-variance
discriminant `4P_bulk^(0)|M_0|^4=-D_hat`.  A robust three-carrier dyadic
guard still gives positive current throughout the certified finite-`L`
amplitude box when its phase is freed.  Thus the open theorem is now
specifically an estimate using the actual `Omega`, common anchor, endpoint,
and correction kernels; neither `D_hat>=0` nor the full signed bound is
claimed.
Section 11.163 proves that the canonical positive-mass triangle transfer is
too large by more than `4.9*10^8` at the required scale, so a signed transfer
or joint endpoint composition is compulsory.  Section 11.164 differentiates
the complete endpoint/Hermitian/transpose current before taking absolute
values and closes the exact flow on `G_0,...,G_5`.  Section 11.165 realifies
that flow as a rank-at-most-eight quadratic, separates the common phase, and
locates the three reciprocal saddle families.  Its Hermitian phase and flow
multiplier are `a^2log(l/k)` and `log(l/k)`, so the reciprocal diagonal is
killed exactly.  Section 11.166 now pairs `(k,l)` with `(l,k)`: the interior
Hermitian main becomes
`-mu[R sin(2pi a^2mu)+I cos(2pi a^2mu)]` on `k<l`.  The ideal real kernel
therefore gains a sine zero, and the offset coordinates `l=k+d` satisfy the
rigorous reciprocal-lattice defect bound
`dist(a^2log(1+d/k),Z)<=d dist(a^2/k,Z)+a^2d^2/(2k^2)`.  The imaginary
correction retains only the single `mu` factor, while exact guards show that
swap symmetry has no sign and integral phase need not kill that correction.
Section 11.167 now supplies the endpoint-complete finite Poisson identity for
all six bare logarithmic moments.  It restores both physical half-endpoints,
keeps the dual sum as a symmetric limit, and derives the negative-curvature
main with phase signature `e(-1/8)`.  Its exact transport law
`M_(j,r)'=iM_(j+1,r)+C_(j,r)` closes first-frequency residual control on the
same six value transforms; the saddle-motion term `C_(j,r)` cannot be
dropped.  The sharp dual main has explicit upper- and lower-cutoff jumps and
is only piecewise `C1`.  A sourced weighted B-process theorem is admissible
after dyadic localization, but its implicit constants and logarithmic loss do
not yet give the required uniform `h^2/400` remainder.  That explicit
remainder, smooth or exactly composed endpoint transitions, the
imaginary-coefficient bound, and signed offset-sum estimate remain open.
Section 11.168 replaces the moving sharp cutoff by one exact global
Morse-Fresnel chart on a fixed dual roster.  The full stationary main and both
half-mains are continuous limits of the same incomplete-Fresnel function,
and the transition residual has an explicit `1/(2pi)` variation majorant.
Two exact integrations by parts reassemble the outside-roster endpoint terms
as stable digamma/Hurwitz-zeta sums before absolute values, leaving an
explicit absolutely convergent tail envelope and no artificial dyadic
boundary.  The six physical transition-variation sums and six tail integrals
are not yet evaluated, so the required `h^2/400` estimate, imaginary-channel
bound, and signed flow theorem remain open.
Section 11.169 performs the first non-enumerative scale audit of that exact
representation.  On a lower-endpoint saddle collar, exact Morse geometry
gives `P<=1`, `Q>=16/25`, and a fixed positive floor `91/96800` for the sum
of the particular `j=0` termwise absolute transition majorants.  At `L>=50`
that chosen majorant exceeds `10^12 h^2`; this is not a lower bound on the
true signed residual, but it retires that absolute-aggregation route at the
required scale.  In contrast, elementary Hurwitz bounds and exact
logarithmic-amplitude moments give six outside-roster integral-tail bounds
`Tail_j<K_j/alpha<=2K_jh^2` with
`K=(6,14,55,336,2738,27936)`.  Both extracted endpoint functionals and the
grouped phase-cancelling roster theorem remain open.
Section 11.170 audits that grouped phase directly.  After the harmless
integer gauge, its first difference approaches zero at the lower endpoint
and its discrete curvature is comparable to `1/alpha`.  An exact Abel
identity gives phase partial sums below `3sqrt(alpha)`, while an actual
terminal subblock of the positive `c_(0,r)(y_1)` coefficients has certified
mass above `1287/[200000sqrt(alpha)]`, so coefficient-blind phase control
naturally stops at `h` scale.  Exact residual integration by parts then
reveals the decisive endpoint coherence: every lower boundary trace has
phase one and every upper trace has the same phase `e(alpha log B)`.  Those
traces must be recomposed with `E_j`, both `calB_j` functionals, and the
terminal recurrence before estimation.  No cancellation of that endpoint
package or endpoint-composed `h^2` theorem is claimed.
Section 11.171 performs that exact recomposition for all six values.  It
collects each bulk moment into an incomplete-Fresnel carrier, stable lower
and upper endpoint packages, a grouped oscillatory `c'` interior, and the
already bounded outside-roster tail.  A one-mode transfer identity proves
that the complete decomposition is invariant under admissible roster
enlargement, while the full-cotangent rewrite is retained only as an audit
because it separates poles that cancel in the stable form.  The real-part
identities prove `|mathcal L_0|>4999/10000` and
`|mathcal U_j|>49A_j(B)/100` for `0<=j<=5`: endpoint composition retains the
hard half-endpoints rather than making them an `h^2` error.  After fixed
normalization, the carrier/remainder split is inserted exactly into the
eight-observation quadratic, with terminal data entering through its
existing linear row.  The grouped `c'` estimate and the resulting signed
current bound remain open.
Section 11.172 now factors that residual through the actual observation image.
Only the eight base and first-derivative observations survive; every residual
component in the common kernel of `U_0^T` and `U_1^T` is exactly invisible to
the flow correction.  Linearity of the Morse-Fresnel representation collapses
all eight linear terms to one complex polynomial amplitude `P_lin` of degree
at most five, so no seventh logarithmic moment is introduced.  The quadratic
part is exactly four signed base/derivative pairings, and the already proved
outside-tail bounds project coefficientwise onto `P_lin`.  A degree-five
witness shows that this compression does not force a factor
`(lambda-log a)` or an endpoint zero.  The physical coefficient expansion,
the grouped `c'` estimate for that single amplitude, and scale control of the
four quadratic pairings remain open.
Section 11.173 expands that single amplitude in the physical terminal-centered
chart.  The exact identity `log N+u_N=log a` removes the large logarithms from
`R` and `R_x`, and every observation row factors through only the correction
channels `C` and `D`.  A recurrence gives all six centered coefficients.  Its
top coefficient is exactly `i rho_2(X_T+V_p)(s_*')^2`, so degree five vanishes
on the total-value fibre; the correction-free model is cubic and becomes
quadratic on simultaneous value/scalar contact.  The remaining Morse interior
is now one explicit off-saddle divided-difference amplitude with a removable
saddle value given by a second-order differential operator on `P_lin` at
`log(alpha/r)`.  No coefficient norm, grouped mode cancellation, quadratic
remainder bound, or signed flow estimate is claimed.
Section 11.174 now certifies the correction scale on the full critical chart.
Keeping the terminal-centered and all-carrier source estimates separate gives
`|rho_2|<h^2/16` and `|rho_(2,x)|<h^4/16` from their exact quotient laws.
Writing the original correction as `frak d` prevents a notation collision and
turns the constant channels into normalized ratios at `ell=log a`.  The six
resulting envelopes are
`|c_0-1|<4379h^2`, `|c_1|<h^2/4`, `|c_2|<h^2/16`,
`|d_0|<16893h^4`, `|d_1|<5h^4/16`, and `|d_2|<h^4/16`.
In particular `2/3<|c_0|<4/3`.  These are physical correction-channel
coefficient bounds; the physical observation rows still have to be propagated
through the exact `p_0,...,p_5` recurrence before a `P_lin` or grouped Morse
norm can be claimed.
Section 11.175 completes the next coefficient propagation without inventing
componentwise observation bounds.  The physical carrier rates obey
`|r_0|<h^3/3000`, `|r_1-i/2|<h^2/6000`, `|t_0|<h^2/16`,
`|t_1|<h^4/16`, and `|delta|<4h^2`; their quadratic combinations remain
within `h^4/15`, `17h^2/8`, and `h^2/3000` of the ideal cubic chart.  Five
row-sensitive envelopes are then applied only after the physical `alpha` and
`beta` rows are assembled.  Combining those joined errors with Section
11.174 gives explicit deviations for `p_0,...,p_3`, an explicit
correction-suppressed bound for `p_4`, and the exact top factor
`p_5=i rho_2(X_T+V_p)(s_*')^2` with
`|p_5|<h^2|X_T+V_p|/63`.  The retained observation gradient is not yet
numerically bounded.  Also, these coefficients multiply
`lambda-log(a)`, so the older `lambda`-monomial tail template requires an
exact basis translation before use.
Section 11.176 performs that translation and separates its real loss from a
coordinate artefact.  The exact binomial map produces centered tail weights
`Khat_k(ell)=sum_(j<=k)binom(k,j)K_j ell^(k-j)`; those weights are sharp for
the old coefficientwise absolute majorant.  On the physical chart, however,
`B_a=sigma-(t/2)log(a)` conjugates the Morse exponent to
`t x^2/4-B_a x`, so no explicit `log(a)` remains in either the first-order
amplitude or the removable second-order saddle operator.  Both the
off-saddle and saddle formulas reduce exactly to two joined contractions of
the physical `alpha` and `beta` rows against explicit four-vector kernels.
If `R=||alpha||_1+||beta||_1`, all departures from the ideal centered cubic
have raw outside-tail contribution below `h^2R/300000000000`.  Thus the
remaining linear theorem is the ideal cubic inside the grouped physical
Morse/endpoint phase sum; no numerical bound on `R`, grouped cubic estimate,
quadratic residual bound, or signed flow estimate is claimed.
Section 11.177 identifies the surviving ideal cubic as a finite nilpotent
jet rather than four unrelated coefficients.  The four-vector field `F_0`
obeys `F_0'=D_0F_0` with `D_0^3=0`, while the joined field
`Y=(F_0,xF_0)` obeys `Y'=Dhat Y` with `Dhat^4=0`.  Hence adjacent Poisson
modes are transported exactly by the terminating exponential
`E(d)=exp(-dDhat)`.  The logarithmic step variation is precisely the
normalized phase curvature,
`d_r-d_(r+1)=Delta^2G/alpha_P`, and on the certified terminal collar this
gives first and second ideal-field variations below
`9||Y_r||_1/(4alpha_P)` and `9||Y_r||_1/alpha_P^2`.  The removable saddle
operator is the explicit matrix polynomial `M(epsilon_g)`; its ideal core
has determinant `12^(-8)`, so nilpotent transport supplies no universal
saddle zero or sign.  The complete physical saddle/Fresnel sum and coherent
endpoint package have not yet been bounded.  The next live step is to place
this exact transport inside endpoint-complete first- or second-order Abel
summation before taking absolute values.

Section 11.178 carries that transport through the centered Morse exponential,
the `1/r` coefficient, the complete off-saddle kernel, and the removable
saddle.  At fixed Fresnel coordinate, one exact block cocycle gives both
first and second mode differences.  Finite Fubini reindexing turns the moving
integration limits into one contiguous active roster.  The resulting
two-variable phase has an exact reciprocal critical point
`(r,u)=(alpha_P/q,q)`, critical value `alpha_P log q`, and Hessian determinant
`-1`; its sequential Fresnel curvatures are exact reciprocals.  Thus the
mode saddle returns the original physical `q`-carrier phase.  Each reciprocal
family has local scalar scale `1/(2alpha_P)`, but a deliberately nonnegative
constant-channel budget over all families exceeds
`200[B^(199/400)-1]/(2587alpha_P)`.  This is a route guard only, not a lower
bound for the signed remainder.  It proves that local second-order Abel gains
must be followed by cancellation across reciprocal families and recomposition
with the retained carrier and endpoint packages; no grouped `h^2` estimate is
yet claimed.

Section 11.179 partitions the reciprocal roster into exact half-open cells
`C_q={r:q-1/2<=alpha_P/r<q+1/2}`, with every half-integer tie transferred to
the upper cell.  Its positive Morse chart has the same analytic Jacobian as
the original negative chart.  For an arbitrary twice-differentiable amplitude
`A`, the inner saddle second jet and the reciprocal leading-amplitude second
jet are both governed by
`mathcal D_qA=q^2A''+2qA'+A/6`, after the reciprocal Jacobian is applied.  The
negative Gaussian moment contributes `+kappa/2` and the positive one
contributes `-kappa/2`, so the first full-line sequential stationary
coefficient cancels exactly.  For the physical weighted polynomial this is
the removable-saddle operator already identified in Section 11.173.  This is
not yet a finite discrete-cell estimate: integer boundary transfers,
incomplete Fresnel tails, physical endpoint packages, and signed `q`-family
summation remain open, and no finite-cell `O(alpha_P^(-2))` bound is claimed.

Section 11.180 restores the real-mode endpoint phases before applying finite
Poisson a second time.  The complete mode
`I_A(r)=integral_1^B A(u)e(alpha_Plog u-ru)du` is exactly the Fourier
transform of the zero-extended physical amplitude; at integer `r` it is the
rejoined `mathcal F`, `mathcal J`, and two coherent roster traces.  Finite
Poisson on the half-open reciprocal cell gives an exact tie term and dual
integrals.  Full-line Fourier inversion returns the starred physical carrier
exactly, so a finite cell differs from it only by its tie, noncentral aliases,
and a symmetric two-sided exterior.  The physical cells `1<=q<=B` tile one
band inside the fixed roster.  Their internal ties, aliases, and cellwise
exteriors recompose, leaving only the exterior of
`(alpha_P/(B+1/2),2alpha_P]`.  Summing that exterior over the dual index gives
exactly the complementary integer-mode package plus the two outer tie halves,
so a second scalar Poisson transform is an involution rather than a source of
smallness.  No signed bound for the joined physical exterior/complement or
the remaining endpoint/terminal package is proved.

Section 11.181 aligns that exact inversion with the terminal decomposition.
The terminal observation is split additively into the genuine endpoint edge
and the finite carriers `B<n<=N`; rotating those carriers with the bulk gives
an all-carrier homotopy whose physical value is exactly `h^2P_ret`.  It still
closes on six moments, uses the same rank-eight flow matrix, and leaves only
the genuine edge in the external linear row.  Extending finite Poisson from
`[1,B]` to `[1,N]` transfers the old half endpoint at `B` into the terminal
interval, places the new half at `N`, and makes the cells `B+1<=q<=N` tile
exactly the terminal reciprocal strip.  They return the omitted terminal
samples with their correct Dirichlet weights.  The two remaining outer mode
families are uniformly nonstationary, but no edge-composed bound for them or
signed full-current theorem is yet proved.

Sections 11.182--11.187 resolve the outer endpoint kernel, specialize the
terminal obstruction, preserve its complex conditional quadrature, complete
the terminal/bulk phases, compress the endpoint-resolved relative lifts, and
rejoin the lower and interior Abel terms.  The resulting exact distinction is
between the reciprocal-band sum `mathscr F_N`, the starred returned carrier
`mathscr M_N`, and their defect `mathscr D_N=mathscr F_N-mathscr M_N`.
The Hermitian relative remainder is `-mathscr D_N`; the transpose remainder
also retains its explicit terminal quadrature.  These are algebraic
reductions, not signed defect estimates.

Section 11.188 replaces the conditional exterior presentation of that defect
by one finite Dirichlet-kernel integral.  Centered cells have zero kernel mean,
leaving an exact endpoint harmonic and source-coupled variation; consequently
the carrier-only join contains two starred carriers.  A second, source-adapted
partition gives a finite physical/reciprocal block matrix.  Its diagonal
contains `u=alpha_P/r`, blocks with `|p-q|>=2` have an explicit phase-gradient
gap, and adjacent blocks have no uniform gap and retain their oriented tie and
terminal recurrence.  Exact Fourier witnesses rule out coefficient-blind
defect smallness.  The signed ideal-cubic diagonal-plus-adjacent composition,
and therefore every completed-current and RH-level conclusion, remains open.

Section 11.189 performs the first ideal-cubic collar audit and finds a precise
global obstruction.  The diagonal and two adjacent strips form the
radius-one collar, but an internal reciprocal tie moves one physical rail out
and the opposite rail in.  The far complement has the opposite jump; only
their sum is the transported band functional.  A simple exact witness proves
that no proper fixed-width collar is coefficientwise invariant.  The far
cells nevertheless compress to two oriented intervals per frequency cell,
so their integration-by-parts boundaries telescope before estimation.  The
diagonal and adjacent phases have exact logarithmic and corner charts, while
finite Fresnel second moments expose the endpoint terms absent from the
full-line first-correction cancellation.  The revised target is the complete
tie-invariant ideal Hermitian/transpose join, including both far rails and the
transpose terminal survivor; no signed aggregate has yet been proved.

Section 11.190 puts the sequential positive and negative Morse transforms in
one exact two-variable chart.  The phase separates as
`alpha_Plog p+(eta^2-y^2)/2`, but the finite physical collar maps to a
curvilinear domain: its limits have the nonzero velocity
`y_U'=(U/p)J(rho)/J(Urho/p)`.  The paired first correction is therefore an
exact boundary flux whose physical rails contain `y_U+eta y_U'`.  The second
term is a moving-endpoint shear omitted by a Cartesian four-endpoint match.
With the nonconstant ideal amplitude, a further interior transport
`eta Gamma_eta-y Gamma_y` survives.  Oriented diagonal/adjacent interfaces,
integer ties, ideal rail values, and the terminal Hermitian null/transpose
survivor are all retained.  The next target is the weighted flux joined to
the rail-compressed far functional and terminal quadrature; no signed ideal
aggregate or RH-level conclusion follows from the geometry alone.

Section 11.191 computes the surviving weighted transport exactly.  Two
removable logarithmic profiles reduce it to `P` and `P'+gP`; both
coefficients vanish at the double saddle, but only to first order.  On the
constant-phase carrier ridge the transport is generally nonzero and reverses
orientation across the saddle, so it cannot be treated as an independent
oscillatory error.  Hyperbolic coordinates turn the phase into `sd` and the
transport into `d partial_s+s partial_d`.  A fixed transverse slice then has
the exact removable kernel `kappa{e(sd_+)-e(sd_-)}/s`.  The ideal derivative
audit also shows that the Hermitian terminal value null leaves derivative
transport.  The next target is to recompose the ridge with both returned
carriers and the transpose terminal block before estimating the transverse
residual; no signed ideal aggregate or RH-level theorem is claimed.

Section 11.192 fixes that ridge coefficient by returning to exact finite
Poisson inversion.  Per cell,
`2M_q-S_q=M_q-tau_q+E_q^(ext)-Lambda_q^(alias)`; after the physical cells
recompose, `2mathscr M_N-mathscr F_N` is exactly
`mathscr M_N-tau_band+E_band`.  Thus one, not two, starred carriers remains,
with the endpoint half-weights unchanged.  In the hyperbolic chart the
finite diagonal splits into a lifted `d=0` trace, its transverse residual,
and outer `s` wings, while the full carrier additionally requires the
two-sided exterior.  An exact piecewise-BV witness at `alpha_P=q=2` has
finite diagonal modulus at most `8/15` but returned carrier `1`, forcing an
exterior of modulus at least `7/15`.  This retires coefficient-blind
finite-diagonal-to-carrier matching and shows that a ridge identity would
require a separate residual-plus-wings-plus-exterior cancellation.  The
complete global exterior/carrier/terminal package still has no signed ideal
bound, so no RH-level theorem is claimed.

## Categories

```text
exact_lemma:
  exact noncircular identities already available for proof development

asymptotic_theorem:
  validated epsilon-dependent theorems beyond a quantified asymptotic threshold

finite_certificate:
  promoted finite interval or manifest-backed certificates

interval_certificate:
  promoted interval-backed theorems on a complete real interval or half-line

diagnostic:
  finite necessary-condition or theorem-search diagnostics

algebraic_reindexing:
  exact algebraic translations such as Toeplitz/Jacobi-Trudi

countermodel_gate:
  executable proof-safety obstructions

theorem_target:
  open bridge theorem needed before any proof promotion

forbidden_promotion:
  invalid proof step rejected by countermodel or logic

hygiene_gate:
  reproducibility, language, status, and reference integrity checks
```

## Current Open Targets

```text
target_direct_coefficient_pf:
  prove all Toeplitz minors of c_k(0)=mu_{2k}(0)/(2k)! are nonnegative;
  this is exactly equivalent to the all-degree/all-shift Jensen-window target
  and to the Stieltjes moment property of a_r=p_(r+1)

target_signed_hankel_jensen_bridge:
  the proposed all-order signed-Hankel/deep-rectangle antecedent is false at
  order ten; identify a weaker Xi/Phi-specific condition that is actually
  satisfied and prove it implies all-degree/all-shift Jensen hyperbolicity,
  Laguerre-Polya membership, or exclusion of positive Newman birth
  specification: outputs/signed_hankel_jensen_bridge_target.md

jensen_window_pf_compound_order11_compact_adaptive_h23_certificate:
  interval theorem on all 129280 quarter blocks over 5700<=t<=38020;
  largest scaled upper 2122.8651669101159 remains below the target 6000

jensen_window_pf_compound_order11_first_summand_curvature_certificate:
  exact four-range composition proving y_1''(t)<=6000/t^2 for every t>=1252;
  zero full-kernel, heat-forward, or RH claims occur in this source theorem

jensen_window_pf_compound_order11_m100_entry_certificate:
  exact fixed-order endpoint composition proving Q_(11,n)(-100)>0 for every
  integer n>=0

jensen_window_pf_compound_order11_lambda0_completion_certificate:
  exact lambda-zero completion proving contiguous and arbitrary-column signed
  Hankel positivity through order eleven only; order twelve and PF-infinity
  remain open

jensen_window_pf_negative_lambda_first_summand_power14_rebalanced_dominance_extension:
  asymptotic theorem proving the complete-to-first-summand defect is below
  2/k^14 for every integer k>=380; not an order-twelve curvature theorem

jensen_window_pf_compound_order12_m100_partial_prefix_certificate:
  interval theorem proving four negative endpoint shifts followed by 1237
  positive shifts through n=1240, with no inconclusive row

jensen_window_pf_compound_order12_m100_endpoint_completion_certificate:
  interval theorem extending coefficient coverage through A_1514, reproducing
  all 1243 inherited Q11 balls by overlap, and proving the remaining 252
  finite endpoint shifts positive through n=1492

jensen_window_pf_compound_order12_lambda0_prefix_certificate:
  interval theorem proving Q_(12,n)(0)>0 at n=0,1,2,3 from fresh 520-digit
  12 by 12 determinant rebuilds

jensen_window_pf_compound_order12_high_cumulant_coarse_corridor:
  exact epsilon-ten vanishing and Cauchy residual bounds supplying normalized
  cumulant orders 21 and 22 on the finite and asymptotic saddle rays

jensen_window_pf_compound_order12_nested_curvature_finite_ray_certificate:
  17,999-block interval theorem proving t^2*v_1''(t)<1109.799<8000 on
  the complete finite saddle ray 2001/1000<=u<=20

jensen_window_pf_compound_order12_nested_curvature_asymptotic_ray_certificate:
  normalized asymptotic theorem proving t^2*v_1''(t)<4573.334<8000 for
  every saddle mode u>=20

jensen_window_pf_compound_order12_sparse_h23_two_unit_propagation:
  exact normalized Taylor transport from H0-H23 anchors through H0-H16 at
  distance at most two under a certified H24 wall; it certifies no source row

jensen_window_pf_compound_order12_compact_refined_h0_h23_cache:
  resumable window-16 H0-H23 compact source with 940/8085 exact rows through
  t=9448, a matching live rebuild at row 939, and next target t=9452

jensen_window_pf_compound_order12_sparse_h23_lower_bridge_prefix:
  rigorous 125-segment, 8000-quarter-block prefix proving the scaled curvature
  target on every real t in 1503..3503; segment 125 resumes at t=3503

jensen_window_pf_compound_order12_curvature_bridge_target:
  open fixed-order target with an exact transfer below 100/k^2 and proved
  finite/asymptotic saddle rays; only the physical lower/compact cover remains

target_jensen_window_pf_bridge:
  prove that every binomially weighted Jensen window
  B^{d,n,0}_j=binom(d,j) A_{n+j}(0) is finite PF-infinity, equivalently prove
  all-order PF-infinity of c_k=A_k/k! or the Stieltjes moment property of
  a_r=p_(r+1), from a weaker Xi/Phi-specific structure that tolerates the
  certified negative order-ten endpoint minors
  specification: outputs/jensen_window_pf_bridge_target.md

target_schur_positive_specialization:
  construct a positive Schur/Edrei-Thoma specialization h_k -> d_k(0)
  specification: outputs/positive_schur_specialization_target.md

target_positive_determinant_integral:
  derive a positive determinant integral formula for every structurally
  nonzero Toeplitz minor

target_edrei_log_power_representation:
  prove that a_r=p_(r+1)=(-1)^r[z^r]H_0'/H_0 is a Stieltjes moment sequence,
  equivalently prove both all-order Hankel columns, a nonnegative S-fraction,
  the global upper-half-plane Phi Pick sign, Suzuki's cofinal all-t Fredholm
  gate, or a Phi-derived positive self-adjoint resolvent
  for the exact logarithmic-derivative ratio; generic backward
  Stieltjes-cone invariance is rejected

jensen_window_pf_newman_strict_laguerre_correlation_target:
  prove strict Fourier positivity, equivalently Wiener translate density,
  for the Xi first-correlation kernel uniformly on 0<t<=1/5

```

Each is explicitly recorded as `open_target`. The checker rejects theorem
targets that lack a current blocker or required upgrade.

## Order-Four Closure

```text
jensen_window_pf_lambda0_first_summand_dominance_transfer:
  covariance monotonicity transfers the lambda=-100 first-summand dominance
  theorem uniformly through lambda=0

jensen_window_pf_uniform_superpolynomial_first_summand_dominance:
  the inherited low/high split suppresses every fixed local higher-theta log
  difference faster than any inverse power, uniformly on -100<=lambda<=0

jensen_window_pf_uniform_first_summand_heat_tilt_asymptotic_theorem:
  O'Sullivan's suitable-multiplier theorem and the Lambert-W recurrence give
  the seven uniform first-summand heat-tilt difference estimates

target_compound_order4_forward_invariance:
  discharged: exact Newton interpolation preserves G_2->1, the universal
  determinant gives one uniform positive tail, and backward cooperative
  variation of constants proves H_(4,n)(lambda)>0 for every n>=0 and
  -100<=lambda<=0

jensen_window_pf_order4_noncontiguous_total_positivity_transfer:
  reversing every finite Hankel column block converts the completed signed
  contiguous layers through order four into positive initial minors; the
  Gasca-Pena criterion makes the block strictly totally positive and proves
  every arbitrary-column order-four sign. The same exact argument transfers
  contiguous signs to arbitrary columns at every fixed order
```

This closes the complete consecutive-row, arbitrary-column signed-Hankel
structure through order four. The next sections record the now-completed
fixed order-five, order-six, and order-seven layers.
Every all-shift compound order at least eight, PF-infinity, RH, and
`Lambda <= 0` remain open.

## Order-Five Closure

```text
jensen_window_pf_compound_order5_uniform_tail_flow_reduction:
  the compact-uniform first-summand expansion through order eleven and exact
  120-permutation determinant algebra give the positive universal term
  294912*G_2^10*h^10 and one eventual H_5 tail on -100<=lambda<=0

  the exact flow is cooperative over the completed H_4 layer, and the uniform
  tail reduces propagation to finite backward induction

jensen_window_pf_compound_order5_m100_prefix_certificate:
  exact factorization H_(5,n)=W_n*J_n with W_n>0 and 325 outward-rounded Arb
  coefficient balls prove J_n(-100)>0, a relative margin above 0.006269, and
  H_(5,n)(-100)>0 for every 0<=n<=316

jensen_window_pf_compound_order5_m100_tail_curvature_reduction:
  with k=n+4, J_n>0 is exactly C_n<-4log(x_k), where
  C_n=Delta^2log(F_n)-Delta^2log(d_(n+3)); the proved defect anchor and a
  coefficient-positive comparison show that C_n<=100/k^2 for k>=321 is a
  sufficient tail theorem

jensen_window_pf_compound_order5_first_summand_curvature_bridge:
  the nested stable coordinates have positive floors away from zero; an
  exact degree-52 coefficient-positive reserve proves
  |C_n-C_n^(1)|<=37/(n+4)^2, reducing the remaining 63-part budget to
  q_1''(t)<=60/t^2 for every real t>=320

jensen_window_pf_compound_order5_nested_curvature_compact_certificate:
  a common three-unit collar and 36 outward-rounded blocks built from the
  hashed 107452-tile cache prove q_1''(t)<=60/t^2 on 320<=t<=V'(2)

jensen_window_pf_compound_order5_nested_curvature_finite_ray_certificate:
  100 collar-extension tiles and 1850 exact-cumulant-corridor blocks prove
  the same bound for every mode 2<=u<=20

jensen_window_pf_compound_order5_nested_curvature_asymptotic_ray_certificate:
  normalized H-derivative boxes and an analytic stable-log interval at
  0<=1/t<=10^(-30) prove t^2q_1''(t)<9.159 for every mode u>=20

jensen_window_pf_compound_order5_m100_entry_certificate:
  the three continuous ranges cover every t>=320; the exact 37+63 transfer
  proves C_n<100/(n+4)^2 on n>=317, and the prior prefix proves
  H_(5,n)(-100)>0 for every n>=0

target_compound_order5_m100_entry:
  discharged as Theorem B6 by the preceding endpoint certificate

jensen_window_pf_compound_order5_uniform_heat_forward_invariance_certificate:
  the endpoint theorem, compact-uniform eventual tail, cooperative flow, and
  fixed-order Gasca-Pena transfer prove contiguous and arbitrary-column
  order five for every shift and every -100<=lambda<=0
```

The exact rational countermodel remains important: the lower layers alone do
not imply order five. The closure uses zeta-specific nested-curvature and heat
flow input. The next section records the corresponding order-six closure; no
all-order sign-regularity or Jensen/PF bridge follows from either fixed-order
theorem.

## Order-Six Closure

```text
jensen_window_pf_compound_order6_uniform_tail_flow_reduction:
  the heat-tilt expansion through order sixteen, superpolynomial higher-theta
  suppression, and exact 720-permutation algebra give the universal signed
  term 1132462080*G_2^15*h^15 and one compact-uniform eventual Q_6 tail

  Desnanot-Jacobi and adjacent Plucker identities make the signed order-six
  flow cooperative over the completed order-five cone

jensen_window_pf_compound_order6_m100_prefix_certificate:
  327 outward-rounded coefficient balls prove K_n>0 and
  Q_(6,n)(-100)>0 for every 0<=n<=316; the weakest relative lower bound is
  above 0.007809 at n=316

jensen_window_pf_compound_order6_m100_tail_curvature_reduction:
  the canonical identity H_(5,n)=A_(n+4)^5*exp(p(n+4)) reduces signed order
  six to D_n=5log(x_k)+P_k<0; the defect anchor shows that P_k<=320/k^2
  for k=n+5>=322 is sufficient

jensen_window_pf_negative_lambda_first_summand_power7_dominance_extension:
  exact monotonicity and endpoint interval gates sharpen the complete-to-first
  moment defect to 0<=delta_k<2/k^7 for every k>=316

jensen_window_pf_compound_order6_first_summand_curvature_bridge:
  the third stable gap satisfies min(S_j,S_j^(1))>=3/(2j); a degree-75
  coefficient-positive audit proves |P_k-P_k^(1)|<100/k^2 and reduces the
  remaining endpoint budget to p_1''(t)<=200/t^2 on t>=321

jensen_window_pf_compound_order6_high_cumulant_coarse_corridor:
  the exact epsilon-ten partition recurrence and unit-disk residual prove
  |kappa_r|q^(r/2-1)/(r-2)!<50000 for r=9,10 and every u>=2

jensen_window_pf_compound_order6_nested_curvature_compact_certificate:
  aligned H-derivative jets through order ten and 38 adaptive Arb blocks prove
  p_1''(t)<=200/t^2 on 321<=t<=V'(2)

jensen_window_pf_compound_order6_nested_curvature_finite_ray_certificate:
  a rigorous mode-two collar and 17,999 rational exact-corridor blocks prove
  the same bound for every saddle mode 2<=u<=20

jensen_window_pf_compound_order6_nested_curvature_asymptotic_ray_certificate:
  normalized H boxes, the high-cumulant corridor, analytic stable logarithms,
  and one dimensionless interval prove t^2p_1''(t)<22.769 for every u>=20

jensen_window_pf_compound_order6_m100_entry_certificate:
  the three continuous ranges cover every t>=321; exact tent integration and
  the full-kernel transfer give P_k<301/k^2<320/k^2 on k>=322, so the
  analytic tail and rigorous prefix prove Q_(6,n)(-100)>0 for every n>=0

target_compound_order6_m100_entry:
  discharged by the preceding endpoint certificate

jensen_window_pf_compound_order6_uniform_heat_forward_invariance_certificate:
  endpoint entry, the compact-uniform eventual tail, cooperative variation of
  constants, and the fixed-order Gasca-Pena transfer prove signed contiguous
  and arbitrary-column order six for every shift and -100<=lambda<=0
```

The order-six countermodel in the flow reduction still matters: completion of
the lower cone alone does not imply the next sign. The zeta-specific proof
uses a new stable logarithm and derivatives through order ten. The next
section records the corresponding order-seven closure. PF-infinity, the
all-order Jensen bridge, RH, and `Lambda <= 0` remain open.

## Order-Seven Closure

```text
jensen_window_pf_graded_kernel_vandermonde_all_order_lemma:
  mixed-kernel coefficient valuations and Cauchy-Binet force the first
  determinant degree D=binom(m,2) at every fixed order m; Vandermonde
  alternants give the positive signed coefficient
  2^D*prod_(j=1)^(m-1)j!*G_2^D

  arbitrary fixed-order saddle truncation, the all-order Lambert recurrence,
  and direct finite-difference control of the remainder prove: for every
  fixed m there exists N_m such that Q_(m,n)(lambda)>0 for n>=N_m uniformly
  on -100<=lambda<=0. The threshold may depend on m and is non-effective

jensen_window_pf_compound_order7_uniform_tail_flow_reduction:
  at m=7 the signed first term is 52183852646400*G_2^21*h^21. Exact
  condensation identifies Q_7>0 with strict log-concavity of the completed
  endpoint Q_6 sequence, and the heat flow is cooperative over order six

  all-shift Q_(7,n)(-100)>0 would therefore propagate through the full heat
  interval and then to arbitrary columns. An exact lower-cone countermodel
  proves that the new endpoint input cannot be inferred from orders <=6

jensen_window_pf_compound_order7_m100_prefix_certificate:
  stable Arb evaluation using A_0,...,A_326, including a twelve-coefficient
  retained-integral repair, proves Q_(7,n)(-100)>0 for 0<=n<=314. The weakest
  relative Q_6 log-concavity margin is above 9/1000 at n=314

jensen_window_pf_compound_order7_m100_tail_curvature_reduction:
  the fourth stable logarithm gives Q_(6,n)=A_(n+5)^6*exp(r(n+5)). For
  k=n+6, order-seven positivity is exactly 6log(x_k)+R_k<0 with
  R_k=r(k-1)-2r(k)+r(k+1)

  the endpoint defect anchor and an exact positive comparison reduce every
  missing n>=315 sign to R_k<=900/k^2 for k>=321

jensen_window_pf_negative_lambda_first_summand_power8_rebalanced_dominance_extension:
  the rebalanced split a(k)=log(k)/10 spends high-region exponential reserve
  to strengthen the low-region tilted probability. Fourteen strict Arb and
  derivative gates prove delta_k<2/k^8 for every k>=300

jensen_window_pf_compound_order7_first_summand_curvature_bridge:
  the completed order-six bounds and two finite endpoint margins prove
  min(T_j,T_j^(1))>=3/(2j) for j>=320. The inverse-eighth-power error survives
  all four stable logarithms, and a degree-102 audit with 103 positive
  coefficients proves |R_k-R_k^(1)|<262/k^2 for k>=321

  consequently r_1''(t)<=600/t^2 on t>=320 would give
  R_k<601/k^2+262/k^2=863/k^2<900/k^2 and close the endpoint tail

jensen_window_pf_compound_order7_shifted_jet_t320_t1000_certificate:
  exact-potential point jets and dimensionless H2-H21 collar remainders prove
  r_1''(t)<=600/t^2 continuously on 320<=t<=1000. All 186 rational blocks
  pass; the worst scaled upper is below 50.911 and the weakest T floor exceeds
  2.069

jensen_window_pf_compound_order7_nested_curvature_compact_certificate:
  aligned H2-H12 interval quadrature, a strict t+-5 collar, and four stable
  logarithms prove r_1''(t)<=600/t^2 continuously on 1000<=t<=V'(2). All 82
  adaptive rational mode blocks pass; the worst scaled upper is below 358.733
  and the weakest T floor exceeds 0.00010247

jensen_window_pf_compound_order7_high_cumulant_coarse_corridor:
  the exact epsilon-ten recurrence and unit-disk residual theorem prove
  normalized kappa_11 and kappa_12 caps below 14001 on 2<=u<=20 and below
  700001 on u>=20, with all 72 formal terms checked exactly

jensen_window_pf_compound_order7_nested_curvature_finite_ray_certificate:
  a mode-two exact collar and 17,999 rational outward-rounded blocks prove
  r_1''(t)<=600/t^2 for every saddle mode 2<=u<=20; the largest scaled upper
  is below 73.543

jensen_window_pf_compound_order7_nested_curvature_asymptotic_ray_certificate:
  normalized H2-H12 boxes, an exact stable-log defect majorant, and one
  interval over 0<=1/t<=10^(-30) prove t^2r_1''(t)<55.541<100 for u>=20

jensen_window_pf_compound_order7_m100_entry_certificate:
  the four continuous ranges cover every t>=320. Exact tent integration and
  the full-kernel transfer give R_k<601/k^2+262/k^2=863/k^2<900/k^2 on
  k>=321, so the analytic tail and rigorous prefix prove
  Q_(7,n)(-100)=-H_(7,n)(-100)>0 for every n>=0

jensen_window_pf_compound_order7_uniform_heat_forward_invariance_certificate:
  endpoint entry, the compact-uniform eventual tail, cooperative variation of
  constants, and the fixed-order Gasca-Pena transfer prove signed contiguous
  and arbitrary-column order seven for every shift and -100<=lambda<=0
```

Thus signed contiguous and arbitrary-column order seven are closed on the
complete heat interval. This does not iterate automatically: the lower-layer
countermodel gate remains valid, so order eight needs a new zeta-specific
endpoint theorem unless a genuinely uniform-in-order mechanism is found.
There is still no common threshold uniform in the order, PF-infinity theorem,
all-degree Jensen bridge, RH proof, or proof of `Lambda <= 0`.

## Fixed Order Eight

```text
jensen_window_pf_compound_order8_uniform_tail_flow_reduction:
  the all-fixed-order Vandermonde theorem gives the positive signed term
  33664847019245568000*G_2^28*h^28 and a compact-uniform eventual Q_8 tail

  signed condensation identifies Q_(8,n)>0 with strict log-concavity of the
  completed Q_7 sequence, while the affine-Hankel/Plucker identity gives a
  cooperative flow with off-diagonal coefficient
  (4n+58)Q_(7,n)/Q_(7,n+1)>0

  the exact rational sequence (1,1,1/2!,...,1/13!,1/87120000000) has every
  available signed contiguous minor through order seven positive but
  Q_(8,0)<0, so lower-cone promotion is impossible

jensen_window_pf_compound_order8_m100_prefix_certificate:
  stable 2048-bit Arb evaluation from 1,257 outward-rounded coefficient balls,
  including 930 retained-integral extension rows, proves every relative Q_7
  log-concavity margin and hence every Q_(8,n)(-100) sign positive for
  0<=n<=1242

  all 1,243 rows are rigorous interval statements; the weakest margin is the
  final row n=1242 and its lower bound exceeds 1/300

jensen_window_pf_compound_order8_m100_tail_curvature_reduction:
  exact fifth-stable-coordinate factorization reduces every n>=1243 endpoint
  sign to W_k<=4300/k^2 for k=n+7>=1250; a coefficient-positive shifted
  quadratic puts that ceiling strictly below the available log buffer

jensen_window_pf_negative_lambda_first_summand_power9_rebalanced_dominance_extension:
  exact retained-summand analysis proves the relative full-kernel moment error
  is below 2/k^9 for every k>=300

jensen_window_pf_compound_order8_first_summand_curvature_bridge:
  exact five-layer gap propagation and a degree-133 coefficient-positive
  inequality prove |W_k-W_k^(1)|<190/k^2 for k>=1250, leaving the continuous
  first-summand target s_1''(t)<=4000/t^2

jensen_window_pf_compound_order8_high_cumulant_coarse_corridor:
  the epsilon-ten formal cumulants 13 and 14 vanish exactly; Cauchy control of
  the existing unit-disk residual proves both normalized exact cumulants below
  one on every saddle mode u>=2

jensen_window_pf_compound_order8_shifted_jet_t699_t999_certificate:
  185 outward-rounded shifted-jet blocks prove s_1''(t)<=2000/t^2 on
  699<=t<=999

jensen_window_pf_compound_order8_nested_curvature_compact_certificate:
  aligned H2-H14 caches and 96 adaptive common-collar blocks prove
  s_1''(t)<=4000/t^2 on 999<=t<=V'(2), with worst scaled upper below 3994.552

jensen_window_pf_compound_order8_nested_curvature_finite_ray_certificate:
  an exact mode-two collar and 17,999 rational interval blocks prove the same
  ceiling for every 2<=u<=20; the worst scaled upper is below 356

jensen_window_pf_compound_order8_nested_curvature_asymptotic_ray_certificate:
  one uniform normalized-H interval over 0<=1/t<=10^(-30) proves
  t^2*s_1''(t)<134.49<200 for every u>=20

jensen_window_pf_compound_order8_m100_entry_certificate:
  continuous curvature, tent integration, and the full-kernel transfer give
  W_k<4001/k^2+190/k^2=4191/k^2<4300/k^2; the analytic tail and prefix prove
  Q_(8,n)(-100)=H_(8,n)(-100)>0 for every n>=0

jensen_window_pf_compound_order8_uniform_heat_forward_invariance_certificate:
  endpoint entry, the compact-uniform eventual tail, cooperative variation of
  constants, and fixed-order Gasca-Pena transfer prove signed contiguous and
  arbitrary-column order eight for every shift and -100<=lambda<=0
```

Thus signed contiguous and arbitrary-column order eight are closed on the
complete heat interval. The next section records the now-completed order-nine
chain and the arbitrary-order propagation reduction.

## Fixed Order Nine And All-Order Heat Reduction

```text
jensen_window_pf_compound_order9_uniform_tail_flow_reduction:
  the all-fixed-order Vandermonde theorem gives the positive leading term
  347485857744891213250560000*G_2^36*h^36 and an eventual Q_9 tail uniform
  on -100<=lambda<=0

  signed condensation identifies Q_(9,n)>0 with strict log-concavity of the
  completed Q_8 sequence; the exact cooperative coefficient is
  (4n+66)Q_(8,n)/Q_(8,n+1)>0

  an exact lower-cone countermodel has every available signed layer through
  order eight positive but Q_(9,0)<0

jensen_window_pf_compound_order9_m100_prefix_certificate:
  1,257 outward-rounded endpoint coefficients, including 38 rigorous
  retained-integral repairs, prove Q_(9,n)(-100)>0 for 0<=n<=1240; the
  weakest relative margin exceeds 1/250

jensen_window_pf_compound_order9_m100_tail_curvature_reduction:
  exact sixth-stable-coordinate algebra reduces the remaining endpoint tail
  to Y_k<=4900/k^2 for k>=1249

jensen_window_pf_compound_order9_first_summand_curvature_bridge:
  complete-to-first-summand propagation through six stable logarithms proves
  |Y_k-Y_k^(1)|<550/k^2 for k>=1251, using a degree-168 numerator with 169
  positive coefficients

jensen_window_pf_compound_order9_high_cumulant_coarse_corridor:
  exact fifteenth- and sixteenth-cumulant corridors supply the highest
  derivative inputs on every saddle mode u>=2

jensen_window_pf_compound_order9_shifted_point_h0_h8_cache:
  validates the hash-bound 8,929-row exact-point jet cache used by the
  localized interval proof

jensen_window_pf_compound_order9_localized_lower_bridge_certificate:
  279 root segments and 874 accepted adaptive blocks prove
  w_1''(t)<=4200/t^2 on 1250<=t<=5700, with largest scaled upper below
  4194.245

jensen_window_pf_compound_order9_nested_curvature_compact_certificate:
  108 adaptive H2-H16 blocks prove the same ceiling on 5700<=t<=V'(2), with
  largest scaled upper below 4199.186

jensen_window_pf_compound_order9_nested_curvature_finite_ray_certificate:
  17,999 rational mode blocks prove the ceiling for 2<=u<=20, with largest
  scaled upper below 2178.822

jensen_window_pf_compound_order9_nested_curvature_asymptotic_ray_certificate:
  one normalized interval proves t^2*w_1''(t)<324.906<500 for u>=20

jensen_window_pf_compound_order9_first_summand_curvature_certificate:
  exact range composition proves w_1''(t)<=4200/t^2 for every real t>=1250

jensen_window_pf_compound_order9_m100_finite_splice_certificate:
  retained-integral balls for A_1257,A_1258 prove the two missing endpoint
  signs n=1241,1242

jensen_window_pf_compound_order9_m100_entry_certificate:
  tent integration gives Y_k^(1)<4201/k^2; together with the transfer this is
  Y_k<4751/k^2<4900/k^2, so the analytic tail, finite splice, and prefix prove
  Q_(9,n)(-100)=H_(9,n)(-100)>0 for every n>=0

jensen_window_pf_compound_order9_uniform_heat_forward_invariance_certificate:
  endpoint entry, the eventual tail, cooperative variation of constants, and
  the initial-minor transfer prove signed contiguous and arbitrary-column
  order nine for every shift and -100<=lambda<=0

jensen_window_pf_all_order_endpoint_heat_reduction:
  arbitrary-order column/row cancellation proves
  Q_(m,n)'=(4n+8m-6)delta(Q_(m,n)); the flag-Plucker identity makes the flow
  cooperative over each completed lower layer

  for each fixed m, its own N_m confines propagation to finitely many shifts;
  ordinary induction in m therefore proves that the full heat-interval
  hierarchy is conditionally equivalent to the static endpoint antecedent
  Q_(m,n)(-100)>0 for every m>=10,n>=0

jensen_window_pf_endpoint_deep_schur_coordinate:
  after the necessary normalization h_k=A_k(-100)/A_0(-100), exact column
  reversal and Jacobi-Trudi give
  Q_(m,n)(-100)=A_0(-100)^m s_((n+m-1)^m)(h)

  arbitrary increasing Hankel columns are in bijection with the deep
  partitions lambda_1>=...>=lambda_m>=m-1. With the fixed-order
  initial-minor theorem, the candidate all-order rectangle hierarchy is
  equivalent to positivity on the whole deep cone

jensen_window_pf_endpoint_order10_counterexample:
  independent 4096-bit Arb Hankel, Jacobi-Trudi, and Toda checks prove
  Q_(10,n)(-100)<0 for n=0,1,2,3, equivalently
  s_((N^10))(h)<0 for N=9,10,11,12; Q_(10,4)>0, and the stable scan is
  positive for every 4<=n<=1240 with no inconclusive rows

jensen_window_pf_endpoint_pf3_boundary_counterexample:
  exact rational propagation of rigorous acb endpoint coefficient balls
  proves s_(1,1,1)(h)<-4.8484864218206096971e-11. Thus the actual endpoint
  sequence is not PF_3 or PF-infinity, but this failed shape is outside D_3
  because its smallest part is 1<2; the separate order-ten artifact supplies
  the counterexample inside the deep cone

jensen_window_pf_endpoint_deep_schur_literature_fit:
  Gasca-Pena supplies the already-used finite initial-minor transfer, and
  Pena's finite reversal result supports the orientation. Edrei is strictly
  stronger and contradicted by the excluded PF_3 minor; published eventual
  total positivity by matrix powers is a different coordinate. No direct
  closing theorem occurs in the bounded audited primary-source set

jensen_window_pf_deep_schur_rectangular_toda_coordinate:
  for tau_(m,N)=s_((N^m))(h), signed Desnanot-Jacobi gives
  tau_(m+1,N)tau_(m-1,N)=tau_(m,N)^2-tau_(m,N-1)tau_(m,N+1). Thus the next
  rectangle order is exactly strict width-log-concavity of the current row;
  the subtractive identity itself supplies no positivity direction

jensen_window_pf_deep_schur_moving_tail_boundary_counterexample:
  resetting b_k^(s)=h_(s+k)/h_s to zero for k<0 changes the Jacobi-Trudi
  boundary. Translation is valid only when the tail shape is itself deep;
  already r=2,s=1,mu=(0,0) has exact defect h_0*h_2/h_1^2>0. Therefore the
  deep cone does not become moving-tail PF by this shortcut

jensen_window_pf_strict_schur_jensen_counterexample:
  H(z)=exp(z/100)/((1-z)(1-2z)) is strictly Schur-positive for every
  partition, but its degree-three shift-zero Jensen polynomial has exact
  discriminant -222484532394597/2000000000000<0. Hence even full unweighted
  Schur/PF positivity cannot replace the Xi/Phi-specific Jensen bridge
```

Thus signed contiguous and arbitrary-column layers are now closed through
fixed order nine. A threshold uniform in `m` is not needed for the dynamical
induction. The equivalent deep rectangular endpoint hierarchy is now
rejected: its first required order has four negative initial shifts. Toda
rewrites this as failure of strict width-log-concavity, and the ordinary
moving-tail PF translation also fails independently at the reset boundary.
The surviving task is a weaker Xi/Phi-specific Jensen-window mechanism;
generic full Schur positivity is excluded as a sufficient bridge.
PF-infinity, RH, and `Lambda <= 0` are not proved.

## New Finite Diagnostic

```text
hankel_sign_consistency_reduction_point_audit:
  exact-rationalized cache point audit for the Grussler-Damm reshaped-Hankel
  finite condition; validates five lambdas, k=2..5, N=18; not an interval
  certificate or all-order sign-consistency proof

hankel_sign_consistency_reduction_finite_certificate:
  Arb/enclosure-backed finite certificate for the same reshaped-Hankel
  frontier; validates 689,795 finite minors for k=2..7, N=20; not an all-order
  sign-consistency theorem

shifted_hankel_sign_consistency_finite_certificate:
  Arb/enclosure-backed finite certificate for shifted reshaped-Hankel blocks;
  staircase manifest validates 3,154,515 finite minors for shifts n=0..20:
  k=2..5 at N=18, k=6 at N=16, k=7 at N=15, and k=8 at N=14; not an
  all-shift or all-order sign-consistency theorem

jensen_hankel_bridge_algebra_gate:
  exact degree-2 signed-Hankel/Jensen identity plus a positive rational
  degree-3 countermodel; blocks promotion from finite low-order reshaped signs
  to Jensen hyperbolicity

jensen_window_pf_obligation_algebra_gate:
  exact low-degree Jensen-window PF obligation algebra; degree 2 matches the
  signed-Hankel threshold, while degree 3 and degree 4 introduce additional
  banded Toeplitz obligations and finite low-order countermodel failures

arb_jensen_window_pf_obligation_diagnostic:
  Arb/enclosure-backed finite diagnostic for selected degree-3/4
  Jensen-window contiguous Toeplitz determinants; validates 1,470/1,470
  finite determinants for the five-lambda grid and shifts n=0..20

arb_jensen_window_sturm_hyperbolicity_diagnostic:
  Arb/Sturm finite diagnostic for selected degree-3/4 Jensen-window
  positive-root counts; validates 210/210 finite root-count rows for the
  five-lambda grid and shifts n=0..20

arb_jensen_window_sturm_d5_hyperbolicity_diagnostic:
  Arb/Sturm finite diagnostic for selected degree-5 Jensen-window
  positive-root counts; validates 105/105 finite root-count rows for the
  five-lambda grid and shifts n=0..20

arb_jensen_window_sturm_d6_d12_hyperbolicity_diagnostic:
  Arb/Sturm finite diagnostic for degree-6 through degree-12 Jensen-window
  positive-root counts; validates 735/735 finite rows on the five-lambda grid
  and shifts n=0..20, bringing the degree-3 through degree-12 total to
  1050/1050 with no failed or inconclusive rows

jensen_window_sturm_pf_finite_consequence:
  finite window-by-window PF consequence of the Arb/Sturm root-count
  manifests; validates 1050/1050 checked Jensen windows across degrees
  d=3..12, five lambdas, and shifts n=0..20

jensen_window_pf_bridge_obligation_ledger:
  theorem-obligation hygiene gate for target_jensen_window_pf_bridge;
  validates 11 exact, finite, open, conditional, rejected, and route-separated
  obligations, with 3 open obligations and finite rows blocked from closing
  the target

jensen_window_pf_theorem_machinery_fit_matrix:
  theorem-search hygiene gate for jwpf_06; validates 7 total-positivity,
  PF, zero-preserver, sign-regularity, downstream, and rejected-shortcut
  theorem-family rows, with 0 ready-to-apply rows and explicit fatal gaps

jensen_window_pf_sign_regular_transfer_gap_matrix:
  finite theorem-search diagnostic combining exact degree-2 signed-Hankel
  contact, degree-3/4 countermodel gates, the rejected all-order antecedent,
  and weaker-Xi/binomial/shift requirements into 9 transfer rows, with 0
  ready-to-apply rows

jensen_window_pf_factorial_multiplier_split_audit:
  finite theorem-search diagnostic separating the valid conditional
  Pólya-Schur factorial multiplier step from the false raw moment-window input
  route; validates 315 raw degree-2 anti-hyperbolic rows and 315 normalized
  degree-2 positive rows, with 0 ready-to-apply rows

jensen_window_pf_reciprocal_gamma_mixture_sign_gate:
  exact fixed-scale theorem and mixture-closure gate; Karlin's reciprocal-gamma
  matrix has the required strict signature in every order, but the exact
  independent-row scale integral has a sign-changing integrand and a positive
  three-atom mixture reverses the order-two sign. For Xi, the completed
  order-two cone is exactly the tilted concentration bound
  CV_n^2<=2/(2n+1); higher compound concentration remains open

jensen_window_pf_reciprocal_defect_compound_order3_gate:
  exact first higher-compound coordinate; the contiguous 3x3 sign is
  equivalent to unit-buffered reciprocal-defect log-concavity C_n>0. A strict
  rational sequence satisfies every currently proved ratio, adaptive-defect,
  and cubic Jensen cone while having C_n=-181/100 and the wrong positive
  Hankel sign. Separate anchored-entry and cooperative-flow theorems now close
  the actual Xi margin at every shift through lambda=0, and planar secant
  transfer closes every increasing column triple

jensen_window_pf_negative_lambda_m100_compound_order3_entry_certificate:
  interval-analytic all-shift contiguous order-three entry theorem; a 318-row
  Arb prefix, the anchor s_319>251/500, and all-k scaled-defect growth give
  C_n(-100)>57613471/66107054971 on the analytic tail and hence the required
  negative contiguous 3x3 Hankel sign for every shift. Forward invariance,
  noncontiguous order-three minors, and higher compounds remain open

jensen_window_pf_compound_order3_forward_invariance_certificate:
  exact cooperative-flow theorem; C_n'/r_(n+2)=alpha_n C_(n+1)+beta_n C_n
  with alpha_n>0, while the proved cubic spatial tail bounds the weighted
  infinite system. The resulting maximum principle gives C_n(lambda)>0 for
  every shift and finite lambda>=-100, hence D_(3,n)(0)<0. Noncontiguous
  order-three minors and compound order four and higher remain open

jensen_window_pf_order3_noncontiguous_secant_transfer_lemma:
  exact planar transfer theorem; positive column scaling turns three-row
  Hankel columns into points whose edge slopes strictly decrease by the
  contiguous theorem. Weighted secant averaging proves
  R_(3,n)(j_1,j_2,j_3)<0 for every shift and increasing column triple at
  lambda=0. Reshaped-Hankel orders two and three are complete; the downstream
  theorem now closes contiguous order four, while arbitrary-column order four
  is the first open compound layer

jensen_window_pf_compound_order4_condensation_gate:
  exact condensation frontier and promotion guard; contiguous order four is
  equivalent to G_(n+1)^2>x_(n+3)^3 G_n G_(n+2). Arb certifies 317 repaired
  lambda=-100 prefix margins and P_n<2/(n+3)^2; the downstream entry theorem
  now proves the sufficient tail P_n<=4/(n+3)^2 for n>=317. A strict rational
  sequence passes the ratio,
  adaptive-defect, cubic, and every available arbitrary-column order-three sign but has
  H_(4,0)<0. At this intermediate gate the flow and arbitrary-column transfer
  remain open; the downstream uniform-tail theorem closes the flow only

jensen_window_pf_compound_order4_first_summand_curvature_bridge:
  exact stable-gap reduction and perturbation theorem; the first-summand gap
  factors as d_1(t)^2*(1-exp(-J_1(t))), with J_1(t)>=1/(7t) for every
  real t>=319. A centered-difference Taylor bound localizes K_1 to same-point
  derivatives and explicit derivative envelopes. The downstream finite and
  analytic exact-corridor theorems now prove the ceiling globally; it leaves
  2/(5k^2) of margin, and the proved full-kernel perturbation fits exactly
  inside that budget

jensen_window_pf_compound_order4_localized_curvature_compact_certificate:
  rigorous compact real-parameter theorem; 107,452 adjacent Arb quadrature
  tiles through standardized moment order eight assemble into 1,073 positive
  localized-curvature blocks and prove K_1(t)<=7/(2t^2) for
  319<=t<=V'(2). The downstream exact-corridor theorems separately close the
  mode tail u>=2

jensen_window_pf_compound_order4_gaussian_cumulant_ray_target:
  exact tilted-Gaussian epsilon^6 algebra for kappa_2 through kappa_8 and the
  alternating factorial leading signature; explicit candidate cumulant
  corridors clear seven conditional full t+-2 collars from u=2 to u=20.
  Its formal and exact corridor obligations are now proved by downstream
  theorems; the continuum corridor-to-U(t) implication remains open

jensen_window_pf_compound_order4_formal_cumulant_corridor_certificate:
  rigorous formal-model interval theorem; 1,800,000 adjacent Arb blocks prove
  all seven epsilon-six formal cumulant corridors on 2<=u<=20, with weakest
  outward-rounded margin above 0.01843. Exact-density remainders remain open

jensen_window_pf_compound_order4_formal_cumulant_asymptotic_ray_certificate:
  exact coefficient-positive formal-model theorem on u>=20; seven buffered
  leading corridor gates, fourteen potential-jet sign gates, and the bound
  |R_r^[6]-F_r|<1/(20u) close the epsilon-six formal model for all u>=2.
  Exact-minus-formal central and tail errors remain open

jensen_window_pf_compound_order4_formal_cumulant_next_parity_certificate:
  exact epsilon-eight formal algebra; 42 epsilon-six coefficients are audited
  term for term and the first omitted scaled hierarchy is q^-3, q^-2, q^-1
  across orders 2-4, 5-6, 7-8. Exact-density errors remain open

jensen_window_pf_compound_order4_formal_cumulant_next_parity_finite_certificate:
  rigorous 1,800-block centered Arb Taylor theorem for all seven first omitted
  coefficient functions on 2<=u<=20, with a full-block seventh-derivative
  remainder preserving the cancellations lost by direct interval substitution

jensen_window_pf_compound_order4_formal_cumulant_next_parity_asymptotic_ray_certificate:
  exact coefficient-positive theorem on u>=20; fourteen leading buffer gates,
  four new order-nine/ten jet gates, and polynomial norm transfer complete the
  signed first omitted coefficient bounds globally on u>=2

jensen_window_pf_compound_order4_exact_cumulant_remainder_budget:
  exact theorem reduction through epsilon ten; scaled exact-minus-epsilon-ten
  errors below 9/1000 on 2<=u<=20 and below 1/(100u) on u>=20 suffice, with
  strict final corridor reserves. Downstream density theorems meet the budgets

jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_certificate:
  exact epsilon-ten subtraction layer; 56 epsilon-eight coefficients are
  audited and the next scaled hierarchy is q^-4, q^-3, q^-2

jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_finite_certificate:
  rigorous 3,600-block centered Arb Taylor theorem for all seven second-next
  coefficient functions and normalized potential jets through order twelve

jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_asymptotic_ray_certificate:
  exact leading-buffer and jet-transfer theorem completing every epsilon-nine/
  ten coefficient bound globally on u>=2

jensen_window_pf_compound_order4_exact_cumulant_complex_disk_contract:
  exact Gaussian-factored partition, logarithm, and Cauchy reduction; 70
  coefficient audits reduce all seven cumulant errors to two unit-disk
  partition residual targets without differentiating away cancellations

jensen_window_pf_compound_order4_exact_cumulant_formal_tail_certificate:
  exact degree-thirty formal-density and Gaussian-hazard theorem closing both
  adaptive formal tails inside the finite and asymptotic partition budgets

jensen_window_pf_compound_order4_exact_cumulant_exact_tail_certificate:
  exact curvature and monotonicity theorem closing both exact-density tails;
  the local standardized curvature ratio stays above 59319/100000

jensen_window_pf_compound_order4_exact_cumulant_partition_extension_finite_certificate:
  rigorous 5,400-block partition and 5,430-block shifted-jet theorem preserving
  the epsilon-eleven through epsilon-fourteen cancellations in 78 functions

jensen_window_pf_compound_order4_exact_cumulant_central_residual_certificate:
  exact Bell-15 and seventeenth-order potential-remainder theorem closing both
  central regimes. With all four tails, zero partition components remain open

jensen_window_pf_compound_order4_exact_cumulant_corridor_theorem:
  exact global theorem proving all seven alternating-factorial cumulant
  corridors for every u>=2. The downstream finite and analytic certificates
  now close the continuum localized-U implication

jensen_window_pf_compound_order4_exact_corridor_localized_curvature_finite_certificate:
  rigorous 20,700-block continuum theorem preserving correlated mode,
  curvature, Hurwitz-zeta, and exact-corridor dependence on 2<=u<=20; 41,400
  shifted-collar gates prove j_0>E_0 and t^2 U(t)<7/2 on every block

jensen_window_pf_compound_order4_exact_corridor_localized_curvature_ray_certificate:
  exact analytic u>=20 theorem; coefficient-positive geometry, midpoint
  Hurwitz-zeta bounds, normalized H-jet boxes, and the logarithmic-defect
  series give t^2 U(t)<3011223637/866377000<7/2. Together with the finite
  theorem this closes corridor-to-U globally on u>=2

target_compound_order4_first_summand_curvature_ceiling:
  discharged exact B4 theorem; compact curvature plus the global
  exact-corridor-to-U theorem prove K_1(t)<=7/(2t^2) for every real t>=319

jensen_window_pf_compound_order4_m100_entry_certificate:
  all-shift contiguous order-four entry theorem at lambda=-100; the repaired
  317-row prefix, global curvature theorem, exact tent transfer, and
  full-kernel perturbation prove H_(4,n)(-100)>0 for every n>=0. At this
  intermediate gate forward invariance and arbitrary-column order four remain
  open; the downstream compact-heat theorem closes forward invariance only

jensen_window_pf_compound_order4_forward_flow_reduction:
  exact affine-Hankel and Plucker flow theorem; the bounded stable margin obeys
  F_n'=alpha_n F_(n+1)+beta_n F_n with alpha_n>0. Uniform tail attainment is
  automatic, and forward propagation is reduced to an upper bound on
  beta_n+alpha_n(n+2)/(n+1) over every compact heat interval

arb_xi_lambda0_order4_prefix_certificate:
  rigorous 24,576-bit direct Arb Xi-series certificate; 507 outward-rounded
  positive coefficient balls prove every raw H_(4,n)(0) determinant and every
  stable margin strictly positive on 0<=n<=500, with zero inconclusive rows

jensen_window_pf_compound_order4_lambda0_eventual_positivity_certificate:
  exact Xi-asymptotic theorem and rigorous prefix; after removing the affine
  ratio factor, the normalized 4x4 determinant vanishes through h^5 and has
  universal first term 768*G_2^6*h^6. Hence H_(4,n)(0)>0 eventually, while
  Arb proves every shift 0<=n<=500. Making the asymptotic threshold explicit
  and closing any intervening finite gap remain open

jensen_window_pf_order_moment_transport_fit_gate:
  primary-source theorem-fit guard; the first summand has the formal
  Gamma-average shape required by the 2026 order-moment transport theorem,
  but Arb proves its transformed kernel increases at the origin and is not
  completely monotone. The theorem also supplies positive Hankel moment
  orientation, not reciprocal sign regularity, so direct promotion is closed

jensen_window_pf_structural_ansatz_matrix:
  structural proof-search hygiene gate for jwpf_06; validates 6 candidate,
  blocked, and rejected ansatz rows against exact degree-2/3/4 hard tests and
  the finite countermodel kill gate, with 0 ready-to-apply rows

jensen_window_pf_cauchy_binet_low_degree_scout:
  symbolic theorem-search diagnostic for the live Cauchy-Binet ansatz;
  validates 15 degree-2/3/4 hard formulas by adjacent-log-concavity ratio
  parametrization with nonnegative Bernstein coefficients, finds 0 kernel
  identities, and keeps the larger-minor countermodel warning active

jensen_window_pf_log_concavity_frontier_scout:
  symbolic frontier diagnostic extending the low-degree scout to larger
  contiguous Jensen-window Toeplitz minors; validates 14 rows and locates first
  Bernstein-certificate failures at degree 3 size 6 and degree 4 size 5, and
  first exact log-concave countermodel negatives at degree 3 size 8 and degree
  4 size 6

jensen_window_pf_ratio_condition_scout:
  ratio-condition theorem-search diagnostic; validates 7 candidate rows,
  rejecting adjacent log-concavity, decreasing ratio contractions,
  second-order log-concavity, and selected low-degree Bernstein positivity by
  exact countermodel, and contraction log-concavity by constructed positive
  extension

jensen_window_pf_contraction_log_concavity_scout:
  ratio-condition rejection diagnostic; validates a positive rational
  extension satisfying x2^2 >= x1*x3 while the degree 3 size 8 and degree 4
  size 6 contiguous Jensen-window minors remain negative

jensen_window_pf_schur_shape_contract:
  exact shape-contract diagnostic; validates the Jensen-window
  Jacobi-Trudi/Schur shape map on a bounded N=8, order<=5 grid, recording
  15,709 finite-band nonzero bounded shapes and 2 hard frontier column-shape
  obligations with mixed-sign h-monomial expansions

jensen_window_pf_column_recurrence_contract:
  column-shape recurrence diagnostic; validates the elementary-symmetric
  recurrence C_m=h_0^m e_m and confirms the degree 3 size 8 and degree 4 size
  6 hard frontier recurrence rows match the exact negative rational
  countermodel values

jensen_window_pf_column_recurrence_finite_coverage:
  finite coverage diagnostic; validates that checked zeta-window grids support
  the column recurrence target with 1,470 direct positive Arb determinant rows,
  210 hard recurrence rows, 315 checked Sturm/PF windows, and 12,600 positive
  recurrence-only stress rows, without promoting this to an all-order
  recurrence theorem

arb_jensen_window_column_recurrence_stress:
  finite Arb stress diagnostic; validates 12,600 positive Jensen-window column
  recurrence rows for degrees d=3..8, sizes m=1..20, five lambda values, and
  shifts n=0..20

jensen_window_pf_reciprocal_coefficient_extended_stress:
  extended reciprocal coefficient diagnostic; validates 72,600 positive
  normalized [t^m]1/H(-t) rows for degrees d=2..12, sizes m=1..40, five
  lambda values, and shifts n=0..32

jensen_window_pf_reciprocal_positivity_route_matrix:
  theorem-search hygiene gate for the column recurrence target
  [t^m]1/H(-t)>=0; validates 9 reciprocal-positivity route rows, with 0
  ready-to-apply rows, 3 live representation candidates, a rejected standard
  positive Stieltjes/Jacobi fraction route, rejected generic ratio shortcuts,
  and finite stress rows kept finite

jensen_window_pf_reciprocal_fraction_scout:
  reciprocal continued-fraction sign diagnostic; validates 3 symbolic rows
  and 735 finite Arb zeta-window sign rows, showing that the standard positive
  S-fraction/J-fraction route for E(t)=1/H(-t) has the wrong first nontrivial
  sign while leaving signed or modified fractions open

jensen_window_pf_reciprocal_signed_j_fraction_scout:
  reciprocal signed J-fraction diagnostic; validates 2 symbolic rows, 3,675
  finite signed reciprocal-Hankel determinant rows, and 2,940 finite
  signed-lambda rows, supporting the signed modified fraction target without
  proving the all-order theorem

jensen_window_pf_reciprocal_signed_jacobi_beta_scout:
  reciprocal signed Jacobi beta diagnostic; validates 3 symbolic rows and
  3,675 finite beta rows, with 2,940 positive rows, 630 negative beta_1 rows,
  and 105 terminal degree-2 zero-containing rows, sharpening the signed
  modified fraction target without proving the all-order theorem

jensen_window_pf_reciprocal_motzkin_path_obstruction_scout:
  reciprocal Motzkin path obstruction diagnostic; validates 3 symbolic rows,
  735 finite mu_2 cancellation rows, and 630 finite beta_1 diagonal obstruction
  rows, rejecting the raw ordinary J-fraction Motzkin path model as manifest
  positivity while leaving modified signed models open

jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_scout:
  reciprocal Motzkin parity-lift obstruction diagnostic; validates 3 symbolic
  rows and 5,145 finite same-length mixed-sign witness rows, rejecting global
  length-parity signs and diagonal sign conjugation as repairs of the raw
  ordinary J-fraction Motzkin path model while leaving state-space modified
  models open

jensen_window_pf_signed_j_fraction_theorem_target:
  theorem-target hygiene gate for the signed modified J-fraction route;
  validates 7 fit/misfit rows and states the missing theorem needed to convert
  all-order signed reciprocal-Hankel/Jacobi signature into coefficientwise
  nonnegativity of E(t)=1/H(-t), with 0 ready-to-apply rows

jensen_window_pf_modified_signed_model_target:
  modified signed-model hygiene gate; validates 9 model rows, rejecting the
  raw ordinary Motzkin/J-fraction model, diagonal sign conjugation, and global
  length-parity sign repairs while leaving 4 modified candidates live only
  conditionally, with 0 ready-to-apply rows

jensen_window_pf_state_space_sign_lift_obstruction_scout:
  derived finite obstruction diagnostic; validates 3 symbolic rows and 735
  mu_2 rows showing that an absolute-value sign-state cover of raw Motzkin
  paths overshoots the actual reciprocal coefficient by 2*kappa_1>0, while
  leaving genuinely modified state-space doubled models open

jensen_window_pf_oscillatory_resolvent_fit_matrix:
  oscillatory/resolvent theorem-search hygiene gate; validates 8 fit/misfit
  rows with 0 ready-to-apply rows, rejects ordinary entrywise Jacobi powers,
  diagonal similarity, absolute-value majorants, classical oscillatory
  spectral conclusions, indefinite moment language, and finite signed
  patterns as standalone coefficient-positivity proofs, and leaves only
  positive spectral-transform and Xi/Phi positive-kernel routes live
  conditionally

jensen_window_pf_positive_readout_theorem_target:
  positive-readout theorem-target hygiene gate; validates 8 candidate rows
  with 0 ready-to-apply rows and 2 live foundational routes, isolating the
  exact positive scalar readout obligation to either a noncircular positive
  spectral transform or an Xi/Phi-specific positive resolvent kernel while
  rejecting wrappers, endpoint factorization, finite quadrature, raw signed
  readouts, and absolute-value majorants as standalone proofs

jensen_window_pf_positive_spectral_moment_obstruction:
  positive spectral moment obstruction diagnostic; validates 3 symbolic rows
  and 735 finite Delta_2 obstruction rows, showing that the raw reciprocal
  coefficients mu_m cannot be ordinary power moments of a positive measure,
  while leaving nonordinary positive transforms and Xi/Phi positive kernels
  open

jensen_window_pf_nonordinary_positive_transform_ansatz_matrix:
  nonordinary positive-transform ansatz diagnostic; validates 8 ansatz rows
  with 0 ready-to-apply rows and 3 live ansatz rows, narrowing the surviving
  positive-readout search to non-power positive functionals, Xi/Phi positive
  kernels, or genuinely modified exact state-space transfer models

jensen_window_pf_nonpower_functional_low_degree_scout:
  nonpower functional low-degree diagnostic; validates 7 scout rows with
  0 ready-to-apply rows and 1 live contract row, recomputing reciprocal
  coefficient formulas through mu_6 and signed composition counts through
  m=8 to show that npt_04 needs a genuine positive cone, basis, and
  functional absorbing signed low-degree cancellations

jensen_window_pf_nonpower_functional_cone_candidate_matrix:
  nonpower functional cone-candidate diagnostic; validates 8 cone rows with
  0 ready-to-apply rows and 2 live cone rows, rejecting raw g-coordinate,
  standalone ratio/log-concavity, tautological, and endpoint PF/LP cones while
  leaving Xi/Phi kernel and Cauchy-Binet/determinant-integral cone routes live

jensen_window_pf_cauchy_binet_cone_frontier_matrix:
  Cauchy-Binet cone frontier diagnostic; validates 8 frontier rows with
  0 ready-to-apply rows and 2 live frontier rows, showing that the route
  cannot rest on selected low-degree Bernstein certificates or adjacent
  log-concavity and must instead construct a positive determinant integral
  for the hard column frontier or an all-shape Cauchy-Binet/Andreief kernel

jensen_window_pf_monotone_contraction_frontier_scout:
  monotone-contraction frontier diagnostic; validates 2 exact hard-row
  certificates with 88 positive Bernstein coefficients and 210 finite zeta
  diagnostic rows, showing that the first hard column-frontier polynomials are
  positive under x1 <= x2 <= x3 while the rational log-concavity countermodel
  violates that sharper condition

jensen_window_pf_monotone_contraction_column_extension_scout:
  bounded exact column-extension diagnostic; validates 25 degree-3/4/5 column
  rows with 3,329 positive Bernstein coefficients under monotone contractions,
  including 3 rows beyond the original first hard frontier and a degree-5 band
  through m=8

jensen_window_pf_monotone_contraction_sparse_degree6_scout:
  bounded exact sparse degree-6 diagnostic; validates 10 degree-6 column rows
  through m=10 with 63,347 strictly positive Bernstein coefficients, 0 zero
  coefficients, and 0 negative coefficients under monotone contractions

jensen_window_pf_monotone_contraction_sparse_degree7_frontier_scout:
  bounded exact sparse degree-7 frontier diagnostic; validates 9 degree-7
  column rows through m=9 with 670,891 strictly positive Bernstein
  coefficients, and records the first one-shot global Bernstein certificate
  obstruction at m=10 with 126 negative Bernstein coefficients and minimum
  -4928

jensen_window_pf_monotone_contraction_sparse_degree7_subdivision_scout:
  bounded exact degree-7 m=10 subdivision diagnostic; repairs the one-shot
  global Bernstein obstruction by splitting s0 into [0,1/2], [1/2,3/4], and
  [3/4,1], certifying 785,400 strictly positive slab Bernstein coefficients

jensen_window_pf_monotone_contraction_all_m_counterexample:
  exact countermodel gate; gives the shift-0 infinite static-cone witness
  x=(19/20,19/20,1,1,1,1,1,...) satisfying every pointwise wall and monotone
  contraction, while its normalized degree-7 m=11 column recurrence is
  negative; the propagated static cone alone is not an all-m theorem

jensen_window_pf_monotone_contraction_theorem_target:
  validated on the needed heat regime; the infinite maximum principle proves
  x_(k+1)>=x_k for every k and finite lambda>=-100. The former all-real-lambda
  statement remains false at lambda=-1156, and contractions alone remain
  insufficient for all-m column recurrence positivity

jensen_window_pf_heat_flow_cone_entry_asymptotic_target:
  validated infinite cone-propagation closure; full entry at lambda=-100 and
  the adjacent-defect maximum principle propagate both pointwise walls and
  x_(k+1)>=x_k to every finite lambda>=-100, including lambda=0

jensen_window_pf_phi_taylor_cone_entry_sign_scout:
  finite Taylor sign certificate; validates 4 Phi coefficient balls and the
  two local sign combinations 2*b-a^2<0 and
  2*(a^3-3*a*b+3*c)>0 needed by the fixed-k cone-entry asymptotic route,
  while leaving the uniform-in-k or collared finite cone-entry theorem open

jensen_window_pf_negative_lambda_cone_entry_prefix_scout:
  finite negative-lambda prefix certificate; validates the canonical 69-row
  ACB prefix through lower/upper ratio-cone walls k<=21 and monotone gaps
  k<=20, the extended k30 93-row prefix through lower/upper walls k<=29
  and monotone gaps k<=28, the extended k50 153-row prefix through
  lower/upper walls k<=49 and monotone gaps k<=48, the extended k60
  183-row prefix through lower/upper walls k<=59 and monotone gaps k<=58,
  the extended k80 243-row prefix through lower/upper walls k<=79 and
  monotone gaps k<=78, the extended k100 303-row prefix through
  lower/upper walls k<=99 and monotone gaps k<=98, the extended k150
  453-row prefix through lower/upper walls k<=149 and monotone gaps k<=148,
  and the extended k200 603-row prefix through lower/upper walls k<=199
  and monotone gaps k<=198 for lambda=-25,-50,-100,
  while leaving the all-k tail theorem or finite-collar flow theorem open

jensen_window_pf_negative_lambda_finite_collar_contract:
  finite-collar accounting diagnostic; extracts from the negative-lambda
  prefix and exact ratio-cone collar rule that the canonical usable active
  finite depth is K=19 with collars x_20,x_21, the extended k30 usable
  active finite depth is K=27 with collars x_28,x_29, the extended k50
  usable active finite depth is K=47 with collars x_48,x_49, the
  extended k60 usable active finite depth is K=57 with collars x_58,x_59,
  the extended k80 usable active finite depth is K=77 with collars
  x_78,x_79, the extended k100 usable active finite depth is K=97 with
  collars x_98,x_99, and the extended k150 usable active finite depth is
  K=147 with collars x_148,x_149, and the extended k200 usable active
  finite depth is K=197 with collars x_198,x_199

jensen_window_pf_negative_lambda_tail_barrier_scout:
  finite theorem-search diagnostic; rewrites the missing negative-lambda tail
  as defect inequalities, certifies 63 one-third-width buffer rows and 60
  defect-monotone rows on the canonical finite prefix, certifies 87
  one-third-width buffer rows and 84 defect-monotone rows on the extended k30
  prefix, certifies 147 cone-buffer rows and 144 defect-monotone rows on the
  extended k50 prefix, certifies 177 cone-buffer rows and 174 defect-monotone
  rows on the extended k60 prefix, certifies 237 cone-buffer rows and 234
  defect-monotone rows on the extended k80 prefix, certifies 297 cone-buffer
  rows and 294 defect-monotone rows on the extended k100 prefix, certifies
  447 cone-buffer rows and 444 defect-monotone rows on the extended k150
  prefix, certifies 597 cone-buffer rows and 594 defect-monotone rows on the
  extended k200 prefix, records the k200 one-third-width buffer frontier at
  179/597 rows, and rejects the
  scaled-defect nonincreasing shortcut without promoting the result to an
  all-k tail theorem

jensen_window_pf_negative_lambda_scaled_defect_frontier_scout:
  finite theorem-search diagnostic; validates 597 exact-cone rows on the k200
  prefix, records that the fixed half-width buffer holds on only 521/597 rows
  with first failures at lambda=-50 and k=191 and lambda=-25 and k=133,
  records that the one-third-width buffer holds on only 179/597 rows with
  first failure at lambda=-25 and k=31, and keeps scaled-defect nonincrease
  rejected without promoting any buffer into an all-k theorem

jensen_window_pf_negative_lambda_defect_recurrence_scout:
  finite theorem-search diagnostic; validates the buffered sufficient condition
  on 63 finite rows, validates defect monotonicity on 60 rows, and rejects the
  direct width-preserving recurrence on all 60 checked adjacent rows without
  promoting the result to an all-k tail theorem

jensen_window_pf_negative_lambda_log_curvature_bridge:
  finite theorem-search diagnostic and exact algebraic reduction; translates
  the buffered defect route into the sufficient conditions
  0<=B_k<=2/(3*(2*k+1)) and B_(k+1)<=B_k, validates 63 simple log-buffer rows
  and 60 curvature-monotone rows on the finite prefix, and identifies the
  monotone part with the open Delta^3 log A theorem target

jensen_window_pf_negative_lambda_bounded_log_curvature_target:
  historical finite diagnostic; records the formerly proposed fixed wall
  B_k=-Delta^2 log A_k<=2/(3*(2*k+1)) and its raw moment-curvature equivalent
  on the old k<=22 prefix, but is no longer a live theorem target after the
  repaired k300 obstruction

jensen_window_pf_negative_lambda_bounded_log_curvature_k300_obstruction:
  finite obstruction gate; rewrites the old wall as C_k=(2*k+1)*B_k<=2/3,
  validates that only 179/897 repaired k300 rows satisfy it while 718/897
  fail it with 0 inconclusive rows, and redirects the live route toward the
  zeta-specific raw corridor or linear curvature barrier

jensen_window_pf_negative_lambda_gaussian_curvature_matrix:
  finite theorem-search diagnostic and exact baseline comparison; identifies
  the bounded log-curvature target as a controlled positive deficit from the
  Gaussian raw-moment baseline, validates 63 positive-deficit and
  bounded-deficit rows, and rejects positive Gaussian scale-mixture arguments
  as an upper-wall proof template

jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix:
  finite theorem-search diagnostic; packages the surviving fixed-k signed
  Gaussian perturbation route, using the certified Phi Taylor signs to record
  positive leading Gaussian deficit and positive monotone correction while
  preserving the uniform-remainder gap

jensen_window_pf_negative_lambda_uniform_remainder_target:
  historical route diagnostic; records why the fixed-k signed-Gaussian route
  needed a two-scale theorem, but every destination of that route is now
  retired or independently closed by the lambda=-100 saddle/corridor chain,
  so it is no longer an active proof-programme blocker

jensen_window_pf_negative_lambda_taylor_moment_budget:
  finite theorem-search diagnostic; derives the exact Gaussian-moment
  normalization for the local q/T route, validates 7 k=22 tail-start samples,
  rejects the low-order Taylor truncation as a finite proof model for
  lambda=-25,-50,-100, and isolates the positivity and log-curvature Taylor
  remainder obligations

jensen_window_pf_negative_lambda_high_order_taylor_scout:
  finite theorem-search diagnostic; automatically generates local Phi Taylor
  polynomials, certifies 8 even coefficient rows through c14, and validates
  35 k=22 truncation rows showing invalid normalizers, upper-wall violations,
  and overbound rows across truncation degrees 6,8,10,12,14

jensen_window_pf_negative_lambda_defect_tail_theorem_target:
  historical simultaneous-target diagnostic. The original all-three-lambda
  statement remains unproved, while the lambda=-100 theorem now proves
  0<=d_k<=2/(2*k+1) and d_(k+1)<=d_k for every k>=1 and discharges the
  one-entry-parameter route

jensen_window_pf_negative_lambda_half_width_tail_target:
  finite diagnostic/rejected-route record; finite-rejects the fixed
  half-width scaled-defect route at k150: s_k<=1/2 holds on only 430/447
  checked rows, first fails at lambda=-25 and k=133, and must be replaced by
  an exact-cone or adaptive scaled-defect target plus the separate monotone
  defect bridge

jensen_window_pf_negative_lambda_adaptive_scaled_defect_target:
  historical simultaneous-target diagnostic. The lambda=-100 composition now
  proves the exact cone 0<=s_k<=1, decreasing defect, and increasing scaled
  defect for every k>=1; no claim is made for the stronger simultaneous
  lambda=-25,-50,-100 formulation

jensen_window_pf_negative_lambda_adaptive_envelope_matrix:
  finite theorem-search diagnostic; validates the k200 monotone-envelope
  pattern supporting the adaptive target: 594/594 adjacent k-increase rows,
  398/398 cross-lambda order rows, maximum checked scaled defect
  0.5376643171065356005 at lambda=-25 and k=199, and 76 finite half-width
  failures. It does not prove k-uniform or continuous-lambda monotonicity

jensen_window_pf_negative_lambda_adaptive_envelope_obligations:
  exact algebraic obligation diagnostic; validates 9 rows separating the
  adaptive target into exact lower-threshold input, the open upper wall
  x_k<=1, the open monotone bridge x_(k+1)>=x_k, the exact scaled
  k-monotone identity 2+(2*k+1)*x_k-(2*k+3)*x_(k+1)>=0, and the rejected
  fixed half-width/one-third buffers

jensen_window_pf_negative_lambda_raw_moment_bridge_matrix:
  exact finite theorem-search diagnostic; translates the adaptive envelope
  route into raw moment ratios R_k=M_(k+1)*M_(k-1)/M_k^2. It validates the
  exact identities x_k=((2*k-1)/(2*k+1))*R_k, 0<=s_k<=1 iff
  1<=R_k<=(2*k+1)/(2*k-1), and the adaptive R_(k+1) corridor, then checks
  597/597 raw-cone rows and 594/594 corridor rows on the k200 prefix while
  preserving the 76 half-width and 418 one-third fixed-buffer failures

jensen_window_pf_negative_lambda_raw_ratio_decrement_corridor_scout:
  exact finite theorem-search diagnostic; rewrites the raw adaptive corridor
  as the decrement recurrence
  2*(R_k-1)/(2*k+1)<=R_k-R_(k+1)<=4*R_k/(2*k+1)^2. It validates 594/594
  decrement-corridor rows, 591/591 theta-k monotonicity rows, and 396/396
  theta lambda-order rows on the k200 prefix, while two exact raw-cone
  monotone counterexamples block using raw decrease as a shortcut

jensen_window_pf_negative_lambda_k300_precision_repair_audit:
  finite precision-repair diagnostic; extends the raw-ratio decrement route
  to a repaired k300 stress. The broad dps160/cutoff6 run produces
  lambda=-100 high-k precision alarms, while local dps220/cutoff7 repairs
  over k220..250 and k245..320 restore 897/897 raw wall rows, 894/894
  decrement-corridor rows, 891/891 theta-k monotonicity rows, and 596/596
  theta lambda-order rows. This remains finite evidence only

jensen_window_pf_negative_lambda_raw_log_decrement_bridge:
  exact finite theorem-search diagnostic; rewrites the raw decrement corridor
  in log-ratio coordinates with p_k=log(R_k) and
  delta_k=p_(k+1)-p_k. It proves the exact equivalence to the two-sided
  delta_k bounds, validates 897/897 raw-log wall rows and 894/894
  log-corridor rows on repaired k300 data, and keeps two exact raw-cone
  counterexamples blocking raw-log decrease as a shortcut

jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge:
  exact finite theorem-search diagnostic; rewrites the raw/log decrement
  route in coefficient-curvature variables
  B_k=-log(((2*k-1)/(2*k+1))*R_k). It proves the exact equivalence to
  log((2*k+3)/(2+(2*k+1)*exp(-B_k)))<=B_(k+1)<=B_k, validates 897/897
  B-wall rows and 894/894 curvature-corridor rows on repaired k300 data,
  and blocks monotone curvature or raw walls as standalone shortcuts

jensen_window_pf_negative_lambda_linear_curvature_barrier_scout:
  exact finite theorem-search diagnostic; proves the exact sufficient lemma
  L_k(B)<=((2*k+1)/(2*k+3))*B for B>=0, reducing the nonlinear lower
  curvature barrier to the stronger linear target
  B_(k+1)>=((2*k+1)/(2*k+3))*B_k. It validates 897/897 B-wall rows and
  894/894 linear-barrier rows on repaired k300 data while keeping the
  analogous defect-width recurrence rejected

jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target:
  interval and analytic all-k theorem at lambda=-100; 16,074 Arb
  interval-Simpson compact blocks, an exact u>=5 ray, the repaired 318-row
  prefix, and the all-k first-summand perturbation transfer prove
  C_(k+1)>=C_k for every integer k>=1. This closes the linear lower
  curvature wall but not PF-infinity or the all-order Jensen bridge

jensen_window_pf_negative_lambda_scaled_curvature_log_ceiling_bridge:
  exact finite theorem-search diagnostic; rewrites C_(k+1)>=C_k as the
  affine log-ratio ceiling
  delta_k<=h_(k+1)-((2*k+1)/(2*k+3))*h_k-2*p_k/(2*k+3). It validates
  894/894 repaired k300 scaled-ceiling rows, 894/894 scaled-log-corridor
  rows, and records that the affine ceiling is sharper than the nonlinear
  raw-log upper wall, with only a tight high-k slack left in the finite data

jensen_window_pf_negative_lambda_relative_gaussian_curvature_bridge:
  exact finite theorem-search diagnostic; rewrites B_k as the negative second
  difference of the relative-Gaussian log moment sequence f_k and rewrites
  C_(k+1)>=C_k as the weighted four-point inequality
  (2*k+1)*f_(k-1)-(6*k+5)*f_k+(6*k+7)*f_(k+1)-(2*k+3)*f_(k+2)>=0. It
  validates repaired k300 B-positive, B-decrease, C-increase, and
  C-lambda-order rows while leaving the weighted all-k theorem open

jensen_window_pf_negative_lambda_relative_gaussian_taylor_stencil_scout:
  finite theorem-search diagnostic; certifies positive fixed-k Taylor leading
  signs for B_k, B_k-B_(k+1), and C_(k+1)-C_k, but validates that finite
  Taylor truncations are unstable proof objects: only 4/35 sampled truncation
  rows are all-positive, so a uniform Taylor-tail remainder theorem is still
  required

jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations:
  exact finite theorem-search diagnostic; decomposes the Taylor-tail log
  error into exact B, companion, and weighted-gap epsilon stencils, records
  4/35 positive baseline rows and 31 blocked baseline rows, and identifies
  the weakest finite half-margin budget 1.166490564421582442E-8. This
  sharpens, but does not prove, the uniform remainder theorem

jensen_window_pf_negative_lambda_relative_gaussian_pointwise_tail_budget:
  exact finite theorem-search diagnostic; converts the exact epsilon-stencil
  margins into pointwise log-tail and multiplicative relative-tail tolerances,
  with weakest half-safety log-tail envelope 1.458113205526978052E-9 and
  relative-tail ratio bound 1.458113204463930993E-9. These are required
  tolerances, not proved analytic tail estimates

jensen_window_pf_negative_lambda_relative_gaussian_next_increment_stencil_stress:
  finite theorem-search diagnostic; shows that the known next Taylor increment
  exceeds the crude pointwise half-safety budget on both tested positive
  baselines while the structured B, companion, and weighted-gap stencil signs
  are preserved. This rejects the current pointwise triangle shortcut and
  sharpens the live direct signed stencil-tail route

jensen_window_pf_negative_lambda_relative_gaussian_degree16_stencil_continuation:
  finite theorem-search diagnostic; computes the degree-16 Taylor coefficient
  and tests all four current positive baselines one Taylor step further. Three
  preserve structured stencil signs, while the degree-14 T=1000 baseline fails
  through the companion stencil and the degree-14 T=2000 baseline survives,
  sharpening the route toward an explicit large-T or q/T collar

jensen_window_pf_negative_lambda_relative_gaussian_degree16_collar_scan:
  finite theorem-search diagnostic; scans every integer T from 900 to 2200 at
  k=22 for the degree-16 continuation. The M=7 baseline first has positive
  stencils at T=918, the M=8 continuation first preserves signs at T=1156, the
  half-safety condition first holds at T=1483, and the pointwise budget fails
  on all 1283 baseline-positive rows

jensen_window_pf_negative_lambda_relative_gaussian_degree16_real_t_collar_scout:
  finite theorem-search diagnostic; certifies, for the rationalized fixed-k
  degree-16 finite surrogate only, that the four normalizers and all three
  structured stencil signs persist on the real half-line T>=1156. The remaining
  upgrade is interval coefficient control plus signed infinite-tail stencil
  bounds on the same collar

jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate:
  finite theorem-search diagnostic; upgrades the fixed-k degree-16 real-T
  surrogate collar to Arb coefficient-ratio balls. Bernstein coefficients
  certify F_21 through F_24, the stripped B product, the stripped companion
  product, and the stripped weighted-gap derivative numerator on
  0<=u<=1/1156, while leaving the infinite residual Taylor-tail theorem open

jensen_window_pf_negative_lambda_relative_gaussian_degree40_arb_collar_ladder_stress:
  finite theorem-search diagnostic; repeats the Arb real-T collar test for
  every even finite surrogate degree 16 through 40 at fixed k=22. All tested
  levels certify positive F_21 through F_24 normalizers and all three
  structured stencils on T>=1156 with zero Bernstein failures, while leaving
  the residual Taylor-tail bound beyond degree 40 and all-k upgrade open

jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget:
  exact finite theorem-search diagnostic; converts the degree-40 Arb collar
  margins into sufficient fixed-k residual targets on 0<=u<=1/1156:
  |R_i(u)|<=5.382819486765314521E-01*u^3 and
  |R_i'(u)|<=9.315354075509573936E-03*u for i=21..24. The finite
  degree-42..80 tail profile consumes less than 0.1% of those budgets, but no
  analytic majorant for the infinite residual tail is proved

jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout:
  finite theorem-search obstruction; extends the formal residual-term profile
  through degree 240 and shows that after decreasing to a least-term region
  near j=103, the formal value and derivative terms grow again from j=103 to
  j=104 for F_21 through F_24. This rejects monotone/geometric infinite
  formal-tail summation and fixed-radius Cauchy termwise summation as proof
  templates, while leaving the actual asymptotic remainder theorem open

jensen_window_pf_negative_lambda_relative_gaussian_asymptotic_remainder_target:
  exact theorem-search diagnostic; converts the degree-40 residual budget and
  formal-tail obstruction into explicit sufficient analytic-remainder
  constants. A common 1000x first-omitted-term theorem would fit inside the
  half-safety value and derivative budgets for F_21 through F_24 on
  0<=u<=1/1156, with common multiplier limit about 1419.939. This is target
  calibration only; the first-omitted-term and least-term remainder theorems
  remain open

jensen_window_pf_negative_lambda_relative_gaussian_actual_endpoint_remainder_scout:
  floating endpoint theorem-search diagnostic; evaluates the actual
  relative-Gaussian multiplier at T=1156 by generalized Gauss-Laguerre
  quadrature. At the selected order N=320, the value and derivative residuals
  for F_21 through F_24 are below one first omitted formal term and below 0.1%
  of the degree-40 half-safety budgets. This is endpoint evidence only; it is
  not interval-certified and does not prove a uniform collar remainder theorem

jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout:
  floating cancellation-reduced theorem-search diagnostic; subtracts the
  degree-40 polynomial inside the Gamma expectation and samples
  T=1156,1500,2000,5000,10000 for F_21 through F_24. Across four quadrature
  orders, all sampled value and derivative residuals are below one first
  omitted formal term, with worst value ratio about 0.9707100590 and worst
  derivative ratio about 0.9693567775. This is finite floating grid evidence
  only, not interval-certified and not a uniform collar theorem

jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target:
  open numerical-certification target; converts the cancellation-reduced
  grid slack into explicit intervalization obligations. A future certificate
  with total ratio error below 1.0e-2 would keep the finite grid below one
  first omitted formal term, but Laguerre node/weight intervals, Phi n-tail
  bounds, quadrature-remainder bounds, coefficient propagation, and a
  grid-to-collar bridge remain open

jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout:
  analytic padded-range theorem-search diagnostic; bounds the n>30 tails in
  Phi, Phi', and Phi(0) on 0<=x<=1 far below the 2.0e-3 per-source
  intervalization cap after normalizing by the diagnostic Phi(0)>=0.44 proxy.
  This only narrows the nlrgit_03 tail source; interval Laguerre-node range
  and Phi(0) lower-bound certificates, plus quadrature, coefficient, rounding,
  and grid-to-collar obligations, remain open

jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate:
  exact/Arb theorem-search diagnostic; certifies the two finite-grid side
  conditions used by the padded Phi-tail scout. Gershgorin/AM-GM bounds every
  recorded generalized Laguerre node by 809<T_min=1156, hence x<=1, and the
  n=1 term in Phi(0) is Arb-certified above 0.44. This does not certify
  individual Laguerre weights, quadrature remainders, coefficient propagation,
  rounding, or the grid-to-collar bridge

jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate:
  finite-grid Phi-tail source certificate; composes the padded n>30 Phi,
  Phi', and Phi(0) tail majorants with the certified side conditions x<=1 and
  Phi(0)>=0.44. The three omitted-n tail sources sit far below the 2.0e-3
  per-source intervalization cap for the recorded finite grid. This does not
  certify finite n<=30 node evaluations, Laguerre weights, quadrature error,
  coefficient propagation, rounding, or the grid-to-collar bridge

jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout:
  high-order floating theorem-search diagnostic; recomputes the worst
  cancellation-reduced grid row T=10000, F_21 at quadrature orders 96 through
  320. All sampled ratios remain below one first omitted formal term, with
  order-spread below 1e-14, calibrating a future rigorous quadrature-radius
  target of 1e-6. This is not a quadrature-remainder theorem or interval
  weight certificate

jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix:
  exact/formula route diagnostic; records the classical N=320 generalized
  Gauss-Laguerre remainder prefactor for alpha=41/2, about 2.79e-159,
  translates the 1e-6 quadrature ratio-radius target into explicit value and
  derivative 640th-derivative supremum caps, and verifies that adding that
  cap to the worst-row finite-plus-tail budget still keeps both ratios below
  one first omitted term. This is not a derivative bound, not interval
  adaptive integration, and not a quadrature-remainder theorem

jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate:
  Arb worst-row far-tail split diagnostic; certifies the finite n<=30
  cancellation-reduced continuum tail y>=200 for T=10000, F_21 using
  monotone Phi majorants and upper incomplete-Gamma moment bounds. The value
  and derivative tail ratios are about 8.66e-26 and 4.50e-27, far below the
  1e-6 quadrature ratio-radius target. Compact integration on 0<=y<=200,
  aggregation, all-row coverage, and grid-to-collar coverage remain open

jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_interval_integration_scout:
  Arb worst-row compact-interval diagnostic; imports the y>=200 far-tail
  split and tests six raw Arb panel hulls on the remaining compact interval
  0<=y<=200 for T=10000, F_21. The raw value and derivative width-to-cap
  ratios are about 1.38e39 and 4.28e37, so plain interval-Riemann hulls are
  rejected. The live route is now a local Taylor/Chebyshev panel model with
  exact Gamma-weighted moments; compact integration remains open

jensen_window_pf_negative_lambda_relative_gaussian_worst_row_chebyshev_panel_moment_scout:
  high-precision floating Chebyshev panel-moment diagnostic; integrates
  Chebyshev interpolants against incomplete-Gamma panel moments on the same
  six compact panels 0<=y<=200. Consecutive degree deltas from 16->20 onward
  are below the unscaled quadrature caps in both value and derivative
  channels, with reference degree-32 estimates about -6.58e-34 and -1.38e-32.
  This calibrates the local-model route but does not provide Arb coefficient
  balls or interpolation-remainder bounds

jensen_window_pf_negative_lambda_relative_gaussian_worst_row_arb_chebyshev_interpolant_moment_scout:
  Arb Chebyshev interpolant-moment diagnostic; promotes the floating
  Chebyshev ladder to Arb-enclosed interpolant coefficients and
  incomplete-Gamma panel moments for degrees 16, 20, 24, and 32. All
  consecutive interpolant Cauchy deltas are below the unscaled quadrature
  caps, but this still certifies only interpolant arithmetic; the
  interpolation remainder between the interpolants and the true compact core
  remains open

jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix:
  interpolation-remainder route diagnostic; quantifies the remaining compact
  panel theorem after Arb interpolant arithmetic. The heaviest Gamma panel is
  20<=y<=50 with mass upper about 0.6021615933248953, and the artifact records
  20 Bernstein degree/rho budgets plus 16 minimal-degree rows for the analytic
  domain and sup-norm certificates needed to bound true-function interpolation
  remainders. This is route sizing only; no analytic-domain, sup-norm,
  endpoint, Taylor-model, or compact-integral certificate is proved

jensen_window_pf_negative_lambda_relative_gaussian_endpoint_parity_repair_matrix:
  endpoint parity-repair diagnostic; makes the y=0 branch obligation explicit
  for the compact interpolation route. Arb series expansion of the finite
  n<=30 Phi truncation records 8 low-order odd Taylor rows through degree 15,
  with first odd coefficient absolute upper about 1.50e-1300. This does not
  prove endpoint analyticity: the admissible repairs are exact infinite-kernel
  evenness plus certified tail charge, or an x-variable endpoint-panel
  certificate

jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_panel_route_matrix:
  endpoint x-panel route diagnostic; quantifies the first-panel repair route
  after the parity matrix. The change y=T*x^2 maps 0<=y<=1 to 0<=x<=0.01,
  the transformed Gamma density has x^42 power for alpha=20.5, and the
  first-panel mass upper is about 1.62e-21. The artifact records 18
  Bernstein degree/rho budgets for a future x-domain and sup-norm certificate,
  but does not prove the exact x moments, x-panel interpolation remainder, or
  compact interval certificate

jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_moment_taylor_certificate:
  finite Arb endpoint certificate for the T=10000, F_21 worst row. It expands
  the finite n<=30 cancellation-reduced value and derivative cores through
  degree 64, integrates all 65 Taylor rows by exact transformed Gamma moments,
  and bounds the true-function tail by a rigorous |z|<=0.2 Cauchy majorant.
  The resulting first-panel integrals are certified negative, with both
  absolute-to-global-cap ratios below 4.10e-47. This closes only 0<=y<=1 for
  the finite core; the separate n>30 tail, five y>=1 panels, all-row coverage,
  and uniform collar theorem remain open

jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_x_moment_taylor_certificate:
  finite Arb compact certificate for the T=10000, F_21 worst row. The complete
  range 0<=y<=200 maps into x<=sqrt(0.02), inside the certified |z|<=0.38
  disk. A degree-128 Taylor model is integrated by 129 exact transformed
  Gamma moments, with true-function remainder radii consuming less than
  4.0e-10 and 1.3e-9 of the value and derivative caps. This replaces the
  six-panel interpolation obligation for the finite n<=30 core. The separate
  n>30 Phi tail, y>200 far tail, all-row coverage, and uniform collar theorem
  remain open

jensen_window_pf_negative_lambda_relative_gaussian_worst_row_full_expectation_certificate:
  complete true-expectation certificate for the T=10000, F_21 row. It
  composes the degree-128 compact finite core, the finite-core y>=200 tail,
  and a new Arb global n>=31 Phi/Phi'/Phi(0) normalization correction. The
  full value and derivative balls are negative, with first-omitted ratio
  uppers about 0.9707101 and 0.9693568 and margins above 0.029. No worst-row
  integral source remains open and generalized Gauss-Laguerre quadrature is
  not used in the final row certificate. The other finite-grid rows,
  all-source aggregation, and finite-grid-to-collar theorem remain open

jensen_window_pf_negative_lambda_relative_gaussian_all_row_direct_expectation_certificate:
  complete recorded-grid direct expectation certificate. A uniform
  degree-384 exact-moment/Cauchy core on x<=8/25, rowwise incomplete-Gamma
  real-tail bounds, and the global n>=31 normalization correction certify all
  20 value and all 20 derivative expectation balls negative and below one
  first omitted term. The worst ratios remain at T=10000, F_21, below
  0.970711 and 0.969357. This retires the recorded-grid integration source,
  but signed finite-degree stencil aggregation, interval coverage in T, the
  real-T collar, sign-regularity bridge, RH, and Lambda <= 0 remain open

jensen_window_pf_negative_lambda_relative_gaussian_recorded_grid_stencil_composition_certificate:
  fixed-k=22 recorded-grid stencil certificate. The 20 expectation balls fit
  the exact rational budgets |R_i|<=(1/2)u^3 and |R_i'|<=(9/1000)u, using
  less than 5.7e-4 of either budget. A fresh Arb perturbation ledger retains
  positive normalizer, B-product, companion-product, and weighted-gap
  derivative margins, certifying all five recorded T systems. The source
  floating thresholds are not proof inputs. Real-T interval coverage,
  remaining k windows, cone entry, sign regularity, RH, and Lambda <= 0
  remain open

jensen_window_pf_negative_lambda_relative_gaussian_finite_collar_segment_stencil_certificate:
  real-interval certificate on 1156<=T<=10000. A two-regime compact Cauchy
  and incomplete-Gamma majorant proves the exact rational residual budgets
  uniformly for F_21 through F_24, then composes them with the degree-40 Arb
  perturbation ledger. The fixed-k=22 stencil system is certified throughout
  the bounded segment. The ray T>10000, remaining k windows, cone entry, sign
  regularity, RH, and Lambda <= 0 remain open

jensen_window_pf_negative_lambda_relative_gaussian_full_kernel_evenness_cauchy_lemma:
  exact analytic lemma. The Jacobi theta transformation and differential
  operator covariance prove the full infinite Riemann kernel even and
  analytic on |z|<pi/8. Subtracting normalized even Taylor terms through
  degree 40 leaves an exact order-42 zero and x^42-factored Cauchy bounds.
  This is not inferred from finite-truncation near-evenness and is not by
  itself a Gamma-expectation, cone-entry, RH, or Lambda <= 0 theorem

jensen_window_pf_negative_lambda_relative_gaussian_full_real_T_fixed_k_stencil_certificate:
  full real-interval certificate for fixed k=22. The order-42 factor, a full
  infinite-kernel disk majorant, full real-tail majorants, geometric n-tail
  bounds, and upper-Gamma hazard monotonicity certify the residual budgets on
  T>=10000. Together with the bounded segment, all four normalizers and the
  three target stencils are certified for every real T>=1156. All-k coverage,
  cone entry, sign regularity, RH, and Lambda <= 0 remain open

jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate:
  Arb denominator-side theorem-search diagnostic; certifies the first omitted
  coefficient r_21 is negative and bounded away from zero, then proves every
  value and derivative first-omitted denominator on the recorded finite grid is
  positive. The worst lows occur at T=10000, F_21, translating a 1e-6
  ratio-radius target into scaled absolute-radius caps about 6.78e-28 and
  1.42e-30. This does not certify the residual numerator, quadrature error,
  rounding, or the grid-to-collar bridge

jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate:
  Arb coefficient-core propagation diagnostic; rebuilds coefficient-ratio
  balls r_0 through r_21 and propagates their radii through exact Gamma
  moments on the recorded finite grid. The worst value and derivative
  coefficient-radius ratios to the first omitted denominators are about
  8.61e-81 and 5.52e-83 at T=10000, F_21, far below the 1e-6 ratio-radius
  target. This retires the coefficient-ball source for finite-grid
  intervalization only; Phi node evaluation, Laguerre node/weight intervals,
  quadrature error, rounding aggregation, and grid-to-collar coverage remain
  open

jensen_window_pf_negative_lambda_relative_gaussian_worst_row_laguerre_root_bracket_certificate:
  Arb worst-row Laguerre node diagnostic; certifies 320 disjoint
  sign-changing brackets for the roots of L_320^(41/2) at the T=10000, F_21
  quadrature row. The widest bracket has width about 4.15e-17. The same row
  records 30 zero floating SciPy weights, so the remaining proof work is a
  non-floating Christoffel-weight interval certificate before this can become
  a quadrature-row interval enclosure

jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_midpoint_scout:
  Arb worst-row Christoffel-weight midpoint diagnostic; evaluates the
  generalized Laguerre weight formula at midpoints of the certified N=320
  root brackets. It gives 320 positive non-floating midpoint weights, repairs
  30 SciPy double-weight underflows, and matches Gamma(43/2) to relative error
  about 1.8e-18. Direct interval denominator evaluation still contains zero on
  all 320 rows, so the midpoint artifact itself remains non-promoted; the
  follow-up interval certificate closes the weight source for this one row

jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_interval_certificate:
  Arb worst-row Christoffel-weight interval diagnostic; evaluates
  L_321^(41/2) on each certified L_320^(41/2) root bracket by a centered
  Taylor enclosure. It separates zero in all 320 denominator intervals,
  certifies 320 positive Christoffel-weight intervals, repairs the 30 SciPy
  double-weight underflows as interval weights, and has a weight-sum interval
  containing Gamma(43/2). Phi/Phi' node evaluation, quadrature remainder,
  rounding aggregation, all-row coverage, and grid-to-collar coverage remain
  open

jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate:
  Arb worst-row finite-part weighted-sum diagnostic; refines the T=10000,
  F_21, N=320 Laguerre node brackets to 120 bisection steps, evaluates the
  finite n<=30 Phi/Phi' terms on those intervals, sums them with certified
  Christoffel weights, and subtracts the polynomial part through exact Gamma
  moments. The value and derivative residual ratios are certified below one
  first omitted term. The n>30 tail composition, quadrature remainder,
  rounding aggregation, all-row coverage, and grid-to-collar coverage remain
  open

jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_plus_tail_budget_certificate:
  Arb worst-row finite-plus-tail budget diagnostic; composes the certified
  finite n<=30 weighted-sum ratio uppers with the full 2.0e-3 Phi-tail
  source cap for T=10000, F_21, N=320. The composed value and derivative
  ratio uppers are 0.9853957992836557769419015895036210773888 and
  0.9714055674762067320093698741711260875260, both below one first omitted
  term. Quadrature remainder, rounding aggregation, all-row coverage, and
  grid-to-collar coverage remain open

jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix:
  exact countermodel gate; validates three positive two-atom Stieltjes moment
  witnesses showing that generic raw moment positivity/log-convexity cannot
  prove the upper raw wall, scaled-upper corridor side, or monotone-bridge
  lower corridor side. It blocks generic moment-sequence promotion while
  leaving the zeta-specific all-k theorem open

jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target:
  exact interval-analytic theorem at lambda=-100. Full ratio-cone entry gives
  B_k>=0 and B_(k+1)<=B_k, the scaled-curvature theorem gives
  B_(k+1)>=((2*k+1)/(2*k+3))*B_k, and the exact linear-to-nonlinear calculus
  lemma proves both raw-ratio decrement-corridor inequalities for every k>=1

jensen_window_pf_negative_lambda_m100_adaptive_defect_certificate:
  exact interval-analytic theorem at lambda=-100. Full cone entry proves the
  defect cone and monotone bridge, while the upper raw-corridor wall proves
  scaled-defect increase; this closes the parameter-specific defect-tail and
  adaptive-defect inputs without claiming the simultaneous all-three-lambda
  theorem

jensen_window_pf_heat_flow_monotone_closure_scout:
  finite heat-flow closure diagnostic; validates 4 exact lambda-flow algebra
  rows, 315 finite Arb threshold rows, and 305 finite Arb flow-bracket rows,
  using the exact boundary-threshold lemma to isolate the remaining global
  flow-invariance and adjacent-log-concavity gaps without promoting them to an
  analytic theorem

jensen_window_pf_heat_flow_ratio_cone_invariance_lemma:
  exact conditional ratio-cone invariance lemma; proves inward-pointing
  boundary algebra for (2*k-1)/(2*k+1)<=x_k<=1 and x_{k+1}>=x_k under the
  heat-flow ratio ODE, while leaving the zeta cone-entry and infinite-flow
  legitimacy theorems open

jensen_window_pf_heat_flow_boundary_threshold_lemma:
  exact boundary-threshold lemma; proves from Phi positivity and
  Cauchy-Schwarz raw-moment log-convexity that
  x_k >= (2*k-1)/(2*k+1), hence the heat-flow boundary threshold
  x_k >= (2*k-1)/(2*k+5), while leaving adjacent log-concavity and global
  monotone contractions open

jensen_window_pf_kernel_mellin_upper_wall_certificate:
  interval-backed theorem certificate; proves y->Phi(sqrt(y)) strictly
  log-concave by a 200-subinterval Arb/Cauchy compact proof plus an analytic
  full-kernel ray. Berwald-Borell then proves x_k(lambda)<=1 for every real
  lambda and every k>=1. The adjacent-k wall x_(k+1)>=x_k remains open

jensen_window_pf_log_concave_mellin_monotone_wall_countermodel:
  exact interval countermodel gate; the log-concave density
  exp(-5y)*1_[0,1](y) has Gamma-normalized Mellin contractions x_1 and x_2
  inside the upper wall but x_2<x_1 by more than 0.027. Generic
  log-concavity therefore cannot prove the remaining adjacent-k wall

jensen_window_pf_negative_lambda_t1156_monotone_wall_counterexample_certificate:
  interval zeta-kernel counterexample gate; rigorous ACB enclosures of
  A_119..A_122 at lambda=-1156 certify x_121-x_120<-1.68e-8 while both
  contractions remain inside the Mellin upper wall. This blocks all-k cone
  entry at T=1156 and any universal promotion of the fixed-k T>=1156 theorem,
  while leaving moderate-parameter finite-collar plus eventual-k routes open

jensen_window_pf_negative_lambda_kernel_summand_shift_lemma:
  exact all-n kernel-shift lemma with a compact interval theorem; proves
  phi_n(u)=n^(-1/2)phi_1(u+(log n)/2) and, at lambda=-100 and k>=300,
  bounds the complete shifted n=2..20 contribution on v<=3/2 below
  2.122e-29 of the n=1 moment. The later all-k dominance theorem discharges
  the shifted far tail; the dominant n=1 adjacent wall remains open

jensen_window_pf_negative_lambda_first_summand_dominance_certificate:
  all-k analytic kernel-tail theorem; exact ratio monotonicity, strict
  first-integrand concavity, and 15 Arb-positive propagation gates prove
  M_k=M_k^(1)(1+delta_k), 0<=delta_k<=2/k^6 for every k>=300. The adjacent
  log-wall perturbation is at most 16/(k-1)^6 for k>=301

jensen_window_pf_negative_lambda_m100_k320_collar_extension_certificate:
  finite Arb collar extension; promotes the repaired A_245..A_320 source into
  76 positive coefficients, 74 ratio-cone rows, and 73 adjacent-wall rows.
  Nineteen newly exposed rows extend the lambda=-100 collar through k=318,
  with minimum new log gap above 3.709e-6

jensen_window_pf_negative_lambda_first_summand_cumulant_bridge:
  exact conditional reduction; writes the cubic moment difference as an
  average of kappa_3,t(2 log U), extracts the exact positive Gamma wall, and
  proves that kappa_3,t(2 log U)>=-37/(50 t^2) for t>=318 is sufficient for
  L_k^(1)>=1/(4 k^2) for every k>=319. The uniform cumulant bound remains open

jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate:
  global interval theorem for the leading cumulant asymptotic; symbolic
  differentiation, 40,740 Arb compact intervals, and an exact analytic ray
  prove caps 13/20, 1/100, and 1/1000 through fifth order for all t>=318.
  The remaining sufficient target is a seventh-order floor -79/1000

jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate:
  large compact all-real interval theorem; exact paired moments, 40,736 Arb
  eighth-derivative envelope intervals, rigorous midpoint errors, and explicit
  two-sided tails prove the normalized remainder floor -79/1000 throughout
  0.9264<=u<=5, reaching t>1.5241916613e10. The same 4,074 Arb blocks certify
  strict negative third log-cumulant, with maximum upper endpoint below
  -9.13e-6. The separate analytic ray certificate completes the half-line

jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate:
  analytic asymptotic-ray theorem; an adaptive sqrt(8 log q) window,
  first-order Gaussian moment comparison, and explicit tails prove
  H_t>=-3/250 on u>=5. Combined with the compact theorem this closes the
  global cumulant hypothesis for every t>=318

jensen_window_pf_negative_lambda_first_summand_saddle_wall_target:
  validated dominant-summand wall closure; exact strict-concavity geometry gives
  one saddle in an all-k bracket, and nine 60-digit samples through k=20000
  satisfy L_k^(1)>1/(4*k^2) with increasing k^2 L_k^(1). The exact cumulant
  bridge and paired compact/ray theorems prove L_k^(1)>=1/(4*k^2) for every
  k>=319. The dominance bound and finite collar then prove the full
  lambda=-100 adjacent wall

jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate:
  interval/analytic theorem composition; a precedence-merged repaired source
  certifies 321 positive coefficients, 319 pointwise cone rows, and 318
  adjacent prefix rows. Exact all-k pointwise walls and the analytic k>=319
  adjacent tail prove full infinite ratio-cone entry at lambda=-100; the
  separate maximum-principle certificate closes forward-flow legitimacy

jensen_window_pf_heat_flow_infinite_cone_invariance_certificate:
  exact infinite-index maximum principle; in adjacent-defect coordinates the
  heat ODE is cooperative with nonnegative source and bounded potential.
  Uniform h_k->0 yields a finite active set for every negative spatial
  minimum, and the Dini-minimum argument propagates the full cone from
  lambda=-100 to every finite lambda>=-100, including lambda=0

jensen_window_pf_defect_complete_monotonicity_scout:
  finite Arb diagnostic with exact guard; certifies 3284 defect and 3288
  negative-log-contraction alternating differences, both complete through
  order 8, while retaining 838 high-order interval inconclusives. The exact
  x=(1/2,1,1,...) Hausdorff-defect witness has cubic Jensen discriminant
  -27/16, so these patterns are not an all-shape bridge

jensen_window_pf_multiplier_complete_monotonicity_frontier_scout:
  finite high-precision Arb certificate; using the dps220 A_0..A_57 sources
  at 250-digit working precision, certifies all 7980 alternating differences
  of y_k=-log(x_k) through order 55 on five nonnegative heat parameters, with
  no inconclusive or negative interval; this is necessary evidence only, not
  a unit-atomic counting-measure construction

jensen_window_pf_multiplier_hausdorff_uniqueness_bridge:
  exact measure-theoretic reduction; proves that the integer log contractions
  determine one finite Hausdorff measure and that any admissible unit-atomic
  multiplier multiset is unique. It characterizes the remaining density
  recovery problem and records the periodic-interpolation guard separating it
  from the continuous Mellin obstruction

jensen_window_pf_multiplier_leading_atom_bound_certificate:
  conditional interval theorem for any putative unit counting measure;
  brackets the order-6 root in (4.863538496,4.863538497), proves
  alpha_min>4.863538496, verifies that all other orders through 55 are weaker,
  and proves N(11/2)<=1 without asserting that such a low atom exists

jensen_window_pf_heat_flow_jensen_hierarchy_lemma:
  exact shift-degree heat hierarchy; identifies partial_lambda J_(d,n) with a
  shifted-window operator and with a first/second z-derivative operator on
  J_(d+1,n). The one-atom Hausdorff defects d_k=(4/9)(9/25)^(k-1) give the
  exact full-cone cubic boundary x=5/9, y=21/25, z=589/625 with
  (partial_lambda F)/r_0=329728/2109375>0; even complete defect monotonicity
  plus the local heat ODE does not provide an invariant higher-minor cone

jensen_window_pf_rank_two_boundary_family_lemma:
  exact all-degree benchmark; the sequence A_k=a^(k-1)(c+kb)/u factors every
  shifted Jensen window into one simple and one repeated negative root. In
  alpha coordinates, finite pointwise products are multiplier sequences and
  give a discrete canonical-product formula for -log x_k. The -937/3456
  mixture and <-27/125 fractional-power guards reject arbitrary positive
  measure weights; a genuine counting-measure representation remains open

jensen_window_pf_cubic_reciprocal_defect_invariance_lemma:
  exact degree-3 coordinate and conditional flow theorem with finite Arb entry;
  q_k=(1-x_k)^(-1/2) turns cubic hyperbolicity into
  q_(k+1)-q_k<=1. At saturation, the heat field points inward exactly when the
  next increment satisfies the same bound. All 318 lambda=-100 prefix margins
  and 310 nonnegative-grid margins are certified. The later cubic tail-entry
  certificate closes the all-k lambda=-100 tail; forward-uniformity and every
  higher-degree minor cone remain open

jensen_window_pf_cubic_m100_tail_entry_certificate:
  all-k degree-3 entry theorem at lambda=-100; 4,074 compact negative-skewness
  blocks plus an analytic ray prove kappa_3<0 for every t>=318. Two-sided
  adjacent log walls imply d_m>=1/(5m+1), q_(k+1)-q_k<1 for every k>=319,
  and q_(k+1)-q_k->0. Together with the 318 prefix margins this proves every
  shifted cubic Jensen polynomial hyperbolic at lambda=-100. The subsequent
  forward-uniform certificate propagates degree 3 through lambda=0; all higher
  degrees remain open

jensen_window_pf_cubic_forward_uniform_tail_certificate:
  exact all-shift degree-3 forward theorem; the q-increment vector field is
  one-sided cooperative and obeys sqrt(k)*g_k'<=7*r_k at the weighted barrier.
  The entry estimate sqrt(k)*g_k<12 and an explicit coercive supersolution give
  sup_[-100,L] g_k<=(12+7*R_L*(L+100))/sqrt(k) for every finite L. This closes
  the infinite first-crossing handoff and proves all shifted cubics hyperbolic
  at lambda=0; degree 4, PF-infinity, and RH remain open

jensen_window_pf_quartic_boundary_flow_obstruction:
  exact countermodel gate; derives Disc(J_4)=256*x^6*y^2*Q and the heat
  derivative of Q. A rational hyperbolic double-root boundary point satisfies
  all local ratio walls and three strict neighboring cubic inequalities, but
  has Q'/r_1=-13108711376416987159336748097/
  20606742971316325673502124987495<0. Therefore the propagated cubic cone is
  not by itself a quartic invariant; this does not show failure of the zeta
  trajectory and instead requires a new coupled degree-4 condition

jensen_window_pf_strong_logconcave_local_quartic_countermodel:
  exact moment-realized local countermodel; the
  `1/5`-strongly log-concave squared-variable density
  `exp(-3*y-y^2/10)1_[1/20,1](y)` has four Gamma-normalized Mellin
  contractions satisfying every local lower/upper ratio wall,
  `x_1<x_2<x_3<x_4`, and three consecutive strict cubic tests. Arb gives
  `Q(x_1,x_2,x_3)<0`, so its shift-zero quartic discriminant is negative.
  Full-support `1/5`-strongly log-concave approximants preserve these strict
  finite signs by dominated convergence. This blocks only a local generic
  promotion: it is not the Xi kernel or a Newman heat trajectory and does not
  challenge the global Xi all-shift cubic theorem

jensen_window_pf_quartic_double_root_threshold_lemma:
  exact local quartic boundary theorem; for
  P=(1+a*w)^2(1+b*w)(1+c*w), p=bc, the heat value at the double root has sign
  u-U with U=(-a^2+2*a+p)(3*a^2-5*a+5*p)/(6*p^2). The complete inward
  condition is (3*a^2-4*a+p)*(u-U)<=0. On the triple-root stratum first-order
  viability requires u=U and the first variation retains (1+a*w)^2. A closed
  contraction-coordinate quartic invariant and higher-time tangency remain open

jensen_window_pf_quartic_signed_hankel_branch_exclusion_lemma:
  exact Xi branch-selection theorem. At a normalized quartic double-root
  boundary, the contiguous order-three Hankel determinant is a positive shift
  factor times `-C^3/216`, where
  `C=3*a^2-4*a+p=(a-b)*(a-c)`. The completed Xi compound-order-three theorem
  gives `D_(3,n)<0`, hence `C>0`: the repeated root lies outside the two simple
  roots, while the middle-root and triple-root strata are impossible. The
  live inward condition therefore reduces from
  `C*(u-U)<=0` to the single unproved outer threshold `u<=U(a,p)`

jensen_window_pf_quartic_outer_threshold_order4_nonpromotion_gate:
  exact finite local countergate. The rational segment `A_1,...,A_10` has an
  outer-double-root shift-one quartic, eight strict pointwise ratio
  coordinates, seven increasing scaled-defect steps, seven reciprocal-defect
  increments below one, seven strict adjacent cubic frontiers, all `120`
  supported order-two signs, all `126` supported arbitrary-column
  order-three signs, and all `56` supported arbitrary-column order-four
  signs across multiple forward shifts. Nevertheless `u-U(a,p)>0`, so its
  quartic heat vector is outward; the adjacent quintic has negative
  discriminant. This blocks finite promotion from the strengthened local
  scalar and compound layers through order four. The segment is not Xi, an
  infinite sign-regular sequence, or a Newman trajectory. The exact next
  `x_10` continuation corridor has positive rational width, and one rational
  interior point preserves the next order-three/order-four signs and scalar
  corridor; the downstream length-13 obstruction proves that this fixed
  prefix nevertheless has no all-length increasing continuation

jensen_window_pf_quartic_outer_branch_length13_obstruction:
  exact fixed-prefix obstruction. Every strict signed-Hankel continuation of
  the outer-contact prefix through `x_12` is parameterized by
  `G_(k-3)=y_k*C_k`, `0<y_k<1`. An increasing `x_13` with the next order-four
  sign would require `Delta_13=C_13-M_13>0`. Exact Bernstein certificates for
  1516 derivative coefficients prove `Delta_13` is coordinatewise
  nondecreasing, but its closed-cube maximum is below
  `-190137/125000000000`. Thus this prefix cannot be an infinite
  countermodel. The alternate-tail gate below proves explicitly that this is
  prefix-specific

jensen_window_pf_quartic_outer_branch_alternate_length13_survivor_gate:
  exact finite survivor and scope-repair gate. Keeping the same outer contact
  only through `x_5`, eight rational corridor parameters construct a distinct
  increasing tail through `x_13`. It clears every stated scalar gate, all ten
  order-three gaps, all eight order-four cap margins, and all `364`, `715`,
  and `792` supported arbitrary-column signed minors of orders two through
  four on `A_1,...,A_14`. Here `Delta_13>0`, whereas this chosen tail has
  `Delta_14<0`. Thus the earlier obstruction is not uniform, and the new
  length-14 obstruction is fixed-tail only within this artifact. Its bounded
  search remains numerical evidence; the separate interval certificate below
  supplies the rigorous all-tail theorem for this one fixed contact

jensen_window_pf_quartic_outer_branch_uniform_length14_interval_certificate:
  rigorous one-contact interval theorem. The first scaled wall is exactly
  `0<y_6<Y_6`, every later signed corridor has `0<y_k<1`, and cancellation-free
  gap monomials place every admissible tail in `[0,Y_6]x[0,1]^7`. A 256-bit,
  hash-chained Arb tree has `74,947` events, `37,473` splits, `37,474`
  certified leaves, maximum depth `23`, and no uncovered box. Independent
  replay proves every mean-value upper enclosure for `Delta_14` is strictly
  negative. Thus every admissible tail from this contact dies at the
  increasing-`x_14` order-four compatibility test. The complete outer-contact
  family, the Xi threshold, PF-infinity, `Lambda<=0`, and RH remain open. The
  contact-level survivor below proves that the one-contact restriction is
  essential

jensen_window_pf_quartic_outer_contact_normal_form_gate:
  exact contact-level reduction. With `delta=1-a`,
  `q=C/(4 delta^2)`, and `1-x_5=delta^2 s`, every contact defect has scale
  `delta^2`; `G_1=delta^6 h_1`, `G_2=delta^4 h_2`, and
  `Delta_14=delta^4 delta14_hat`. The opposite outer-root branch violates
  `x_2<x_3`. The first extension has two exact regimes according as the
  scaled wall or order-four cap is smaller. This normal form removes the
  singular contact scale. It also factors the adjacent-quintic discriminant
  as a positive prefactor times `epsilon^2 R(epsilon)`, with
  `epsilon=T-s`; for `0<q<15/16`, one exact positive root `epsilon_+`
  bounds a negative-discriminant outward collar. It does not decide the
  complete tail-cell sign or establish that Xi enters this collar

jensen_window_pf_quartic_outer_contact_length14_survivor_gate:
  exact finite countermodel. The rational contact
  `delta=13/200`, `q=1/2`, `s=T-1/100000` has `x_5-U>0`. Setting every
  `y_6,...,y_14` to `99/100` constructs a strict tail through `x_14`.
  All `4,043` supported arbitrary-column signed minors of every possible
  order two through eight on `A_1,...,A_15` have positive 512-bit Arb lower
  endpoints, while `Delta_14>0` and `Delta_15>0` exactly. This rejects
  contact-uniform length-14 promotion from the finite recorded hypotheses.
  Yet the adjacent quintic has `P_5(-1/a)<0` and negative exact
  discriminant, hence three real roots and one nonreal conjugate pair. This
  separates the finite signed-Hankel conditions from degree-five
  hyperbolicity, but it is not Xi, an infinite sequence, or a Newman
  trajectory

jensen_window_pf_quartic_quintic_polar_contact_lemma:
  exact adjacent-degree theorem; normalized windows obey
  P_d=P_(d+1)-w*P_(d+1)'/(d+1). For a positive-root hyperbolic quintic, a
  double root of the quartic polar derivative forces a quintic triple root.
  In the quartic boundary coordinates P_5(-1/a) is a nonzero prefactor times
  -(u-U), so the heat threshold is exactly quintic triple contact. This points
  to a general adjacent-degree hierarchy but does not close its infinite
  top-down dependence

jensen_window_pf_cofinal_degree_polar_closure_lemma:
  exact all-degree reduction; repeated zero-pole polar descent proves that one
  hyperbolic terminal degree D closes every lower degree at the same shift.
  Therefore an unbounded sequence D_j of hyperbolic terminal degrees is enough
  for all finite degrees, and cofinal terminal sequences at every shift suffice
  for the complete Jensen family. Current Sturm evidence reaches only degrees
  3 through 12 on a finite grid; no cofinal zeta terminal theorem is proved

jensen_window_pf_cofinal_scaling_limit_equivalence_gate:
  exact noncircularity theorem; P_(D,n)(z/D) converges locally uniformly to the
  fixed-shift entire function F_n. Hence an unbounded hyperbolic degree
  subsequence forces F_n into the Laguerre-Polya class, and Jensen's theorem
  gives the converse. At shift zero the cofinal terminal theorem is already an
  endpoint-equivalent statement, while the bounded 1050-row Sturm ladder does
  not enter the limit

jensen_window_pf_coefficient_pf_equivalence_gate:
  exact route correction; c_k=A_k/k! has ordinary generating function equal
  to the exponential generating function of A, and differentiation gives each
  shifted tail. ASW/Edrei, Polya-Schur, derivative closure, and finite ASW make
  all-order coefficient PF-infinity exactly equivalent to PF-infinity of every
  shifted Jensen window. Finite evidence and signed-Hankel data still do not
  prove the common endpoint theorem

jensen_window_pf_edrei_stieltjes_equivalence_gate:
  exact route correction; for the entire normalized H, the shifted signed
  Edrei-log sequence a_r=p_(r+1) is exactly the signed Taylor sequence of
  H'/H. Sokal's criterion identifies its all-order Stieltjes moment property
  with H in LP+, hence with coefficient PF-infinity and every shifted Jensen
  window. Strict positivity of the two Hankel columns s=1,2 at every order is
  sufficient; the 4205 finite rows remain diagnostics only

jensen_window_pf_polar_heat_collision_cascade_lemma:
  exact collision-cascade theorem; at a multiplicity-m root the adjacent
  degree controls the full low heat jet. If every higher degree is hyperbolic,
  polar lifting raises the common-root multiplicity once per degree. Such a
  cascade forces the entire coefficient function to be exponential times a
  bounded-degree polynomial. Therefore a non-exponential-polynomial LP
  boundary is strict in every fixed Jensen degree, and failure can approach it
  only through degrees tending to infinity

jensen_window_pf_scaled_double_zero_boundary_layer_lemma:
  exact degree-scaled boundary theorem; the first two Jensen corrections are
  -z^2*F''/(2D) and (z^3*F'''/3+z^4*F''''/8)/D^2. Near a nondegenerate double
  zero rho<0, the joint layer is eta^2+8*rho*tau-rho^2. The finite pair has an
  unscaled D^(-3/2) gap, and its D^(-2) collision correction contains
  F'''/F''=3U'/U, the regularized global root external field. The nested
  finite-degree eventual thresholds exhaust Lambda, but the required degree-
  and zero-height-uniform remainder theorem remains open

jensen_window_pf_newman_root_external_field_lemma:
  exact root-field translation; at a squared double zero the D^(-2) Jensen
  correction is governed by E_x=sum_(j!=*)m_j/(x_j^2-x^2), while its signed
  Newman stiffness is the strictly positive sum K_x=sum_y 1/(x-y)^2. The
  squared pair gap satisfies q'=8-4qS and has correction -16*K_x*(t-t_*)^2.
  Two finite LP products realize opposite field signs, so generic LP or
  coefficient positivity cannot replace the open Xi-specific balance theorem

jensen_window_pf_newman_classical_field_balance_gate:
  exact reference-field and compactness gate; the Riemann-von Mangoldt
  continuum field tends to -pi/8, matching the -pi/4 high-positive-time zero
  drift. Published asymptotics make every sufficiently high fixed-time zero
  simple, so a positive boundary collision lies below exp(C/Lambda). Exact
  even-lattice perturbations show that bounded classical-location error still
  allows either unbounded field sign; a lambda-uniform reciprocal-gap estimate
  or compact-region exclusion remains open

jensen_window_pf_newman_local_odd_count_reduction_lemma:
  exact Stieltjes localization theorem; the published uniform zero count
  controls every field contribution outside radius H=log(4c)^2, leaving
  B(c)=-pi/8+S_H+O(1/log c), where S_H is the inverse-square weighted odd
  local counting discrepancy. An exact even backward-heat polynomial has the
  same -pi/8 field and -pi/4 drift while still producing positive square-root
  birth. Thus both a lambda-uniform Xi odd-count theorem and a separate global
  collision obstruction remain necessary

jensen_window_pf_newman_boundary_energy_direction_gate:
  exact collision-energy and directionality gate; every double-zero birth has
  renormalized pair energy asymptotic to 1/(8*(t-t_*)), so finite integrated
  energy down to the boundary would exclude collision. The exact classical-
  field model has positive stiffness and positive cubic gap jet while its
  energy decreases after birth from a nonintegrable trace. The published
  Rodgers-Tao theorem assumes Lambda<0 and starts at Lambda/4, away from the
  boundary, so an Xi-specific boundary-trace estimate remains open

jensen_window_pf_newman_positive_boundary_attainment_lemma:
  exact compactness composition; under Lambda>0 the absolute strip bound and
  uniform positive-time high-zero theorem trap nonreal zeros below the boundary
  in one compact rectangle, forcing a finite real multiple zero of H_Lambda.
  With the published Lambda<=1/5 bound, Lambda<=0 is equivalent to simplicity
  for every 0<t<=1/5. The universal
  multiplicity-m Hermite split has ordered cluster energy
  m(m-1)/(8*(t-Lambda))+O((t-Lambda)^(-1/2)), so every multiplicity has a
  nonintegrable endpoint trace; Xi simplicity or endpoint control remains open

jensen_window_pf_newman_positive_boundary_delta_localization_gate:
  exact cofinal compact-strip reduction. For every fixed delta>0, the dominant
  ray confines possible contacts with t>=delta to
  |x|<4*pi*exp(25/delta); delta_j=1/(5j) gives a countable family equivalent
  to Lambda<=0. The exact model G_t=x^2-2t has positive-time simplicity and a
  double endpoint zero with vanishing first-jet margins. Thus a uniform floor
  down to t=0 is a valid stronger route, not a logically necessary burden;
  none of the Xi compact strips is certified here

jensen_window_pf_newman_positive_boundary_diagonal_exhaustion_gate:
  exact independent compact-exhaustion theorem. Because a positive boundary
  supplies one fixed finite collision, any delta_j->0 and R_j->infinity
  detect it eventually, with no required relation between their rates.
  Composing with the certified |x|<=38 core gives the linear shell family
  1/(5j)<=t<=1/5, 38<|x|<=38+j, exactly equivalent to Lambda<=0. An
  arbitrary-height even polynomial heat flow realizes the classical field
  -pi/8 and stiffness tending to pi^2/64 at any prescribed collision height
  and shifted boundary time. Thus neither the exponential radius nor a local
  field replacement is logically mandatory; the unbounded Xi shell family
  remains open

jensen_window_pf_newman_first_jet_winding_gate:
  exact all-multiplicity contact-degree and boundary-flux reduction. The
  universal Hermite heat polynomial shows that an m-fold real contact creates
  floor(m/2) pairs and has local first-jet index +floor(m/2) in the standard
  (x,t) orientation. Therefore all
  interior contacts contribute with one sign and zero boundary winding is
  equivalent to no contact. For Z=H+iH_x, its phase one-form is exactly
  -(L/Q)dx+(L_x/Q)dt. On the half rectangles
  [1/(5j),1/4]x[0,38+j], axis positivity and top-time simplicity close two
  edges, reducing each 2D stage to bottom/right first-jet separation and one
  integer phase-flux condition. Those unbounded Xi edge estimates remain open

jensen_window_pf_newman_theta_first_diagonal_shell_interval_certificate:
  rigorous first-stage interval theorem. Five 160-bit Arb/Taylor boxes with
  analytic 10^-800 tails prove first-jet separation at t=1/5 on 38<=x<=39,
  with no subdivision, no unresolved box, and normalized margin above 6/5.
  Evenness, the certified |x|<=38 core, and simplicity for every t>1/5 close
  Q_1=[1/5,1/4]x[0,39] and give zero first-jet winding on its boundary. This
  certifies j=1 only; Q_2 and the cofinal family remain open

jensen_window_pf_newman_theta_second_diagonal_shell_first_block_route_guard:
  rigorous pointwise route guard for the second stage. At t=1/10, Arb upper
  bounds prove that both first-block sufficient ratios are below one at
  x=397/10,199/5,399/10,40; the largest failed upper ratio is below 24/25.
  The neighboring x=198/5 still passes the derivative branch. Therefore
  subdivision alone cannot promote the unchanged first-block bars across the
  bottom edge of Q_2. This does not indicate a contact. At the guard stage
  Q_2 required sharper tails, more theta blocks, or a different winding-edge
  certificate; the succeeding two-block theorem supplies that repair

jensen_window_pf_newman_theta_second_diagonal_shell_two_block_interval_certificate:
  rigorous second-stage interval theorem. Retaining the actual n=1,2 theta
  contributions gives the required cancellation, while an exact positive-tail
  estimate proves M_0<10^-10 and M_1<10^-12 for n>=3. Twenty-five 160-bit
  Arb/Taylor boxes, with no subdivision or unresolved box, prove H_t<0 on
  [1/10,1/5]x[38,40]; the minimum normalized value-sign ratio is above 39000.
  Composition with the |x|<=38 core and simplicity above 1/5 closes
  Q_2=[1/10,1/4]x[0,40] and its winding. This certificate certifies j=2 only;
  the separate successor below closes Q_3, while the cofinal family remains
  open

jensen_window_pf_newman_theta_third_diagonal_shell_two_block_interval_certificate:
  rigorous third-stage interval theorem. The exact retained n=1,2 kernel and
  analytic n>=3 tail theorem apply unchanged. Forty-eight 160-bit Arb/Taylor
  boxes, with no subdivision or unresolved box, prove H_t<0 on
  [1/15,1/5]x[38,41]; the minimum normalized value-sign ratio is above 12000.
  Composition with the |x|<=38 core and simplicity above 1/5 closes
  Q_3=[1/15,1/4]x[0,41] and its winding. This certificate certifies j=3 only;
  the separate successor below closes Q_4, while the cofinal family remains
  open

jensen_window_pf_newman_theta_fourth_diagonal_shell_two_block_interval_certificate:
  rigorous fourth-stage mixed-branch interval theorem. One hundred twenty
  160-bit Arb/Taylor boxes, with no subdivision or unresolved box, prove full
  first-jet separation on [1/20,1/5]x[38,42], with minimum normalized
  disjunction ratio above 88000. Of these, 105 prove H_t<0; the remaining 15
  cover [1/20,1/5]x[83/2,42] and require strict derivative separation.
  Composition with the core and simplicity above 1/5 closes
  Q_4=[1/20,1/4]x[0,42] and its winding. This certifies j=4 only; Q_5 and the
  cofinal family remain open

jensen_window_pf_newman_theta_fixed_block_cofinal_obstruction_gate:
  exact route guard for fixed arithmetic truncations with static absolute-tail
  bars. Uniform double-exponential decay and two integrations by parts give
  retained bounds A_N/x^2 and B_N/x^2. Therefore the normalized value and
  derivative margins against fixed positive moment bars both fall below one
  beyond an explicit finite threshold, so no fixed N can certify every
  diagonal shell. For N=2, the omitted n>=3 tail exactly cancels the positive
  retained endpoint derivative a_2=8.2652795777...e-8; taking absolute moments
  erases that cancellation. The hash-chained j=5..30 scout remains a finite
  diagnostic and promotes no new Q_j. This rejects only the static-tail proof
  mechanism, not Xi transversality; a growing-block or cancellation-aware
  remainder theorem remains open

jensen_window_pf_newman_theta_modular_blend_gate:
  exact cancellation-preserving theta partition. The entire switch produces
  positive even entire blocks summing ordinarily to Phi, and every bounded
  positive Newman deformation remains uniformly Schwartz. Transform
  derivatives and the coupled Laguerre matrix therefore converge normally,
  while repeated Fourier integration by parts gives arbitrary-power tail
  decay. Useful derivative constants and a retained lower profile remain
  open. Negative block-transform and nonreal Jensen witnesses block termwise
  spectral positivity and individual-block Laguerre-Polya promotion

jensen_window_pf_newman_theta_modular_blend_high_frequency_scout:
  ten high-precision point diagnostics for the first three blended blocks and
  full theta Laguerre expression at two times and five frequencies. The
  arithmetic tail cancels more than five digits at x=120, eleven at x=150,
  and twenty at x=200. This blocks fixed three-block absolute-tail dominance
  but proves no interval sign or adaptive remainder theorem

jensen_window_pf_newman_theta_modular_blend_adaptive_saddle_gate:
  exact modular-transition saddle geometry with finite adaptive diagnostics.
  The a=5,9 component saddles cross the modular switch at
  n_*(x,t)=sqrt(x/(4*pi))(1+O(x^-1)) for bounded positive Newman time,
  identifying the Riemann-Siegel scale as the natural growing block count.
  Ten high-precision rows track the first scale-level count from three to six
  and fit a two-index empirical collar. The saddle law proves no contour
  deformation or remainder sign, and the sampled global -L_t' branch is
  rejected by the separate Xi Lehmer certificate. The live target is a
  cancellation-aware direct first-jet remainder on the adaptive scale

jensen_window_pf_newman_theta_adaptive_modular_c1_remainder_contract:
  exact composition of the modular partition, arbitrary-power tail budgets,
  direct J/J' error bars, and preliminary square-root saddle scale. Either one
  of two explicit retained inequalities excludes full-Xi contact. The later
  derivative-envelope theorem supplies the corrected absolute-tail scale, but
  no retained first-jet lower profile or transition-cell theorem. The rejected
  global -L_t' condition is not assumed

jensen_window_pf_newman_theta_modular_tail_derivative_envelope_gate:
  exact heat-factor polynomials and positive-half-line derivative envelopes
  for d_(N,0,m),d_(N,1,m). Minimizing the reflected modular phase proves
  effective C*N^A*exp(-c*N^(4/3)) bounds for every fixed weighted derivative
  norm. Thus N_K(x)=ceil(K*(1+x)^(3/4)) yields exponential-in-x absolute tail
  control, and the corrected count is max(N_sad,N_K). This is an asymptotic
  tail theorem, not a retained first-jet separation theorem. Its formerly
  effective constants are instantiated at T=1/5 by the later direct-tail
  explicit-constant gate

jensen_window_pf_newman_theta_modular_tail_derivative_budget_scout:
  finite T=1/5 stress test of m=5..9 envelopes and 64 retained first jets.
  The spectral ladder is stable to about 3.4e-31 relatively, while the
  absolute-value quadrature has 1.93e-2 worst node drift and remains
  diagnostic. At x=200 the kappa=2 rows first pass at m=9 and m=8. At x=300,
  kappa=2,3,4 fail every sampled order and kappa=5 passes at both times. This
  calibrates the growing-count route but proves no collar obstruction or
  interval separation

jensen_window_pf_newman_theta_modular_tail_arb_quadratic_pilot:
  exact weighted Cauchy-Schwarz reduction of every positive heat-envelope
  L1 term to a weight mass and an analytic quadratic integral. A 96-bit Arb
  pilot rigorously encloses the order-nine n=7..16 block on u in [0,11/5].
  It does not include n>=17, the outer-u tail, the remaining derivative
  matrix, or retained first-jet lower balls, so it is not an exact r_6 bound

jensen_window_pf_newman_theta_arbitrary_n_stable_remainder_gate:
  exact all-N successor to the stable forward-remainder gate. For every
  N>=1 and q<=9, the cap-free identity writes r_N^(q) as the finite
  switch-defect derivative sum plus T_(N+1,q). Every K>=2 tail has a closed
  coefficient-norm majorant, with uniform endpoint checks
  16*pi-45>0 and 11-5*pi<-4. This removes the fixed forward-tail cap but
  does not by itself bound the arbitrary-N switch defect, assemble full
  d0/d1 budgets, or prove retained first-jet separation

jensen_window_pf_newman_theta_switch_defect_explicit_constant_gate:
  exact all-N direct-tail successor. Kernel polynomials and switch Laurent
  recurrences through order nine reduce every derivative to forward,
  compact-reflected, and large-reflected positive templates. The large
  phase retains c_safe=3*3^(2/3)*pi^(2/3)/16, and directed upper incomplete
  gamma integrals plus a unimodal sum lemma give explicit bounds for every
  required I_(p,k). Exact Hermite composition yields cap-free m=9 d0/d1
  upper budgets at T=1/5 for every N>=2. A 256-bit Arb compiler and an
  independent implementation agree on eight witnesses. This closes the
  omitted-tail upper-budget obligation, not the adaptive retained J_N/J_N'
  lower margin or terminating transition cover

jensen_window_pf_newman_theta_adaptive_gamma_scale_tail_gate:
  exact cofinal successor for N(x)=ceil((1+x)^(3/4)). Gamma-normalized
  errors increase inside each adaptive cell, while every forward,
  compact-reflected, and large-reflected endpoint template decreases from
  N=63 onward. Directed N=63 witnesses are below 0.214 and 0.217 of the
  natural exp(-pi*x/8) scale. Therefore both omitted value and derivative
  errors are below one quarter of that scale for every 0<=t<=1/5 and
  x>=245. The retained quarter-gamma lower margin remains open

jensen_window_pf_newman_theta_forward_six_term_bridge_tail_gate:
  exact bounded-bridge direct-tail theorem. On u>=0 the ordinary theta
  summands are positive, so retaining n=1,...,6 and bounding K>=7 gives
  explicit E0/E1 Gaussian arithmetic tails without using the invalid
  termwise-even integration-by-parts shortcut. At the worst endpoint x=245,
  the directed J and J' errors are below 9.646e-13 and 9.804e-13 of
  exp(-pi*x/8), hence below 10^-11 of that scale throughout
  38<=x<=245. Five retained terms are an explicit failure guard. This is an
  omitted-tail theorem only, not retained first-jet separation

jensen_window_pf_newman_theta_forward_adaptive_sqrt_tail_gate:
  exact tunable cofinal direct-tail theorem. For x>=245 and h>=0, retain
  ordinary theta terms through
  N_(F,h)=ceil(sqrt(x/8+2log(1+x)+h))-1. The Gaussian exponent, ceiling
  geometry, and two decreasing directed envelopes prove both omitted
  first-jet errors below exp(-3h-pi*x/8)/50. For fixed h the count is
  asymptotic to sqrt(pi/2) times the Riemann-Siegel saddle count, replacing
  the earlier x^(3/4) modular count on the preferred error side. Increasing
  h follows any prescribed positive tolerance, so the approximation does
  not impose a t-independent floor near t=0. Retained arithmetic
  transversality remains open

jensen_window_pf_newman_theta_forward_sqrt_to_corrected_rs_C1_transfer_gate:
  exact critical-layer C1 transfer. The closed normalizer amplitude gives
  A_t>=exp(-pi*x/8)x^(7/4)/32 and converts the adaptive ordinary tail into
  normalized errors E_F0=exp(-3h)/(25x^(23/4)) and
  E_F1<0.53E_F0. Thus
  T_L[O_h]>2exp(-6h)/(625x^(23/2)) is a direct sufficient no-contact
  target. The exact dual split with the corrected Riemann-Siegel main gives
  scaled first-jet distance below 5600exp(-3L/4). This algebraic transfer
  does not preserve the tiny ordinary threshold: the two currently proved
  envelopes differ by a certified cost O(exp(5L+3h)). The bulk A+B
  envelope itself is already only at the exp(-3L/4) scale, so higher
  endpoint corrections alone cannot close the gap. Retained arithmetic
  separation remains open on either quantitative route

jensen_window_pf_newman_theta_forward_six_term_finite_bridge_interval_certificate:
  rigorous 352-bit finite Q207 theorem. One shared order-(30,24) bivariate
  Taylor model per half-unit x panel encloses the first six ordinary theta
  terms; absolute moments through order 86 include an analytic u>11/5 tail.
  The append-only hash chain certifies all 414 panels and all 7,102 leaves on
  [1/1035,1/155]x[38,69] union [1/1035,1/5]x[69,245], with no subdivision or
  unresolved leaf and minimum retained-to-tail ratio above 2.2234e25.
  Composition with the compact and Q31 theorems proves no contact and zero
  winding on Q_207=[1/1035,1/4]x[0,245], hence closes Q_1 through Q_207.
  This finite theorem is not a fixed-count cofinal promotion

jensen_window_pf_newman_theta_q31_margin_geometry_audit:
  rigorous derived finite diagnostic from a fresh replay of all 1,240 Q31
  cache rows and 1,566 terminal leaves. The weakest tail-relative ratio box
  contains a retained zero and has strictly positive derivative, so the
  1.45e20 ratio reflects the roughly 10^-23 tail bar rather than an intrinsic
  cofinal margin. In s=sqrt(x/(4*pi)), partial_s=4sqrt(pi*x)partial_x, the
  fixed Q31 cover has a positive finite saddle-jet floor. Q31 keeps N=7,
  while even the kappa=1 x^(3/4) count runs from 16 to 25 there, so neither
  the raw ratio nor the finite floor is extrapolated

jensen_window_pf_newman_strict_laguerre_correlation_target:
  open RH-equivalent theorem target; the first Laguerre expression is exactly
  the Fourier transform of K_(1,t)(v)=integral phi_t(s+v)phi_t(s-v)s^2 ds.
  Positive-boundary attainment, Wiener, and the published Lambda<=1/5 bound
  imply Lambda<=0 iff these kernels' translates are dense in L1 for every
  0<t<=1/5. The kernels and every correlation are now proved uniformly
  strongly log-concave and admissible on that interval, but K_(0,0) gives an
  exact Xi-specific shape-and-tail counterexample with double Fourier zeros.
  The theta-curvature reduction supplies exact component formulas and rational
  B_0/B_1 pointwise contact bars through characteristic order three. A rational
  origin estimate and 1,900-box Arb/Taylor certificate now promote those bars
  to a uniform no-contact theorem for 0<=t<=1/5 and |x|<=38. Their polynomial
  growth still blocks global promotion. The exact fixed-block obstruction now
  proves that no static finite truncation can be cofinal. A positive even
  modular partition supplies arbitrary-power Fourier tail budgets. The saddle
  law fixes a sqrt(x) spectral transition, while the exact derivative-envelope
  theorem proves the larger N=max(N_sad,ceil(K*(1+x)^(3/4))) absolute-tail
  architecture. The direct-tail explicit-constant gate supplies referee-ready
  arbitrary-N m=9 upper budgets, and its adaptive successor proves quarter-
  gamma omitted-tail errors for all x>=245. A direct six-term theorem and
  complete 414-panel, 7,102-leaf interval cover close the finite bridge and
  Q_207.  The later closed-boundary winding theorem separately closes Q_208,
  making Q_209 the first unresolved linear stage with two exact open shell
  regions. The ordinary tail now also has a tunable cofinal square-root
  successor: its fixed-h count is only sqrt(pi/2) times the Riemann-Siegel
  saddle count and both errors are below exp(-3h-pi*x/8)/50. Thus the
  x^(3/4) count is no longer needed on the preferred error side. The
  high-frequency route now lacks only retained arithmetic transversality.
  On the critical overlap, exact normalizer conversion gives the direct
  ordinary target T_L[O_h]>2exp(-6h)/(625x^(23/2)) and an algebraic
  transfer to the corrected main. The current corrected remainder is
  certified only at the much coarser exp(-3L/4) amplitude scale, so the
  corrected joint small-ball or phase-critical-value theorem is an
  alternative quantitative target, not a free inheritance of the smaller
  ordinary threshold.
  Weighted Cauchy-Schwarz and the finite compact Arb certificate separately
  supply sharper finite-N budgets through Q31; only the retained cofinal
  first-jet lower profile remains open on this branch.
  The stronger subroute -L_t'>0 is now rigorously rejected at the Lehmer
  stress point and, by continuity, at sufficiently small positive time. The
  target still needs direct s^2-weighted or hierarchy-coupled structure, or
  corrected C1 double-zero transversality

jensen_window_pf_newman_theta_curvature_probability_operator_gate:
  exact theta-primitive probability and operator theorem; every summand of
  S_t(u)=exp(tu^2)R(u) is strictly decreasing and convex for
  0<=t<=1/5. Thus dmu_t=2S_t''du is a probability and the primitive cosine
  transform C_t=(1-A_t)/(2x^2) is strictly positive. Its theta weights are
  independent of time, with w_1>2799/2800 and explicit uniform C0-C3 tail
  budgets for all n>=2. Modular odd-jet identities nevertheless make that
  tiny tail cancel a first-block fifth endpoint jet above 300, blocking global
  promotion from mass dominance. The endpoint-subtracted operator factors as
  D_t=-(2t*d_x-x)^2-1, and H_t becomes an explicit
  second-order characteristic expression J_t/(16x^4). The ratio
  G_t=8H_t/C_t obeys a weighted Doob diffusion with nonnegative instantaneous
  Sturm-Liouville generator, but a quadratic guard shows that its reverse-time
  heat semigroup can lose two nodes at a double contact. Modular endpoint
  asymptotics make C_t eventually log-convex, blocking global
  contracting-drift and positive-curvature shortcuts. The fixed theta split
  gives exact formulas for J_(n,t), J_(n,t)' and rational B_0/B_1 bars, so
  either strict first-block inequality excludes contact pointwise. Their raw
  polynomial growth prevents global promotion. Boundary collision is exactly
  J_t=J_t'=0; at t=0 this is tangential contact with the unit
  exponential characteristic function. A smooth Neumann-lift model preserves
  all generic normalizer properties while retaining a double Fourier zero, so
  arithmetic C1 separation beyond the certified compact window remains open

jensen_window_pf_newman_theta_compact_transversality_interval_certificate:
  rigorous compact interval theorem. Positivity of the full Xi kernel and
  exact rational moment bounds give H_t(x)>248/371925 for |x|<=1/4.
  On 1/4<=x<=38, 160-bit Arb retained integrals, an analytic u>2 tail below
  10^-800, and second-order-time/third-order-frequency Taylor enclosures
  certify |J_1|>B_0 or |J_1'|>B_1 on 1,900 rational leaf boxes. Only ten
  one-level subdivisions are required, no box is unresolved, and the weakest
  rigorous normalized ratio exceeds 1.7521852560914668. Hence
  (H_t,H_t')!=(0,0) for every 0<=t<=1/5 and |x|<=38. The bounded-L band
  38<|x|<4*pi*exp(50), finite existential-threshold shoulders, and the
  residual L>=50 corrected phase layer 0<tL<=c_*+o(1) remain open, as do
  Lambda<=0 and RH

jensen_window_pf_newman_backward_pick_collision_bridge_audit:
  exact countermodel and literature gate. The even Cartwright heat flow
  E_t(z)=z^2+a-2t has top-time real zeros and Pick-positive logarithmic
  derivative, but a double collision births nonreal zeros and a negative-Pick
  half-disc backward in time. The zero speed diverges at the collision, and a
  fixed lower cutoff y>=rho hides the negative bubble for a time rho^2/2.
  Thus cutoff exhaustion and collision bridging do not commute. This
  identifies the missing separation/cutoff uniformity in the cited 2025
  backward-Pick preprint without ruling out a genuinely Xi-specific repair

jensen_window_pf_newman_correlation_hierarchy_gaussian_mixture_gate:
  exact hierarchy and route-elimination gate; all Xi correlations evolve by
  partial_t K_n=2v^2 K_n+2K_(n+1), and a multiple boundary root has a rigid
  hierarchy contact, with F_2=3F_1'' at a double root. Complete monotonicity in
  v^2 would force strict Fourier positivity, but the exact double-exponential
  Phi tail makes K_(1,t) decay faster than every Gaussian and therefore rules
  out Gaussian mixtures and direct PF-infinity membership. A separate
  55-digit scout detects the corresponding negative log-convexity minor

jensen_window_pf_newman_positive_time_strong_logconcavity_gate:
  exact published-theorem composition; strict concavity of
  log(Phi(sqrt(r))) gives (log Phi)''<=-kappa globally, and an Arb origin
  certificate gives kappa=74.9076... with kappa/2>37.45. Hence phi_t and every
  K_(n,t) are uniformly strongly log-concave throughout 0<=t<=1/5, the latter
  by Prekopa marginalization. The exact identity Fourier[K_(0,t)]=2H_t^2 and
  Hardy's theorem show that this full Xi-specific shape package still permits
  double Fourier zeros; a K_(1,t) proof must exploit its s^2 weight or hierarchy

jensen_window_pf_newman_weighted_strong_logconcavity_countermodel_gate:
  explicit theta-tail, root-concave admissible weighted-correlation countermodel;
  exp(-u^2-delta*u^4-epsilon(cosh(4u)-1))*(1+u^4/10) has curvature at most
  -(2-6/sqrt(10))-12delta*u^2-16epsilon*cosh(4u). At delta=1/10 and
  epsilon=1/1000 it has strict root-variable concavity and theta-type
  double-exponential decay. A 256-bit Acb/Arb certificate gives
  L_1[F](21/5)=-0.0002134784805...<0, so its first s^2-weighted correlation is
  not positive definite. Xi arithmetic or modular coupling, not even this
  Xi-style shape-and-tail package plus generic weighting, is needed

jensen_window_pf_newman_strict_laguerre_monotonicity_scout:
  rigorous rejection of the stronger monotonicity route. The implication
  M_t=-L_t'>0 => L_t>0 is exact, as are its sine-transform and theta-primitive
  forms, but the premise is false for Xi. At x=1401016343/100000, a 256-bit
  Arb third jet certifies M_0/A_0^2=-3.2418674697...e-6<0 and
  M_0/L_0=-0.0001463088540...<0. Continuity preserves the negative sign for
  sufficiently small positive heat times. The close-pair identity
  M[F](m)=-4a^2G(m)G'(m)+a^4(G(m)G'''(m)-G'(m)G''(m)) explains why 12000
  moderate rows and 20 selected rows through x=200 missed it. This retires
  global strict decrease only; direct L_t positivity and corrected C1
  double-zero transversality remain live

jensen_window_pf_newman_theta_summand_spectral_square_gate:
  exact theta/Mellin and finite-block obstruction; Phi=(R''-R)/8 with
  R'(0)=-1/2 gives H_t=1/16+D_t[C_t]/8 for an explicit second-order D_t and
  reduces L_t to one curvature expression in A_t=D_t[C_t]. The bilateral
  shifted-profile transform reconstructs xi itself. Every finite theta
  truncation retains a positive odd endpoint jet A_N and consequently has
  first Laguerre tail -2A_N^2/x^6+O(x^-8)<0. The infinite modular cancellation
  is essential; finite self/cross spectral squares cannot prove the target

jensen_window_pf_newman_gasper_fake_xi_remainder_gate:
  exact Gasper comparison and route-elimination gate; the scaled fake-Xi block
  is a genuine real-zero benchmark. Arb midpoint certificates with explicit
  tails prove that the best scalar absolute-remainder budget at x=25 exceeds
  one at t=0 and t=1/2; 80-digit values are 2.5852... and 2.5692.... Exactly,
  the kernel ratio Phi/Psi is greater than one at the origin and tends to one
  at infinity, so it cannot be a nonconstant positive cosh-mixture. The
  one-block triangle estimate and direct Cardon convolution transfer are
  closed; sign-aware and multi-block routes remain

jensen_window_pf_newman_gasper_residual_two_block_gate:
  exact and interval-certified two-block obstruction; Phi-Psi_9 is strictly
  positive, but the transform of the residual fails the first Laguerre test.
  More generally, tail positivity confines the 9/5 coefficient beta to
  0<=beta<=pi-3/(2pi). Acb Cauchy derivative bounds and Arb convexity cover
  that full interval with L[R_beta](66)<0 or L[R_beta](50)<0. Larger beta has
  a negative residual tail, and the tail-matched multiplier also fails the
  standard universal-factor imaginary-zero test. Signed coupled blocks remain

jensen_window_pf_newman_classical_three_block_residual_gate:
  exact and interval-certified classical three-block obstruction; the Polya
  P2 and de Bruijn real-zero kernels both leave strictly positive residuals,
  but Arb certifies negative first Laguerre values at x=86 and x=52. Tail and
  origin inequalities confine every nonnegative 9/5/1 block with a globally
  nonnegative residual to a rational triangle. Three Acb spectral certificates
  cover all 64,908 closed parameter boxes in that triangle. Gasper's published
  single-shift square does not sign the mixed products of a multi-shift sum;
  signed higher blocks or a new coupled square remain open

jensen_window_pf_newman_signed_universal_factor_residual_gate:
  exact and interval-certified signed universal-factor obstruction; the full
  signed 9/5/1 Polya multiplier condition is equivalent to one quartic having
  all roots in [0,1]. Exact endpoint and critical-point constraints plus the
  Xi origin bound place every globally positive residual candidate in one
  rational rectangle. Acb residual tests at x=86 and x=122 and exact quartic
  and critical discriminants classify 4,094 adaptive leaves grown from 3,416
  base boxes, with maximum depth six and no unresolved box. Independent LP
  residuals are closed for the whole standard signed three-shift multiplier
  cone; higher shifts and genuinely coupled mixed-term squares remain open

jensen_window_pf_newman_theta_bessel_higher_shift_regularization_gate:
  exact higher-shift and non-Fubini gate; every symmetrized theta summand has
  an absolutely transformable fixed-index Bessel-K expansion with coefficient
  c_(n,m)=pi^m*n^(2m)*(2*pi^2*n^4-3m)/m!. The n=1 coefficients turn negative
  at m=7. Yet the zero-frequency fixed-block transforms sum as a negative
  multiple of sum n^(-1/2), while the complete Xi transform is finite. Naive
  termwise higher-shift summation is closed; modularly grouped mixed-term
  matrix identities or sign-preserving renormalization remain open

jensen_window_pf_newman_theta_cell_renormalization_gate:
  exact endpoint renormalization; subtracting each continuous theta-index cell
  leaves Phi unchanged and produces normal convergence in every polynomially
  weighted L1 norm. The resulting transforms give Euler's convergent zeta
  sum, positive contributions at x=0, and an absolutely convergent coupled
  Laguerre matrix. Each block nevertheless retains a nonzero exp(-5u) tail,
  so positive Newman weighting makes it nonintegrable; a t-compatible modular
  grouping remains open

jensen_window_pf_newman_polymath15_oscillatory_zeta_handoff_theorem:
  asymptotic theorem with exact threshold
  c_*=4911678521/1933561194=2.540223984760008...; an eleven-exponent-pair
  envelope and the Polymath-15 remainder transfer prove exact first-Laguerre
  positivity for every fixed tL>=c_*+epsilon once L is sufficiently large.
  It supplies no practical threshold and does not cover tL<=c_*+o(1)

jensen_window_pf_newman_polymath15_cancellation_zero_free_wall_gate:
  exact weighted-block frontier and route-classification gate; it reproduces
  c_*=4911678521/1933561194, maps c=2 to the zeta 1-line, and gives a strict
  radius-proportional cancellation condition sufficient for a conditional
  c=2 zeta handoff. Fixed c<2 enters Re(s_*)<1, where current zero-free
  inputs do not supply the nonzero phase-amplitude floor; the inner
  Wronskian theorem remains open

jensen_window_pf_newman_polymath15_gaussian_legendre_duality_gate:
  exact finite Gaussian-shift identity and Legendre-equivalence theorem; using
  the current pointwise partial-sum profile in the Gaussian heat average gives
  exactly the existing weighted dyadic frontier. The equality point is
  q_*=4800718975/7734244776 at c_*, and c=2 retains the exact exponent deficit
  3133668399/48144906818. A semigroup rewrite alone therefore cannot improve
  the threshold; a genuinely stronger profile or non-pointwise cancellation
  theorem is required

jensen_window_pf_newman_polymath15_antedb_beta_frontier_audit:
  exact current-source audit at ANTEDB commit
  99668603896af86e6cda90ed6755cf3116aab0ac. Four Tao-Trudgian-Yang and two
  Cushing post-2023 pairs improve intermediate radii, while the exact current
  hull and one-pass direct-beta envelope retain alpha_*=62831/155153,
  beta_*=220633/620612, and c_*=4911678521/1933561194. Twelve finite beta-only
  van der Corput passes preserve the same contact. This closes the stale
  2023-table concern for the pinned finite audit but supplies no stronger beta
  bound, infinite transform closure, c=2 theorem, or inner Wronskian separation

jensen_window_pf_newman_polymath15_lambda01965_provenance_audit:
  provenance and nonpromotion gate separating the peer-reviewed interval
  0<=Lambda<=1/5 from a June 2026 preprint claiming Lambda<=393/2000. The
  deposit MD5, all 505 internal manifest entries, exact row arithmetic, and
  22 portable certificate invocations were checked locally. Three FLINT/Arb
  producer packages, independent theorem-to-code review, and peer review were
  not supplied, so 0.1965 remains a reproduced unrefereed candidate rather
  than established literature. Even acceptance would only shrink the positive
  interval by 7/2000 and would not prove Lambda<=0 or RH

jensen_window_pf_newman_polymath15_critical_scaled_coercivity_target:
  open corrected finite Riemann-Siegel curvature target on the remaining
  asymptotic strip 0<=tL<=c_*+o(1); the exact correction and C2 remainder
  transfer are available, but the required phase-sensitive sign is not

jensen_window_pf_newman_polymath15_critical_component_wronskian_gate:
  exact corrected-component Wronskian reduction. Its Hermitian rate matrix has
  rank at most two and, for unequal component rates, inertia (1,1,m-2).
  An exact ordered-negative-speed countermodel still has a nonzero double
  real-part crossing, while four corrected t=0 diagnostics show both signs
  and cancellation above 700. This blocks an unrestricted positive-
  semidefinite rate-matrix proof but leaves the joint arithmetic small-ball
  theorem open

jensen_window_pf_newman_one_sided_phase_moment_bridge_gate:
  exact zero-free one-sided lift and score-probability bridge. It turns double
  contact into the joint equations E[sin(xU)]=E[U cos(xU)]=0 with no
  exceptional complex-lift zero, identifies the signed-Hankel coefficients as
  normalized odd score moments, and gives every shifted Jensen window as an
  explicit positive generalized-Laguerre scale mixture. An independent
  uniform/Beta pushforward turns that mixture into one positive Abel measure
  with moments `M_k=k!A_k`, Bessel transform `2H_t(sqrt(-z))`, and closed
  tail-source flow. Its quadratic Jensen condition is the sharp lower
  concentration ratio `(n+1)/(n+2)`; measure positivity supplies only the
  opposite-side Cauchy-Schwarz upper bound. Full support also rules out a
  common interlacer for the complete fixed-scale Laguerre family in every
  degree at least two, while leaving a weighted Xi-specific total-positive
  connection open. Two exact minors of the fractional Abel kernel have
  opposite signs, so the bare Abel operator is neither TP2 nor sign-regular
  of order two; the Xi weight must participate in any surviving composition.
  The sine observable obeys a closed radial backward-heat equation. Restoring
  the exact lift amplitude makes its first jet uniformly equivalent to the
  normalized Polymath first jet on the corrected overlap. The raw phase is
  rapidly flat and an exact smoothed-triangle score countermodel has a double
  characteristic zero, so the Xi-specific joint-avoidance theorem remains
  open

jensen_window_pf_newman_score_abel_radial_dimension_lift_gate:
  exact all-shift radial geometry of the score/Beta Abel measure, with the
  proved Xi bounded-degree layers reconciled explicitly. For each
  shift `n`, the density
  `r_t(|y|^2/4)/((4*pi)^(n+1)A_n)` has mass one on `R^(2n+2)` and Fourier
  profile `mathcal F_t^(n)(-|xi|^2)/A_n`; at `n=0` its coordinate marginal
  is the normalized Newman kernel. Successive profiles obey the classical
  radial dimension walk, so every shifted Jensen window belongs to one
  explicit dimension of the ladder. This does not satisfy Schoenberg's
  fixed-profile all-dimensions premise. A two-Gaussian Newman flow has the
  entire ladder and arbitrarily strong root-variable log-concavity, but its
  generating function has an explicit nonreal zero lattice and every shifted
  quadratic Jensen discriminant is negative; its squared-variable profile is
  strictly log-convex. The Xi profile is different: strict log-concavity of
  `Phi_t(sqrt(y))` plus Berwald-Borell closes every shifted quadratic, and the
  reciprocal-defect heat theorem closes every shifted cubic throughout
  `0<=t<=1/5`. The independent real-zero-band/sector theorem closes every
  shifted Xi degree through 361 on that interval. Thus generic radial positivity
  remains decisively insufficient, but the earlier quartic and quintic
  frontier is no longer open for Xi. The strong-log-concave Mellin witness and
  length-14 outer-contact survivor remain exact generic nonpromotion guards.
  The first unproved direct degree is 362; the terminal need is an unbounded
  cofinal degree sequence or an all-degree theorem.

jensen_window_pf_newman_zero_slab_degree71_sector_certificate:
  rigorous bounded-degree continuum theorem. A 192-bit directed
  near-minus-tail estimate proves `Re H_t(x+iy)>0` for
  `0<=t<=1/5`, `|x|<=84/5`, and `|y|<=1`. De Bruijn's strip theorem and the
  map `F_t(s)=2H_t(i*sqrt(s))` put every zero of `F_t` in the negative-axis
  sector with `sin(delta)=840/7081`. Chasse's sector theorem and derivative
  closure then prove every shifted Xi Jensen polynomial through degree 71
  hyperbolic throughout the target heat interval. This excludes the
  quartic/quintic survivor for Xi but proves neither degree 72 nor any
  unbounded/all-degree conclusion.

jensen_window_pf_newman_real_zero_band_degree361_sector_certificate:
  rigorous bounded-degree continuum theorem. A 192-bit Arb/Taylor enclosure
  proves both vertical sides of `|Re z|<=38`, `|Im z|<=1` zero-free. The
  de Bruijn strip, the classical endpoint horizontal boundary, published
  `Lambda<=1/5` top-time reality, and the independent no-real-collision
  certificate then prove by homotopy that every central-band zero is real.
  The map `s=-z^2` gives `sin(delta)=76/1445`, and Chasse's sector theorem
  proves every shifted Xi Jensen polynomial through degree 361 hyperbolic
  on `0<=t<=1/5`. Degree 362 and every unbounded/all-degree conclusion
  remain open.

jensen_window_pf_laguerre_scale_mixture_gate:
  exact kernel theorem and preservation guard; each unshifted Jensen window is
  a positive scale integral of L_D^(-1/2), and each fixed scale is hyperbolic.
  Exact two-atom and exponential-density countermodels show that positive
  mixing and log-concavity do not preserve this property. Half-integer Gamma
  mixing is nevertheless hyperbolic in every degree by Euler-Jacobi
  factorization. Full support makes global common interlacing of all fixed
  scales impossible in degree at least two, isolating a weighted
  Xi-specific total-positive connection as the surviving version

jensen_window_pf_multiplier_counting_measure_target:
  rejected sufficient subclass; the order-6 moment magnitude forces every
  unit atom above 4.863538496, while the next-moment ratio exceeds the maximum
  atom ratio allowed above that cutoff. The resulting interval contradiction
  retires the unit-atomic product without rejecting general multiplier
  sequences or the other all-order Jensen/PF routes

jensen_window_pf_multiplier_unit_atomic_obstruction_certificate:
  interval-certified route closure; combines the atom cutoff with exact
  monotonicity of g_(m+1)(alpha)/g_m(alpha) and an Arb ratio gap above
  7.68e-4 to rule out the convergent unit-atomic elementary multiplier product
  for the normalized lambda-zero zeta coefficients

jensen_window_pf_mellin_multiplier_power_sum_obstruction:
  interval-certified continuous-interpolation guard; eleven Arb-enclosed Phi
  log moments produce nine candidate power sums and six shifted Hankel tests.
  Three determinants are strictly negative, including the shift-2 size-4 value
  near -2.588644974358276e-19. This rejects the natural Gamma-normalized Mellin
  product, but not a product identity asserted only at integer indices

jensen_window_pf_monotone_contraction_stress:
  finite monotone-contraction stress diagnostic; validates 2875 finite
  Arb-classified zeta-window rows across degrees d=3..12 and the k<=64 cache,
  all satisfying adjacent log-concavity and increasing ratio contractions
  with no zero-containing required positivity gaps

jensen_window_pf_edrei_heat_flow_boundary_gate:
  exact radial Burgers/Edrei-moment heat hierarchy and exact repeated-zero LP+
  boundary model; every shifted rank-one 2x2 wall points inward under forward
  heat, while the quadratic orbit has nonreal zeros and negative shifted
  Stieltjes minors immediately under backward heat

jensen_window_pf_edrei_hankel_boundary_flux_gate:
  exact polynomial divided-difference flux for every shifted Edrei Hankel
  form; integer residues make finite-rank walls forward-nonnegative and
  repeated factors strict. Infinite support forces canonical boundary
  witnesses to unbounded order. An exact orthogonal-polynomial Schur
  quotient then proves that a repeated atom gives an exponential high-shift
  determinant-flux signature. The required Phi-specific subexponential
  bound remains open

jensen_window_pf_edrei_raw_moment_collision_resolution_gate:
  exact Phi raw-moment to Taylor to Edrei-moment transfer and exact symmetric
  LP+ split model. At a repeated leading atom the shifted determinant has
  condition number of order (beta/q)^s, while a gap epsilon is not resolved
  until shifts of order 2*log(beta/epsilon)/log(beta/q). Any fixed string of
  larger simple atoms preserves that scale at arbitrary determinant rank.
  Every fixed-rank collision indicator has noncommuting split and high-shift
  limits, so independent algebraic-accuracy raw-moment estimates cannot prove
  the open Phi determinant-flux bound

jensen_window_pf_phi_pick_kernel_target:
  exact Krein/Sokal endpoint equivalence between LP+ and positivity of
  P_Phi(z)=-Im(F_0'(z)conj(F_0(z))) throughout the upper half-plane;
  polarization gives an explicit Phi double integral, while the first-wall
  variance identity and R'(0)=7/19200 two-scale witness block positive-mixture
  promotion

jensen_window_pf_xi_pick_suzuki_hankel_bridge:
  exact normalization M(w)=xi((1+w)/2)/4 and hyperbolic Pick direction;
  Suzuki's published criterion replaces an unspecified positive resolvent by
  explicit all-t Fredholm nonvanishing and a terminal canonical-kernel
  condition. The determinant-only reduction makes the latter cofinal-sequence
  redundant. The continuum Fredholm expansion gives signed-Hankel contact,
  while exact horizontal-growth, rank-one TN, and contour-residue
  countermodels block three tempting shortcuts

jensen_window_pf_suzuki_spectral_frontier:
  exact fixed-space truncation and Hilbert-Schmidt formula; continuity from
  zero and norm monotonicity make all-history determinant nonvanishing
  equivalent to pointwise strict contraction. Hilbert-Schmidt domination and
  Suzuki's sup-contour certificate have finite reach, while the desired HB
  case forces the finite-t norms to approach one

jensen_window_pf_suzuki_determinant_only_reduction:
  compatible all-time finite Hankel forms extend to a bounded causal
  convolution multiplier. Its H-infinity continuation forces any off-line
  zero to generate distinct shifted zeros accumulating at itself along a
  cofinal omega sequence. Thus the determinant family alone is equivalent
  to RH, but the all-time determinant premise remains open

jensen_window_pf_suzuki_fixed_omega_phase_diagram:
  exact fixed-shift equivalence between all-time determinant nonvanishing,
  meromorphic innerness, and horizontal shifted-zero multiplicity
  cancellation. Outside a countable exceptional set this is the zero-free
  half-plane condition; failure of RH creates a whole small-omega failure
  interval. Suzuki's Jordan-totient L2 residual and eventual-sign condition
  give a scalar arithmetic form of the same open gate

jensen_window_pf_suzuki_jordan_totient_sign_scout:
  finite double-precision diagnostic with 6000 positive sampled values for
  omega=1/2 through 1/32 and x through 5000, plus a signed-weight guard; no
  eventual-sign, L2-tail, all-time determinant, RH, or Lambda<=0 promotion

jensen_window_pf_suzuki_cofinal_monotonicity_hierarchy:
  Suzuki's published logarithmic smoothing hierarchy and Landau implication,
  composed with the fixed-omega phase theorem, reduce RH to eventual
  one-sign behavior along any one sequence omega_j->0 at arbitrary selected
  smoothing orders k_j>=1. Every smoothing weight remains signed, and the
  cofinal equivalence is an internally audited theorem candidate rather than
  a proof of its arithmetic antecedent

jensen_window_pf_suzuki_jordan_error_kernel_reduction:
  exact Abel-summation reduction showing that Suzuki's weight annihilates
  the positive Jordan-totient residue main term and leaves only a signed
  Mobius summatory-error convolution. Elementary absolute bounds reach only
  O(x^(1/2-omega)log x), at every finite smoothing order, so the missing
  step is a genuine cancellation or order theorem

jensen_window_pf_suzuki_cofinal_l2_hierarchy:
  central Taylor subtraction turns every Suzuki smoothing order into an
  exact L2 residual criterion for fixed-shift innerness. Along any sequence
  omega_j->0, finite residual energy at arbitrary selected orders k_j>=1 is
  RH-equivalent. The corresponding Jordan-error energy estimate remains an
  open arithmetic antecedent, not a proof

jensen_window_pf_jordan_muntz_causal_energy_bridge:
  exact generalized Muntz expansion and finite Mellin-Plancherel formula for
  the unsmoothed Jordan-totient error. Its positive-time energy gives an
  internally audited fixed-shift/cofinal criterion, while finite natural
  Mobius partial sums expose the remaining all-height uniform mollifier
  inequality. The inequality remains open and RH-strength

jensen_window_pf_jordan_muntz_burnol_hardy_intertwiner:
  exact finite Hardy-operator equivalence between the Jordan-Muntz partial
  sums and Burnol's weighted natural fractional-part approximants at
  epsilon=2omega. Explicit two-sided norm constants source-back the cofinal
  criterion through Burnol and Balazard-Saias, while leaving the same
  uniform natural-mollifier estimate open

jensen_window_pf_burnol_cell_energy_tail_obstruction:
  exact reciprocal-cell decomposition of the weighted Burnol norm. Uniform
  boundedness is equivalent to a discrete discrepancy energy and forces the
  reciprocal-zeta partial sums at 1+2omega to have the RH-strength rate
  N^(-1/2-omega). The rate and the remaining short-interval cancellation
  are necessary consequences of the still-open uniform energy target

jensen_window_pf_burnol_tail_discrepancy_dyadic_reduction:
  exact split of the finite discrepancy into its limiting error, stable
  reciprocal-zeta tail, and post-prefix omitted-divisor fluctuation. The
  latter is equivalently a summable dyadic square function of weighted
  Mobius increments on (k/(q+1),k/q]. The square-function estimate remains
  open and RH-strength

jensen_window_pf_weighted_fractional_autocorrelation_gram_bridge:
  exact stationary Gram and spectral representation of the same weighted
  Burnol norm on the logarithmic divisor axis. Its diagonal is uniformly
  bounded; every possible growth lies in a signed Mobius off-diagonal.
  Positive definiteness supplies no uniform upper bound, so the weighted
  off-diagonal cancellation estimate remains open and RH-strength

jensen_window_pf_weighted_autocorrelation_ou_tail_energy_reduction:
  explicit period averaging splits the stationary kernel into a leading
  Ornstein-Uhlenbeck kernel plus an exponentially smaller pointwise
  remainder. The leading Gram is exactly a weighted reciprocal-Mobius tail
  energy, and boundedness along every member of one cofinal alpha sequence
  is RH-equivalent. The remainder density changes sign, so this does not
  compare the leading energy with the full Burnol norm or prove RH

jensen_window_pf_ou_mertens_mean_square_reduction:
  weighted Hardy/Copson identities identify the limiting OU tail energy
  with a weighted Mertens mean square. Dyadic decomposition makes its
  cofinal RH criterion equivalent to subpower normalized block energy and
  to one signed origin-anchored integrated Mobius-correlation estimate.
  Generic Hardy loses N^(1-alpha), while averaged-Chowla and almost-all
  short-interval theorems do not control the anchored cumulative slice

jensen_window_pf_mertens_weighted_prefix_affine_defect_reduction:
  the weighted Mobius prefix is an exact Volterra-Hardy coordinate for the
  same reciprocal-zeta energy. The first q=floor(k/d) block observes its
  detrended increments but misses a constant anchor mode, as an exact
  constant-tail countermodel proves. Applying Vaughan after correlation
  collapse is endogenous and recycles the Mertens energy; applying it
  before collapse gives the concrete signed handoff B_1(X)-TI_X+TII_X

jensen_window_pf_mertens_anchor_logarithmic_tail_energy_reduction:
  adjacent reciprocal tails reconstruct the weighted prefix exactly across
  all cutoffs. Weighted Hardy/Copson bounds identify its energy with
  sum_N N^alpha|r_N|^2 and hence with the logarithmic stable-prefix series
  sum_N P_(alpha/2,N)/N. A critical scalar countermodel proves that
  sup_N P_(omega,N)<infinity is one logarithm short; the Mobius-specific
  logarithmic gain remains open and RH-equivalent on a cofinal sequence

jensen_window_pf_mertens_dyadic_cosine_mode_bottleneck:
  exact Neumann cosine diagonalization of each dyadic reciprocal-tail
  block proves every mode above sqrt(K) summable from coefficient size
  alone. The cofinal RH-equivalent target is reduced to the mean and
  O(sqrt(K)) smooth modes. Davenport's arbitrary logarithmic savings lose
  the power K^(1-alpha); square-root additive-twist cancellation is only
  a sufficient conditional calibration and is not proved

jensen_window_pf_mertens_local_path_cosine_transference:
  the surviving modes are exactly the DCT of the anchored
  reciprocal-weighted local Mertens path, with an affine incoming-tail
  defect in the mean. Finite Abel transforms compare the full weighted
  and ordinary local Mertens paths; that full-path comparison alone does
  not control low projection. An exact non-Mobius scalar model has zero dyadic endpoint
  tails and zero means but a divergent first nonconstant mode, ruling out
  generic coefficient-size or large-sieve closure while leaving the
  Mobius-specific two-component gate open

jensen_window_pf_mertens_centered_bridge_vaughan_handoff:
  the complete coefficient Gram has summable diagonal trace, leaving a
  signed off-diagonal Mobius form. The Abel map preserves constants and
  is uniformly invertible on the quotient by constants, so the
  nonconstant low modes are equivalent to centered local Mertens energy
  with Brownian-bridge kernel min(i,j)-ij/K. Its off-diagonal correlation
  is exactly C_K=-TI_K+TII_K before rowwise absolute values. Averaged
  Chowla and separate Type I/II bounds remain one power short. The affine
  mean is not controlled here; the next exact lemma localizes and unifies
  it without supplying the missing cancellation

jensen_window_pf_mertens_affine_tent_bridge_handoff:
  the scale filter u_K-(1/2)u_(2K) cancels the affine mean's infinite
  plateau and is invertible after normalization. Its compact tent square
  and the Brownian-bridge energy combine into one local nonnegative kernel
  of size O(K^(-1-alpha)) on K<n<n+h<4K, with a summable diagonal and an
  exact two-interval Vaughan identity. Before ordinary-Mertens transfer,
  the affine rank-one term exactly fills the bridge's constant direction,
  yielding an equivalent positive compact tail-lattice square. A second
  finite Abel isomorphism removes the current-block power weights and
  produces an equivalent ordinary-Mobius suffix-square criterion with one
  explicit bounded-coefficient future anchor on [2K,4K). The resulting
  signed, weighted-positive, and unweighted-anchored criteria are all
  equivalent to the reciprocal-tail energy, but no arithmetic gain is
  established. The sufficient scale target is H_K=O(K^(2+epsilon)) for
  epsilon<alpha, whereas coefficientwise control gives only O(K^3).
  All-interval results lack the needed power, almost-all results leave
  the fixed dyadic bases exceptional, and no alpha=0 limit is permitted

jensen_window_pf_mertens_ordinary_vector_vaughan_handoff:
  the ordinary anchored suffix energy is exactly the squared norm of a
  K-dimensional Mobius feature sum with deterministic coefficients in
  [0,1]. Its explicit Gram kernel has current-current block min(i,j),
  current-future block i*b_K(n), and rank-one future block
  K*b_K(n)b_K(m). The complete diagonal is at most (5/4)K^2 and is
  summable after normalization. Applying finite Vaughan componentwise on
  the two standard intervals before squaring gives
  B_K=-V_(I,K)+V_(II,K), making the live target the weighted square norm
  of this signed vector difference. The identity is exact, but no
  vector Type I/II power gain is proved; coefficientwise control remains
  one full power short. Even the all-interval H*log^(-A)X theorem,
  applied coordinatewise, remains at K^3*log^(-2A)K

jensen_window_pf_mertens_mixed_boundary_sine_vaughan_handoff:
  the min(i,j) suffix Gram has an exact mixed Dirichlet-Neumann sine
  spectrum at theta_(K,r)=(2r-1)pi/(2K+1). Unfolding the future anchor
  turns each mode into one compact ordinary-Mobius sine test with
  O(K^(-1/2)) coefficients. The RH-equivalent energy is uniformly
  comparable to the half-odd square function
  sum_K K^(-alpha)sum_r |X_(K,r)|^2/(2r-1)^2, and finite Vaughan gives
  X_(K,r)=-T_(I,K,r)+T_(II,K,r) before squaring. Davenport, Parseval,
  and scalar all-interval logarithmic savings remain one power short;
  square-root-plus-epsilon cancellation would suffice but is unproved

jensen_window_pf_mertens_spectral_anchor_high_mode_reduction:
  the half-odd coefficient splits exactly into a current Mobius mode and
  one endpoint anchor. Parseval gives a K/R^2 current tail, while the
  endpoint envelope gives a beta_K^2/(KR) anchor tail. With only the
  trivial bound |beta_K|<=2K, the cutoff
  R_(alpha,K)=ceil(K^(1-alpha/2)) makes every high-mode contribution
  dyadically summable. The reciprocal-tail criterion is therefore
  equivalent, for each fixed positive alpha, to the first
  O(K^(1-alpha/2)) modes. This is not uniform at alpha=0, and the
  reduced low-mode square remains open

jensen_window_pf_mertens_truncated_half_odd_kernel_handoff:
  the surviving low-mode square is exactly the quadratic form of a
  truncated half-odd kernel. The kernel is entrywise positive and
  positive semidefinite, nests by rank-one increments, and completes
  at infinite odd rank to pi^2(2K+1)^(-2)min(i,j). Pullback to ordinary
  Mobius support gives a uniformly bounded per-scale diagonal, so the
  RH-equivalent reduced criterion is exactly a weighted signed
  off-diagonal bound. Its endpoint column is strictly increasing and
  gives an exact Abel average of local Mobius suffixes; finite Vaughan
  decomposition preserves the joint signed Gram form. No required
  weighted off-diagonal cancellation is proved

jensen_window_pf_mertens_shift_kernel_variation_handoff:
  the signed off-diagonal is exactly a sum over shifts with three
  consecutive current/current, current/future, and future/future
  blocks. Every kernel weight is O(1/K), and its discrete variation
  is O(1/K) along either a fixed shift or a fixed base point. Exact
  Abel formulas expose maximal two-point correlations and local
  Mertens intervals, while total and per-shift kernel masses are
  explicitly bounded. Absolute summation along either coordinate
  remains power-short even under conditional square-root input, so
  the live estimate must retain cancellation jointly across base
  point, shift, Vaughan type, and possibly scale

jensen_window_pf_mertens_planar_abel_handoff:
  zero extension and two-dimensional Abel summation turn the signed
  off-diagonal into a pairing of joint base/shift prefixes with mixed
  kernel curvature. The total mixed variation is bounded by
  3*pi^2*(1+log(2R))/K: the current interior costs the L1 norm of the
  finite odd Dirichlet kernel, while all joins and future blocks
  telescope through convex decrements. Weighted Cauchy-Schwarz isolates
  a curvature-energy target of size K^(1+epsilon). That target remains
  open; the stronger planar maximum condition already contains
  RH-scale block-Mertens input, and separate axiswise square-root
  prefixes do not imply it

jensen_window_pf_mertens_planar_edge_gram_vaughan_handoff:
  lifting each triangular-band pair to an edge indicator turns the
  planar curvature energy into an exact Gram quadratic form. Its
  diagonal is below (27/2)*pi^2*K*(1+log(2R)), so that entire part is
  already at the required K^(1+epsilon) scale. The remaining energy
  gate is therefore exactly the signed off-diagonal edge-pair bound,
  containing both one-vertex collisions and four-distinct-endpoint
  terms. Supported finite Vaughan decomposition and symmetrization
  preserve their joint cancellation; the U=V=1 guard shows why the
  four large directed pieces cannot be bounded separately. Existing
  large-sieve, Gowers-uniformity, and averaged-Chowla theorems do not
  supply this quadratic nested-edge estimate

jensen_window_pf_mertens_planar_incidence_antidiagonal_handoff:
  the complete signed edge-pair form has exact sorted-endpoint
  incidence formulas C=C_3+C_4, while a bounded-sign witness has C=0
  through perfect cancellation C_3=210 and C_4=-210. Random cuts give
  an exact bilinear representation of the strong energy gate but no
  gain. Before absolute mixed curvature, the current-current interior
  collapses to a signed odd-frequency anti-diagonal projection with a
  concrete square-function target for Y_(K,r). That target and the
  transition/current-future/future-future cells remain open; recent
  almost-all short-interval results are logarithmic and have the wrong
  averaging quantifiers

jensen_window_pf_mertens_planar_boundary_flux_anchor_handoff:
  one exact modewise discrete-flux law covers every planar curvature cell.
  The endpoint jump carries a q_r^2 smoothing factor, both artificial
  interfaces join before absolute values, and the top current-interior
  stencil contains a mandatory current-future shoulder A_CF. The complete
  signed decomposition is
  O=O_(CC,int)+B_CC-A_CF+beta*gamma+(c/2)(beta^2-Q), and its remaining
  terms form one endpoint-spectral projection sum_r Z_(K,r)/q_r^2.
  Orthogonal completion gives L=E_perp+c(beta+gamma/c)^2. These identities
  are exact; the Y, joined-Z, and lossless low-mode Mobius estimates remain
  open

jensen_window_pf_mertens_planar_boundary_suffix_energy_handoff:
  the complete joined boundary/future remainder is exactly
  R_B=(E_B-D_B)/2+C_delta. Here E_B is a positive anchored suffix energy,
  D_B is uniformly bounded, and the endpoint-smoothed defect is O(log R).
  Thus subpower control of the formerly separate transition/future blocks is
  equivalent to one explicit E_B weighted Mertens-energy estimate. Uniform
  odd-sine control also gives E_B<=pi(1/2+2pi)K^(-2)H_K for the earlier
  affine-tent anchored energy H_K. The reduction is exact, but the H_K/E_B
  and inherited Y_(K,r) estimates remain open

jensen_window_pf_mertens_planar_boundary_suffix_localization_gate:
  the odd-sine square-wave tail gives a uniform lower bound for every
  boundary increment outside an O(K^(alpha/2)) terminal collar. Bounded
  adjacent Mobius increments recover that collar, proving that all-epsilon
  subpower control of E_B along one fixed cofinal positive-alpha sequence is
  RH-equivalent. A bounded terminal-anchor witness has K^2 E_B/H_K tending
  to zero, so there is no uniform reverse norm comparison. The exact result
  classifies this boundary gate as RH-strength but supplies no E_B or
  inherited Y_(K,r) estimate

jensen_window_pf_mertens_planar_joined_energy_equivalence_gate:
  the candidate join J_B=O_(CC,int)+E_B/2 is exactly
  L/2+(D_B-D_L)/2-C_delta, where L is the inherited lossless positive
  low-mode energy, both diagonals are bounded, and C_delta=O(log R).
  Therefore J_B, the full signed off-diagonal, and L have equivalent
  all-epsilon and fixed-alpha weighted dyadic criteria. The joined route
  returns the old lossless gate and supplies no L, E_B, or Y_(K,r) estimate

jensen_window_pf_mertens_planar_curvature_energy_scout:
  finite float64 diagnostic on 28 rows with dyadic K through 1024 and
  alpha in {1/8,1/4,1/2,3/4}. Direct and double-Abel evaluations agree
  within 6.7e-16; on this grid max K*V_K is below 6.571, max E_K/K is
  below 0.831, and max M_K/K is below 2.096. This is a checked
  falsification scout only and proves no uniform or asymptotic estimate

jensen_window_pf_fixed_shift_inner_reciprocal_boundary_separation_gate:
  a real even four-zero toy completed function has a rational-inner
  fixed-shift quotient with finite exact causal H2 energy, while reciprocal
  square energy on the same shifted line diverges at the canceled boundary
  zeros. This blocks a generic same-shift promotion from limiting
  Jordan/Burnol innerness to the logarithmic stable-prefix energy, without
  blocking the cofinal full-Burnol implication

jensen_window_pf_newman_convex_phase_cell_unwrapping_lemma:
  a continuous closed first-jet path covered arcwise by cyclic convex
  origin-free cells is homotopic to an exact rational intersection-witness
  polygon. Its winding is the polygon's exact signed half-open ray-crossing
  count, with no floating-point angle or branch-cut decision. The globally
  regular proxy F=16(1+x^4)H is joined to the ordinary first jet by a
  positive-determinant triangular homotopy. All 20 stored Q208 right-edge
  cells transfer to F_t'(246)>0. The complete bottom/top cell chains,
  cyclic Q208 winding, Q208, the cofinal family, Lambda<=0, and RH remain
  open

jensen_window_pf_newman_q208_bottom_phase_cell_certificate:
  492 fixed-time cells rigorously cover t=1/1040 and 0<=x<=246 for
  F=16(1+x^4)H. All cells exclude the origin without subdivision; 313 use
  value separation and 179 derivative separation. Exact dyadic witnesses
  give a 493-vertex open chain with signed crossing contribution -19

jensen_window_pf_newman_q208_top_phase_cell_certificate:
  492 fixed-time cells rigorously cover t=1/5 and 0<=x<=246 for the same
  proxy. All cells exclude the origin without subdivision; 319 use value
  separation and 173 derivative separation. Its forward open-chain crossing
  contribution is also -19, and the Q208 boundary traverses it in reverse

jensen_window_pf_newman_q208_closed_boundary_winding_certificate:
  492 bottom, 20 transformed-right, 492 reversed-top, and one positive-axis
  cell form a cyclic 1005-cell cover with 1005 exact dyadic witnesses. The
  exact signed crossing count is zero. Sign-definite heat-contact index and
  the independently proved Lambda<=1/5 bound therefore certify Q_208 and,
  by containment, Q_1 through Q_208. No stage j>=209 or cofinal theorem is
  proved

jensen_window_pf_newman_cofinal_phase_cell_scaling_diagnostics:
  exact shell identities place the cofinal bottom in the tL->0 layer, outside
  the proved dominant-saddle ray. Finite source audits give 408/414 Q207/Q208
  and 470/492 Q208 bottom/top branch agreements, no half-unit witness turn
  reaching pi/2, and at least 5.85e24 tail domination. The scale comparisons
  are diagnostics only and prove no all-j first-jet theorem

jensen_window_pf_newman_adiabatic_phase_cell_successor_lemma:
  P_(j+1) is exactly P_j plus a global bottom collar of time width
  1/(5j(j+1)) and a one-unit right strip. If every bottom phase cell absorbs
  its heat-jet displacement and the new strip lies in one strict first-jet
  half-plane, contact-freeness and zero degree pass from P_j to P_(j+1).
  The induction is exact and conditional; the Xi all-j antecedent is open

jensen_window_pf_newman_q207_q208_adiabatic_bottom_collar_refined_tail_certificate:
  379 half-unit prefix panels plus 222 order-six quarter-unit tail panels
  rigorously certify [1/1040,1/1035]x[0,245] by heat-jet transport from the
  new-edge Q208 bottom phase cells. The largest transport/cell-distance ratio is
  0.7201430404231473655, and the largest refined-tail ratio is
  6.9552399397749245e-5. This proves the collar but, by direction, is a
  reverse-endpoint calibration rather than the old-edge induction hypothesis

jensen_window_pf_newman_q207_q208_forward_adiabatic_successor_certificate:
  all 414 stored Q207 old-bottom cells transform exactly to the regular proxy.
  Their 303 half-unit and 222 quarter-unit transport panels all pass, with
  maximum ratio 0.7630253717885428289. The compact core and derivative-positive
  new right strip then realize one genuine forward Q207-to-Q208 successor.
  No Q209 or all-j theorem is inferred

jensen_window_pf_newman_single_carrier_adiabatic_benchmark:
  an exact hypothetical carrier with a phase-speed floor and uniform
  logarithmic jets has heat-transport condition number O(ell^2), while
  delta_j*log((j+38)/(4pi))^2 tends to zero. A four-carrier exact double-zero
  model proves that positive amplitudes and ordered component speeds do not
  supply the required floor for a sum. The corrected Xi crossing small-ball
  estimate remains open

jensen_window_pf_newman_polymath15_crossing_slope_gap_reduction:
  the scaled corrected-component first jet is a rank-at-most-two Gram form
  that collapses to rank one on the crossing hyperplane. The exact crossing
  observable is a positive/negative real-mass weighted slope gap plus the
  zero-real-component term. A stronger four-carrier model, with no singular
  component slopes, makes that gap vanish and forces the remaining theorem
  to use the actual zeta phases, amplitudes, and endpoint relation

jensen_window_pf_newman_multiplicity_compatible_adaptive_jet_benchmark:
  the existing arbitrary-multiplicity Hermite split becomes quantitatively
  adiabatic in the jet (H,sqrt(2t)H_x): each fixed multiplicity has relative
  condition number C_m/(2t), and a cofinal successor costs only
  (C_m/2)log(1+1/j). This permits endpoint clearance to vanish without a
  positive-time contact. The Xi relative estimate and entry-strip theorem
  remain open

jensen_window_pf_newman_time_dependent_scaled_successor_lemma:
  every continuous positive derivative scale preserves first-jet contacts
  and degree. The full scaled heat derivative supports both the original
  additive cell bound and a multiplicative relative transport theorem. The
  latter carries inherited clearance through every positive collar without
  an endpoint floor and permits a rigorous interface between separate
  arithmetic-entry and multiplicity-descendant scales

jensen_window_pf_newman_ray_aligned_parabolic_frequency_reduction:
  the cofinal schedule t_j=25/(100+j) and
  log(R_j/(4pi))=101+j puts every genuinely new strip exactly in the proved
  tL>=25 dominant-saddle region. The smooth reciprocal-square scale enters
  in the 1/L coordinate and tends to sqrt(2t) on fixed descendants, with an
  exact alpha^2+beta^2=1 chart partition. Equivalently, the one remaining
  target is a relative theorem below
  tau(x)=min(1/4,25/log(x/(4pi))). The separate new-strip arithmetic
  hypothesis is removed. Bare existence of the relative bound would be
  equivalent to noncontact, so the open content is an a priori Xi majorant
  that never divides by the unknown first-jet norm

jensen_window_pf_newman_q207_q208_parabolic_frequency_relative_diagnostic:
  all 525 stored first-successor panels give rigorous finite relative costs
  in four coordinates. The blend agrees with the parabolic cost within a
  factor 1.007053, but both are much worse than the existing unit coordinate
  on this finite collar. The asymptotic blend is therefore not promoted
  backward into the already certified finite transition

jensen_window_pf_newman_parabolic_frequency_contact_normal_hierarchy_gate:
  every contact lies in a sharp corrected-main C1 remainder box, so the live
  target needs derivative separation only on its thin value band. A relative
  heat-jet proof exposes normalized second and third jets and receives no help
  from the scale-time term at collision. The exact chart split puts q<=1 only
  in t<=1/(2L^2), and two endpoint deficits control the optional smooth q>=1
  surrogate. The frequency-layer box reduces further to one robust weighted
  slope gap on the corrected value band, without a strip-wide cone. A
  universal correlation signature and quadratic Newman-flow countermodel
  reject generic hierarchy log-convexity. The Xi crossing-band theorem and
  its three domain pieces remain open

jensen_window_pf_newman_polymath15_critical_RS_C1_endpoint_peeling_contract:
  the first omitted Riemann-Siegel coefficient is
  C_1(p,sigma)=F'''(p)/(12pi^2)+i(2sigma-1)F'(p)/(4pi). On the critical line
  only the F''' term remains. Retaining it gives an explicit signed endpoint
  peel, removes the a^-1 remainder, and exposes an a^-2 endpoint target.
  Critical C_1 is odd, so its signed value matches at adjacent cutoffs; the
  core derivative jump is already O(T^-1)

jensen_window_pf_newman_polymath15_critical_dirichlet_first_correction_gate:
  the exact shifted finite-block ratio is a centered scaled-gamma expectation.
  Its first signed correction is
  d_(t,n)=1/(6s)+alpha'(s)(t/4+t^2 alpha_n^2/8), and the remaining logarithmic
  error splits exactly into M_0 Taylor, scaled-gamma, and rational-shift
  remainders

jensen_window_pf_newman_polymath15_critical_dirichlet_second_order_remainder_certificate:
  t|alpha_n|<=27 on the critical collar. A T^(1/3) central/tail split gives the
  explicit per-term bound 400000/T^2. Summation and Cauchy yield fixed-cell
  finite-Dirichlet bounds 8000000exp(-7L/4) in value and
  8000000Lexp(-7L/4) in the first derivative. The heat-integrated endpoint
  a^-2 estimate and adjacent-cutoff signed lift remain open

jensen_window_pf_newman_polymath15_critical_first_order_global_remainder_certificate:
  the heat-integrated C_0+C_1/a endpoint, variable-order negative-sigma tail,
  holomorphic lift, finite Dirichlet residual, and adjacent-cutoff comparison
  give |r_[1]|<100000exp(-5L/4) and
  |r_[1],x|<200000Lexp(-5L/4) on L>=50 and 0<tL<=25. The strict Xi
  contact inequality remains open

jensen_window_pf_newman_polymath15_critical_first_order_signed_contact_reduction:
  the retained finite d_n and endpoint C_1/a terms form one complex main
  E_[1]. A true contact lies in
  |Re E_[1]|<50000exp(-5L/4),
  |Re E_[1],x|<100000Lexp(-5L/4). The saddle-centered logarithmic moment
  gives the exact frequency-layer observable. The q<1 layer remains a
  separate multiplicity-compatible problem

jensen_window_pf_newman_polymath15_first_order_cofinal_boundary_reduction:
  the first-order squared boundary error is
  50000000000exp(-5L/2), a factor (3125/2)exp(-L) relative to the old
  budget. Ray-aligned right edges are already closed, bottom edges split
  exactly at q=1, and the same margin automatically joins adjacent cutoff
  lifts because their vector-square difference is only 13/500 of the budget.
  The two bottom margins, finite shoulders, and one-sided proxy phase bound
  remain open

jensen_window_pf_newman_polymath15_first_order_oriented_successor_winding_reduction:
  the standard (x,t) orientation gives positive contact degree. Exact
  oriented-chain cancellation confines every successor increment to the
  descendant slab, where it is a nonnegative integer. The scaled Pruefer form
  and positive-imaginary-ray intersections give signed bottom and connector
  crossing formulas. After first-order boundary transfer, any rigorous
  winding upper bound below one forces the increment to vanish; exact zero
  winding is stronger than necessary. The q>=1 and q<1 margins, finite phase
  shoulders, and one-sided all-stage bound remain open

jensen_window_pf_newman_polymath15_first_order_wronskian_crossing_reduction:
  the improved first-order remainder converts a hypothetical contact into an
  exp(-5L/4) Wronskian small-ball. At a nonzero complex-main crossing, the
  vector margin is exactly a 50000sqrt(5) normalized Wronskian margin, and the
  sign of W_[1]Im(E_[1]) identifies upward crossings. Zeros of E_[1] itself
  form a second explicit crossing class and cannot be divided away. The
  q>=1 signed Wronskian budget, q<1 chart, and finite phase cells remain open

jensen_window_pf_newman_polymath15_first_order_centered_real_residual_reduction:
  saddle centering gives S_a=Re(R_a)-v_a Im(E_[1])=U-u_a X, so S_a
  classifies every horizontal crossing, including E_[1]=0. Exact
  critical-line cancellation gives |u_a|<exp(-L), |v_a|<3/x^2, and the
  existing mass and endpoint bounds reduce U on the contact band to an
  explicit centered-moment-plus-endpoint scalar A_a with error below
  10^-6 exp(-5L/4). The q>=1 lower bound and A_a-positive crossing budget
  remain open; q<1 and finite phase cells remain separate

jensen_window_pf_newman_polymath15_first_order_centered_complex_zero_scout:
  a selected N=6 solve at L=3.7216 has |E_[1]| below 10^-80 while the
  real crossing slope is 0.5561 and the two-variable complex-zero Jacobian
  is 0.3954. The Wronskian vanishes, but S_a and A_a retain the simple
  upward crossing. Seventy- and ninety-five-digit reruns agree beyond
  fifty places. This is a moderate-height diagnostic, not an interval
  certificate or evidence for the L>=50 scalar theorem

jensen_window_pf_newman_polymath15_first_order_centered_adjacent_saddle_recurrence:
  differentiating the exact C_0 recurrence through order four gives the
  corrected adjacent endpoint and its x derivative. The adjacent A_a jump
  is one exact entering-moment/endpoint expression. In the reflected sharp
  ratio, the complete a^-1 coefficient vanishes before absolute values;
  selected L>=50 rows show exp(-3L/4) pieces cancelling at the
  exp(-7L/4) scalar scale

jensen_window_pf_newman_polymath15_first_order_centered_adjacent_chart_stability_certificate:
  on L>=50 and 0<=tL<=25 the corrected adjacent ratio satisfies
  |rho-1|<11/a^2 and |rho_x|<3/a^3. The centered endpoint rate and imported
  nuisance budgets then give
  |A_(a,N+1)-A_(a,N)|<5000exp(-7L/4), below 10^-7 of the contact scale at
  L=50. This closes adjacent endpoint bookkeeping, not the chart-invariant
  bulk scalar lower bound or signed crossing theorem

jensen_window_pf_newman_polymath15_first_order_centered_bulk_pair_transfer_gate:
  the chart-invariant pair has exact reflected-sharp and finite-Abel bulk
  coordinates. The full corrected ratio lock makes adjacent amplitudes
  decrease as rho_n<exp(-2h_n/5), and every interior locked pair has
  det(B_n)>=h_n^2/16 and a first-jet gap
  h_n^2|f_n|/(48L). Exact component and isolated locked-pair countermodels
  show why this local theorem does not control the aggregate, while an exact
  common-phase turning guard rejects relative ratio phases alone. The next
  input must use the ratios between successive pairs together with the
  normalizer/endpoint phase anchor, or an equivalent anchored Xi-specific
  signed turning/winding theorem

jensen_window_pf_newman_polymath15_first_order_centered_absolute_phase_anchor_reduction:
  the nonzero first coefficient gives a canonical branch-free phase anchor
  and relative coefficient chain with q_1=1. Exact endpoint cancellation
  reduces its direction to T_0+i and puts the finite sum and endpoint in two
  relative shapes Z_0 and Z_A. Their real anchored projections are exactly
  X and A_a. The natural aggregate determinant collapses to the centered
  Wronskian plus the already isolated frame and d_(n,x) terms, loses the
  E_[1]=0 branch, and has no certified complex adjacent-chart bound. The
  live theorem is therefore still the direct anchored real-projection lower
  bound and one-sided crossing count, with the complex-main-zero branch kept

jensen_window_pf_newman_polymath15_first_order_centered_direct_projection_regime_reduction:
  the anchored chain now has exact positive saddle amplitudes and unit
  carriers, explicit endpoint projections, and a derivative identity whose
  complete frame error is exponentially small. With tau_L=L exp(-L), the
  q>=1 theorem splits into an ordinary logarithmic phase/radial-current
  inequality on |Z_0|>=tau_L and directional real transversality on
  |Z_0|<tau_L. At L=50 the target-to-tau ratio is below 0.746, so the
  ordinary input is genuinely order one. An exact backward-heat solution
  shows that simple complex zeros and a nonzero two-variable Jacobian do not
  imply the exceptional real-direction slope. Both Xi-specific estimates
  and their oriented crossing count remain open

jensen_window_pf_newman_polymath15_first_order_centered_carrier_kernel_abel_prefix_reduction:
  pi is traced from the completed-zeta gamma factor through the
  Riemann-Siegel saddle, cutoff-cell width, and complex-exponential period;
  the prefix polygon defines none of it. The endpoint-complete radial and
  tangential currents now have exact
  pairwise and O(N) finite-Abel prefix forms. One continuous scalar
  mathcal_C_N equals the physical centered slope at every real crossing,
  including W_0=0; the Wronskian loses that class only by multiplying it
  by -Im(W_0). Across the L>=50 contact band, the remaining terminal real
  correction is below 1.34e-35 of the target. Exact positive-c and
  zero-fiber carrier models reject generic current and directional-slope
  floors. The actual Xi relative carriers rotate at speed greater than
  1/3 and complete a full relative turn in every cutoff cell, so the live
  theorem is an aggregate anchored prefix lower bound and signed crossing
  count, not a fixed termwise sector

jensen_window_pf_newman_polymath15_first_order_centered_normalized_prefix_phase_flux_reduction:
  removing the first-coefficient phase makes every prefix exactly
  G_k=sum_(n<=k)q_n with G_1=1 while preserving the external projector and
  endpoint. On q>=1 every nonfirst coefficient moves inward, all angular
  currents are ordered, and every adjacent pair advances by more than 4pi
  in a complete cutoff cell. An exact two-carrier model nevertheless has
  six simple zeros, three upward crossings, zero net signed count, and
  first-jet winding -3 in one phase cell. Restoring the required n=4
  carrier produces a completed dyadic cubic of discriminant -1696 with
  only one upward crossing. This revives complete multiplicative blocks as
  a structured arithmetic route, but does not prove their cross-block or
  endpoint composition. The viable integer target is
  therefore the O(N) endpoint-complete first-jet argument flux after full
  successor boundary composition, not an individual-carrier turn count

jensen_window_pf_newman_polymath15_first_order_centered_prime_power_heat_block_composition_guard:
  the absolute carrier rate is
  nu_n=v_a+b*u_n+Im(d_(n,x)/(1+d_n)). Every nonterminal n<=N-1 carrier
  rotates strictly in the same direction; only n=N can be near stationary,
  inside a fixed-t collar narrower than 50701/x. Complete p-power amplitudes
  decrease with ratio below p^(-12/25). After the moving p-free external
  phase is factored out, the correction-free normalized block is zero-free
  in the unit disk, and its one-step shifted pure-power block has winding
  one. The heat quadratic is an exact Gaussian mixture, but its terminal
  tilted expanding-ratio sector is not uniformly negligible. Most
  importantly, an exact length-three/length-two dyadic model has two
  origin-free shifted blocks of winding one whose sum has winding three.
  The live theorem must therefore keep the joined p-free prefixes, terminal
  recurrence, endpoint, and first-jet flux together; no actual Xi
  cross-block lower bound or successor theorem is claimed

jensen_window_pf_newman_polymath15_first_order_centered_joined_dyadic_odd_prefix_first_jet_reduction:
  unique dyadic valuation and an exact deterministic heat shift reorganize
  the complete corrected prefix into nested odd-prefix layers. Its value and
  first jet close on H_0,H_1,H_2,D_0,D_1 plus endpoint data, and one cutoff
  change enters exactly one valuation layer together with the recurrent
  endpoint. A polynomial-Mellin phase lift separates the dyadic and odd
  phases. In the stipulated correction-free synchronized-phase hypothetical,
  the joined z-polynomial has geometrically decreasing positive coefficients,
  so Enestrom-Kakeya proves unit-disk zero-freeness and an explicit margin.
  The physical odd-phase/correction transport, joined Xi lower bound, and
  strict successor flux theorem remain open

jensen_window_pf_newman_polymath15_first_order_centered_phase_cylinder_jacobian_transport_guard:
  the joined phase lift has an exact oriented z/xi Jacobian with an explicit
  unordered-pair kernel. Along the physical phase linkage its tangent is
  iH_1, and the cylinder degree theorem transports unit-circle winding by
  the signed Jacobian zeros. Odd-phase, d_n-correction, and endpoint-strength
  restorations therefore form three exact transport cylinders. A genuine
  n=2,3 kernel changes sign, a zero-free dyadic block can have nonmonotone
  linked argument, and endpoint addition can create a boundary zero. The
  Jacobian survives as a local degree current, while the Xi signed crossing
  budget remains open

jensen_window_pf_newman_polymath15_first_order_centered_endpoint_schur_cohn_first_jet_guard:
  the endpoint-augmented joined polynomial has an exact outside-disk
  Schur-Cohn recursion and inverse. Strict recursive reflection pivots are
  equivalent to closed-unit-disk zero exclusion and give the conditional
  boundary margin |b_0| product_d(1-|alpha_d|). Fixed-z x tangents, the
  centered Z_A companion, and their mixed tangent propagate through the same
  recursion while retaining H_2,D_0,D_1 and endpoint derivatives. The first
  physical inequality is |r_0+B_(0,0)|>|B_(K,0)|. A strict first pivot,
  decreasing coefficient moduli with free phases, and prefix zero-freeness
  before endpoint addition are each insufficient; every actual Xi pivot
  remains open

jensen_window_pf_newman_polymath15_first_order_centered_endpoint_first_pivot_odd_small_ball_guard:
  the terminal dyadic layer is the singleton m=1, and the shared
  normalization denominator cancels exactly, while the m=1 correction
  remains inside O_N. The first pivot is therefore the endpoint-shifted odd
  small-ball condition
  O_N notin Dbar(-R_N,rho_K), with one explicit endpoint correlation as the
  missing term. In a correction-free N=5 model, Kronecker density and the
  triangle formed by lengths 1, 1/sqrt(3), and 1/sqrt(5) prove that even the
  genuinely linked logarithmic phases can make the pivot negative; a whole
  neighborhood of small fixed endpoints has the same obstruction. This
  rejects amplitude, linked-phase, and endpoint-smallness promotion, not an
  actual Xi pivot. The live alternatives are an Xi-specific odd-disk
  avoidance theorem or the weaker linked-point/signed-degree route

jensen_window_pf_newman_polymath15_first_order_centered_endpoint_odd_fibre_correlation_feasibility_guard:
  after restoring the physical unit phase, write the actual generally
  complex endpoint as g_0=A_N+iB_N.  The missing first-pivot correlation is
  A_NX_odd+B_NY_odd.  Rotating into that actual endpoint direction gives
  |S_odd+g_0|^2=(P_odd+D_N)^2+Q_odd^2 for
  D_N=sqrt(A_N^2+B_N^2)>0, with D_N=0 retained without division. A four-term dyadic
  finite difference on q_3,q_6,q_12,q_24 fixes H_0,H_1,H_2,D_0,D_1, the
  endpoint, and the terminal carrier in the unrestricted correction-free
  joined class, while changing the odd fibre arbitrarily. Separately,
  F_lambda(w)=1+lambda(w-1)^3 fixes its linked value and first two derivatives
  but changes the first Schur-pivot sign. These are exact current-input
  insufficiency guards, not Xi counterexamples. Full Schur disk stability is
  retained as a falsification coordinate; the endpoint-complete signed
  linked-point theorem is again the primary open route

jensen_window_pf_newman_polymath15_first_order_centered_abel_scalar_shear_flux_reduction:
  on each fixed-N q>=1 chart the physical slope is the real shear
  mathsf_A=mathcal_C_N+c u_N mathsf_X. Conditional on the proposed
  contact-band gap, a determinant-one homotopy removes that shear and every
  positive jet scale, reducing the signed coordinate to the division-free
  O(N) proxy Psi_ell=mathsf_X+i mathcal_C_N/ell. The established nuisance
  budgets make the true horizontal derivative differ from mathcal_C_N by
  less than 2.03e-14 A_L on the band, so the gap fixes crossing orientation
  even at W_0=0. The exact backward-heat family
  exp(m^2t)sin(mx+pi/4) nevertheless has an arbitrarily large local gap,
  m upward crossings, and winding -m per period. Therefore the Xi
  Abel-scalar gap and the completely composed one-turn phase budget are
  independent open theorems; neither full Schur stability nor generic
  transversality supplies the latter

jensen_window_pf_newman_polymath15_first_order_centered_dominant_ray_connector_phase_cap:
  on every ray-aligned successor connector, the exact normalized Xi first
  jet remains within a uniform 5/9 relative cone around the one-saddle
  ellipse. The exact heat-normalizer identity moves the saddle phase by less
  than pi/64, and the certified ellipse eccentricity varies by at most 1/50.
  Combining these bounds confines the complete connector to angular range
  9719/6720<3/2<pi/2, uniformly through every cutoff transition. After the
  two proxy-switch tracks are assigned once to the complementary path, the
  one-turn theorem reduces from an unknown 2pi full-boundary budget to the
  signed horizontal-composite target H_j<3pi/2. The pointwise Abel gap,
  q<1 arcs, finite shoulders, endpoint tracks, and that horizontal
  inequality remain open

jensen_window_pf_newman_polymath15_first_order_centered_uniform_q_ge_1_degree_excision_reduction:
  the moving interface L_*(t)=max(50,(2t)^(-1/2)) cuts each successor into
  the exact q>=1, L>=50 outer collar and the q<1/bounded-L inner region.
  The outer collar is covered by the open Abel regime tL<=25 and the already
  proved dominant-saddle noncontact regime tL>=25. Conditional on proving
  the Abel gap uniformly on the whole closed first regime, the exact outer
  jet is nonzero and its degree vanishes, so the successor integer localizes
  to the inner region. A paired many-turn guard validates horizontal
  cancellation under uniform interior noncontact, while a quadratic heat
  contact shows that boundary-only nonvanishing cannot justify excision.
  The uniform Abel gap and localized inner successor theorem remain open

jensen_window_pf_newman_polymath15_first_order_centered_oscillatory_spliced_outer_collar_reduction:
  for every fixed positive epsilon below 25-c_*, raising the outer cutoff
  to B_epsilon=max(50,L_epsilon) lets the proved oscillatory-zeta theorem
  replace the open Abel obligation above c_*+epsilon.  Together with the
  dominant theorem above 25, this gives a closed three-regime collar cover.
  Conditional outer degree excision therefore needs the Abel gap only on
  q>=1 and tL<=c_*+epsilon, while q<1 and the fixed bounded-L shoulder stay
  in the inner theorem.  The current July-2026 ANTEDB frontier leaves
  c_*=4911678521/1933561194 unchanged.  L_epsilon remains existential, so
  no finite-height completion or RH conclusion is promoted

jensen_window_pf_newman_polymath15_first_order_centered_reciprocal_saddle_self_duality_gate:
  the Poisson saddle u_nu=a_omega^2/nu conjugates the logarithmic phase,
  while the leading heat amplitude is exactly reciprocal-self-dual up to
  the scalar defect Delta=2sigma-1-tlog(a_omega).  For the physical Xi
  shift that defect is O(t/x), with a uniform exp(+-13/(2x)) modulus
  enclosure.  The active pointwise block r_*=125662/155153 maps to the
  outside-cutoff tail radius 184644/155153, exposing an endpoint-complete
  primal/tail composition as a new cancellation target.  Equal reciprocal
  moduli can reinforce as well as cancel, so the global Xi phase, corrected
  endpoint, adjacent recurrence, and first derivative remain mandatory.
  The discrete paired theorem and any improvement of c_* remain open

jensen_window_pf_newman_polymath15_first_order_centered_reciprocal_normalizer_phase_reinforcement_guard:
  the exact Polymath-15 normalizer locks the reciprocal stationary phase
  to the conjugate physical carrier.  The heat-shifted angular error is
  below 1/x, and the combined amplitude/phase partner differs from the
  conjugate carrier by less than 8/x.  A real carrier is therefore
  asymptotically doubled, while a small nearly imaginary projection is
  only the original oscillation.  Raw reciprocal pairing cannot supply
  the N^(-d_2-eta) gain needed at c=2.  The surviving reciprocal target
  is the signed weighted-minus-unweighted tail symbol, with endpoint and
  first derivative retained; no improved c_* or RH claim is promoted

jensen_window_pf_newman_polymath15_first_order_centered_signed_handoff_reciprocal_symbol_guard:
  the signed zeta-handoff difference has exact continuous reciprocal
  symbol Sigma_N=1-exp(-t log(u_nu)^2/4) after division by the dominant
  weighted saddle.  At the active radius r_*=125662/155153 this is
  1-N^(-c r_*^2/8+o(1)), so it tends to one rather than zero.  At c=2
  the exact weighted-minus-ordinary exponent gap is
  r_*^2/4=3947734561/24072453409.  The fixed-chart first derivative of
  the signed factor is power-small, the active dual block is separated
  from the endpoint seam by N^(1-r_*), and the checked correction and
  endpoint scales cannot reverse the leading symbol uniformly.  Signed
  primal/tail pairwise cancellation is therefore retired.  A reciprocal
  improvement would now require genuinely new arithmetic cancellation
  inside the dual exponential sum; the direct Xi Abel-phase theorem
  remains the primary route

jensen_window_pf_newman_polymath15_first_order_centered_contact_signed_transport_reduction:
  after terminal centering, the endpoint and every carrier contribute one
  real value mass c_j and one real slope contribution d_j, with
  sum c_j=mathsf_X and sum d_j=mathcal_C_N.  On the contact surface the
  positive and negative real masses balance and the Abel scalar is exactly
  M(h_+-h_-)+D_perp.  Ordering nonzero-real components by effective slope
  gives a second exact form with nonnegative slope gaps multiplying signed
  cumulative real masses.  This isolates a concrete Xi-specific transport
  inequality.  For H_a=H_R+iH_I and J_a=J_R+iJ_I, the corrected endpoint
  atom is c_0=kappa(T_0H_R-H_I) and
  d_0=kappa[T_0J_R-J_I-c u_N(T_0H_R-H_I)].  Its effective slope is divided
  only when T_0H_R-H_I!=0; the full zero-real-projection fibre is retained
  without division.  The associated two-feature
  Gram is PSD but drops to rank at
  most one on contact, and an endpoint-shaped synthetic null model confirms
  that symmetry, decreasing amplitudes, and ordered distances do not create
  coercivity.  The complete H_2,D_0,D_1 first jet and adjacent cutoff law
  remain attached.  No Xi cumulative-mass bound or Abel-scalar gap is proved

jensen_window_pf_newman_polymath15_first_order_centered_interior_projective_current_gate:
  for every interior carrier n<=N-1, the division-free current
  J_n=c_n d_(n,x)-d_n c_(n,x) cancels the radial amplitude current and is
  bounded above by -(5/64)u_n^2|z_n|^2.  Thus every interior projective
  atom turns strictly clockwise, including through c_n=0 and at q=1.  The
  terminal carrier is exceptional: at u_N=0 its current is
  -b u_x X_NY_N and has no source-level fixed sign.  It is therefore
  packaged with the complex recurrent endpoint into one exact edge block.
  The within-chart J_a^(der) and adjacent-cutoff J_a^(adj) are explicitly
  distinguished.  The Xi cumulative-mass estimate and edge sign remain open

jensen_window_pf_newman_polymath15_first_order_target_reconciliation_gate:
  the cutoff-uniform first-order main J_[1] and its exp(-5L/4) C1 remainder
  supersede the historical zeroth-order exp(-3L/4) handoff.  A contact lies
  in the exact rectangle |X|<50000exp(-5L/4),
  |U|<100000Lexp(-5L/4).  The radial energy threshold
  50000000000exp(-5L/2) excludes that rectangle but is strictly stronger
  than the box-optimal signed band.  The endpoint-complete, division-free
  Abel gap implies the box target through the certified terminal-shear and
  real-residual bounds, including W_0=0.  The remaining outer arithmetic
  obligation is restricted to one fixed-epsilon low-c, q>=1 wedge.  Exact
  component cancellation forbids promotion from diagonal positivity, and
  no Xi Abel gap or pointwise contact exclusion is claimed

jensen_window_pf_newman_polymath15_first_order_rectangular_boundary_degree_reduction:
  the matched first-order value and derivative errors define a rectangle,
  so their weighted sup gauge gives a componentwise boundary-Rouche
  homotopy strictly weaker than the Euclidean radial margin.  Positive
  diagonal error whitening, the normalizer shear, and the adjacent-cutoff
  sub-box all preserve winding.  Since every real heat contact has positive
  local index, a transferred successor winding below one excludes contacts.
  On bottom j the remaining fixed-epsilon low-c q>=1 Abel domain is the
  explicit one-dimensional interval
  max(B_epsilon,sqrt((100+j)/50))<=L
  <=min(L_j,(c_*+epsilon)(100+j)/25).  One may prove the Abel gap uniformly
  on the full outer wedge and excise it, or prove it only on the joined
  boundary and separately establish H_j<3*pi/2.  Neither open inequality is
  promoted

jensen_window_pf_newman_polymath15_first_order_centered_ray_bottom_logarithmic_flow_reduction:
  each nonempty hard ray is partitioned exactly into canonical fixed-cutoff
  cells, with cutoff equality owned by the new chart.  On each cell
  D_j=d/dL=x partial_x, so the normalized carriers, every prefix, the
  recurrent endpoint, and the moving absolute projector have exact
  endpoint-complete logarithmic derivatives.  For the Abel proxy
  Psi=mathsf_X+i mathcal_C_N/L, its phase numerator is the negative square
  -Lx mathcal_C_N^2 plus the endpoint-complete signed defect
  mathsf_X(LD_jmathcal_C_N-mathcal_C_N) and the already tiny orientation
  error.  The pointwise Abel gap would therefore preserve crossing
  orientation in L, including W_0=0, but it does not bound recrossings.
  Abel proxies join through the physical slope coordinate at cutoffs, and
  opposite rays must be compared only inside the fully joined H_j ledger.
  No Abel gap, aggregate defect sign, or horizontal phase bound is promoted

jensen_window_pf_newman_polymath15_first_order_centered_geometric_prime_power_phase_monotonicity_gate:
  a correction-free geometric prime-power block has an exact phase-speed
  formula and the lower margin
  mu_M(a)=1/(1+a)-Ma^M/(1-a^M).  This proves strict one-turn phase for
  every block over p>=5, for ternary blocks of length at least four, and
  for dyadic blocks of length at least eight.  Multiplicative completion is
  not enough by itself: a complete two-level dyadic or ternary block at
  external phase pi/2 has four simple real crossings, two upward crossings,
  and first-jet winding -2 even though the complex block winds once.
  Therefore short p=2,3 chains remain explicit and the positive long-block
  margins become C1 perturbation budgets for the heat and coefficient
  corrections.  P-free bases, singleton chains, cutoffs, and the recurrent
  endpoint still have to be rejoined before estimating the signed ray defect

jensen_window_pf_newman_polymath15_first_order_centered_prime_power_logarithmic_phase_flow_gate:
  the actual fixed-ray chain derivative separates into an angular moment,
  radial drift, and normalized coefficient-current term.  This gives one
  exact division-free current for the shifted internal block and separately
  exposes the moving p-free external phase.  The normalized heat exponent
  can be order one across a long chain.  A complete actual dyadic M=8
  cutoff-cell witness proves that the former absolute geometric-C1 target
  fails by more than 1/2 even though its heated phase speed stays positive.
  Direct one-sided heat-current positivity replaces that invalid norm.
  No uniform long-chain, joined p-free, endpoint, Abel-gap, or successor
  theorem is promoted

jensen_window_pf_newman_polymath15_first_order_centered_prime_power_heat_starlikeness_base_certificate:
  exact tensor-Bernstein arithmetic proves correction-free current margins
  above 1/300 for the dyadic M=8 base block and above 1/25 for the ternary
  M=4 base block.  A direct two-level formula proves a margin above 1/20
  for every p>=5, with an independent exact p=5 Bernstein cross-check.
  Restoring the saddle offset, d_n coefficients, radial drift, and epsilon_n
  current costs less than 1/1000.  The resulting actual fixed-ray margins
  are 7/3000, 39/1000, and 49/1000.  Length propagation, short and singleton
  chains, p-free external phases, the recurrent endpoint, and the joined
  ray defect remain open

jensen_window_pf_newman_polymath15_first_order_centered_prime_power_heat_starlikeness_length_propagation_gate:
  the exact append-one recurrence isolates the contracted old current and
  the complete new cross current.  A 32-leaf exact tensor-Bernstein
  certificate proves J_0>1/20 for the dyadic M=9 box.  For every p>=5,
  the heat-coefficient ratios make the Fourier current sequence decreasing
  and convex at every length.  Its exact Fejer-kernel decomposition gives
  J_0>1/20, and length-uniform coefficient sums retain
  J_ray>49/1000 after the actual saddle, d_n, radial, and epsilon_n terms.
  The remaining complete-chain length families are p=2,M>=10 and
  p=3,M>=5; short, singleton, p-free, endpoint, and joined ray obligations
  remain open

jensen_window_pf_newman_polymath15_first_order_centered_prime_power_heat_starlikeness_ternary_length_gate:
  for every complete ternary chain p=3,M>=5, exact Fourier grouping leaves
  only one possibly negative second difference.  Five exact rational
  small-d cluster certificates, a uniform positive D_0 reserve, twelve
  exact finite defect-absorption certificates, and an analytic M>=17 tail
  prove J_0>17/1000.  Length-uniform geometric coefficient sums retain
  J_ray>2/125 after the actual saddle, d_n, radial, and epsilon_n terms.
  The remaining all-length complete-chain family is p=2,M>=10; short,
  singleton, p-free, endpoint, Abel-gap, and joined winding obligations
  remain open

jensen_window_pf_newman_polymath15_first_order_centered_prime_power_heat_starlikeness_dyadic_length_gate:
  for every complete dyadic chain p=2,M>=10, exact Fourier grouping
  localizes possible defects to D_0 and D_(M-2).  A completed-square
  terminal-pair identity, five finite seven-kernel certificates, one
  length-free four-kernel certificate, exact finite terminal guards, and
  an analytic decreasing terminal tail prove J_0>3/1000.  Length-uniform
  geometric coefficient sums retain J_ray>1/500 after the actual saddle,
  d_n, radial, and epsilon_n terms.  Together with the M=8 and M=9 gates,
  complete dyadic chains are covered for M>=8.  Short, singleton, p-free,
  endpoint, cutoff, adjacent-chart, Abel-gap, and joined winding obligations
  remain open

jensen_window_pf_newman_polymath15_first_order_centered_prime_power_heat_starlikeness_short_family_counter_gate:
  at the common linked interior point q=s=9999/10000, eight exact
  Q(sqrt(p)) witnesses prove J_0<-1/100 for every short dyadic length
  p=2,M=2..7 and both short ternary lengths p=3,M=2,3.  The witnesses
  use simple rational angular coordinates, exact component or square-sign
  guards, and no floating-point sign decision.  Thus M=8 and M=4 are
  sharp thresholds for the present blockwise small-prime theorem.  The
  desired universal short-block positivity statement is false; short
  chains must be retained in the joined p-free, singleton, endpoint,
  cutoff, Abel-gap, and winding problem

jensen_window_pf_newman_polymath15_first_order_centered_joined_phase_current_polarization_gate:
  the endpoint-complete phase current polarizes exactly into diagonal
  component currents and all pair cross currents.  At a real-projection
  contact it is -mathsf_Y(x mathcal_C_N+E_j), while each physical
  prime-power theorem enters only through one explicit diagonal term.
  The exact family exp(i theta)+a exp(3i theta) has positive self currents
  but joined current (1-a)(1-3a) at theta=pi/2.  At a=1/3 its real
  projection is (4/3)cos(theta)^3 and has a cubic contact; at a=1/2 the
  joined current is -1/4.  Thus the long-chain theorems are valid diagonal
  input, but an Xi-specific cross-current aggregate or the direct
  division-free Abel gap remains necessary

jensen_window_pf_newman_polymath15_first_order_centered_joined_pair_kernel_pfree_rejoin_gate:
  every physical carrier pair has an exact division-free current kernel
  whose ordered logarithmic rates multiply unrestricted complex
  correlations.  Summing all diagonal and pair rows collapses exactly to
  the H_0,H_1,D_0 moment current.  The heat-shifted p-free decomposition
  telescopes linearly, but its quadratic phase current rejoins with every
  cross term and reconstructs that original full-prefix moment.  An exact
  two-layer jet has diagonal reserve 4/3, cross current -4/3, and joined
  current zero.  Thus the prime symmetry remains a valid coordinate
  organization, while an Xi-specific correlation estimate or the direct
  endpoint-complete division-free Abel gap remains necessary

jensen_window_pf_newman_polymath15_first_order_centered_mangoldt_abel_contact_gate:
  the retained first Dirichlet correction is exactly quadratic in log(n),
  so the endpoint-complete Abel bulk closes on four correction-free
  logarithmic moments H_0 through H_3.  The three nonconstant moments are
  exact transpose-symmetric von Mangoldt bilinear forms on dm<=N, with the
  heat quadratic isolated as exp[(t/2)log(d)log(m)].  This gives a genuine
  arithmetic symmetry and the exact contact coordinate
  mathcal_C_N=Re[g+s_*'U]-c u_N mathsf_X, including W_0=0.  A prime edge
  sqrt(N)<p<=N gives the strict principal-minor determinant
  -log(p)^(2j)/4 for each j=1,2,3.  The symmetry is therefore useful for a
  cancellation-preserving Type-I/II decomposition, but it is not a Gram
  positivity proof and supplies no Abel gap by itself

jensen_window_pf_newman_polymath15_first_order_centered_mangoldt_contact_centering_gate:
  the corrected nonconstant moments combine exactly into the one physical
  logarithmic derivative Z_Lambda=sum_(n<=N)log(n)z_n, which is also the
  symmetric Mangoldt hyperbola sum.  Contact centering gives
  mathfrak B_N=E_N-s_*'Z_Lambda+s_*'log(a)W_0=g+s_*'U and therefore
  retains W_0=0 without division.  The identity
  Z_Lambda=log(N)S-sum_k log((k+1)/k)F_k proves that the Mangoldt and Abel
  prefix expressions are dual coordinates of the same jet, not
  independent evidence.  A floor(sqrt(N)) split gives one complete square
  and one exact hyperbolic wing.  The signed centered-jet estimate and Abel
  gap remain open

jensen_window_pf_newman_polymath15_first_order_centered_mangoldt_adjacent_cutoff_transport_gate:
  at fixed (t,x), the anchored endpoint-value jump cancels exactly between
  Delta E_N and log(a)Delta W_0.  The remaining centered jump is
  Delta g+s_*'log(a/n)z_n, exactly the established physical adjacent jet,
  so its certified real-projection bound transfers without a new complex
  observable.  At nonsquare cutoffs only the hyperbolic wing changes.  At
  n=r^2, the new complete-square row moves out of the old wing and cancels
  internally; the diagonal plus new divisor boundary equals log(n)z_n.
  An independent formal prime-log audit passes 79 consecutive transitions,
  including all eight square transitions through n=81, with no residual
  floor term.  The signed centered-jet estimate and Abel gap remain open

jensen_window_pf_newman_polymath15_first_order_centered_mangoldt_endpoint_composed_vaughan_gate:
  because the physical moment is sum log(n)z_n, the Vaughan decomposition
  of Lambda retains its final convolution by 1.  The resulting exact
  coefficient identity log=C_0+C_I+C_II composes with the recurrent
  endpoint as
  mathfrak B_N=R_(0;U,V)-s_*'Z_I-s_*'Z_II, without dividing on W_0=0 or
  mathsf_X=0.  At the canonical cube-root cutoff, all historical
  perfect-cube transfers cancel coefficientwise, including the square/cube
  overlap at n=64.  Exact 1,5,25 and 1,11 route guards show that separate
  Type-I/II absolute values and endpoint-free estimates can both erase the
  required cancellation.  The missing input is an Xi-specific pointwise
  signed endpoint/Type-I/Type-II correlation theorem; no such lower bound
  or Abel gap is proved

jensen_window_pf_newman_polymath15_first_order_centered_complex_endpoint_source_normalization_gate:
  the published endpoint function C_0(p)=F(p), not the Vaughan coefficient
  C_0(n), is generally complex.  Exact parity gives C_0'''(0)=0 while
  Im C_0(0)<0, so H_a(0) is nonreal at every half-integer saddle.  The
  corrected endpoint has Re e=kappa(T_0H_R-H_I),
  Re g=kappa(T_0J_R-J_I), and
  gconj(e)=kappa^2(T_0^2+1)J_aconj(H_a).  The projective exceptional fibre
  is T_0H_R-H_I=0, the exact endpoint/carrier contact minors remain
  unsigned, and the endpoint-composed Vaughan square has rank one.  Seven
  historical real-H_a endpoint rows are quarantined; the unaffected
  branch-free, Abel, Mangoldt, cutoff, and Vaughan identities remain exact.
  No Xi phase law, Abel gap, winding cap, contact exclusion, or RH result is
  proved

jensen_window_pf_newman_polymath15_first_order_centered_endpoint_relative_phase_current_recurrence_gate:
  the corrected endpoint self-current is
  Im(e_xconj(e))=kappa^2[(T_0^2+1)Im(H_(a,x)conj(H_a))-|H_a|^2/2].
  The common frame phase cancels from the division-free current against
  each physical carrier.  On the physical q=1 cofinal family its terminal
  scaled limit is -1/4 at p=0 but
  [1-sqrt(4+2sqrt(2))/4]/4>0 at p=1, so no uniform terminal sign exists.
  Terminal tails telescope before projection into an older endpoint plus
  exact recurrence defects, while the phase current retains an uncontrolled
  mixed defect.  At contact the exact cumulative identity
  sum_n K_(0n)=c_0 mathcal C_N recovers the full scalar when c_0!=0 and
  becomes blind when c_0=0.  The terminal-composed real-edge sign,
  signed cumulative-minor estimate, Abel gap, winding cap, contact
  exclusion, and RH result remain open

jensen_window_pf_newman_polymath15_first_order_centered_real_edge_projective_current_gate:
  composing the corrected endpoint and terminal carrier before real
  projection gives the exact division-free current
  J_edge=C_edge Re(G_x)-Re(G)C_(edge,x)-alpha_x C_edge^2.  On the physical
  q=1 fixed-p cofinal family its second-jet symbol collapses to
  K_edge=[A A''-(A')^2]/(16pi^2), with A=-Re C_0(p-2).  A 192-bit Arb cover
  of 3072 rational intervals, using analytic sinc charts at both removable
  points, proves (A')^2-AA''>3/50 and therefore K_edge<-3/8000 throughout
  -1<=p<=1.  This is not yet a uniform finite-height sign; the explicit
  remainder, q>1 extension, adjacent-cutoff splice, Abel gap, winding cap,
  contact exclusion, and RH result remain open

jensen_window_pf_newman_polymath15_first_order_centered_q1_finite_height_real_edge_remainder_gate:
  on q=1 and L>=50, the terminal carrier has a cancellation-stable logarithm
  Z=-Q exp(W) with |W|<6/a and |W_x|<2/a^2.  Two factored complex Cauchy
  charts and a 4096-box Arb cover bound C_0 through order five.  The four
  normalized real-edge jets have defects bounded by 100, 250, 250, and 900,
  giving |a^2J_edge/(kappa T_0)^2-K_edge|<5500/a<1/10000000.  Combined with
  the strict cofinal margin, the retained first-order model current is below
  -3749/10000000 throughout -1<=p<=1.  This does not bound the omitted
  higher-order Xi error, prove a q>1 or Xi-level edge sign, splice adjacent
  cutoffs, or prove the Abel gap, winding cap, contact exclusion, or RH

jensen_window_pf_newman_polymath15_first_order_centered_critical_ray_finite_height_real_edge_gate:
  replacing t<=1/5000 by the complete critical bound t<=1/2 and recomputing
  27 normalized rational majorants gives |W|<4/a, |W_x|<1/a^2,
  |mu|<1/a^2, and |mu_x|<1/a^4.  The same four edge-jet defect constants
  100, 250, 250, and 900 remain valid, so the retained first-order model
  current stays below -3749/10000000 for L>=50, 0<tL<=25, and -1<=p<=1.
  This contains every q>=1 critical point.  The adjacent-cutoff signed
  splice, omitted Xi/source current remainder, cumulative Abel gap, winding
  cap, contact exclusion, and RH remain open

jensen_window_pf_newman_polymath15_first_order_centered_adjacent_cutoff_real_edge_projective_splice_gate:
  at an integer cutoff, the aligned old p=-1 and new p=1 edge rows have
  exact source values A_-=A_+=-sqrt(2+sqrt(2))/4,
  B_-=-3sqrt(2-sqrt(2))/16, and B_+=-sqrt(2-sqrt(2))/16.  Their leading
  oriented wedge is therefore -sqrt(2)/32.  The certified finite-height
  value and slope defects give
  Delta_cut=-sqrt(2)h/32+O(500h^2+50000h^3)<-h/25.  Both aligned real
  values remain below -1/4, so the affine projective join never vanishes
  and has strict clockwise current.  This closes every adjacent cutoff for
  the retained first-order q>=1 critical-ray edge model without a second
  adjacent x-jet.  The omitted Xi/source remainder, Xi-level edge sign,
  cumulative-minor estimate, Abel gap, winding cap, contact exclusion, and
  RH remain open

jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_full_support_ideal_cubic_symmetric_outer_pairing_gate:
  the global exterior/tie combination is exactly the common symmetric
  outer-mode complement.  Every finite common cutoff splits into one finite
  near sum with contiguous mean-one Dirichlet kernel and opposite remote
  pairs.  The remote pairs have an explicit absolutely convergent C^2
  bound.  The ideal Hermitian endpoint jump is nonzero and gives opposite
  harmonic terms in the separate positive and negative tails, so those
  half-series diverge and cannot receive independent budgets.  The cubic
  reflection moves u to a^2/u and preserves neither the physical interval
  nor the Fourier phase; conjugation also reverses alpha_P.  The signed
  carrier-plus-near estimate, physical remote derivative bound, transpose
  terminal completion, and every RH-level conclusion remain open

jensen_window_pf_fixed_root_shift_tangency_rigidity_gate:
  a fixed nonzero multiplicity-m root in every degree-d Jensen shift forces
  the repeated-root recurrence (1+rE)^(d-m+1)A_(n+m-1)=0.  In the corpus
  exponential-generating-function normalization this is exactly
  F(z)=exp(qz)Q(z), q=-1/r, with deg Q<=d-m; its converse also holds.  The
  old ordinary-series pole is therefore not an entireness contradiction.
  More usefully, p+1 consecutive recurrence defects form one Hankel
  matrix-vector product, so p+1 exact rows force its determinant to vanish.
  Equivalently, p consecutive hyperbolic upper-degree extensions propagate
  one multiplicity-m root through exactly the p+1 rows needed for a finite
  first-boundary collision obstruction.
  Four fixed-root quartic double tangencies would force H_(4,n+1)=0 and are
  excluded throughout -100<=lambda<=0 by the certified positive all-shift
  H_4 forward-invariance theorem.
  Exact Newton and singular-value bounds quantify near-rigidity for one
  reference root.  Heat-uniform Hankel margins, a compact nonzero root
  annulus, adjacent-root drift, and every RH-level conclusion remain open

jensen_window_pf_first_boundary_eventual_hankel_escape_gate:
  fix a finite degree d, shift n, and multiplicity m before setting the
  determinant order k=d-m+2.  If all degree-(d+1) shifts from n onward are
  hyperbolic, polar lifting and equal-degree transfer propagate the same
  nonzero root arbitrarily far in shift.  Choosing the fixed-order eventual
  threshold N_k and then moving the chain beyond it forces one H_k both zero
  and nonzero.  Thus no attained finite-degree, finite-shift collision can
  witness a global first boundary.  The degree-20 matrix audits all 190
  degree/multiplicity pairs, but the theorem is symbolic for every fixed
  finite degree.  No threshold uniform in k is obtained; unbounded degree,
  unbounded shift, non-attainment, all-degree hyperbolicity, and RH remain
  open

jensen_window_pf_newman_arbitrary_multiplicity_jensen_boundary_layer_gate:
  if an entire coefficient function has a finite negative multiplicity-m
  zero rho, its scaled degree-D Jensen polynomial on the parabolic heat layer
  converges to exp((4rho*tau-rho^2/2)d_eta^2)eta^m.  This is a real-rooted
  scaled Hermite polynomial on one side and has a nonreal pair on the other,
  so a degree sequence has collisions with D_j(lambda_j-lambda_*) tending to
  rho/8.  Under Lambda>0, positive-boundary attainment gives rho=-c^2 at the
  single shift n=0 and therefore t_j=Lambda-c^2/(8D_j)+o(D_j^-1).  The
  theorem removes the double-zero and shift-escape assumptions but supplies
  no estimate uniform in degree, no collision exclusion, and no RH-level
  conclusion

jensen_window_pf_newman_arbitrary_multiplicity_secondary_cubic_layer_gate:
  at multiplicity m>=3, the critical tau=rho/8 layer has the second zoom
  epsilon=D^(-1/6).  The normalized cluster tends to the universal polynomial
  exp(q*d_y^2+d_y^3/12)y^m.  Its real-rooted parameter set is exactly a lower
  ray ending at a finite q_m<0, which refines the collision parameter by
  (rho*q_m/4)D^(-4/3) and the root by rho*y_m D^(-5/3).  For m=3,
  q_3=-2^(-7/3) and y_3=2^(-2/3).  Neither the fixed shift nor the local Xi
  unit enters this layer, so it gives no Xi-specific sign, degree-uniform
  exclusion, first-jet boundary theorem, or RH-level conclusion

jensen_window_pf_newman_arbitrary_multiplicity_tertiary_universal_layer_gate:
  extending the fractional zoom by one graded step gives the complete
  relative-epsilon^2 correction
  Q_(m,n)=((r-y/2)d_y^2+(2n+1)d_y/4)P_m.  It is still universal: every
  independent fourth-derivative contribution cancels or recombines.  At a
  nondegenerate threshold contact, r_m=y_m/2 and s_m=(1-2n)/4.  For m=3,
  r_3=2^(-5/3), giving the explicit D^(-5/3) heat correction and a D^(-2)
  root-center term matching the double-zero layer.  The first linear local
  Xi-unit coefficient is delayed to relative epsilon^4 and equals
  rho[U'(rho)/U(rho)]P_(m+1).  The rest of that epsilon^4 layer, its sign,
  a degree-uniform remainder, collision exclusion, and every RH-level
  conclusion remain open

jensen_window_pf_newman_arbitrary_multiplicity_first_xi_jet_layer_gate:
  the complete next layer is
  P_m+epsilon^2 A P_m+epsilon^4{R P_m+ell P_(m+1)}+O(epsilon^6), with
  ell=rho U'(rho)/U(rho) and an explicit universal fifth-order operator R
  whose final coefficient is 1/80.  At the exact m=3 contact this gives
  v_3=-(2ell+3n+2)/8 and next root coefficient
  2^(-7/3)(2ell+n+4), hence the first source-dependent D^(-2) collision-time
  refinement.  The genus-zero product identifies ell as a signed
  neighboring-zero field.  Exact all-negative-root, positive-coefficient
  models realize ell=-1/3 and ell=4/3, so local Laguerre--Polya geometry
  does not fix its sign.  An actual Xi field bound, degree-uniform remainder,
  collision exclusion, and every RH-level conclusion remain open

jensen_window_pf_newman_regular_field_first_loss_nonclosure_gate:
  for every real ell, every c>0, and every multiplicity m>=2, the exact
  product (s+c^2)^m(s+alpha*c^2)(s+beta*c^2), with
  alpha=1-1/(ell^2+2) and beta=1+1/(ell^2-ell+2), has only negative roots,
  strictly positive coefficients, and regular field rho U'(rho)/U(rho)=ell
  at rho=-c^2.  Its even lift has a Newman-style heat threshold exactly at
  zero: every nonnegative-time heat image is real-rooted, while the universal
  backward cluster exp(d_z^2)z^m immediately has a nonreal pair.  The field
  changes only the common forward-cluster drift through
  V'(c)/V(c)=(m+4ell)/(2c).  Hence first-loss dynamics and the shared local
  Laguerre--Polya geometry impose no sign or bound on ell.  An actual Xi field
  estimate, degree-uniform Jensen remainder, collision exclusion, and every
  RH-level conclusion remain open

jensen_window_pf_newman_fourier_moment_regular_field_gate:
  at an exact multiplicity-m Xi heat-flow zero, the regular field is the
  parity-correct adjacent Fourier-moment ratio
  Re[i^(m+1)M_(m+1)]/((m+1)Re[i^mM_m]), equivalently an adjacent signed
  score-moment ratio.  For every m>=2 and every real ell, the exact sinc
  product sinc(pi*x/c)^m sinc(y*x/c), pi<y<2pi, is the characteristic
  function of a positive compact log-concave density, has a Newman-style
  threshold exactly at zero, and realizes that field.  Explicit Gaussian
  smoothing makes the density positive and uniformly strongly log-concave
  while preserving arbitrary-field surjectivity.  Hence generic Fourier
  positivity, score positivity, strong log-concavity, and local first-loss
  dynamics do not bound ell.  An actual Xi adjacent-moment estimate using
  simultaneous theta/modular and contact structure, a degree-uniform Jensen
  remainder, collision exclusion, and every RH-level conclusion remain open

jensen_window_pf_newman_modular_primitive_contact_nonclosure_gate:
  the endpoint-subtracted theta coordinate A_t=D_t[C_t]=8H_t-1/2 transfers
  every multiplicity-m contact and its regular field to an exact adjacent
  primitive-jet ratio.  A sinc-Bessel Laguerre--Polya family combines an
  arbitrary-order contact, every real ell, a smooth positive log-concave
  double-exponential-tail frequency density, and an exact first-loss heat
  boundary.  Pre-tilting places the contact at any t_*>0, while exact cosh
  normalization constructs a positive decreasing convex primitive with
  R''-R=8Phi, R'(0)=-1/2, R(-u)-R(u)=sinh(u), and a unit-mass triangular
  curvature law.  Thus the modular reflection consequence and positive
  primitive package do not bound ell.  The discrete n^(-1/2) Jacobi-theta
  translate roster, fixed arithmetic weights, coupled infinite contact
  cancellations, degree-uniform Jensen remainder, Xi collision exclusion,
  and every RH-level conclusion remain open

jensen_window_pf_newman_discrete_theta_contact_ladder_scout:
  at multiplicity two, normalization, the two contact equations, and a
  prescribed regular field form four division-free linear equations in the
  component weights.  Contact vertices therefore have support at most three,
  and Caratheodory gives prescribed-field witnesses with support at most four.
  A finite scout compares continuous, geometric, and actual integer shift
  locations.  Independent 80-decimal-place recomputation validates 24 main
  witnesses and two targeted integer refinements; all four tested rosters
  realize ell in {-20,-5,-1,0,1,5,20} when their weights are freed.  This is
  high-precision finite evidence that the tested shift locations alone do not
  supply the missing rigidity.  It is not an interval certificate or an
  all-real-field theorem.  The exact unit theta coefficients, complete
  infinite cancellation, fixed-weight Xi contact exclusion, degree-uniform
  Jensen remainder, and every RH-level conclusion remain open

jensen_window_pf_newman_fixed_theta_mass_profile_dual_scout:
  the genuine unit-coefficient profile p_(n,t)=M_(n,t)/sum_k M_(k,t)
  defines exact finite total-variation and Fisher projections onto the
  multiplicity-two contact slice.  A double-precision depth scout locates a
  sensitive region near t=1/5 and x=133.675 but fails its own resolution
  gate.  Independent 80-decimal-place integration gives selected-point
  closure distance 1.23179830387e-16 and an L1 KKT slack above 6.92, while
  exposing relative discrepancies about 0.9904 in the locator distance and
  0.999992 in its Fisher dual.  The apparent arithmetic dual pattern is
  therefore rejected.  No interval continuum optimum, infinite arithmetic
  separator, Xi contact exclusion, degree-uniform Jensen remainder, or
  RH-level conclusion is obtained

jensen_window_pf_newman_theta_fourth_summand_tangency_homotopy_scout:
  for F_lambda=H_1+H_2+H_3+lambda H_4, the exact heat identities give the
  contact Jacobian [[-F_xx,0],[-F_xxx,F_xx]] and the implicit fourth-coefficient
  response.  A 50-decimal-place continuation begins at the three-summand
  contact t=0.4320045242609588054..., x=135.5616320426760547..., remains
  regular on lambda=0,0.02,...,0.22, and crosses t=0 transversely at
  lambda=0.2016916799221276872..., x=135.8619892565584601..., with
  dt/dlambda=-2.4064083881484346....  Independent 75-decimal-place phase-split
  integration validates four anchors with residual below 1.5e-55 and
  relative jet/response drift below 5.2e-35.  This is a concrete drifting-root
  mechanism outside the fixed-root theorem, not interval continuation or a
  count of every branch.  Full-theta contact exclusion, a degree-uniform
  Jensen remainder, and every RH-level conclusion remain open

jensen_window_pf_newman_theta_fourth_summand_tangency_interval_continuation_certificate:
  eleven parameter-uniform ACB/Krawczyk charts rigorously cover
  0<=lambda<=0.22 for the tracked four-term contact.  Every chart proves one
  unique root in its local box, H_4>0, F_xx<0, and dt/dlambda<0.  Ten strict
  endpoint containments identify one connected branch, and a direct
  (lambda,x) Krawczyk solve certifies its unique t=0 crossing near
  lambda=0.2016916799221276872405642 and x=135.8619892565584600777999052.
  A 224-bit stronger-order replay independently validates three anchor charts
  and the crossing.  This does not count contacts outside the chart union,
  certify a larger compact boundary degree, add n>=5, control the complete Xi
  transform, or prove Lambda<=0, RH, or any prize-level conclusion

jensen_window_pf_newman_theta_fourth_summand_compact_contact_degree_certificate:
  352 parameter-boundary ACB/Taylor cells certify that (F,F_x) is nonzero on
  the boundary of Omega=[-0.06,0.45]x[135.5,136] uniformly for
  0<=lambda<=0.22.  At lambda=0, 32 convex image segments have winding -1 in
  (t,x), independently reproduced with three rational rays.  Boundary
  nonvanishing preserves this degree.  The interval branch contributes one
  regular contact of index -1, while every possible multiplicity-m contact
  has index -floor(m/2); therefore the tracked branch is the unique contact
  in the cylinder.  This finite four-term compact count does not continue the
  fourth coefficient to one, add n>=5, control complete Xi, or prove RH

jensen_window_pf_newman_theta_fourth_summand_positive_time_degree_certificate:
  1248 ACB/Taylor cells certify that (F,F_x) is nonzero on the boundary of
  Omega_+=[0,0.45]x[135.5,136] uniformly for 0.22<=lambda<=1.  At lambda=0.22,
  32 convex image segments have two oppositely oriented positive-ray
  crossings and winding zero, independently reproduced with three rational
  rays.  Uniform nonvanishing preserves degree zero.  Every possible contact
  has index -floor(m/2), so the rectangle contains no contact throughout this
  coefficient range.  Combined with the prior branch count, the unique local
  four-term contact exits at lambda_*=0.2016916799... and never re-enters
  before lambda=1.  Other frequency rectangles, n>=5, complete Xi, and every
  RH-level conclusion remain open

jensen_window_pf_newman_theta_fifth_summand_positive_time_degree_certificate:
  one exact coefficient slab and 32 ACB/Taylor cells certify that (F,F_x) is
  nonzero on the boundary of Omega_+=[0,0.45]x[135.5,136] uniformly for
  F=H_1+H_2+H_3+H_4+mu H_5 and 0<=mu<=1.  The minimum separation exceeds
  1.94658100827e-21, while the selected H_5 box is at most 1.317636672e-12
  of its margin.  At mu=0, 32 convex image segments have winding zero,
  independently reproduced with three rational rays.  Uniform boundary
  nonvanishing preserves degree zero, and every contact would have index
  -floor(m/2), so the finite five-term family has no contact in this local
  rectangle.  The complete theta tail and every global RH-level conclusion
  are outside this finite theorem

jensen_window_pf_newman_theta_complete_tail_positive_time_degree_certificate:
  the canonical Jacobi-theta normalization gives the pointwise component
  estimate sum_(n>=6)phi_n(u)<10^-9 phi_5(u) for every u>=0.  Rigorous
  positive moments, enlarged by an exact safety factor two, bound both the
  complete tail and its first x derivative by 6.08549613537e-42.  This is
  less than 4e-21 of the five-term boundary margin.  The homotopy from the
  five-term endpoint to the complete theta sum therefore preserves winding
  zero, and the same-sign all-multiplicity index theorem excludes every
  contact in Omega_+=[0,0.45]x[135.5,136].  This is a complete-theta theorem
  only on that local rectangle; other frequencies, negative time, the global
  de Bruijn--Newman constant, Lambda<=0, RH, and every prize-level conclusion
  remain open

jensen_window_pf_newman_theta_complete_spatial_tile_degree_certificate:
  two adjacent complete-theta tiles, [0,0.45]x[135,135.5] and
  [0,0.45]x[136,136.5], are covered by 64 ACB/Taylor boundary cells after
  adding the rigorous n>=6 derivative radii to the finite five-term boxes.
  Every initial cell succeeds, the minimum selected component separation is
  above 1.56095210922e-21, and the maximum tail-to-margin ratio is below
  3.899e-21.  Both complete boundary images are homotopic inside convex
  origin-excluding boxes to finite endpoint polygons of winding zero.  Every
  possible multiplicity-m contact has index -floor(m/2), so neither tile can
  contain one.  Composing with the central complete-theta theorem gives no
  contact on [0,0.45]x[135,136.5].  Other frequencies, negative time, all
  real x, the global de Bruijn--Newman constant, Lambda<=0, RH, and every
  prize-level conclusion remain open

jensen_window_pf_newman_outer_frequency_frontier_consolidation_gate:
  exact dependency audit showing that the 1,900-box compact theorem already
  closes 0<=t<=1/5 and |x|<=38, while the rigorous Q208 winding theorem is the
  current finite cofinal base and Q209 is the first unresolved linear stage
  with exactly two open shell regions.  The complete-theta band near x=135 is
  a valid outer calibration, not a global theorem.  The remaining direct
  domain is x>38 and 0<t<=1/5; a source-level cofinal collar-transport theorem
  or Xi-specific phase-critical joint-avoidance theorem is still required.
  This is route hygiene and proves no new outer contact exclusion or RH

jensen_window_pf_newman_c1_remainder_supersession_frontier_gate:
  target-supersession audit proving that the old endpoint-peeling remainder
  obligation is already closed uniformly on L>=50 and 0<tL<=25, including
  adjacent cutoffs.  The certified C1 bounds confine a contact to
  |X_[1]|<50000exp(-5L/4) and |U_[1]|<100000Lexp(-5L/4), and the whole-jet
  boundary transfer needs no C2 residual bound.  The live corrected-main
  target is instead the signed carrier-plus-near package of Formal Core
  Section 11.193, with endpoint half-weights and symmetric remote pairing
  preserved.  The ordinary saddle-scale target remains quantitatively
  independent.  This is dependency hygiene and proves no signed reserve,
  outer contact exclusion, Lambda<=0, or RH

jensen_window_pf_newman_c1_carrier_near_odd_harmonic_anchor_gate:
  exact carrier-near half-cell decomposition.  If D_N is the odd reciprocal
  sum over the positive reciprocal band, the near kernel has half-cell
  integrals 1/2+iD_N/pi and 1/2-iD_N/pi.  This yields complete ideal splits
  mathcal J_H^0=-A_H+R_H and mathcal J_T^0=-A_T+R_T+Q_T after the symmetric
  remote pair is composed.  On the fixed physical q=2tL^2=1 chart,
  D_N>L/4 and the raw anchors satisfy A_H>14900 and A_T>14300 for L>=50.
  The stationary carrier-variation terms R_H,R_T remain open.  The terminal
  phase rotates the negative raw centres and odd roster transfers move D_N,
  so no physical-current sign, other-q theorem, contact exclusion, or RH
  conclusion is promoted

jensen_window_pf_newman_c1_reciprocal_stationary_disk_geometry_gate:
  exact anchored reciprocal decomposition and route guard.  The Hermitian
  residual is W_H=R_H=A_H+mathcal J_H^0, while the transpose correction can
  be recombined as W_T=R_T+Q_T=A_T+mathcal J_T^0.  This improves the
  transpose handoff to |R_T+Q_T|<=A_T-delta.  Per reciprocal cell,
  2M_q-S_q=M_q-tau_q+E_q^(ext)-Lambda_q^(alias), so the anchored residuals
  retain every tie, alias, exterior, and terminal term.  The exact identity
  |A+J|^2-A^2=|J|^2+2A Re J proves that the disk is tangent at zero and
  requires a negative raw real part.  A source-specific interior block shows
  that both fully termwise carrier budgets already exceed their anchors for
  every L>=50 on q=1, retiring separated carrier/variation triangle budgets.
  No signed complete join, disk inclusion, physical-current sign, all-q
  theorem, contact exclusion, Lambda<=0, or RH conclusion is proved

jensen_window_pf_newman_c1_terminal_phase_winding_covariant_two_channel_gate:
  exact complete-cell phase theorem and two-channel current compression.  On
  fixed physical q=1, a(L)^2=exp(L)+1/(32L^2) is strictly increasing, and the
  terminal phase phi_N=Omega log N gains more than 4pi N log N>2pi across
  every complete fixed-N cell.  Hence exp(i phi_N) covers the full unit
  circle and no fixed proper terminal sector supplies a cell-uniform raw
  anchor sign.  If h=Re(conjugate(tau_0)mathcal J_H^0),
  t=tau_0mathcal J_T^0, and z is the complex source, the joined ideal current
  is (|z|^2h+Re(z^2t))/2.  Its real symmetric source-phase matrix has
  eigenvalues (h-|t|)/2 and (h+|t|)/2.  Thus the exact actual-phase target is
  h+Re(exp(2i gamma_phys)t)<0, while the stronger phase-uniform target is
  h+|t|<0.  Bare static centres fail the robust target somewhere in every
  complete cell.  No signed residual theorem, physical-current sign, all-q
  transport, contact exclusion, Lambda<=0, or RH conclusion is proved

jensen_window_pf_newman_c1_actual_source_saddle_phase_lock_center_symmetry_gate:
  exact actual-source phase theorem and bare-centre symmetry.  At xi=Omega,
  nu=eta/S_a=omega_c/(B_0T_0), omega_c=(-1)^Neta, and the transpose source
  phase is eta^2.  If chi_N=phi_N+arg eta, the reciprocal normalizer lock
  gives 2chi_N=F_N(Omega)+pi/4+psi+2arg(1+d_1), with
  F_N(v)=v[1+log(2pi N^2/v)] and |psi|<1/x.  An exact alternating-log
  estimate and a 1947/N^2 physical correction budget prove
  Delta chi_N<-2pi on every complete q=1 cell, tending to -2pi as N grows.
  Thus the absolute terminal, source, and transpose phases all cover the
  unit circle.  The actual-source bare centre factors as
  -rho|nu|^2A_T[cos(theta_eta)sin(phi_N+theta_eta)
  +(u_N/log N)sin(phi_N)], with 0<=u_N/log N<2/N.  The joint phase-product
  sign and the complete tie-invariant residual current remain open.  No
  physical-current sign, all-q transport, contact exclusion, Lambda<=0, or
  RH conclusion is proved

jensen_window_pf_newman_c1_actual_source_raw_center_sign_obstruction_gate:
  exact two-scale actual-source obstruction.  For
  chi_0=F_N(Omega)/2+pi/8, the post-saddle part of every complete q=1 cell
  has Delta chi_0<-2pi and |d chi_0/d phi_N|<2/N.  A safe positive-sine arc
  therefore contains more than one full source-phase turn.  The physical
  phase error is below 487/N^2<1/1000 and the unequal-anchor term is below
  2/N<1/1000, producing one roster-interior point with normalized centre
  bracket B>49/100 and another with B<-49/100.  Thus the actual-source bare
  centre current takes both signs on every complete q=1 cell.  At the
  positive-centre witness, negativity of the complete ideal current requires
  a grouped residual below -(98/100)rho A_T, hence
  |R_H|+|R_T+Q_T|>(98/100)A_T.  This retires a uniformly favourable raw
  centre and proves the compensation must be macroscopic and correctly
  signed on the adverse arc; it may be a uniform negative bias or a
  phase-correlated response.  It gives no signed residual
  compensation, complete current sign, all-q transport, contact exclusion,
  Lambda<=0, or RH conclusion

jensen_window_pf_newman_c1_actual_source_complete_compensation_ledger_gate:
  exact normalization, two-tone identity, and adverse-threshold gate.  The
  bracket has the form
  B=sin(phi_N+2theta_eta)/2+(1/2+u_N/log N)sin(phi_N).  This gate recorded a
  provisional eight-package non-centre split.  The later nonterminal
  determinant-completion gate proves that split omitted E_nr, so its
  package-completeness and seven-secondary statements are superseded by the
  corrected nine-package ledger.  At an adverse B<-49/100 witness the sharp
  necessary threshold remains
  E_tot<-(98/100)rho A_T; the reverse inequality at one such witness proves
  a complete-current sign change.  Coefficient and atomic phase smallness
  still lack the required multiplying-row bounds.  Uniform negative bias and
  phase-correlated response remain logically admissible proof mechanisms.
  No signed main theorem, complete-current sign, all-q transport, contact
  exclusion, Lambda<=0, or RH conclusion is proved

jensen_window_pf_newman_c1_actual_source_terminal_secondary_budget_gate:
  exact common-unit terminal budget.  On q=1, the physical exponent gives
  |E_N|<1.  The logarithmic kernel at alpha_P/N obeys S_(1,N)<L+9, so
  |kappa|<1/6 gives rho<(L+9)/6.  Since (L+9)h^3<1, the ideal transpose
  terminal survivor satisfies |E_term|<rho A_T/68640.  Using w_N=nu E_N in
  the pure terminal envelope gives |E_aa|<rho A_T/22880.  Hence the two
  distinct terminal secondary packages, each attached exactly once, obey
  the joint bound |E_term+E_aa|<rho A_T/17160.  The adverse-witness finite-package
  threshold consequently sharpens to
  E_F<-[98/100-delta_*-1/17160]rho A_T once the remaining five-package
  secondary sum has modulus at most delta_*rho A_T.  No paired-remote,
  correction, affine, moving-tail, quadratic, signed finite-package,
  complete-current, all-q, contact-exclusion, Lambda<=0, or RH theorem is
  proved

jensen_window_pf_newman_c1_actual_source_full_phase_remote_budget_gate:
  exact common-unit remote budget.  For every k>floor(2alpha), the two full
  phases alpha log u-ku and alpha log u+ku are uniformly nonstationary.  The
  paired first endpoint denominator is
  2(alpha/u)/((alpha/u)^2-k^2), so its k^(-1) terms cancel before moduli and
  the common-cutoff tail is O(k^(-2)).  Two exact integrations by parts with
  the nonoscillatory ideal-cubic amplitude give a sub-unit combined
  second-boundary/interior remainder.  Converting the two endpoint amplitudes
  through D_N>L/4 and A_T>14300 proves
  |mathscr C_N[P_H^0]|+|mathscr C_N[P_T^0]|<A_T/6, hence
  |E_R|<rho A_T/6.  Together with the terminal budget, only the physical
  correction, fixed affine, moving-tail, and nonterminal quadratic packages
  remain without common-unit bounds.  No four-package budget, signed finite
  carrier-near theorem, complete-current sign, all-q transport, contact
  exclusion, Lambda<=0, or RH theorem is proved

jensen_window_pf_newman_c1_actual_source_physical_correction_budget_gate:
  exact common-unit physical-correction budget.  The explicit terminal
  factorizations and the certified C, D, R, R_x, and delta envelopes give
  |Delta B_H|,|Delta B_T|<3h^2(L+70)^2((L+3)/2)^2 and corresponding lifted
  polynomial bounds.  The complete finite-plus-outer functional is returned
  exactly to 2mathscr M_N-mathscr F_N.  The reciprocal roster is contiguous,
  so its Dirichlet kernel has unit-period L1 norm at most 1+log m rather than
  its pointwise height.  With m<=2a^2, N<=a, and S<=0, both nonterminal
  correction channels are below 679/100 and the transpose terminal correction
  is below 1/1000.  Hence |E_Delta|<7rho<rho A_T/2000.  Joining this with the
  remote and terminal budgets leaves only the fixed affine, moving-tail, and
  nonterminal quadratic secondary packages unbounded, with normalized reserve
  697361/858000 before their budget.  No three-package budget, signed finite
  carrier-near theorem, complete-current sign, all-q transport, contact
  exclusion, Lambda<=0, or RH theorem is proved

jensen_window_pf_newman_c1_actual_source_edge_affine_ownership_gate:
  exact edge-affine ownership and signed-main reclassification.  After the
  terminal carriers are moved into the full-support moment vector, the
  external row is t_E=2J omega_E; retaining t_T repeats the moved-carrier
  interaction.  The unique affine package is
  R_aff=Re{nu c_xi mathscr O_N[P_E]} with
  P_E=i(lambda-log a)F_(mathcal N_E,V_E,-Q_E,-A_E).  At the physical source,
  E_aff=2B_0T_0 Re{omega_c mathscr O_N[P_E]}.  Its correction-free polynomial
  has leading term -iV_E(lambda-log a)^3/4, and on every p=0 midpoint
  V_E>3/16.  Thus no h^2 coefficient argument makes the complete affine
  functional perturbative.  The affine ownership and non-suppression remain
  valid; the later determinant-completion gate supersedes this gate's
  provisional two-secondary target by restoring E_nr to the signed main.  No
  absolute affine bound, signed finite-affine theorem, complete-current sign,
  all-q transport, contact exclusion, Lambda<=0, or RH theorem is proved

jensen_window_pf_newman_c1_actual_source_edge_affine_remote_budget_gate:
  exact source-normalized affine paired-remote budget.  The corrected
  endpoint normalizer and terminal multiplier give
  B_0T_0/rho<198/329<61/100 on fixed physical q=1.  The genuine edge alone
  satisfies |V_E|<3/5, |A_E|,|Q_E|<4h, and |mathcal N_E|<4h^2.  Keeping the
  two full remote phases and one common cutoff through two integrations by
  parts bounds the ideal affine remote functional by 13A_T/100.  The physical
  field drift is controlled with the contiguous-kernel L1 estimate and costs
  less than A_T/1000, so |mathscr C_N[P_E]|<131A_T/1000 and
  |E_aff^C|<4rho A_T/25.  The near affine piece remains signed in
  E_FAN=E_F+E_aff^N.  Together with the previously closed remote, terminal,
  and physical-correction packages, the absolute-secondary fraction is
  280759/858000.  This bound remains valid; the later determinant-completion
  gate supersedes the provisional remaining-package target by restoring E_nr
  to the signed main.  No affine-near sign, completed nonterminal-current
  estimate, complete-current sign, all-q transport, contact exclusion,
  Lambda<=0, or RH theorem is proved

jensen_window_pf_newman_c1_actual_source_moving_tail_budget_gate:
  exact source-free integrated moving-tail budget.  The all-carrier transfer fixes
  R_move^tr as minus the physical gap times the integral of the complete
  observation-level derivative defect, attached exactly once.  On the fixed
  q=1 chart, the correction-free weights obey
  q_n<(8/7)n^(-1/2), while all four carrier observation polynomials are below
  P_*=((L+3)/2)^2/2.  This gives all-carrier mass and first-moment envelopes
  of orders h^(-1/2) and R_*h^(-1/2), and a terminal first-moment envelope
  below 4h^(3/2)K^2.  The exact J bilinear then reduces the integrated defect
  to three monotone terms.  Source normalization, K^18<=a,
  B_0T_0/rho<61/100, rho>(329/3300)Lexp(-L/4), and A_T>14300 give
  |E_move^tr|<rho A_T/1000000.  The later ownership gate supersedes the
  former insertion of this transfer-only number into the pointwise-current
  ledger.  No pointwise moving-current or quadratic budget, signed
  finite-near-affine theorem, complete-current sign, all-q transport, contact
  exclusion, Lambda<=0, or RH theorem is proved

jensen_window_pf_newman_c1_actual_source_current_transfer_quadratic_primitive_gate:
  exact current/transfer ownership correction and quadratic determinant
  primitive.  At xi=Omega the moving-carrier displacement vanishes, and the
  established source-free observation envelopes prove the pointwise bound
  |E_move^cur|<rho A_T/100.  The corrected six-package pointwise absolute
  fraction is 289339/858000.  For the nonterminal remainder observation e,
  the four quadratic constituents are exactly the xi derivative of
  Q_quad=(2/|nu|^2)(e_Ve_mathcalN-e_Ae_Q).  Its complex polarization retains
  both Hermitian phase-difference and transpose phase-sum channels, while its
  transfer telescopes to Q_quad(Omega)-Q_quad(T_0).  The six-package budget
  gives the necessary adverse threshold 551501/858000 and triangle-safe
  sufficient threshold 1130179/858000.  The later nonterminal
  determinant-completion gate preserves those constants but supersedes the
  provisional signed-main left side by restoring E_nr.  No signed main
  inequality, endpoint determinant bound, complete-current sign, all-q
  transport, contact exclusion, Lambda<=0, or RH theorem is proved

jensen_window_pf_newman_c1_actual_source_quadratic_residual_two_carrier_kernel_gate:
  exact nonterminal quadratic two-carrier closure.  On each fixed reciprocal
  roster cell, the normalized remainder functional turns Q_quad into one
  symmetric Hermitian kernel functional plus one symmetric transpose kernel
  functional.  Their physical coefficient factorizations retain R, R_x,
  delta, C, and D exactly.  Differentiation gives the phase-difference
  multiplier i(lambda-mu) and phase-sum multiplier
  i(lambda+mu-2log a), so the Hermitian current diagonal vanishes but the
  transpose diagonal remains.  In the correction-free chart the
  half-normalized kernels are exactly the Section 11.160 squared-sum and
  squared-difference kernels, including -i*u_(N,x)/2.  The quadratic residual
  therefore creates no third phase family.  The later ledger correction
  identifies E_nr, not E_FAN, as its exact determinant partner; these kernels
  remain valid inside E_nr+E_quad^cur before the completed channel is joined
  to E_FAN.  No signed joint arithmetic estimate, endpoint determinant bound,
  complete-current sign, all-q transport, contact exclusion, Lambda<=0, or RH
  theorem is proved

jensen_window_pf_newman_c1_actual_source_nonterminal_determinant_completion_gate:
  exact nonterminal determinant-completion and compensation-ledger correction.
  The raw residual current has five blocks: retained/nonterminal linear,
  nonterminal quadratic, terminal mixed, pure terminal, and genuine-edge
  affine.  The old eight-package ledger represented the last four but omitted
  R_nr=2p^TJdot(e)+2dot(p)^TJe.  The original full-support polynomial already
  carries this term as Re{nu c_xi mathscr R_N[P_nr]}, with
  P_nr=F_(alpha_N)+i(lambda-log a)F_(beta_N-beta_E).  Restoring
  E_nr=2R_nr/|nu|^2 gives the exact completed primitive
  partial_xi[(2/|nu|^2){(p+e)^TJ(p+e)-p^TJp}]=E_nr+E_quad^cur.  The corrected
  pointwise ledger has nine packages and signed main
  E_FAN+E_nr+E_quad^cur.  The necessary and triangle-safe sufficient constants
  remain 551501/858000 and 1130179/858000 because no absolute E_nr allowance is
  introduced.  The later complete outer-determinant gate preserves this
  correction and joins it with the terminal and genuine-edge blocks.  No
  estimate or sign for the completed nonterminal current, E_FAN, complete
  current, all-q transport, contact exclusion, Lambda<=0, or RH theorem is
  proved

jensen_window_pf_newman_c1_actual_source_complete_outer_determinant_recombination_gate:
  exact complete outer-determinant recombination.  On a fixed reciprocal-roster
  cell, put o=a+e, b=omega_E+p, and Q(x)=x^TJx.  The retained/nonterminal,
  nonterminal quadratic, terminal mixed, pure terminal, and genuine-edge
  affine currents are exactly
  partial_xi[(2/|nu|^2){Q(b+o)-Q(b)}].  Eight pointwise packages therefore
  collapse to E_geom^cur, and the corrected ledger is
  E_tot^corr=E_geom^cur+E_move^cur.  The canonical outer split o=n+c gives a
  near-first determinant increment followed by a common-cutoff paired-remote
  increment, which owns the near/remote cross.  Equivalently, O_N=-D_N and
  D_N=F_N-M_N give the finite-cell form
  Q_geom=(2/|nu|^2){Q(omega_E+2p-f)-Q(omega_E+p)}, with no conditional infinite
  series.  At the calibrated level point B=-49/100, the pointwise moving-current
  budget makes the geometric thresholds 97/100 necessary and 99/100 sufficient
  independently of the moving-current sign.  The later reciprocal-block
  symmetry gate preserves this determinant and supplies its invariant matrix
  and diagonal-first arithmetic charts.  No sign or bound for E_geom^cur,
  either determinant increment, the finite-cell determinant, complete current,
  all-q transport, contact exclusion, Lambda<=0, or RH theorem is proved

jensen_window_pf_newman_c1_actual_source_finite_determinant_reciprocal_block_symmetry_gate:
  exact determinant symmetry and reciprocal-block telescope.  The observation
  vector x=(V,mathcal N,A,Q) is the matrix [[V,A],[Q,mathcal N]], whose
  determinant is x^TJx.  Thus the finite primitive is division-free:
  det(B-D)-det(B)=det(D)-tr(adj(B)D).  Simultaneous left-right actions with
  det(L)det(R)=1 preserve it, including SL(2,R)xSL(2,R), so arbitrary
  coordinate norms are not intrinsic proof objects.  A nonsingular base has
  a frozen pointwise scalar or split-diagonal normal form according to the sign
  of det(B)=h^2P_ret; rank-zero and rank-one bases remain separate null-cone
  strata unless P_ret is proved nonzero.  Splitting the finite band into
  diagonal, adjacent, and far blocks gives d_0=f^(0)-p, d_1=f^(1), d_2=f^(2)
  and an exact diagonal-first primitive and current telescope.  The order owns
  every cross; using the original base for all three increments would omit
  2d_0^TJd_1+2d_0^TJd_2+2d_1^TJd_2.  No physical block sign or bound,
  P_ret nonsingularity, geometric-current inequality, complete-current sign,
  all-q transport, contact exclusion, Lambda<=0, or RH theorem is proved

jensen_window_pf_newman_c1_calibrated_roster_interior_phase_root_scout:
  high-precision physical phase diagnostic.  In the first three complete q=1
  cells above L=50, N=72004899338, 72004899339, and 72004899340, the exact
  M_t(s), d_1, eta, Omega, and u_N formulas locate roster-interior points with
  B=-49/100, exercising both cutoff parities.  The first point has
  a-N=0.901387818864118... and N+1-a=0.0986121811358815..., and independent
  90-, 130-, and 170-digit reconstructions converge to the same row.  A
  separate checker evaluates the normalizer through log M_t and reproduces
  the phase identities.  This is not interval arithmetic.  The physical
  omega_E, retained observation, diagonal/adjacent/far block observations,
  tangents, det(B_mat)=h^2P_ret, and near/paired-remote cross-check remain
  unevaluated.  No root uniqueness, determinant sign or bound, complete-current
  sign, all-q transport, contact exclusion, Lambda<=0, or RH theorem is proved

jensen_window_pf_newman_c1_fixed_cell_low_n_determinant_evaluator_gate:
  finite fixed-cell evaluator regression certificate.  On tie-safe fixtures
  N=4,5,7, a reusable degree-six monomial basis evaluates the starred carrier,
  reciprocal diagonal, adjacent and far sectors, contiguous near kernel, and
  common-cutoff paired remote tail.  The 184 finite modes reproduce the band
  independently through the reciprocal block table, direct Dirichlet kernel,
  and endpoint-harmonic cell-variation formula.  Two integrations by parts
  and the periodic Bernoulli B_2 Fourier identity make the entire remote tail
  finite without truncating or separating its conditional halves.  Synthetic
  four-row projections then validate the diagonal-first determinant primitive
  and current telescopes and the exact i(lambda-log a) tangent lift.  An
  adaptive reciprocal-cell checker, distinct from the builder's fixed Gauss
  mode matrix, reproduces the basis and determinant current.  This is a
  floating-point low-N implementation contract, not interval arithmetic or
  physical Xi data.  No physical-N block value or sign, retained determinant
  bound, complete-current inequality, all-q transport, contact exclusion,
  Lambda<=0, or RH theorem is proved

jensen_window_pf_newman_c1_physical_observation_polynomial_compiler_gate:
  exact physical-observation coefficient compiler.  In observation order
  (V,mathcal N,A,Q), the carrier laws compile the four rows as
  (C,RQ+R_xC,RC,Q), with degrees (2,4,3,3), and agree coefficient by
  coefficient with the independent moment-vector rows (v_0,n_0,a_0,q_0).
  Translation ell=log a+x reproduces the centered cancellation before norms,
  and dot(P)=i(ell-log a)P gives the exact lifted degrees (3,5,4,4).  The
  physical alpha, beta_p, and beta_E rows assemble into one degree-at-most-five
  P_lin with the genuine edge owned exactly once.  A rational-complex
  degree-five fixture and three synthetic low-N basis replays validate the
  executable interface independently.  The complex coefficient rate is named
  chi_rate to prevent confusion with the later phase angle.  Its physical
  value, the genuine endpoint jets, retained observations, and compressed
  physical block rows remain unevaluated.  No physical determinant value,
  interval enclosure, signed-current theorem, all-q transport, contact
  exclusion, Lambda<=0, PF-infinity, or RH theorem is proved

jensen_window_pf_newman_c1_physical_scalar_edge_adapter_gate:
  high-precision physical scalar and genuine-edge interface diagnostic.  At
  the three immutable calibrated q=1 roots N=72004899338, 72004899339, and
  72004899340, the exact fixed-t, fixed-N formulas evaluate rho_1,rho_2 and
  their derivatives, s_*',s_*'', chi_rate, u_(N,x), and the genuine endpoint
  row omega_E=(V_E,mathcal N_E,A_E,Q_E).  The carrier-normalizer and stable
  terminal-W_x routes independently agree for chi_rate, and the
  90-, 130-, and 170-digit ladder is stable.  A separate checker interpolates
  the exact correction quotient, directly differentiates the carrier and
  endpoint fields, and reconstructs twelve base rows, twelve tangent rows,
  and three single-owner degree-five edge-affine polynomials.  This is not
  interval arithmetic.  No retained alpha or beta_p observations,
  compressed physical blocks, determinant sign or bound, complete-current
  inequality, all-q transport, contact exclusion, Lambda<=0, PF-infinity,
  or RH theorem is proved

jensen_window_pf_newman_c1_shifted_hardy_diagnostic_bridge_gate:
  shifted-Hardy retained-carrier diagnostic.  At the actual carrier height
  T_0=Omega=2pi*a^2, not the unshifted x/2, the parity-adjusted source phase
  aligns numerically with the stationary Riemann-Siegel endpoint direction
  at all three calibrated roots.  Fifteen real samples at delta_j=0.01j
  implement fitted complex kernels for all four base and four tangent rows.
  Across three roots and twenty-four endpoint-inclusive held-out grids, the
  maximum relative interpolation error is below 3.81e-11.  The pinned
  GPL-3.0 external evaluator is repaired by four minimal OpenMPI portability
  patches and independently calibrated at T=10^10 and T=10^12.  Wider
  spacings 0.02 and 0.04 fail direct-main tests and are rejected.  The largest
  coefficient l1 norm is about 2.312e6, no validated asymptotic error constant
  or interval enclosure is supplied, and the non-resumable physical job is
  not launched.  No physical alpha or beta_p observation, determinant sign or
  bound, complete-current inequality, all-q transport, contact exclusion,
  Lambda<=0, PF-infinity, or RH theorem is proved

jensen_window_pf_newman_c1_hardy_resumable_evaluator_equivalence_gate:
  finite computational reproducibility certificate for the accepted h=0.01
  fifteen-value Hardy evaluator.  A GPL-3.0 derivative commits all direct,
  Gaussian, and disjoint Riemann--Siegel tail units through a full-precision
  atomic snapshot and fsynced append-only journal.  At T=10^10 the exact
  uninterrupted and four-times-parked paths each contain 40 states; at
  T=10^12 the uninterrupted and twice-parked paths each contain 56 states.
  The independent checker compares 192 arbitrary-precision records, exact
  final snapshots and displayed values, six parked invocations, and two
  fail-closed provenance/chunk mismatches.  Both heights reproduce their
  accepted ten-decimal fixtures with zero displayed difference.  This proves
  restart equivalence only: no uniform external error constant, joint
  correlated-error bound, physical carrier value, interval enclosure,
  retained observation, determinant sign or bound, complete-current
  inequality, all-q transport, contact exclusion, Lambda<=0, PF-infinity,
  or RH theorem is proved

jensen_window_pf_newman_c1_hardy_joint_correlated_error_ledger_gate:
  finite two-height correlated-error diagnostic.  The fifteen accepted
  shifted-Hardy residuals are decomposed exactly in a discrete orthonormal
  polynomial basis before contraction with all twenty-four saved physical
  coefficient vectors.  Twenty-one vectors annihilate a constant error to
  relative sensitivity below 1e-12, and the observed correlations reduce
  independent L1 bounds by factors up to about 1.77e6.  However, while the raw
  maximum residual falls by a factor 0.398752 from T=10^10 to T=10^12,
  twenty-one of twenty-four projected errors increase, by as much as 3.48541.
  Six pseudoinverse cutoffs from 1e-8 through 1e-16 preserve the same 21/24
  obstruction and nearly the same maximum projections.  Scalar two-height
  extrapolation is therefore retired in favor of explicit physical-height
  bounds for modes zero through five and the rough L2 remainder.  No
  asymptotic mode theorem, physical carrier value, interval enclosure,
  determinant sign or bound, complete-current inequality, all-q transport,
  contact exclusion, Lambda<=0, PF-infinity, or RH theorem is proved

jensen_window_pf_newman_c1_hardy_modewise_derivative_handoff_gate:
  exact conditional derivative-to-mode reduction and primary-source audit.  On
  the fifteen-point grid delta_j=j/100, exact integer Gram polynomials and
  Taylor's theorem turn uniform bounds on the external error derivatives
  through order six into the six low-mode bounds and the orthogonal L2 tail
  required by Section 11.236.  The rough-tail multiplier has exact square
  4069880729/64800000000000000000000000000.  One fixed-branch holomorphic
  disk bound is a sufficient alternative by Cauchy's estimate.  The audited
  paper equations (120)--(123) reduce the MGS error to the largest local
  iteration error, but the paper explicitly leaves the precise higher-order
  bound unformulated; its final accuracy study is sample evidence.  The
  pinned fifteen-value source shares a central hierarchy but supplies no
  derivative contract.  No evaluator derivative or disk bound, transition
  exclusion, external absolute error constant, physical carrier value,
  interval enclosure, determinant sign or bound, complete-current
  inequality, all-q transport, contact exclusion, Lambda<=0, PF-infinity,
  or RH theorem is proved

jensen_window_pf_newman_c1_hardy_exact_low_mode_projected_kernel_gate:
  exact rational low-mode cleanup with finite kernel replay.  Orthogonal
  subtraction from the twenty-four saved decimal coefficient vectors gives
  nearby rational vectors satisfying eighty-four selected Gram-mode moments
  exactly.  Their largest L2 change is about 8.26410e-8, while a 180-digit
  endpoint-inclusive replay keeps every target-normalized held-out error below
  9.322e-11.  The resulting largest normalized derivative multipliers through
  order six are about 0.484617, 0.0193847, 0.000928849, 3.10155e-5,
  6.49848e-7, 2.03324e-8, and 5.20556e-8.  The 21/24 two-height
  scalar-extrapolation obstruction remains.  No external derivative theorem,
  interval kernel approximation, physical carrier value, retained
  observation, determinant sign or bound, complete-current inequality,
  all-q transport, contact exclusion, Lambda<=0, PF-infinity, or RH theorem
  is proved

jensen_window_pf_newman_c1_hardy_local_defect_obligation_ledger_gate:
  exact local-defect reduction and componentwise source audit.  Paper W1,
  W2, W3/W4, W5, and the endpoint half-sum map respectively to Fortran t5,
  t2, t4, t1, and t3.  Thirteen analytic, model, numerical,
  outer-representation, multi-shift, and transition obligations are kept
  separate.  The executable fixes ip=3 even though its q-routine comment
  says ip near twenty was found adequate; it also contains finite
  special-function approximations and explicitly disabled or omitted
  corrections that require bounds.  Independently of those internals,
  equation (120) defines each local defect as a finite difference.  Exact
  affine expansion proves that directly evaluating a parent shell removes
  every skipped descendant defect, including the first post-kernel defect
  identified by the paper as usually largest.  No W1 contour or tail bound,
  transformed-level adapter, evaluator C6 or disk theorem, transition
  exclusion, outer Hardy representation bound, physical carrier value,
  interval kernel approximation, determinant sign or bound,
  complete-current inequality, all-q transport, contact exclusion,
  Lambda<=0, PF-infinity, or RH theorem is proved

jensen_window_pf_newman_c1_hardy_block20_native_q_required_compensation_gate:
  finite exact-arithmetic inversion of the emitted block-20 recurrence.  For
  each of 374 replayed parent calls, the unique additive q correction needed
  to make the logged parent equal its replayed children is reconstructed at
  high precision.  Every required correction is nonzero, with magnitude from
  about 4.67e-3 to 2.06e-1, and the inversion and orientation identities close
  far below the reported scale.  This localizes the missing contribution; it
  does not identify its analytic source, prove a contour-remainder bound,
  validate the external Hardy value, or prove Lambda<=0 or RH

jensen_window_pf_newman_c1_hardy_block20_ecor_redundancy_gate:
  exact symbolic source audit of the commented quartic endpoint correction.
  Expanding the cubic stationary-point expression used by the executable
  reproduces the displayed ecor coefficient exactly, while source inspection
  finds no executable phase use of ecor.  Adding it separately would therefore
  double count that quartic term in the current source model.  This does not
  reconstruct an unavailable Maple version, bound omitted higher saddle
  terms, validate equation-(69) contour deformation, or prove Lambda<=0 or RH

jensen_window_pf_newman_c1_hardy_block20_w1_finite_ip_completion_gate:
  finite interval completion of every index displayed in the paper's W1
  endpoint formula at source setting ip=3.  Of 374 calls, 204 already cover all
  displayed indices and the other 170 omit 175 terms.  Their full correlated
  contribution is at most about 4.56e-4, under 3.56 percent of the native-q
  requirement, and all 374 completed residuals remain nonzero.  This excludes
  displayed-index clipping as the explanation, but it is not an enclosure of
  the exact equation-(69) contour remainder or of unprinted asymptotic terms

jensen_window_pf_newman_c1_hardy_block20_t5_sqrt2_literal_kind_gate:
  finite literal-kind audit of the two unsuffixed Fortran sqrt(2.0) calls.
  The source first rounds sqrt(2) to binary32, shifting individual arguments by
  at most about 5.25e-7 and correlated q values by at most about 3.75e-8.
  Correcting this alone, and stacking it with the complete displayed W1 index
  set, leaves all 374 native-q residuals nonzero.  The gate isolates source
  literal semantics only; it does not alter pi normalization, bound the exact
  contour remainder, validate the external Hardy value, or prove Lambda<=0 or
  RH

jensen_window_pf_newman_c1_hardy_block20_exact_eq69_saddle_family_gate:
  rigorous finite replacement of all 544 equation-(69) entire-function
  integrals across the 374 saved recurrences.  Independent 70/110-digit Arb
  evaluations overlap, the two cancellation constructions agree, and every
  post-W1 residual remains nonzero between about 5.29e-2 and 8.05e-2.  This
  clears the complete finite W1 saddle family as the sole explanation but does
  not bound nonsaddle W2--W5, outer Hardy, or height-uniform errors, or prove
  Lambda<=0 or RH

jensen_window_pf_newman_c1_hardy_block20_t2_lplus1_mismatch_gate:
  exact source/paper audit of the source-only `-1/(L+1)` inside `t2`.  Its
  induced q term is fixed algebraically and absent from published equation
  (89).  Removing it after exact W1 improves all 374 residuals and reduces the
  worst magnitude ratio below 0.02, while residuals from about 1.32e-6 to
  1.54e-3 remain nonzero.  This is a finite source-correction certificate, not
  a uniform W2--W5 estimate, evaluator patch authorization, or proof of
  Lambda<=0 or RH

jensen_window_pf_newman_c1_hardy_block20_w2_w5_correspondence_gate:
  rigorous source/paper orientation audit of the endpoint half-sum and W2--W5.
  The final source conjugation is retained explicitly, exact dyadic xi and
  endpoint transports are recorded, and 90/150-digit Arb evaluations cover
  165 W3 and 209 W4 calls.  Both complete residual constructions overlap;
  all 374 residuals remain nonzero between about 1.32e-6 and 1.53e-3, with
  52 improved and 322 worsened.  This is finite published-formula evidence,
  not an exact nonsaddle contour or uniform recurrence theorem

jensen_window_pf_newman_c1_hardy_block20_cubic_vertical_contour_admissibility_gate:
  exact cubic far-field audit of equations (82)--(91).  The missing cubic
  term makes the stated (67b) vertical ray grow on all 187 negative-Phi3
  calls and the stated (67c) ray grow on all 187 positive-Phi3 calls.  No
  saved cubic admits both printed infinite vertical contours.  The initial
  pi/6 cubic-decay candidate is sharpened by the following sign-aware gate;
  the height-uniform recurrence remains open

jensen_window_pf_newman_c1_hardy_block20_cubic_sector_contour_repair_gate:
  exact finite-roster repair of the failed vertical contours.  The negative-
  Phi3 (67b) branch uses 5*pi/6, the positive-Phi3 (67c) branch uses pi/6,
  and the opposite signs retain pi/2.  Strict convexity and positive linear,
  quadratic, and cubic decay are certified on all 374 calls, making the far
  connector vanish and the outer 1/n sums absolutely interchangeable with
  the rays.  The quadratic models rotate exactly to the paper's vertical
  error-function contours.  The remaining target is to enclose each full
  cubic endpoint ray against that quadratic W2--W4 model; no uniform-height
  recurrence, outer Hardy control, Lambda<=0, or RH theorem follows

jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_reduction_gate:
  exact summation of the repaired outer 1/n tails into the principal-logarithm
  kernel `K_m`.  Integration by parts and explicit source-orientation
  bookkeeping reduce `P^- - I69^-` to four logarithmic endpoint rays, one
  optional positive-Phi1 zero mode, and endpoint terms, with no infinite
  numerical n cutoff.  Four sign-complete mpmath comparisons close between
  about 2.1e-20 and 5.1e-49 but are labelled non-rigorous diagnostics.  The
  following two interval gates supply the finite Arb enclosure; no uniform
  recurrence, Lambda<=0, or RH theorem follows

jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_arb_pilot_gate:
  rigorous four-branch Arb enclosure of the logarithmic endpoint rays.  An
  analytic majorant covers the omitted logarithmic origin, an explicit
  polynomial-times-exponential bound covers infinity, and a cancellation-free
  hypergeometric kernel covers the compact pieces.  At 70 and 110 decimal
  digits both exact nonsaddle identities enclose zero for chains 1, 2, 3, and
  36.  This is a sign-complete finite pilot, not an all-height theorem

jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_arb_full_gate:
  rigorous extension of the Arb construction to all 374 saved calls.  Exact
  adaptive tail bounds select cutoffs from 64 to 1024 and stay below 1e-30;
  both independently assembled identities enclose zero on every call with
  maximum gap below 4.156e-14.  Every exact full-cubic correction excludes
  zero between 1.317477937491e-6 and 1.530207801447e-3.  This closes the saved
  finite contour calculation, but not a uniform-height recurrence, outer
  Hardy control, Lambda<=0, PF-infinity, or RH

jensen_window_pf_newman_c1_hardy_block20_endpoint_linearization_decomposition_pilot_gate:
  exact identification of the paper's hybrid endpoint model.  The generic
  endpoint-linear rays sum by digamma identities to the endpoint half-sum and
  W5; W2 and W3 are exceptional quadratic-mode replacements and W4 is the
  positive-Phi1 zero-mode replacement.  Thus the exact correction equals
  `R_b+R_c+zero_mode-W2-(W3 or W4)`.  Arb validates this and the W5 identity
  on all four sign branches, but does not yet give an all-height bound

jensen_window_pf_newman_c1_hardy_block20_endpoint_linearization_decomposition_full_gate:
  rigorous 70/110-digit extension of both decomposition identities to all 374
  calls.  The maximum identity gaps are below 1.824e-106 and 5.237e-14.  A
  fully componentwise triangle bound can exceed the exact correction by about
  138940, exposing essential cancellation rather than a usable uniform error
  budget

jensen_window_pf_newman_c1_hardy_block20_paired_exceptional_cancellation_gate:
  rigorous pairing of each exact endpoint remainder with the exceptional
  W2--W4 approximation it replaces.  The worst triangle/correction ratio
  falls to 210.274137, more than a 660-fold improvement, but chain 260 still
  requires cancellation between the upper and lower endpoint families.  The
  ratios are finite diagnostics, not fitted constants or a recurrence theorem

jensen_window_pf_newman_c1_hardy_block20_exceptional_mode_recombination_pilot_gate:
  exact identification of W2 and W3 as linear-minus-quadratic exceptional
  mode corrections and W4 as a quadratic zero-mode plus endpoint boundary.
  The paper models use global Phi2.  Direct Arb mode integrals validate the
  resulting reciprocal-free identities on all four sign branches; this is a
  finite pilot rather than a height-uniform bound

jensen_window_pf_newman_c1_hardy_block20_exceptional_mode_recombination_full_gate:
  rigorous 70/110-digit extension of the reciprocal-free recombination to all
  374 calls, including 165 W3 and 209 W4 branches.  Model gaps stay below
  1.857e-106 and complete correction gaps below 5.237e-14.  A four-piece
  triangle can still lose a factor about 6484, so endpoint-family coupling
  remains mandatory

jensen_window_pf_newman_c1_hardy_block20_shared_contour_exceptional_majorant_gate:
  exact common-ray homotopy bounds for the W2 and W3 exact-minus-quadratic
  modes.  The Gaussian majorants contain neither delta nor eta and remain
  finite at hypothetical exact tangency.  Arb verifies positive slack on all
  374 W2 and 165 W3 modes.  The W4 zero-mode pair, generic coupled tails,
  recurrence propagation, outer Hardy control, Lambda<=0, and RH remain open

jensen_window_pf_newman_c1_hardy_block20_w4_contour_endpoint_hurwitz_majorant_gate:
  exact Cauchy deformation of the positive-Phi1 zero mode into endpoint
  half-rays, identifying the W4 endpoint boundary as the linear n=0 ray at N.
  A common-ray homotopy sums all generic exact-minus-linear modes by explicit
  Hurwitz zeta moments.  Together with the exceptional bounds this gives a
  selector-gap-uniform majorant for the complete local endpoint correction,
  certified on all 374 calls.  Global recurrence and outer Hardy budgets
  remain open

hardy_block20_accumulation_weight_fixture:
  exact finite source-equivalence fixture containing all 3180 binary128
  block-20 amplitude/phase rows.  The generated observer leaves the checkpoint
  state and all fifteen displayed Hardy values unchanged.  This inventories
  weights only; it proves no approximation or outer Hardy theorem

jensen_window_pf_newman_c1_hardy_block20_transported_endpoint_budget_gate:
  rigorous transport of the first complete local endpoint majorants through
  the exact 374 recursive and 50 direct source weights.  The absolute endpoint
  budget is about 0.00858014 on every output and exceeds 0.005 by a factor at
  most 1.716028, while the separate direct-kernel budget is below 1.886e-32.
  The finite signed correction is diagnostic and is not used as a proof bound

jensen_window_pf_newman_c1_hardy_block20_retained_decay_endpoint_majorant_gate:
  exact common-ray theorem retaining the positive homotopy decay through
  `(1-exp(-x))/x`, with selector-correct W2/W3/W4 linear gaps that remain
  finite at their boundaries.  Arb certifies all 374 local applications at
  60/90 digits.  Exact transport gives a maximum endpoint-model budget below
  0.002189215, so all fifteen block-20 outputs lie below the nominal 0.005
  local scale without signed cancellation.  Coefficient neighborhoods,
  Legendre tails, other blocks, and outer Hardy control remain open

jensen_window_pf_newman_c1_hardy_block20_coefficient_neighborhood_transport_gate:
  rigorous transport of the retained-decay endpoint functional over all 374
  occupied coefficient boxes, preserving all 165 W3 and 209 W4 selector
  choices.  The maximum exact-weight output budget is below 0.002196524, so
  all fifteen outputs remain below 0.005.  Separate phase, child, multiplier,
  saddle, and fracL arithmetic columns show that the broad selector cells do
  not close the complete recurrence; this is an endpoint-only finite theorem

jensen_window_pf_newman_c1_hardy_block20_complete_cubic_legendre_tail_budget_gate:
  complete saved-point cubic Legendre phase-and-amplitude remainder from the
  exact stationary radical, covering all 753 child indices in the 374
  recursive calls.  Exact weight transport contributes less than 0.000150954;
  endpoint plus tail remains below 0.002347477 on all fifteen outputs without
  signed cancellation.  Coefficient-cell tail transport, source arithmetic,
  other blocks, and outer Hardy control remain open

jensen_window_pf_newman_c1_hardy_block20_coefficient_cell_complete_cubic_legendre_tail_budget_gate:
  exact z-squared factorization of the complete cubic Legendre phase remainder
  and common-phase removal over all 374 selector cells and 753 child indices.
  The transported cell-tail column is below 0.000169449; endpoint plus cell
  tail remains below 0.002365972 on all fifteen outputs without signed
  cancellation.  Parent/saddle arithmetic, other blocks, and outer Hardy
  control remain open

jensen_window_pf_newman_c1_hardy_block20_coefficient_cell_joint_recurrence_identity_gate:
  exact correlated identity W69(M,C)=I69-M*C on the physical slices of all 374
  occupied coefficient cells.  Child and multiplier changes cancel inside W69,
  reducing the corrected recurrence error exactly to Q_exact-Q_paper.  The
  transported endpoint majorant is below 0.002196524 on all fifteen outputs.
  Source exact-W1 replacement, t2, special-function and floating arithmetic,
  other blocks, and outer Hardy control remain open

jensen_window_pf_newman_c1_hardy_block20_source_to_corrected_model_bridge_gate:
  signed finite composition of exact W1 replacement, source-only t2 deletion,
  2*pi normalization, and complete W2--W5 replacement on all 374 calls.  Every
  inter-gate interval handoff overlaps the independent endpoint correction.
  PSI/ERF replacement accounts for the complete nonsaddle transport up to
  1.582e-9, and the final output majorant is below 0.002196524 on all fifteen
  outputs.  Source recompilation, full-cell arithmetic, other blocks, and outer
  Hardy control remain open

jensen_window_pf_newman_c1_hardy_t1e10_crossblock_structural_atlas_gate:
  exact-binary finite route atlas for all 6784 cubic calls in blocks 20--35 at
  t=10^10.  It records 1414 recursive and 5370 direct calls, proves strict
  selector, dual-floor, curvature, and stationary-radical margins on every
  recursive parent, and locates the final recursive calls in block 28.  Later
  endpoint errors, direct kernels, outer weights, and height uniformity remain
  open

jensen_window_pf_newman_c1_hardy_t1e10_crossblock_retained_endpoint_gate:
  retained common-ray endpoint theorem on every one of the 1414 recursive
  calls in blocks 20--28.  The roster-wide maximum remains the block-20 value
  below 0.004879658, and every later block maximum is at most 0.856694 of it.
  These are local call majorants; cross-block output weights, coefficient-cell
  transport outside block 20, and all direct-kernel errors remain open

hardy_crossblock_accumulation_weight_fixture:
  source-equivalent exact-binary128 observer fixture for all 50880 outer
  accumulation rows in blocks 20--35 at t=10^10.  Its deterministic roster
  carries 407040 payloads and leaves the evaluator checkpoint and fifteen
  displayed values exactly unchanged.  This records coefficients but proves
  no error bound

jensen_window_pf_newman_c1_hardy_t1e10_crossblock_recursive_output_budget_gate:
  absolute outer-weight transport of all 1414 recursive endpoint majorants.
  The block-20 subtotal reproduces the admitted output values exactly; blocks
  21--28 add at most 0.003616586, and the complete recursive column is below
  0.005805801, or 1.161161 times the nominal 0.005 scale.  Direct kernels,
  later coefficient cells, and outer Hardy control remain open

jensen_window_pf_newman_c1_hardy_t1e10_crossblock_direct_kernel_output_gate:
  two-precision exact-source finite-sum enclosures and analytic source-to-2*pi
  phase transport for all 5370 direct calls in blocks 20--35.  Their complete
  outer contribution is below 3.851e-29.  Adding it to the recursive column
  gives a complete 6784-call exact-point bound below 0.005805801, so the direct
  route is closed and negligible while the recursive triangle bound remains
  above the nominal scale

jensen_window_pf_newman_c1_hardy_t1e10_crossblock_hybrid_exact_output_gate:
  unsigned hybrid transport using the rigorous full-contour correction
  magnitudes for all 374 block-20 calls, retained majorants for the 1040 later
  recursive calls, and every certified direct kernel.  The corrected-model
  endpoint/direct column is below 0.003807296 on all fifteen outputs, leaving
  more than 0.001192704 below 0.005.  Later full-cell, source-arithmetic,
  outer-Hardy, and height-uniform bridges remain open

jensen_window_pf_newman_c1_hardy_t1e10_blocks21_28_coefficient_cell_signed_bridge_pilot_gate:
  two nonzero-radius selector-stable cells at the worst block-21 and later
  exact/retained-ratio witnesses.  Correlated 70/110-digit enclosures preserve
  the signed Q_exact-Q_paper correction with less than 2.79 percent inflation
  and retain the finite outer margin.  Source rounding and complete recurrence
  transport remain open

jensen_window_pf_newman_c1_hardy_t1e10_blocks21_28_coefficient_cell_signed_bridge_full_gate:
  complete analytic corrected-model cell transport for all 1040 recursive
  calls in blocks 21--28.  Every nonzero cell overlaps its exact-point
  correction and remains inside the retained local majorant.  Exact outer
  transport keeps all fifteen complete outputs below 0.000448799 without
  cross-call cancellation.  Source binary128 cell rounding, complete later
  recurrence tails, outer-Hardy control, and height uniformity remain open

hardy_binary128_rounding_mode_probe:
  pinned source-environment probe verifying radix-2, 113-bit binary128 and
  startup round-to-nearest mode under the accepted Docker image and compiler
  flags; source special-function and full operation accuracy remain separate

jensen_window_pf_newman_c1_hardy_t1e10_later_complete_recurrence_stress_cell_pilot_gate:
  four-cell complete-tail and source-rounding pilot on chains 626, 653, 3406,
  and 3423; exact rounding preimages, complete factored Legendre tails, and
  source complex multiply-add envelopes preserve all fifteen outer margins

jensen_window_pf_newman_c1_hardy_t1e10_later_complete_recurrence_cell_full_gate:
  complete later-cell transport of the factored cubic Legendre tail and pinned
  source multiply-add rounding on all 1040 calls.  All fifteen partial complete
  outputs stay below 0.001185344 without cross-call cancellation.  Source q
  and special-function cell variation remains open

jensen_window_pf_newman_c1_hardy_lower_fold_airy_fresnel_detuning_ode_gate:
  exact finite-Y Airy-Fresnel detuning equation with both endpoint sources
  retained; all 84 saved coefficients are nondegenerate and three independent
  moment residuals contain zero.  All-detuning propagation and the outer join
  remain open

jensen_window_pf_newman_c1_hardy_lower_fold_detuning_lattice_quadrature_gate:
  complete 84-mode interval lattice at the saved height.  Its grouped sum
  overlaps the independent Airy-first Dirichlet-kernel evaluation componentwise
  to better than 1.5e-27 in midpoint discrepancy.  Height uniformity remains
  open

jensen_window_pf_newman_c1_hardy_lower_fold_height_transport_identity_gate:
  exact fixed-selector law lambda_t=-1/beta and complete grouped point
  derivative at t=10^10.  This identity alone is not a nonzero-radius theorem

jensen_window_pf_newman_c1_hardy_lower_fold_nonzero_height_cell_gate:
  two independent interval partitions prove |Ai|<0.54 on the full local
  argument range and hence |S''|<0.146 on |t-10^10|<=0.01.  The Taylor
  remainder is below 7.3e-6 and total movement below 0.001313.  Selector-fold
  joins and source-height coverage remain open

jensen_window_pf_newman_c1_hardy_lower_fold_displayed_height_window_gate:
  the source roster contains exactly the 15 offsets -0.07 through +0.07.
  The exact inverse-root margins preserve the same 84 turning modes and the
  same 42+42 split throughout.  Two independent interval partitions prove
  |S''|<0.146, with Taylor remainder below 0.000358 and total canonical
  lower-fold movement below 0.009494.  Printed Hardy values and the complete
  source formula are not promoted

jensen_window_pf_newman_c1_hardy_lower_fold_turning_event_atlas_handoff_gate:
  exact triangular-number ordering of all 399 fixed-selector turning events,
  with minimum spacing pi.  Every internal ownership change preserves the
  exact finite mode sum by a one-mode handoff.  The C-to-C+2 selector jump is
  correctly separated as a 399-mode removed-lower-strip problem, not a
  one-mode event

jensen_window_pf_newman_c1_hardy_lower_fold_selector_strip_poisson_reassembly_gate:
  exact fixed-upper-endpoint splice when C advances to C+2.  The endpoint
  half-current plus the complete symmetric one-cell Poisson strip reconstructs
  the one removed source lattice term.  The 399 positive turning modes alone
  are explicitly insufficient; zero, negative, and endpoint pieces remain
  mandatory

jensen_window_pf_newman_c1_hardy_upper_endpoint_selector_boundary_stability_gate:
  the pinned 35-block endpoint recurrence reproduces saved B=5122421.  The
  source transition advances the effective start about 8640.59 before the
  selector boundary, after which B=5122423 is interval-stable through the
  boundary and for 200000 height units above it.  There is no simultaneous
  upper strip at C-to-C+2; reassembly of the earlier coupled transition term,
  direct-roster shift, and added upper lattice point remains open

jensen_window_pf_newman_c1_hardy_source_transition_endpoint_scheduling_gate:
  exact finite-roster decomposition of the earlier source transition.  The
  moving rosters differ by removal of alpha=159577 and addition of
  alpha=5122423, while the RS cutoff remains 621.  The source change is exactly
  transition-defect plus added-upper-endpoint; the papers classify the
  transition and hybrid as approximate.  The theorem chart must therefore
  hold B=5122423 fixed and keep source approximation error in a separate
  ledger.  The quadrature value is diagnostic, and no complete source-error or
  quantitative 399-mode bound is proved

jensen_window_pf_newman_c1_hardy_fixed_B_selector_boundary_399_mode_fourier_completion_gate:
  exact one-period geometry for the fixed-B selector cell and rigorous interval
  quadrature of its 399 joint-saddle canonical modes.  The modes are a
  symmetric Fourier block; Dirichlet-Jordan completes all zero, negative, and
  outer-positive modes jointly as the endpoint midpoint minus that block.
  Adding the endpoint half-current reconstructs the canonical removed source
  term exactly.  The exact finite-t Kummer-to-canonical error, complete
  T_upper theorem, Lambda<=0, and RH remain open

jensen_window_pf_newman_c1_hardy_fixed_B_selector_boundary_finite_t_completed_strip_error_gate:
  exact all-mode collapse of the finite-t and canonical selector strips to
  one y=0 fold integral at the boundary.  A rigorous pi/6 contour enclosure
  bounds the equation-(9) normalized difference below 1.2e-16.  The
  integrated beta^-2 correction cancels and the first formal term is of order
  beta^-4.  A nonzero height cell and ordinary-Morse join remain separate

jensen_window_pf_newman_c1_hardy_fixed_B_selector_boundary_nonzero_height_completed_strip_cell_gate:
  event-free extension to t*-pi/16 through t*+pi/16.  The exact Fourier
  completion persists, and a center derivative plus completed-difference
  second-derivative majorant bounds the normalized strip error below 8.9e-13
  throughout.  The Airy-to-ordinary-Morse remainder join and complete T_upper
  theorem remain open

jensen_window_pf_newman_c1_hardy_lower_fold_turning_event_fresnel_retention_gate:
  exact conversion between the Airy defect, turning height, logistic
  characteristic, and reduced-saddle Fresnel endpoint coordinate.  Every one
  of the 399 event-isolation envelopes remains within 0.004976 Fresnel units
  of the half-saddle, so replacing the transition by its full-interior value
  has relative defect above 0.496.  This rules out that replacement, not the
  complete grouped mode; a Fresnel-retaining join remains open

jensen_window_pf_newman_c1_hardy_lower_fold_fresnel_retaining_logistic_transform_gate:
  exact positive-mode logistic-Morse transform with the endpoint exponential
  and Fresnel difference kept in one amplitude.  The phase is globally
  Gaussian in the Morse coordinate, q_C=0 at a turning event, the Jacobian is
  regular, and no characteristic inverse is used.  Amplitude transport, the
  Airy/logistic remainder, and 399-event propagation remain open

jensen_window_pf_newman_c1_hardy_lower_fold_hyperbolic_morse_fresnel_overlap_core_gate:
  exact common hyperbolic phase chart for the retained-Fresnel logistic and
  Airy descriptions of mode 39894.  Independent interval covers at the event
  and both pi/16 faces certify |R_t|<0.082 and |A_t-1|<0.001 on |u|<=200,
  with carrier mismatch below 2e-12.  Integration before absolute values,
  the |u|>200 tails, the complete one-mode join, and 399-event propagation
  remain open

jensen_window_pf_newman_c1_hardy_lower_fold_hyperbolic_morse_fresnel_compact_integral_gate:
  closed-form Fresnel evaluation of the inner quadratic phase and direct
  interval integration of the signed outer difference for mode 39894.
  Independent partitions certify normalized compact error below 0.00049 and
  physical error below 0.000221 at the event and both pi/16 faces.  The
  complementary real line is deliberately a separate claim

jensen_window_pf_newman_c1_hardy_lower_fold_hyperbolic_morse_fresnel_full_line_join_gate:
  common physical (z,y) contour completion of the retained-Fresnel logistic
  and Airy charts for prototype mode 39894.  Arb Taylor/Cauchy enclosures,
  vertical connectors, signed lifted segments, and analytic far tails prove
  normalized full-line error below 0.000019 and physical error below
  0.0000086 at three heights.  A continuous-height theorem, propagation
  through the other 398 events, and complete T_upper remain open

jensen_window_pf_newman_c1_hardy_lower_fold_hyperbolic_morse_fresnel_continuous_height_cell_gate:
  exact common transport D_t=exp(i(t-tau)z/beta)D_tau and certified contour
  moments through degree 12 upgrade the three-height prototype join to every
  |t-tau|<=pi/16.  The uniform normalized bound is below 0.000017246 and the
  physical bound below 0.000007764; a finer degree-14 reconstruction passes.
  Event-index uniformity, the other 398 cells, and complete T_upper remain open

jensen_window_pf_newman_c1_hardy_lower_fold_event_parameter_contour_geometry_gate:
  exact 399-event detuning roster and saddle formulas show that R=9 misses the
  outer saddle for 246 events.  A common R=14, Im z=1 contour has positive
  saddle, Cauchy-disk, hyperbolic, and far-tail margins throughout the atlas.
  No uniform signed remainder estimate follows from geometry alone

jensen_window_pf_newman_c1_hardy_lower_fold_event_parameter_extreme_full_line_pilot:
  rigorous signed endpoint diagnostic on the common contour.  The bare Airy
  defects at the two extreme detunings are 0.0123535 and 0.0127987, so the
  prototype 0.000019 budget cannot be propagated unchanged.  Intermediate
  detunings and a corrected normal form are separate claims

jensen_window_pf_newman_c1_hardy_lower_fold_event_parameter_extreme_corrected_normal_form_gate:
  the exact beta^-2 amplitude and phase correction removes over 99.8 percent
  of both extreme defects, but independently certified residuals 2.5067e-5
  and 2.5456e-5 remain above target.  This rejects stopping at first order,
  not the complete finite-height normal form

jensen_window_pf_newman_c1_hardy_lower_fold_event_parameter_extreme_second_order_normal_form_gate:
  the complete beta^-4 multiplier, integrated with H_0 through H_4 on the
  common contour, gives extreme residuals below 3.571e-8 and physical
  residuals below 1.608e-8.  An independent denser reconstruction passes.
  The intermediate-detuning and all-cell height theorems remain open

jensen_window_pf_newman_c1_hardy_lower_fold_event_parameter_second_order_all_event_atlas_gate:
  event-centered rational real-contour partitions and factored exact and
  canonical phase derivatives certify the complete beta^-4 residual at all
  399 exact event centers.  The maximum production upper bound is 2.150e-6
  normalized and 9.677e-7 physical; eight higher-order nonmatching
  reconstructions pass.  Continuous height in all event cells remains open

jensen_window_pf_newman_c1_hardy_lower_fold_event_parameter_second_order_all_event_continuous_height_gate:
  exact transport by exp(i h z/beta), seven signed contour moments, an
  explicit degree-six exponential remainder, and transported analytic tails
  certify every |h|<=pi/16 cell at all 399 exact lower-fold events.  Maximum
  bounds are 2.153e-6 normalized and 9.690e-7 physical, both at event 2;
  eight higher-order moment and direct-endpoint reconstructions pass.  The
  ordinary-mode join and complete T_upper remain open

jensen_window_pf_newman_c1_hardy_lower_fold_ordinary_mode_coverage_ledger_gate:
  exact target-mode and selector-height ownership ledger.  The source target
  is 622..39894, the event-atlas overlap is 39695..39852, and a disjoint
  projection is ordinary 622..39694 plus fold-owned 39695..39894.  The 399
  event buffers leave 400 open corridors, so this is not a quantitative join

jensen_window_pf_newman_c1_hardy_ordinary_full_saddle_interface_barrier_gate:
  exact grouped-current witnesses reject the bare full-line carrier as a
  termwise overlap object at both the proposed and structural interfaces.
  The result is not a lower bound for the complete outer integral

jensen_window_pf_newman_c1_hardy_ordinary_frozen_endpoint_current_aggregate_barrier_gate:
  independent 125-digit summation over all 39273 source modes leaves
  two-real frozen-current discrepancy -0.0962793319571....  This rejects
  freezing the retained endpoint current, not the non-frozen Morse method

jensen_window_pf_newman_c1_hardy_ordinary_global_morse_gamma_bulk_gate:
  exact Abel-regularized Gamma evaluation of the full-line ordinary bulk.
  Its common factor supplies the saved-height Riemann-Siegel phase correction
  and leaves aggregate two-real error below 2.186e-19 over all target modes.
  Finite endpoint tails remain separate

jensen_window_pf_newman_c1_hardy_ordinary_finite_endpoint_tail_decomposition_gate:
  exact split into full-line bulk, lower paired endpoint/Fresnel tail, and
  upper paired endpoint/Fresnel tail.  Odd endpoint parity proves that bare
  endpoint exponentials are not independently summable; complete symmetric
  Poisson reassembly remains open

jensen_window_pf_newman_c1_hardy_ordinary_symmetric_endpoint_tail_reassembly_gate:
  exact Dirichlet-kernel reassembly of the bulk-extracted residual.  The
  endpoint half-current and every zero, negative, low, target, and outer mode
  remain in one symmetric limit, where the one-over-m current cancels inside
  each signed pair.  A quantitative residual bound remains open

jensen_window_pf_newman_c1_hardy_half_kummer_reflection_branch_reduction_gate:
  exact odd-alpha reflection reduces the full source to twice the real
  half-domain integral.  At the saved height only modes 1..39894 have a
  positive half-domain stationary branch; the upper branch is already its
  conjugate.  The half-domain complement remains open

jensen_window_pf_newman_c1_hardy_half_kummer_B_crossing_tangent_profile_barrier_gate:
  the complete grouped tangent B profile has two-real aggregate
  2.17818022827e-5, above target.  Omitting its q-density creates a false
  pass.  This rejects the isolated tangent model, not the exact nonlinear
  pair residual

jensen_window_pf_newman_c1_hardy_half_boundary_linear_roster_pair_anchor_gate:
  exact midpoint audit: odd-square parity simplifies the direct integer
  roster but not its canonical real-u interpolation.  Explicit nonzero pair
  coefficients reject midpoint pair cancellation; the parity identity is
  retained only as boundary data

jensen_window_pf_newman_c1_hardy_half_paired_target_residual_normal_form_gate:
  exact half-domain source-minus-target partition into endpoint-half, zero,
  low, ordinary, fold-owned, and outer signed-pair blocks.  Every negative
  mode and target carrier occurs once, without a false modewise half-Gamma
  interpretation.  Quantitative pair-block bounds remain open

jensen_window_pf_newman_c1_hardy_half_signed_pair_logistic_morse_transform_gate:
  exact one-phase half-domain Morse integral for each +m/-m pair.  Signed
  Fresnel coordinates differ but their outer m-squared phase is common; the
  half-boundary limit crosses the saddle exactly at 39894|39895.  Pair
  amplitude bounds remain open

jensen_window_pf_newman_c1_hardy_half_pair_endpoint_transport_volterra_gate:
  exact roster PDE and induced first-order x transport ODE.  Its zero initial
  pair gives an endpoint-driven Volterra solution and triangular half-Kummer
  representation that retains endpoint/Fresnel cancellation.  A quantitative
  Volterra operator theorem remains open

jensen_window_pf_newman_c1_hardy_half_volterra_incomplete_gamma_outer_kernel_gate:
  exact path-defined incomplete-Gamma kernel.  Monotone nonstationary phase
  proves the complete m>=39895 kernel l1 norm is below 0.000103 at the saved
  height.  The oscillatory endpoint-driver convolution is not yet bounded

jensen_window_pf_newman_c1_hardy_half_pair_separable_triangle_geometry_gate:
  exact separable pair integral over 0<z<x<1/2.  The B, A, and half-boundary
  events are joint-saddle crossings of one triangle's faces, unifying the
  former endpoint and fold walls.  A uniform interior/face/corner theorem
  remains open

jensen_window_pf_newman_c1_hardy_half_pair_exact_bimorse_face_fold_gate:
  exact global endpoint and outer Morse coordinates reduce the complete
  phase to (P^2-S^2)/2.  The B face is noncharacteristic, while the A event
  and half-boundary corner meet the same exact null-face condition.  Uniform
  transformed-amplitude bounds remain open

jensen_window_pf_newman_c1_hardy_half_pair_face_trace_lattice_deghosting_gate:
  exact face trace and saved-height interval roster.  Only 621/622 are local
  to the B face; mode 257 is separated by an exact normal derivative above
  1.70e13 in magnitude.  The A lower local pair is 39852/39853

jensen_window_pf_newman_c1_hardy_half_positive_B_crossing_step_tail_normal_form_gate:
  exact sign-adapted B tail plus bulk-step identity.  Its artificial jumps
  cancel, and the first paired boundary current has a cotangent lattice sum
  after local poles are removed.  The pole-subtracted quantitative bound
  remains open

jensen_window_pf_newman_c1_hardy_half_pair_abel_theta_modular_dual_roster_gate:
  exact Abel-theta sum-first triangle and Jacobi dual phase.  The B dual
  indices recover all 2,481,423 source alphas in reverse.  This removes
  modewise face poles but is involutive and does not itself prove a gain

jensen_window_pf_newman_c1_hardy_half_B_face_gaussian_window_fresnel_remainder_gate:
  the natural |xi|<=70 B window contains exactly crossings 621/622.  Every
  other positive mode has |q_B|>75, and the complete nonlocal third-order
  Fresnel remainder, including its analytic infinite tail, is below 7.48e-6.
  The two exact modes and three rational lattice currents remain open

jensen_window_pf_newman_c1_hardy_half_B_face_pole_subtracted_cotangent_current_gate:
  exact removal of the 621/622 poles from the first paired cotangent current.
  The analytic background is strictly decreasing and has one zero at
  xi approximately -6.567.  Its raw absolute integral is too costly to
  discard, so it must remain signed with the exact two-mode replacement

jensen_window_pf_newman_c1_hardy_half_B_face_discrete_joint_gradient_gap_gate:
  the B-face tangential and normal coordinates are monotone, and the integer
  lattice misses their simultaneous zero by more than 22.40 natural units.
  This supplies a denominator above 501 for a smooth two-direction
  integration-by-parts partition; its derivative constants remain open

jensen_window_pf_newman_c1_hardy_half_B_local_exact_quadratic_domain_gap_gate:
  the exact local phase is -y^2+pi*u^2/2.  The curved residual domains for
  modes 621 and 622 exclude the joint Gaussian saddle and have minimum
  largest gradient components above 57.17 and 45.63.  Transformed-amplitude
  derivatives and the resulting local integral bounds remain open

jensen_window_pf_newman_c1_hardy_half_pair_abel_theta_dual_continuation_stationary_cancellation_gate:
  matched B/A continuation labels share the same reduced phase, and their
  universal leading z-stationary channel cancels exactly.  Independent
  rational pi bounds prove A-2<sqrt(8t/pi)<A, leaving k=0 as the sole
  face/fold exception; a uniform full stationary remainder is not claimed

jensen_window_pf_newman_c1_hardy_half_pair_abel_theta_dual_weber_completion_cancellation_gate:
  an exact Gaussian coordinate evaluates the completed kernel by a Weber
  integral and makes it independent of endpoint D for every fixed matched
  label.  The completed pair cancels, while honest common-cutoff bookkeeping
  exposes an L-term A shift defect instead of silently reindexing infinity

jensen_window_pf_newman_c1_hardy_half_pair_abel_theta_dual_cutoff_shift_defect_vanishing_gate:
  scaling z=xq and retaining the x factor gives an explicit high-index
  integration-by-parts bound.  The weighted mode norm is O(n^-1/2) plus
  O(A^2n^-2), so the fixed-width L-term shift defect vanishes at the saved
  height.  No height-uniform version is claimed

jensen_window_pf_newman_c1_hardy_half_pair_abel_theta_dual_weber_source_roster_tail_reassembly_gate:
  the finite completed B block has exactly the equation-(9) normalization
  and recovers labels B-2 through A.  All remaining positive dual terms form
  one same-index endpoint-difference theta tail whose driver vanishes on its
  new face.  Its separate zero/direct sector is not yet closed in this gate

jensen_window_pf_newman_c1_hardy_half_pair_abel_theta_dual_double_weber_exact_source_reconstruction_gate:
  the second Jacobi-Weber transform cancels every nonzero endpoint difference
  by odd endpoint parity and its zero mode by an Abel primitive.  The theta
  tail, dual zero, direct subtraction, Poisson zero, and endpoint halves then
  reconstruct the complete finite source roster.  This is an exact
  involution guard, not a small source-minus-target estimate

jensen_window_pf_newman_c1_hardy_half_B_local_normal_safe_fresnel_remainder_gate:
  exact rational split points and monotone distance panels bound the
  normal-safe positive B-tail remainders for modes 621 and 622 by a combined
  1.97e-10 after the equation-(9) projection.  The complementary local
  currents and completion sectors remain coupled and open

jensen_window_pf_newman_c1_hardy_half_B_face_pole_subtracted_cotangent_signed_window_gate:
  stable removable-pole charts and signed integration place the first
  cotangent B-window current below 4.792e-6 in the physical real projection.
  Its complex magnitude is about 2.07024e-5, so a modulus-first budget would
  erase the observed cancellation.  This is one current, not a complete B sum

jensen_window_pf_newman_c1_hardy_half_B_face_higher_boundary_current_absolute_window_gate:
  the next two pole-subtracted boundary currents reduce to lattice sums R_2
  through R_5.  Finite interval sums and analytic tails bound their physical
  absolute window contribution by 1.491e-6.  The post-third remainder, local
  replacements, and outside-window tails remain open

jensen_window_pf_newman_c1_hardy_equation9_endpoint_residual_source_hybrid_nonidentification_guard_gate:
  the exact endpoint object Q_K, the published hybrid Q_P, and implemented
  source output Q_S are distinct.  The saved -0.006983... discrepancy cannot
  be transferred to Q_K-T until both bridge defects are bounded.  The exact
  Gamma bulk changes Q_K-T by less than 2.186e-19, but does not close it

jensen_window_pf_newman_c1_hardy_lower_fold_beta4_canonical_airy_derivative_reduction_gate:
  the complete beta-minus-four corrected fold model reduces exactly, in the
  symmetric Abel sense, to one y integral with only Ai(-lambda-y) and its X
  derivative.  Four explicit rational weights replace every z moment through
  degree ten.  The Airy-to-ordinary carrier remainder theorem remains open

jensen_window_pf_newman_c1_hardy_lower_fold_beta4_canonical_exact_hankel_branch_factorization_gate:
  for positive X the two corrected Airy channels split exactly into conjugate
  Hankel branches of orders 1/3 and 2/3, with forced pi/6 connection phases.
  The branch basis gives the ordinary saddle phases; its prefactored branches
  have finite limits at X=0, where the written Hankel products require the
  regular Airy basis or explicit limits.  The quantitative splice is open

jensen_window_pf_newman_c1_hardy_lower_fold_beta4_hankel_carrier_event_selection_gate:
  every certified event cell has X>=pi/(16beta)>9e-5.  The exact extracted
  cubic carrier selects H^(1) for positive d_m and H^(2) for negative d_m,
  crosses at the event center, and leaves the opposite carrier with endpoint
  slope above 0.023.  Scaled-Hankel amplitude derivatives remain open

jensen_window_pf_newman_c1_hardy_lower_fold_beta4_scaled_airy_branch_amplitude_envelope_gate:
  after exact extraction of the cubic carrier, the complete beta-minus-four
  Airy branch amplitude and its first y derivative satisfy |S|<0.356 and
  |S_y|<0.258 on all 399 event cells.  Arb atlases and the signed DLMF
  derivative-modulus remainder certify the envelopes; the logistic bridge is open

jensen_window_pf_newman_c1_hardy_lower_fold_beta4_opposite_branch_roster_reflection_pairing_gate:
  the old lower-edge, shared-boundary Fourier, and adjacent lower-edge
  399-mode rosters each reduce exactly to 199 conjugate reflection pairs plus
  one edge mode.  The shared-boundary block uses old-strip labels -793..799,
  not adjacent-chart labels -795..797.  This organizes the already-completed
  Poisson object and does not replace its complement or endpoint half-current

signed_hankel_jensen_dependency_graph:
  dependency hygiene gate connecting finite evidence, countermodel gates, open
  theorem targets, and the `lambda_le_0_goal` node; validates that finite and
  diagnostic nodes only support open targets and have no direct proving edge
  to the not-proved conclusion

target_jensen_window_pf_bridge:
  theorem target that reformulates all-degree/all-shift Jensen hyperbolicity
  as finite PF-infinity of every binomially weighted window
  B^{d,n,0}_j=binom(d,j) A_{n+j}(0); open and not proved

countermodel_gates:
  validates 11 proof-safety examples, including local heat birth, finite
  coefficient-prefix promotion, finite Jensen-window rectangle promotion,
  finite Schur-prefix promotion, finite signed-Hankel grid promotion, finite
  moment-recurrence promotion, and Stieltjes/Hankel-to-Toeplitz promotion traps
```

## Boundary

Passing the ledger checker means the corpus has not confused finite
certificates, diagnostics, or algebraic translations with the missing bridge
theorem. It does not prove any open target.
