# Newman C1 Calibrated Roster-Interior Phase-Root Scout

Date: 2026-08-05

Status: high-precision physical phase diagnostic validated; determinant observations remain open; not a proof of RH.

## Physical Row

On q=1 use t=1/(2L^2), x=4pi exp(L), s=(1-ix)/2, s_*=s+t alpha(s)/2, Omega=-Im(s_*), and a^2=exp(L)+1/(32L^2). The first complete cell above L=50 has N=floor(a)=72004899338.

Use the exact Polymath-15 M_t(s)=exp[t alpha(s)^2/4]M_0(s), d_1=1/(6s)+alpha'(s){t/4+t^2alpha(s)^2/8}, and eta=(M_t/|M_t|)(1+d_1)/|1+d_1|.

With phi_N=Omega log N and u_N=log(a/N), evaluate B=Re(eta)Im(eta exp(i phi_N))+(u_N/log N)Im(exp(i phi_N)). This equals the trigonometric centre factor without choosing an argument branch.

```text
N = 72004899338
L = 50.000000000042094782678209870093208056456958815668741152203319508061561578709709317330579871713691733251987056925971604559559956654841405091942530844118747807
a-N = 0.90138781886411846440343380712520520260367977241598784887543317475284747281750624060629837925669046552385399717862497678439477029376219689262798156474148437918
N+1-a = 0.098612181135881535596566192874794797396320227584012151124566825247152527182493759393701620743309534476146002821375023215605229706237803107372018435258515620822
B = -0.49000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000020539978679775180588245933659129132484283998
B+49/100 = -2.0539978679775180588245933659129132484283997717320806050591529471538846016217254238032975208077484206026890161846391655528961642666237668055360627938137922752e-115
```

## Construction

After Omega crosses K_N=2pi N^2, solve chi_0=F_N(Omega)/2+pi/8 at the congruence chi_0=pi/2 mod 2pi. In that safe arc solve theta_0=chi_0-phi_N first at 3pi/2 and then at pi modulo 2pi. The resulting endpoint values straddle B=-49/100, so bisection constructs a roster-interior numerical level point.

The early bracket value is `-0.00000000000000000000000075689513757685148634278087966409207415296766657913327485838672797460970945268485072968504607553487624610202857700628026341732544141617324256853248708745520895` and the negative bracket value is `-1.0000000000005007369371501879980698875244732538557328260609215514096062091361916072169478185720372133559245992155820181284594427828232445281277845111236531263`. The final L-bracket width is `7.8320878922458355204535684960564447882269904304994213968161746092885713710259200135748906576129497677623552130514498423839931054412963622537630607257372103738e-139`.

## Precision Ladder

The row was rebuilt from scratch at 90, 130, and 170 decimal digits. Consecutive L differences are:

- `90 -> 130`: `1.8722728359138705093173305798717136917308580354613e-60`
- `130 -> 170`: `2.3939515956583657745595599566548414050919425308441e-99`

## Adjacent Cells

The same branch-free construction was independently run at 170 digits in the next two complete cells:

- `N=72004899339`: `L=50.000000000060361016824637911076260756838903417107508224425985557003096079054273978075384048957710631930867688084609208006783182504833206670773191348960981097`, `B=-0.48999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999983271877824287526765308474116863016706078195`, `a-N=0.55901699437429876862945623262022762749958643814103732642249847582007786933901414228772315226971255162544664949156473311690011138591664790900305665463571071091`
- `N=72004899340`: `L=50.000000000097646558136123683070915607814660709703558232238377963326914568657663272825228675580727227053603474662147107709606008575250622231201525379002199738`, `B=-0.49000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000011900577431958983894545472282702980430408987`, `a-N=0.90138781886406267845907558608875852465757554596039599314082921486580284674663054228071986266727746402456940621851147994249871846881962835741319573184070977495`

This is a convergence diagnostic, not an interval enclosure.

## Determinant Input Gap

The following fields are still missing:

- `omega_E and dot(omega_E)`
- `retained observation p and dot(p)`
- `finite block observations f^(0),f^(1),f^(2) and their tangents`
- `independent near observation n and common-cutoff paired-remote observation c`
- `det(B_mat)=h^2P_ret and every determinant current increment`

The calibrated cell has N=72004899338. Direct enumeration of N carriers or an N by N reciprocal block table is not an admissible next algorithm. The missing evaluator must use the certified stationary diagonal, adjacent recurrence, and far paired/compressed representations.

## Next Action

Use this saved root as the immutable input row for a determinant-evaluator contract. First implement omega_E, p, and the diagonal-plus-adjacent observation through the exact physical coefficient and Morse/adjacent formulas; cross-check it against the near-first chart. Add the far sector only through a resumable paired or Type-I/II compression. Do not infer a sign theorem from this phase scout.

## Pi Provenance

Every pi here is inherited from the completed-zeta normalizer, T_0=2pi a^2, Q(v)=exp(i[v log(v/(2pi))-v-pi/4]), and the period of the complex exponential. No fitted circle or polygon supplies pi.

## Proof Boundary

This is reproducible high-precision diagnostic evidence for calibrated physical phase levels in three consecutive fixed q=1 roster interiors. It is not interval arithmetic and proves no uniqueness of any B root, no all-cell statement, determinant-row value, P_ret nonsingularity, geometric-current sign or bound, complete-current sign, all-q transport, Lambda<=0, PF-infinity, RH, or prize-level conclusion.

built Newman C1 calibrated roster-interior phase-root scout: 12 rows, 0 issues, 6 source audits, 3 precision ladders, 3 physical q=1 cells, 3 calibrated B=-49/100 roots, 3 fixed-roster margins, 0 determinant rows, 1 open evaluator row
