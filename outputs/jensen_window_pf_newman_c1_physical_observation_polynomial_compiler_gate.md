# Newman C1 Physical Observation Polynomial Compiler Gate

Date: 2026-08-05

Status: exact coefficient compiler and low-N interface validated; physical scalar calibration and signed determinant remain open; not a proof of RH.

## Compiler Contract

The compiler maps the twelve explicit scalar inputs into ascending coefficient arrays for P_V=C, P_N=mathcal N, P_A=RC, and P_Q=Q. It derives u_N=log_a-L_N internally, preventing an inconsistent duplicate input.

Each base row is constructed both from the five moment vectors v_0,v_1,v_2,e_0,e_1 of Section 11.159 and from the carrier laws C,D,R,R_x,Q,mathcal N of Section 11.160. All twenty symbolic coefficient differences vanish.

Translation ell=log_a+x exactly reproduces R=-s_prime*x-c*u_N, R_x=-s_second*x-c_x*u_N+i*b*u_N_x, Q=(R+delta)C+D, and mathcal N=RQ+R_xC with delta=chi_rate-i*b*u_N.

The tangent compiler applies dot(P)(ell)=i(ell-log_a)P(ell), giving the exact ascending recurrence dot(p)_j=i[p_(j-1)-log_a*p_j].

## Physical Ownership

For alpha=(N_(p,xi),V_(p,xi),-Q_(p,xi),-A_(p,xi)), beta_p=(N_p,V_p,-Q_p,-A_p), and beta_E=(N_E,V_E,-Q_E,-A_E), P_lin=F_alpha+i(ell-log_a)F_(beta_p+beta_E) splits exactly into one bulk term plus the single edge lift. The edge is not inserted into alpha or counted twice.

## Validation

- Symbolic row identities: `4`.
- Centered translation identities: `4`.
- Tangent identities: `4`.
- Exact fixture degree: `5`.
- Low-N fixtures: `3`.
- Largest low-N lift projection discrepancy: `3.340e-15`.
- Largest low-N joined projection discrepancy: `3.553e-15`.

## Next Input Layer

A physical numerical row now requires explicit certified values of chi_rate, rho_1,rho_2 and their x-derivatives, c,b,c_x,b_x,u_(N,x), plus the genuine endpoint edge observations and retained carrier observations used in alpha and beta. The compiler deliberately does not invent them.

The compiler calls the complex coefficient rate `chi_rate`. This is intentionally distinct from the later phase angle that reuses the printed symbol `chi_N` elsewhere in the corpus.

## Pi Provenance

No new pi is introduced by the coefficient algebra. Any pi inside u_(N,x), calibrated rates, or the low-N Fourier basis is inherited from the completed-zeta/Fourier normalizations already audited by the source gates.

## Proof Boundary

This gate proves exact symbolic compatibility, coefficient ownership, degree ceilings, and a low-N software interface. Its low-N rates are synthetic regression data. It does not evaluate the physical chi_rate or endpoint jets, enumerate the physical carrier, prove interval enclosures, a determinant sign, a flow inequality, Lambda<=0, RH, or a prize-level conclusion.

built C1 physical-observation polynomial compiler gate: 4 base rows from 2 exact derivations, 4 centered translations, 4 tangent lifts, single edge ownership, 3 low-N interface fixtures, 0 physical numerical rows and 0 signed bounds
