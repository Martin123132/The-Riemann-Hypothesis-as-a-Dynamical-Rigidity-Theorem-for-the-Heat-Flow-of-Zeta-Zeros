# Short Small-Prime Heat-Starlikeness Counter-Gate

Date: 2026-07-30

Status: exact countermodel gate. This is not a proof of a joined
Xi current theorem, `Lambda<=0`, PF-infinity, or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_prime_power_heat_starlikeness_short_family_counter_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_prime_power_heat_starlikeness_short_family_counter_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_prime_power_heat_starlikeness_short_family_counter_gate.py
```

## Outcome

The single linked interior point

```text
q=s=9999/10000
u=log(p)/2
t=-4log(q)/(log(p))^2>0
```

lies inside both certified small-prime boxes. At the rational
angular coordinates below, exact arithmetic proves `J_0<-1/100`.

| p | M | x=cos(theta) | exact J_0 at q=s=1 | J_0 at q=s=9999/10000 |
|---:|---:|---:|---|---:|
| 2 | 2 | -99/100 | `2/1 + (-297/200)sqrt(2)` | -1.00087079701092639481475128605388921e-1 |
| 2 | 3 | -1/2 | `7/4 + (-11/8)sqrt(2)` | -1.94335337195652692370441320700010075e-1 |
| 2 | 4 | 0/1 | `-1/4 + (0/1)sqrt(2)` | -2.49199730622557028090992396830896887e-1 |
| 2 | 5 | 1/3 | `163/432 + (-179/432)sqrt(2)` | -2.06718954934978299055695752829033371e-1 |
| 2 | 6 | -1/2 | `-1/16 + (-1/64)sqrt(2)` | -8.34297779656366846629422844380184757e-2 |
| 2 | 7 | -1/4 | `167/256 + (-481/1024)sqrt(2)` | -1.01636271767902379230144020337109703e-2 |
| 3 | 2 | -99/100 | `5/3 + (-99/100)sqrt(3)` | -4.79873305832594778199386496369142615e-2 |
| 3 | 3 | -1/2 | `4/3 + (-7/9)sqrt(3)` | -1.35230545606210568936882293540893020e-2 |

Every displayed decimal is secondary. The checker reconstructs
the exact pair in `Q(sqrt(p))`, shifts it by `1/100`, and
proves negativity by rational signs or an exact square
comparison.

## Consequence

```text
p=2,M=2..7: universal blockwise J_0>0 is false
p=3,M=2..3: universal blockwise J_0>0 is false
```

Together with the positive base and all-length gates, this
shows that `M=8` for p=2 and `M=4` for p=3 are sharp
thresholds for the present blockwise heat-current theorem.
The short rows are no longer an unfinished Bernstein search:
the desired uniform statement is false.

The strict witnesses are interior in q, s, and x, so
polynomial continuity gives nonempty negative neighborhoods.
They are also on the exact q-s heat linkage. This still does
not assert that a physical Xi ray attains any listed angular
phase.

## Route

Short chains must be retained inside the fully rejoined
endpoint current with their p-free base phases, singleton
chains, recurrent endpoint, and cutoff terms. A sum of
independent blockwise positivity or winding claims cannot
close that joined object.

The angular variable uses the ordinary `2*pi` period of
`exp(i theta)`. Any `pi` in the saddle scale remains inherited
from completed-zeta normalization.

## Boundary

The exact witnesses reject universal blockwise heat-starlikeness for p=2,M=2..7 and p=3,M=2..3. They do not assert that an Xi ray attains a witness phase, and they do not prove a joined p-free, singleton, endpoint, cutoff, Abel-gap, or winding theorem, contact exclusion, Lambda<=0, PF-infinity, RH, or a prize-level result.
