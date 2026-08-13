# Pole-subtracted cotangent current on the B window

Date: 2026-08-13

Status: exact meromorphic reduction plus interval certificate; not a proof of
the complete B-face estimate

For `c=Bx/2`, the paired first boundary current from Section 11.362 contains

```text
S_all(c)=sum_(m>=1)1/[4(c^2-m^2)]
        =pi*cot(pi*c)/(8c)-1/(8c^2).                 (PC1)
```

Remove the two exact local modes before estimating:

```text
S_hat(c)=S_all(c)-1/[4(c^2-621^2)]
                  -1/[4(c^2-622^2)].                (PC2)
```

The apparent singularities at 621 and 622 are removable.  In general,

```text
lim_(c->k){S_all(c)-1/[4(c^2-k^2)]}=-3/(16k^2),    (PC3)
```

so the second local term can simply be subtracted at each limit.  The
certified values have opposite signs:

```text
S_hat(621)=[0.0002006401038312234182812424536759239998873228301863990432429055236415183477691725940309100313040884758 +/- 3.79e-104],
S_hat(622)=[-0.0002016109487122156260000326223830321420388205859014927302117361509978153067897651190917023588004155869 +/- 3.80e-104].         (PC4)
```

Every retained summand satisfies

```text
d/dc {1/[4(c^2-m^2)]}=-c/[2(c^2-m^2)^2]<0.         (PC5)
```

Thus `S_hat` is strictly decreasing on the entire `|xi|<=70` B window.  Its
endpoint and trace values are

```text
S_hat(c_low)=[0.0002287566947123639952553857976628996750018483712014892160637408938285838743517519866931881390746364378 +/- 1.60e-96],
S_hat(c_B)  =[-2.172069043349191876854041778990298768995550684588141779948687807995543444418218036693214559184788890e-5 +/- 5.15e-98],
S_hat(c_high)=[-0.0002876512597187688412166760086876603807571715240745292081569167684594467730083398761789170711908660455 +/- 2.02e-97],                (PC6)
```

and `|S_hat|<0.000288` throughout.  There is exactly one zero, with

```text
c_* in [621.49827876663, 621.49827876665],
xi_* in ['[-6.567007143028825808100663896987657910485049276244088511671553400001965446920898900527238225416705001 +/- 1.57e-93]', '[-6.567007140753546488949481559310762051760043866683290092645188876896354812530400227851542161297611015 +/- 1.57e-93]'].                 (PC7)
```

The corresponding analytic first boundary current is

```text
C_hat(x)=-2(2+i*pi*B^2*x)S_hat(Bx/2)/pi^2.           (PC8)
```

Its absolute value at the trace saddle exceeds `88000`.  Since `S_hat` is
negative and decreasing for `0<=xi<=1`, even that one-unit strip has raw
absolute integral above

```text
[0.0003021963971352309997844190771000589891374309201815937335534283217595626059962724763752441306153849558 +/- 7.16e-97] > 0.000302. (PC9)
```

Equation (PC9) is a cancellation guard, not a lower bound for the signed
residual.  It proves that the analytic cotangent background cannot be
triangled away: its oscillatory integral must be combined with the exact
621/622 currents and the endpoint/negative/outer completion.

Pi provenance: `pi` comes from the equation-(9) face current and the integer
cotangent partial fraction.  No fitted constant is used.

Proof boundary: exact pole removal, monotonicity, unique-zero enclosure, and
saved-height size bounds for the first analytic boundary current only.  The
exact local currents, the other rational terms, outside-window estimate, and
A-fold splice remain open.  This result does not establish the complete
paired residual, `T_upper`, `Lambda<=0`, PF-infinity, RH, or a prize-level
conclusion.
