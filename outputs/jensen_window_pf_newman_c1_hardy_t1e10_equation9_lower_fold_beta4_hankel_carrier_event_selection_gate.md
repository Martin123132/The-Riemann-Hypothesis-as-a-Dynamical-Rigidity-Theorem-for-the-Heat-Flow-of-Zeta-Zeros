# Beta^-4 Hankel carrier: exact event selection

Date: 2026-08-13
Status: exact carrier-selection gate; not a proof of the amplitude splice

Factor the positive-`X` Hankel branches from Section 11.378 as

```text
C_+(X)=e^(i[xi-pi/4]) B_+(X),
C_-(X)=e^(-i[xi-pi/4])B_-(X),       xi=2X^(3/2)/3,    (ES1)
```

where `B_+` contains the exact scaled Hankel amplitudes of orders `1/3` and
`2/3`, and `B_-=conj(B_+)` for real data.  This factorization is exact: no
large-`X` replacement is made.  Restoring the normal carrier
`exp(i[y^2/(4beta)-d_m*y])` gives extracted phase derivatives

```text
F_+(y)=y/(2beta)-d_m+sqrt(lambda+y),
F_-(y)=y/(2beta)-d_m-sqrt(lambda+y).                  (ES2)
```

At event `tau_m`, write `t=tau_m+delta`.  The exact event identity is

```text
lambda=d_m^2-delta/beta.                              (ES3)
```

Set `sigma=sign(d_m)`, choosing `C_+` for positive `d_m` and `C_-` for
negative `d_m`.  Rationalization gives the exact endpoint orientation

```text
F_sigma(0)
 =-d_m+sigma*sqrt(lambda)
 =-sigma*delta/[beta(sqrt(lambda)+|d_m|)],             (ES4)

F_-sigma(0)=-sigma(|d_m|+sqrt(lambda)).                (ES5)
```

Thus the selected carrier crosses the lower endpoint exactly once, at the
event center, while the opposite carrier does not approach that endpoint
saddle.

For the 399 cells, `|4m-C|` runs through the nonzero odd labels of absolute
size `1,...,797`.  Since `beta^3=pi*C^2/8` and `|delta|<=pi/16`,

```text
lambda >= (beta/C)^2-pi/(16beta)
        = pi/(16beta)
        = [9.113735018465280815095579752470019118589176704071891157128083836080253086767956786676652574538943211e-5 +/- 1.18e-103] >9e-5.                (ES6)
```

Hence every point of every certified event cell and every `0<=y<=64` lies
strictly in `X=lambda+y>0`.  The entire fold strip also obeys

```text
X<[179.7826212842364358783291133355320921861894327921030789194750313894923904623866089297296628709318970 +/- 8.67e-98]<180<beta^2.             (ES7)
```

Consequently the opposite extracted carrier has fixed derivative sign on
the fold strip, with endpoint magnitude greater than
`[0.02304750320303096570934008097507117659999729076109874125176871511415862386522358812650639225263348796 +/- 1.10e-101]`.  The selected carrier is
strictly monotone in `y` and has at most one carrier saddle.  Its endpoint
slope at either cell face has magnitude below
`[0.003954326391965416039225027083227013015504189783863254265614367474068863631079466840803494222975977140 +/- 3.07e-102]`.

Equations (ES1)--(ES7) select the correct cubic carrier and orient every
event cell.  They do not bound derivatives of the exact scaled Hankel
amplitudes `B_+/-`; therefore they do not yet license integration by parts on
the opposite branch or identify the selected branch with the exact
logistic/Gamma carrier.  No corridor remainder, `Q_K-T`, complete `T_upper`,
height-uniform theorem, `Lambda<=0`, PF-infinity, RH, or prize-level result is
proved.
