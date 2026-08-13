# Uniform completed cutoff-family Dirichlet bound

Date: 2026-08-13

Status: rigorous ordinary-top-corridor bound for every beta-minus-four
completed cutoff; exact finite-integral remainder and final splice remain open

For `L=39696`, `U=40094`, and `L<=k<=U`, Section 11.393 defines

```text
S_k=sum_(m=L)^k G4_m,
R_k=H4+C4+sum_(m=k+1)^U G4_m
   =g4_lambda(0)-S_k.                                (CD1)
```

The coefficient formula from the completed projection gives

```text
G4_(39895+n)=Integral_0^1 g4_lambda(s)e^(-2pi i n s)ds.
```

Hence a prefix of length `N=k-L+1` is one contiguous Fourier packet,

```text
S_k=Integral_0^1 g4_lambda(s)e^(2pi i 199s)D_N(s)ds,
D_N(s)=sum_(j=0)^(N-1)e^(-2pi i j s).                (CD2)
```

No coefficient is bounded separately.  For `0<=s<=1/2`,
`sin(pi s)>=2s` and `|sin(pi Ns)|<=min(pi Ns,1)`.  Splitting at
`s=1/(pi N)` and using symmetry proves

```text
Integral_0^1 |D_N(s)|ds <=1+log(pi N/2)
 <=1+log(399pi/2)
 =[7.440544122179318321235820642148441115093749227258795776124987781963111005950006887553192580496414532 +/- 3.49e-100]<7.441.              (CD3)
```

On the ordinary half of the selector top cell,

```text
0<=lambda<=pi/(16beta),
0<=y<=Y=pi A/(2beta),
0<=X=lambda+y<117.                                   (CD4)
```

The saved Airy modulus atlas covers the range above its first positive
point; one direct interval closes the remaining interval down to `X=0`.
Triangle evaluation of the exact beta-minus-four polynomials gives

```text
|U|<[1.000002490295157223407932902122253432662501450145145550773236756972637430163650659991613289369000432 +/- 8.06e-101],
|V|<[0.0001458203042186848072503435467396129693445750274697412422140751583552388872052681861127070812782193710 +/- 3.43e-104],
sup_s |g4_lambda(s)|
 <[82.64569549258764647268742038077922092391217238921981714578486483549692205462395822010653975912092037 +/- 6.44e-99]<82.66.                         (CD5)
```

Combining (CD1)--(CD5) proves simultaneously for all 399 cutoffs

```text
max(|g4_lambda(0)|,max_k|R_k|)
 <[697.5746393133824413779534853310194797586451439776906343413574545001966529628299427057473584106925355 +/- 6.25e-98]<698.              (CD6)
```

Inserting the certified endpoint-augmented weight variation from Section
11.393 gives a cancellation-preserving uniform bound for the auxiliary
leading-defect-weighted full beta-minus-four roster:

```text
|W_full|=|sum_(m=L)^U w_m G4_m|
 <[0.02157723351885976688022908813615363079367348881460337244572615573875037000450061266826820973865474117 +/- 5.87e-102]<0.02159.                 (CD7)
```

This closes the full-roster partial-completion information deficit at a
conservative scale.  It is not yet a selected-branch estimate: exact branch
algebra gives `W_sel=W_full-W_opp`, and `W_opp=sum w_mG_opp,m` remains to be
bounded.  It also does not identify the exact finite-integral amplitude
remainder with `w_m` or splice (CD7) into the branch initial-data and source
coordinates of Section 11.392.  Sharpening (CD3) or exploiting the profile's
variation may reduce the constant, but is not needed for validity.

Pi provenance: `Y`, the height interval, and the Fourier kernel inherit the
exact identities `beta^3=pi A^2/8` and `hY=2pi`; (CD3) uses the standard
Fourier exponential `e^(-2pi i n s)`.  No fitted circle constant enters.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_completed_cutoff_family_uniform_dirichlet_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_completed_cutoff_family_uniform_dirichlet_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_completed_cutoff_family_uniform_dirichlet_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_completed_cutoff_family_uniform_dirichlet_gate.py
```

No weighted opposite-branch estimate, selected-branch weighted bound, exact
finite-integral amplitude remainder, completed branch initial-data or source
splice, complete `Q_K-T` or `T_upper`, all-corridor theorem,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
