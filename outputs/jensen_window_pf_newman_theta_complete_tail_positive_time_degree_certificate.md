# Complete-Theta Local Positive-Time Degree Certificate

Date: 2026-08-04

Status: rigorous complete-theta local contact exclusion. This is not a proof of RH.

```text
work/rh_compute/results/jensen_window_pf_newman_theta_complete_tail_positive_time_degree_certificate.json
python work/rh_compute/scripts/jensen_window_pf_newman_theta_complete_tail_positive_time_degree_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_theta_complete_tail_positive_time_degree_certificate.py
```

## Theorem

The complete theta-kernel heat transform and its first x derivative have no
common zero in

```text
0<=t<=0.45,  135.5<=x<=136.
```

The certified finite boundary separation is `[1.946581008278572725404437041561276623575684854292304814661743878009124e-21 +/- 4.28e-91]`.
The complete n>=6 tail budget is `[6.085496135366523828070687011210104660166492461570447146370589647700371e-42 +/- 1.54e-112]`.
Their ratio is at most `[3.126248591497424204626508401912865241041842156376065885300607362926303e-21 +/- 3.87e-91]`.

The boundary homotopy from the five-term endpoint to the complete theta sum
therefore stays nonzero and preserves winding zero. The all-multiplicity
same-sign index theorem then excludes every interior contact.

## Proof Boundary

This is a rigorous complete-theta contact exclusion only in the declared local positive-time rectangle. It does not control other frequency rectangles, negative heat time, all real x, the global de Bruijn-Newman constant, a degree-uniform Jensen remainder, Lambda<=0, RH, or a prize-level conclusion.
