# hka_pert.py — Stage-B perturbation operator + the GAUGE-MODE VALIDATION test (the honest gate).
#
# STATUS (see RESULTS_hka_beta.md): Stage B is NOT landed. The perturbation operator L(x;kappa) as
# assembled from the transcribed HKA coefficients (5.5)-(5.13) FAILS the rigorous gauge-mode exactness
# test, so no eigenvalue (and no beta) is reported from it. This module ships the operator, the
# validated pieces, and the falsification test, so the wall is reproducible and precisely located.
#
# THE VALIDATION GATE (a correct perturbation operator MUST pass it):
#   The pure-gauge mode (HKA 5.20) Psi_g(x;kbar) = ((lnA)'_x, (lnN)'_x + kbar, (lnom)'_x, V'_x) is an
#   EXACT solution of Psi' = L(x;kbar) Psi for EVERY kbar (it is pure coordinate freedom). We test
#   |d/dx Psi_g - L(x;kbar) Psi_g|.  For the correct L this is ~0.
#
# WHAT IS VERIFIED:
#   - The auxiliary identity (5.14) holds EXACTLY on Psi_g (residual ~1e-11) -> the Abar row and the
#     gauge-mode form are correct.  [hka_check_514.py]
#   - The reduced background flow Jacobian J3 (constraint surface) reproduces the kappa=0 gauge mode
#     (flow tangent) to ~1e-9.  [hka_pert_reduced.py]
# WHAT IS NOT CLOSED:
#   - The kappa-coupling of the Nbar and (ombar,V) rows: assembling them from (5.5)-(5.10)+(5.13) as
#     written does NOT make Psi_g exact for kbar != (special values). The gauge-mode constraints
#     (K e_N = 0 and K D = -L0 e_N) UNDER-determine the s-coupling; the transcribed coefficients do
#     not supply the missing piece. So the operator is not trusted, and no kappa is extracted.  See
#     hka_pert_symbols2.py / hka_pert_jac2.py / hka_pert_red3.py for the systematic falsification.
#
# Eq #s per HKA_beta_equations.md.
import numpy as np, math
import hka_pert_core as PC       # the (5.5)-(5.13) assembled operator L(x;kappa)=P+kappa Q
import hka_ec as E

G=4.0/3.0

def gauge_mode(x, kbar, bg_state):
    """HKA (5.20) pure-gauge perturbation vector at x, eigenvalue kbar."""
    A,N,om,V,ombar_x,Vx=bg_state(x); oV2=1-V*V
    Ax=A*(1-A+2*om*(1+(1/3.0)*V*V)/oV2); Nx=N*(-2+A-(2-G)*om)
    return np.array([Ax/A, Nx/N+kbar, ombar_x, Vx],complex)

def gauge_mode_residual(kbar, bg_state, xs, xlist=None):
    """|d/dx Psi_g - L Psi_g| (relative) at several x. ~0 for a correct operator."""
    if xlist is None: xlist=[-8.0,-4.0,-1.0,xs-0.05]
    out=[]
    for x in xlist:
        Psi=gauge_mode(x,kbar,bg_state)
        dP=(gauge_mode(x+1e-6,kbar,bg_state)-gauge_mode(x-1e-6,kbar,bg_state))/2e-6
        L=PC.Lnum(bg_state(x),complex(kbar))
        res=dP-L.dot(Psi)
        out.append((x, np.linalg.norm(res)/max(np.linalg.norm(dP),1e-12)))
    return out

if __name__=="__main__":
    import hka_beta4 as B; B.bg(); B.bg_path()
    xs=B.bg()['xs']
    print("GAUGE-MODE VALIDATION GATE for the (5.5)-(5.13)-assembled operator (hka_pert_core.Lnum):")
    print("(a correct perturbation operator gives residual ~0 for ALL kbar; large residual => operator wrong)")
    for kbar in [0.35699, 1.0, 2.81055255]:
        res=gauge_mode_residual(kbar, B.bg_state, xs)
        print(f"  kbar={kbar}: rel residuals = {[f'{r:.1e}' for _,r in res]}")
    print("\n=> operator FAILS the gate (residuals O(1)-O(1e4)); Stage B blocked, no beta reported.")
    print("   (The identity 5.14 DOES hold on the gauge mode -> Abar row + gauge-mode form verified;")
    print("    the Nbar/(ombar,V) kappa-coupling is the unresolved piece. See RESULTS_hka_beta.md.)")
