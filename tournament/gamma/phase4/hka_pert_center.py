# hka_pert_center.py — center behaviour of the perturbation; build the regular center subspace.
#
# Near the center (x->-inf), L(x;kappa) -> L_c (constant, evaluated at the center fixed point). Its
# eigenvalues are the indicial exponents: {0,0,-1,-2} (seen numerically: modes 1,1,e^{-x},e^{-2x}).
# As x->-inf: e^{-x},e^{-2x} -> +inf (GROWING), the two lambda=0 modes are bounded. Regularity/BC
# (md line 102-105): the perturbation must be BOUNDED at the center, and the e^{-x} mode is excluded
# by identity (5.15). So a regular center perturbation lives in the span of the bounded modes; the
# eigenvalue kappa is fixed by ALSO requiring analyticity at the sonic point.
#
# We build the center launch: at a deep center x0, seed the perturbation in the bounded (regular)
# eigenspace of L_c, integrate outward. The dimension of the regular launch space determines the
# shoot. We examine it here.
import numpy as np, math
import hka_pert_core as PC
import hka_ec as E

def center_L(kappa, xc=-10.0, N_inf=1.0, oi=0.375):
    z=math.exp(xc); N=N_inf/z; om=oi*z*z; V=E.M_CENTER/N_inf*z
    fld=PC.bg_fields(N,om,V)
    return PC.Lnum(fld,complex(kappa)), fld

if __name__=="__main__":
    for kappa in [2.81055255, 1.0, 0.5+1j]:
        print(f"\nkappa={kappa}:")
        for xc in (-8,-10,-12):
            L,fld=center_L(kappa,xc)
            ev,V=np.linalg.eig(L)
            idx=np.argsort(ev.real)
            print(f"  xc={xc}: eig(L_c)={ev[idx].round(5)}")
            # bounded modes: Re(lambda)>=0 as x->-inf means e^{lambda x}->0 for Re<0? No:
            # as x->-inf, e^{lambda x} -> 0 if Re(lambda)>0, -> inf if Re(lambda)<0, const if 0.
            # So BOUNDED (regular) = Re(lambda) >= 0. The negative-Re modes blow up at the center.
            bounded=[i for i in range(4) if ev[i].real>-1e-6]
            print(f"        bounded (Re lambda>=0, regular at center) modes: {len(bounded)}  eigs={ev[bounded].round(4)}")
