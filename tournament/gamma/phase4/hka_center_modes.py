# Determine the center regular subspace precisely: eigenvectors of L_c and which satisfy Abar_p=0
# and which is the e^{-x} mode excluded by (5.15). Also confirm the e^{-2x} growing mode.
import numpy as np, math
import hka_pert_core as PC
import hka_ec as E

def L_center(kappa, xc=-14.0, N_inf=1.0, oi=0.375):
    z=math.exp(xc); N=N_inf/z; om=oi*z*z; V=E.M_CENTER/N_inf*z
    fld=PC.bg_fields(N,om,V)
    return PC.Lnum(fld,complex(kappa)), fld

if __name__=="__main__":
    for kappa in [2.81055255, complex(2.0,1.0)]:
        L,fld=L_center(kappa)
        ev,Vc=np.linalg.eig(L)
        idx=np.argsort(ev.real)
        print(f"\nkappa={kappa}:")
        for i in idx:
            v=Vc[:,i]; v=v/ (v[np.argmax(np.abs(v))])   # normalize by largest comp
            print(f"  lambda={ev[i].real:+.4f}{ev[i].imag:+.4f}j  eigvec(Abar,Nbar,ombar,V)={np.round(v,4)}")
        # md center modes: {1,1,e^{-2x},e^{-x}}. lambda=0 (x2), -2, -1.
        # Abar_p=0 BC picks the lambda=0 mode with zero Abar component.
        # e^{-x} mode (lambda=-1) excluded by (5.15).
        # e^{-2x} (lambda=-2) is the growing mode to kill.
