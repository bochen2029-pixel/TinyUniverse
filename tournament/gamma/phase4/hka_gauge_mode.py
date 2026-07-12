# Identify the GAUGE mode (5.20) among the center lambda=0 modes, so the physical center-regular
# subspace is 1D. Gauge mode (md line 108): h_gauge = {d/dx(ln A_ss), d/dx(ln N_ss)+kbar, d/dx(ln om_ss),
# ... } -- more precisely, in the log variables (Abar,Nbar,ombar,V)_p the gauge mode is
#   Abar_p = A_ss'(x)/A_ss * s-factor? The md gives h_gauge = e^{kbar s}{ h'_ss for Abar,ombar,V ;
#   h'_ss + kbar for Nbar }, where h'_ss = d/dx of the background log field. So as a perturbation VECTOR
# at fixed x:  Psi_gauge = ( (ln A)'_x , (ln N)'_x + kbar? , (ln om)'_x , V'_x )  -- but the mode has
# its own kbar eigenvalue. In our L(x;kappa) framework the gauge mode is the solution with kappa=kbar
# whose vector is the x-derivative of the background logs. We build it and check it's a solution.
import numpy as np, math
import hka_pert_core as PC
import hka_ec as E

def bg_log_derivs(x, N_inf=1.0, oi=0.375):
    """(ln A)'_x, (ln N)'_x, (ln om)'_x, V'_x along the background (using center asymptotics valid deep)."""
    z=math.exp(x); N=N_inf/z; om=oi*z*z; V=E.M_CENTER/N_inf*z; A=E.A_of(N,om,V)
    omx,Vx=E._fluid_slopes(A,N,om,V)      # om', V'
    Ax=A*(1-A+2*(1/3.0)*V*V/ (1-V*V) *0 + 2*om*(1+(1/3.0)*V*V)/(1-V*V))  # A' (use full 4.1a)
    Ax=A*(1-A+2*om*(1+(1/3.0)*V*V)/(1-V*V))
    Nx=N*(-2+A-(2-4/3.0)*om)              # N'
    return np.array([Ax/A, Nx/N, omx/om, Vx])

if __name__=="__main__":
    # The gauge mode has some eigenvalue kbar. At the sonic gauge Nbar_p(0)=0 the md says kbar=0.357;
    # at origin gauge Nbar_p(-inf)=0, kbar=1. So the gauge mode as a center-regular (bounded) mode is
    # the lambda=0 mode that equals the background-log-derivative vector. Let's compare the center
    # lambda=0 eigenvectors to bg_log_derivs.
    for xc in [-13.0]:
        z=math.exp(xc); N=1.0/z; om=0.375*z*z; V=E.M_CENTER*z
        for kappa in [1.0, 0.35699, 2.81055255]:
            L=PC.Lnum(PC.bg_fields(N,om,V),complex(kappa))
            ev,Vr=np.linalg.eig(L); idx=np.argsort(np.abs(ev.real))
            m0=Vr[:,idx[0]]; m1=Vr[:,idx[1]]     # two lambda~0 modes
            g=bg_log_derivs(xc)
            # normalize
            g=g/np.linalg.norm(g)
            print(f"xc={xc} kappa={kappa}: bg-log-deriv vec (Abar',Nbar',ombar',V')={np.round(g,4)}")
            print(f"   center lambda=0 mode0 = {np.round(m0/np.linalg.norm(m0),4)}")
            print(f"   center lambda=0 mode1 = {np.round(m1/np.linalg.norm(m1),4)}")
            # overlap of g with the 2D lambda=0 span
            B=np.column_stack([m0,m1]); Q,_=np.linalg.qr(B)
            proj=Q.dot(Q.conj().T.dot(g)); rej=np.linalg.norm(g-proj)
            print(f"   |g - proj_(lambda=0 span)| = {rej:.3e}  (small => gauge mode is IN the regular span)")
            print()
