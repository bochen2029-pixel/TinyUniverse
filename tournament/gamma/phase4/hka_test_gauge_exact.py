# Test whether the GAUGE mode (5.20) is an EXACT solution of Psi'=L(x;kappa)Psi, which validates the
# perturbation operator L. Gauge mode (md line 108):
#   Psi_gauge(x) = ( (lnA)'_x , (lnN)'_x + kbar , (lnom)'_x , V'_x )
# with eigenvalue kappa=kbar. (The +kbar in the Nbar component is the gauge shift.) For ANY kbar this
# should solve the perturbation equations (it's pure gauge). We verify d/dx(Psi_gauge) = L(x;kbar).Psi_gauge
# along the EC background. If residual ~0, L is correct AND we have the true gauge-mode vector.
import numpy as np, math
import hka_pert_core as PC
import hka_ec as E

G=4.0/3.0

def bg_full(x, N_inf=1.0, oi=0.375):
    """Background (A,N,om,V) and log-derivatives along the EC solution via the true integrated bg."""
    # use the reduced background integrated (accurate) — reuse hka_beta4.bg_state which integrates.
    import hka_beta4 as B; B.bg(); B.bg_path()
    fld=B.bg_state(x)   # (A,N,om,V,ombar_x,V_x)
    return fld

def gauge_vec(x, kbar):
    A,N,om,V,ombar_x,Vx=bg_full(x)
    # (lnA)'_x = A'/A ; A' from 4.1a
    oV2=1-V*V
    Ax=A*(1-A+2*om*(1+(1/3.0)*V*V)/oV2)
    Nx=N*(-2+A-(2-G)*om)
    lnAx=Ax/A; lnNx=Nx/N; lnomx=ombar_x   # ombar_x IS om'/om
    return np.array([lnAx, lnNx+kbar, lnomx, Vx],complex)

def gauge_vec_deriv(x, kbar, dx=1e-5):
    return (gauge_vec(x+dx,kbar)-gauge_vec(x-dx,kbar))/(2*dx)

if __name__=="__main__":
    import hka_beta4 as B; B.bg(); B.bg_path()
    xs=B.bg()['xs']
    print("Test: is the gauge mode Psi_gauge an exact solution of Psi'=L Psi?  (residual d/dx Psi - L Psi)")
    for kbar in [0.35699, 1.0, 2.81]:
        print(f"\n kbar={kbar}:")
        for x in [-8.0, -4.0, -1.0, xs-0.05]:
            Psi=gauge_vec(x,kbar)
            dPsi=gauge_vec_deriv(x,kbar)
            L=PC.Lnum(B.bg_state(x),complex(kbar))
            res=dPsi-L.dot(Psi)
            print(f"   x={x:+.3f}: |Psi|={np.linalg.norm(Psi):.3f} |d/dx Psi - L Psi|={np.linalg.norm(res):.3e}  rel={np.linalg.norm(res)/np.linalg.norm(dPsi):.2e}")
