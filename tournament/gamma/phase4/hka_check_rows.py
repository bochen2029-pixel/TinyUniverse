# Check each candidate perturbation ROW on the gauge mode independently by finite-differencing the
# gauge vector. This isolates which row(s) are wrong.
import numpy as np, math
import hka_beta4 as B
B.bg(); B.bg_path(); G=4.0/3.0

def gv(x,kbar):
    A,N,om,V,ox,Vx=B.bg_state(x); oV2=1-V*V
    Ax=A*(1-A+2*om*(1+(1/3.0)*V*V)/oV2); Nx=N*(-2+A-(2-G)*om)
    return np.array([Ax/A, Nx/N+kbar, ox, Vx],complex),(A,N,om,V,ox,Vx)

def gvp(x,kbar,dx=1e-6):
    return (gv(x+dx,kbar)[0]-gv(x-dx,kbar)[0])/(2*dx)

for kbar in [1.0, 2.81]:
    print(f"\nkbar={kbar}:")
    for x in [-3.0,-1.0,-0.3]:
        Psi,(A,N,om,V,ox,Vx)=gv(x,kbar); dP=gvp(x,kbar); oV2=1-V*V
        Ab,Nb,ob,Vp=Psi
        # ROW Abar (5.14): Abar' = -kbar Ab - (2gNVom/oV2)(Nb+ob) - (2gNom(1+V^2)/oV2^2)Vp
        rA=-kbar*Ab-(2*G*N*V*om/oV2)*(Nb+ob)-(2*G*N*om*(1+V*V)/oV2**2)*Vp
        # ROW Nbar (lin 4.1b): Nbar'=-kbar Nb + A Ab + (g-2)om ob
        rN=-kbar*Nb+A*Ab+(G-2)*om*ob
        print(f"  x={x:+.2f}: dAbar'={dP[0]:+.4e} pred={rA:+.4e} diff={abs(dP[0]-rA):.2e}")
        print(f"           dNbar'={dP[1]:+.4e} pred={rN:+.4e} diff={abs(dP[1]-rN):.2e}")
        # ROW ombar: try lin of (lnom)'=om'/om. Background (lnom)' depends on all. Perturbation should be
        # d/dx ombar_p = [Jacobian of (lnom)' wrt log vars].Psi_(A,N,om) + (V-part) MINUS kappa*(s-coupling).
        # We'll just report dP[2], dP[3] for now.
        print(f"           dombar'={dP[2]:+.4e}  dV'={dP[3]:+.4e}")
