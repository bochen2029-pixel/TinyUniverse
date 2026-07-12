# Compare my derived Euler ODEs vs KHA's fluid ODEs by solving each 2x2 for (om_x,V_x) at
# random (N,A,om,V) points. If dV/dx, dom/dx differ, KHA (as transcribed) has an error and
# my covariant derivation is the ground truth. Also re-check the CENTER limit dV/dx with MINE.
import sympy as sp, pickle, numpy as np
en_s,mo_s=pickle.load(open("fluid_manual.pkl","rb"))
en=sp.sympify(en_s); mo=sp.sympify(mo_s)
A0,N0,o0,V0,Ax,Nx,ox,Vx,Es=sp.symbols('A0 N0 om0 V0 A_x N_x om_x V_x E')
# my equations are en=0, mo=0. Substitute KHA's exact A_x,N_x (metric slopes, verified) then
# solve the 2x2 {en,mo}=0 for (ox,Vx).
Axv=A0*(1-A0+(2*o0/(1-V0**2))*(1+V0**2/3))
Nxv=N0*(-2+A0-sp.Rational(2,3)*o0)
en2=en.subs({Ax:Axv,Nx:Nxv})
mo2=mo.subs({Ax:Axv,Nx:Nxv})
sol=sp.solve([en2,mo2],[ox,Vx],dict=True)
print("my-derivation solutions:",len(sol))
omx_mine=sp.simplify(sol[0][ox]); Vx_mine=sp.simplify(sol[0][Vx])
f_omx_mine=sp.lambdify((A0,N0,o0,V0,Es),omx_mine,'numpy')
f_Vx_mine =sp.lambdify((A0,N0,o0,V0,Es),Vx_mine,'numpy')

# KHA's fluid eqs -> solve for (u=omega_x/omega, Vx):
def kha_slopes(N,A,om,V):
    Ap=A*(1-A+(2*om/(1-V*V))*(1+V*V/3.0)); Np=N*(-2+A-2.0*om/3.0)
    omV2=1-V*V
    m11=(1+N*V);m12=4*(N+V)/(3*omV2);m21=(4*V+N+3*N*V*V);m22=4*(1+V*V+2*N*V)/omV2
    r1=-(-N*V*Ap/(3*A)+4*V*Np/3.0+2*N*(1+4*om/(9*omV2)))
    r2=-( N*omV2*Ap/A+4*(1+V*V)*Np+2*N*(1+3*V*V))
    det=m11*m22-m12*m21
    u=(r1*m22-m12*r2)/det; Vp=(m11*r2-r1*m21)/det
    return u*om, Vp   # dom/dx, dV/dx

# E=e^x cancels? my ODEs may still contain E. Check if omx_mine depends on E:
print("omx_mine depends on E?", omx_mine.has(Es), " Vx_mine depends on E?", Vx_mine.has(Es))
print()
print(" test points:  (dom/dx, dV/dx) MINE  vs  KHA")
np.random.seed(0)
for _ in range(6):
    N=np.random.uniform(1.2,3); A=np.random.uniform(0.8,2.5); om=np.random.uniform(0.05,0.8); V=np.random.uniform(-0.5,0.5)
    E=np.exp(np.random.uniform(-1,1))
    om_mine=float(f_omx_mine(A,N,om,V,E)); v_mine=float(f_Vx_mine(A,N,om,V,E))
    om_k,v_k=kha_slopes(N,A,om,V)
    print(f"  N={N:.2f} A={A:.2f} om={om:.2f} V={V:+.2f} E={E:.2f}:")
    print(f"     dom/dx MINE={om_mine:+.4f}  KHA={om_k:+.4f}   {'MATCH' if abs(om_mine-om_k)<1e-6 else 'DIFF'}")
    print(f"     dV/dx  MINE={v_mine:+.4f}  KHA={v_k:+.4f}   {'MATCH' if abs(v_mine-v_k)<1e-6 else 'DIFF'}")

# CENTER limit with MINE: V=0,A=1,om->0,N large
print("\n CENTER limit (A=1,V=0,om->0,N large) dV/dx MINE:")
for N in [1e3,1e5]:
    print(f"   N={N:.0e}: dV/dx={float(f_Vx_mine(1.0,N,1e-9,0.0,1.0)):+.5f}  dom/dx={float(f_omx_mine(1.0,N,1e-9,0.0,1.0)):+.3e}")
