# hka_pert_reduced.py — perturbation on the CONSTRAINT SURFACE (3D), first-principles.
#
# Background reduced flow, Y3=(lnN,lnom,V):  dY3/dx = F3(Y3), with A=A_of(N,om,V) (constraint 4.2).
# The perturbation delta Y3 = psi(x) e^{kappa s} satisfies psi' = [J3(x) - kappa Ks] psi, where
#   J3 = dF3/dY3   (flow Jacobian, first-principles, exact for kappa=0),
#   Ks = the s-coupling. For the fluid (lnom,V) block the s-coupling is Minv Ms (HKA 5.5,5.6); the
#   lnN (metric) equation's s-coupling we DETERMINE by requiring the reduced gauge mode to be exact.
#
# Reduced gauge mode: the gauge shift acts as delta(lnN)= (lnN)'*? ... In the reduced variables the
# gauge mode is psi_g = d/dx(Y3_bg) + kappa * e_gauge, where e_gauge is the gauge direction. We find
# e_gauge and Ks by validation.
#
# We build J3 symbolically and test.
import numpy as np, sympy as sp, math, pickle

def build():
    g=sp.Rational(4,3)
    lN,lo,V=sp.symbols('lN lo V',real=True)
    N=sp.exp(lN); om=sp.exp(lo); oV2=1-V**2
    A=1+2*om*(1+(g-1)*V**2)/oV2+2*g*N*V*om/oV2         # A_of via constraint (4.2)
    Nx=N*(-2+A-(2-g)*om)
    RHS_d=(3*(2-g)/2)*N*V-((2+g)/2)*A*N*V+(2-g)*N*V*om
    RHS_e=(2-g)*(g-1)*N*om+((7*g-6)/2)*N+((2-3*g)/2)*A*N
    a=(1+N*V)/om; b=g*(N+V)/oV2; c=(g-1)*(N+V)/om; d=g*(1+N*V)/oV2
    det=a*d-b*c
    omx=(d*RHS_d-b*RHS_e)/det; Vx=(a*RHS_e-c*RHS_d)/det
    F=sp.Matrix([Nx/N, omx/om, Vx])          # (lnN)', (lnom)', V'
    J=F.jacobian([lN,lo,V])
    J=sp.simplify(J)
    pickle.dump(dict(J=sp.srepr(J)),open("hka_pert_reduced.pkl","wb"))
    return J,(lN,lo,V)

if __name__=="__main__":
    import time
    t0=time.time(); J,syms=build(); lN,lo,V=syms
    print(f"built reduced Jacobian J3 in {time.time()-t0:.1f}s")
    fJ=[[sp.lambdify((lN,lo,V),J[i,j],'numpy') for j in range(3)] for i in range(3)]
    import hka_beta4 as B; B.bg(); B.bg_path(); G=4.0/3.0
    def y3(x):
        A,N,om,V,ox,Vx=B.bg_state(x); return (math.log(N),math.log(om),V)
    def F3(x):
        A,N,om,V,ox,Vx=B.bg_state(x)
        return np.array([N*(-2+A-(2-G)*om)/N, ox, Vx])   # (lnN)',(lnom)',V'
    # reduced gauge mode: delta from an infinitesimal x-shift is dY3/dx = F3. That's the kappa=0 mode.
    # Test kappa=0: psi=F3 should satisfy psi'=J3 psi (variational eq along the flow). Check.
    print("\nkappa=0 check: does psi=F3 (flow tangent) satisfy psi'=J3 psi?")
    for x in [-3.0,-1.0,-0.3]:
        psi=F3(x); dpsi=(F3(x+1e-6)-F3(x-1e-6))/2e-6
        Jn=np.array([[fJ[i][j](*y3(x)) for j in range(3)] for i in range(3)])
        res=dpsi-Jn.dot(psi)
        print(f"  x={x:+.2f}: |psi'-J3 psi|={np.linalg.norm(res):.2e} rel={np.linalg.norm(res)/np.linalg.norm(dpsi):.2e}")
    # eigenvalues of J3 at center and sonic (indicial exponents of the reduced perturbation)
    print("\nJ3 eigenvalues (reduced indicial exponents):")
    xs=B.bg()['xs']
    for x in [-12.0, -6.0, xs-0.05]:
        Jn=np.array([[fJ[i][j](*y3(x)) for j in range(3)] for i in range(3)])
        print(f"  x={x:+.3f}: eig(J3)={np.sort_complex(np.linalg.eigvals(Jn)).round(4)}")
