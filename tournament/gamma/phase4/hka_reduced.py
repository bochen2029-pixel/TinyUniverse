# hka_reduced.py — HKA background on the CONSTRAINT SURFACE (4.2), a stable 3D system.
#
# KEY FACT (hka_constraint.py): the constraint (4.2) obeys dC/dx = -A*C, so C=0 is invariant but
# only marginally stable -> integrating the full 4D resolved flow toward the center (-x) amplifies
# tiny C-violations exponentially (the chaos the prior effort hit). CURE: use (4.2) to eliminate A
#   A = 1 + 2 om(1+(g-1)V^2)/(1-V^2) + 2 g N V om/(1-V^2)                         (from 4.2)
# and integrate only (N, om, V) via (4.1b,d,e). The trajectory then stays EXACTLY on C=0 by
# construction -> no runaway. The regular center (4.11): A->1, V->0, om->0, N->+inf (NV->-2/(3g)).
#
# Eq numbers per HKA_beta_equations.md.

import numpy as np, math

G=4.0/3.0; GM1=G-1.0; CS=math.sqrt(GM1)
M_CENTER=-2.0/(3.0*G)

def A_of(N,om,V):
    """A from the constraint (4.2)."""
    oV2=1-V*V
    return 1 + 2*om*(1+GM1*V*V)/oV2 + 2*G*N*V*om/oV2

def sonic_data(V0):
    N0=(1-CS*V0)/(CS-V0)
    A0=(G**2+4*G-4+8*GM1**1.5*V0-(3*G-2)*(2-G)*V0**2)/(G**2*(1-V0**2))
    om0=2*CS*(CS-V0)*(1+CS*V0)/(G**2*(1-V0**2))
    return A0,N0,om0

def Dson(N,V):
    return (1+N*V)**2 - GM1*(N+V)**2

def rhs3(Y3):
    """dY3/dx for Y3=(N,om,V) on the constraint surface. A eliminated via (4.2). Uses (4.1b) for N',
    and the fluid pair (4.1d,e) for (om',V'). Singular at the sonic point (det=0)."""
    N,om,V=Y3; g=G; oV2=1-V*V
    A=A_of(N,om,V)
    Nx=N*(-2+A-(2-g)*om)                                     # (4.1b)
    RHS_d=(3*(2-g)/2)*N*V-((2+g)/2)*A*N*V+(2-g)*N*V*om
    RHS_e=(2-g)*(g-1)*N*om+((7*g-6)/2)*N+((2-3*g)/2)*A*N
    M=np.array([[(1+N*V)/om, g*(N+V)/oV2],[(g-1)*(N+V)/om, g*(1+N*V)/oV2]])
    b=np.array([RHS_d,RHS_e])
    omx,Vx=np.linalg.solve(M,b)
    return np.array([Nx,omx,Vx])

def desing3(Y3):
    """Desingularized 3D flow dY3/dxi = d*(N',om',V') with d = Cramer denominator (∝Dson). Smooth at
    the sonic point. Uniform scaling => Jacobian eigenvector is the analytic direction."""
    N,om,V=Y3; g=G; oV2=1-V*V
    A=A_of(N,om,V)
    Nx=N*(-2+A-(2-g)*om)
    RHS_d=(3*(2-g)/2)*N*V-((2+g)/2)*A*N*V+(2-g)*N*V*om
    RHS_e=(2-g)*(g-1)*N*om+((7*g-6)/2)*N+((2-3*g)/2)*A*N
    # fluid pair written as [[1+NV, g(N+V)/oV2 * om],[(g-1)(N+V), g(1+NV)/oV2 * om]] . (om_x/om? no)
    # Cleanest: M.(om_x,V_x)=b with M=[[(1+NV)/om, g(N+V)/oV2],[(g-1)(N+V)/om, g(1+NV)/oV2]], b=(RHS_d,RHS_e).
    # Multiply row-eqs by om: M'=[[1+NV, g(N+V)om/oV2],[(g-1)(N+V), g(1+NV)om/oV2]], b'=(om RHS_d, om RHS_e).
    a11=(1+N*V); a12=g*(N+V)/oV2*om; a21=(g-1)*(N+V); a22=g*(1+N*V)/oV2*om
    d=a11*a22-a12*a21                    # ∝ om*Dson ; det of M'
    numOm=(om*RHS_d)*a22 - a12*(om*RHS_e)      # = om_x * d
    numV = a11*(om*RHS_e) - (om*RHS_d)*a21      # = V_x  * d
    return np.array([Nx*d, numOm, numV])

def eigdir3(V0, eps=1e-7):
    """Leading eigen-direction of the desingularized 3D flow at the sonic point."""
    A0,N0,om0=sonic_data(V0); Y0=np.array([N0,om0,V0])
    f0=desing3(Y0); J=np.zeros((3,3))
    for j in range(3):
        Yp=Y0.copy(); Yp[j]+=eps; J[:,j]=(desing3(Yp)-f0)/eps
    w,Vc=np.linalg.eig(J)
    return Y0,J,w,Vc

def rk4_3(Y3,h):
    k1=rhs3(Y3);k2=rhs3(Y3+0.5*h*k1);k3=rhs3(Y3+0.5*h*k2);k4=rhs3(Y3+h*k3)
    return Y3+(h/6)*(k1+2*k2+2*k3+k4)

def inward(V0, es=+1, h=5e-5, nmax=600000, launch=1e-6, xstop=-20.0):
    """Launch off the sonic point along the real eigendirection, integrate 3D flow in -x toward
    the center. Return (status, x, Y3, A, diag)."""
    Y0,J,w,Vc=eigdir3(V0)
    order=np.argsort(-w.real); vec=Vc[:,order[0]]; vr=np.real(vec)
    if np.linalg.norm(vr)<1e-12: vr=np.real(Vc[:,order[1]])
    vr=vr/np.linalg.norm(vr)
    Y=Y0+es*launch*vr; x=0.0
    Vprev=Y[2]; vz=0; Amin=A_of(*Y); Nmax=Y[0]
    for i in range(nmax):
        Y=rk4_3(Y,-abs(h)); x-=abs(h)
        if not np.all(np.isfinite(Y)): return dict(status='nan',x=x,Y=Y,vz=vz,Nmax=Nmax,Amin=Amin)
        N,om,V=Y; A=A_of(N,om,V)
        Amin=min(Amin,A); Nmax=max(Nmax,N)
        if V*Vprev<0: vz+=1
        Vprev=V
        if N>1e6 or abs(A)>1e4 or abs(V)>10 or om>1e3:
            break
        if N>200 and abs(V)<1e-3 and abs(A-1)<1e-2 and om<1e-2:
            return dict(status='center',x=x,Y=Y,A=A,vz=vz,Nmax=Nmax,Amin=Amin,M=N*V)
        if x<xstop: break
    N,om,V=Y; A=A_of(N,om,V)
    st='center' if (N>1e3 and abs(A-1)<0.2 and abs(V)<0.1) else 'diverge'
    return dict(status=st,x=x,Y=Y,A=A,vz=vz,Nmax=Nmax,Amin=Amin,M=N*V)

if __name__=="__main__":
    import sys,time
    V0=-0.25
    if len(sys.argv)>1: V0=float(sys.argv[1])
    A0,N0,om0=sonic_data(V0)
    print(f"V0={V0}: sonic A0={A0:.5f} N0={N0:.5f} om0={om0:.5f}  A_of(N0,om0,V0)={A_of(N0,om0,V0):.5f} (check=A0)")
    Y0,J,w,Vc=eigdir3(V0)
    print("desing3 Jacobian eigenvalues at sonic:")
    for i in range(3):
        print(f"  lam={w[i].real:+.5f}{w[i].imag:+.5f}j  vec={np.real_if_close(Vc[:,i]).round(4)}")
    print()
    t0=time.time()
    for es in (+1,-1):
        r=inward(V0,es)
        print(f" es={es:+d}: {r['status']:8s} x={r['x']:+.3f} A={r.get('A',float('nan')):.4f} "
              f"N={r['Y'][0]:.4g} om={r['Y'][1]:.5f} V={r['Y'][2]:.4f} M=NV={r.get('M',float('nan')):.4f} "
              f"Vz={r['vz']} Amin={r['Amin']:.4f}")
    print(f"# {time.time()-t0:.1f}s")
