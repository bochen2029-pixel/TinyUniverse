# hka_center.py — Evans-Coleman background by shooting FROM the regular center OUTWARD.
#
# The regular center (4.11) has EXACTLY ONE unstable direction (md line 52). So the natural shoot is
# center -> outward (+x): launch along the center asymptotics (4.12)-(4.13), integrate the resolved
# flow in +x, and shoot the single center parameter so the trajectory reaches the SONIC point
# ANALYTICALLY (passes through Dson=0 smoothly instead of blowing up). This is the standard EC/HKA
# construction (md §IV E) and avoids the unstable inward integration.
#
# Center asymptotics (4.12)-(4.13), z=e^x->0:
#   N = N_inf / z ,  A = 1 + Ai z^2 ,  om = oi z^2 ,  V = Vi z
#   with  Ai = (2/3) oi ,  N_inf * Vi = -2/(3g) = -1/2 .
# Gauge: fix N_inf (x-translation). Free physical parameter: oi (central density) OR equivalently Vi.
# Given N_inf and oi: Vi = -1/(2 N_inf) is FIXED by (4.13)?? No -- (4.13) says N_inf Vi = -1/2, so Vi
# is fixed by N_inf alone. Then oi is the independent shooting parameter (Ai=(2/3)oi follows).
# We verify the seed against the ODE and shoot oi (with N_inf as gauge) to hit the sonic point.
#
# On the constraint surface: A is redundant (A_of(N,om,V) from 4.2). We can either integrate 4D and
# monitor the constraint, or integrate reduced 3D (N,om,V). Reduced is stable -> use it.
# Eq #s per HKA_beta_equations.md.

import numpy as np, math

G=4.0/3.0; GM1=G-1.0; CS=math.sqrt(GM1)
M_CENTER=-2.0/(3.0*G)     # = -1/2

def A_of(N,om,V):
    oV2=1-V*V
    return 1 + 2*om*(1+GM1*V*V)/oV2 + 2*G*N*V*om/oV2

def Dson(N,V):
    return (1+N*V)**2 - GM1*(N+V)**2

def _fluid_slopes(A,N,om,V):
    """(om_x, V_x) via analytic 2x2 inverse (fast; no np.linalg.solve). M.(om_x,V_x)=b,
    M=[[(1+NV)/om, g(N+V)/oV2],[(g-1)(N+V)/om, g(1+NV)/oV2]], b=(RHS_d,RHS_e)."""
    g=G; oV2=1-V*V
    RHS_d=(3*(2-g)/2)*N*V-((2+g)/2)*A*N*V+(2-g)*N*V*om
    RHS_e=(2-g)*(g-1)*N*om+((7*g-6)/2)*N+((2-3*g)/2)*A*N
    a=(1+N*V)/om; b=g*(N+V)/oV2; c=(g-1)*(N+V)/om; d=g*(1+N*V)/oV2
    det=a*d-b*c
    omx=(d*RHS_d-b*RHS_e)/det
    Vx=(a*RHS_e-c*RHS_d)/det
    return omx,Vx

def rhs3(Y3):
    """dY3/dx, Y3=(N,om,V), on the constraint surface (A eliminated via 4.2)."""
    N,om,V=Y3; g=G
    A=A_of(N,om,V)
    Nx=N*(-2+A-(2-g)*om)
    omx,Vx=_fluid_slopes(A,N,om,V)
    return np.array([Nx,omx,Vx])

def rhs4(Y):
    """Full 4D resolved dY/dx (A independent) — to monitor constraint drift as a check."""
    A,N,om,V=Y; g=G; oV2=1-V*V
    Ax=A*(1-A+2*om*(1+GM1*V*V)/oV2)
    Nx=N*(-2+A-(2-g)*om)
    omx,Vx=_fluid_slopes(A,N,om,V)
    return np.array([Ax,Nx,omx,Vx])

def center_seed(N_inf, oi, z):
    """Center asymptotic state (4.12)-(4.13) at z=e^x. Vi fixed by N_inf*Vi=-1/2; Ai=(2/3)oi."""
    Vi = M_CENTER/N_inf            # = -1/(2 N_inf)
    Ai = (2.0/3.0)*oi
    N = N_inf/z
    A = 1 + Ai*z*z
    om = oi*z*z
    V = Vi*z
    return np.array([A,N,om,V])

def rk4_3(Y3,h):
    k1=rhs3(Y3);k2=rhs3(Y3+0.5*h*k1);k3=rhs3(Y3+0.5*h*k2);k4=rhs3(Y3+h*k3)
    return Y3+(h/6)*(k1+2*k2+2*k3+k4)

def shoot_out(N_inf, oi, x0=-14.0, h=2e-4, nmax=400000):
    """Launch from the center at x0 with (N_inf,oi); integrate 3D outward (+x) until the sonic point
    (Dson=0) is reached or the solution breaks. Return the state and Dson history for shooting."""
    z0=math.exp(x0)
    Y4=center_seed(N_inf,oi,z0); Y3=Y4[[1,2,3]]  # (N,om,V)
    x=x0
    Dprev=Dson(Y3[0],Y3[2])
    traj=[(x,Y3.copy())]
    for i in range(nmax):
        Yn=rk4_3(Y3,abs(h)); xn=x+abs(h)
        if not np.all(np.isfinite(Yn)):
            return dict(status='nan',x=x,Y3=Y3,Dprev=Dprev,traj=traj)
        Dn=Dson(Yn[0],Yn[2])
        if Dprev>0 and Dn<=0:
            # crossed the sonic locus (Dson: + inside center side -> 0 at sonic). Interpolate.
            frac=Dprev/(Dprev-Dn)
            Ycross=Y3+frac*(Yn-Y3); xcross=x+frac*abs(h)
            return dict(status='sonic',x=xcross,Y3=Ycross,Dprev=Dprev,Dn=Dn,traj=traj,
                        V=Ycross[2],N=Ycross[0],om=Ycross[1],A=A_of(*Ycross))
        Y3=Yn; x=xn; Dprev=Dn
        if i%50==0: traj.append((x,Y3.copy()))
        if Y3[0]<1e-3 or abs(Y3[2])>5 or Y3[1]<0 or Y3[1]>50:
            return dict(status='break',x=x,Y3=Y3,Dprev=Dprev,traj=traj)
    return dict(status='maxstep',x=x,Y3=Y3,Dprev=Dprev,traj=traj)

def verify_seed(N_inf, oi, x0=-14.0):
    """Check the center seed satisfies the ODE: compare analytic dY/dx (from series) to rhs4(Y)."""
    z0=math.exp(x0); Y=center_seed(N_inf,oi,z0)
    # analytic derivative wrt x: dN/dx=-N_inf/z*? d/dx(N_inf e^{-x})=-N; d/dx(1+Ai e^{2x})=2 Ai e^{2x};
    Vi=M_CENTER/N_inf; Ai=(2.0/3.0)*oi
    dA=2*Ai*z0*z0; dN=-N_inf/z0; dom=2*oi*z0*z0; dV=Vi*z0
    f=rhs4(Y)
    return Y, np.array([dA,dN,dom,dV]), f

if __name__=="__main__":
    import sys
    N_inf=1.0
    print("verify center seed vs ODE (should match to O(z^2)):")
    for oi in [0.5,1.0,2.0]:
        Y,dana,f=verify_seed(N_inf,oi)
        print(f" oi={oi}: A,N,om,V={Y.round(6)}")
        print(f"   analytic dY/dx={dana}")
        print(f"   rhs4(Y)      ={f}")
        print(f"   reldiff N,om,V: {abs((dana[1]-f[1])/f[1]):.2e} {abs((dana[2]-f[2])/f[2]):.2e} {abs((dana[3]-f[3])/f[3]):.2e}")
