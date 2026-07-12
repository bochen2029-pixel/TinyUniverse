# hka_survey.py — FAST coarse survey of the inward shoot (uses hardcoded numpy RHS, no sympy loop).
# Launch off the sonic point along the desingularized real eigenvector (computed once per V0 from a
# small numeric Jacobian of the direct-numpy desing flow), integrate the RESOLVED flow inward (-x),
# classify: reaches regular center (A->1,V->0,N large) vs diverges. Locate the EC critical V0 bracket.

import numpy as np, math

G=4.0/3.0; GM1=G-1.0; CS=math.sqrt(GM1)

def sonic_point(V0):
    N0=(1-CS*V0)/(CS-V0)
    A0=(G**2+4*G-4+8*GM1**1.5*V0-(3*G-2)*(2-G)*V0**2)/(G**2*(1-V0**2))
    om0=2*CS*(CS-V0)*(1+CS*V0)/(G**2*(1-V0**2))
    return np.array([A0,N0,om0,V0])

def rhs(Y):
    """Resolved dY/dx (direct numpy). Finite away from sonic. Y=(A,N,om,V)."""
    A,N,om,V=Y; g=G; oV2=1-V*V
    Ax=A*(1-A+2*om*(1+GM1*V*V)/oV2)
    Nx=N*(-2+A-(2-g)*om)
    RHS_d=(3*(2-g)/2)*N*V-((2+g)/2)*A*N*V+(2-g)*N*V*om
    RHS_e=(2-g)*(g-1)*N*om+((7*g-6)/2)*N+((2-3*g)/2)*A*N
    M=np.array([[(1+N*V)/om, g*(N+V)/oV2],[(g-1)*(N+V)/om, g*(1+N*V)/oV2]])
    b=np.array([RHS_d,RHS_e])
    omx,Vx=np.linalg.solve(M,b)
    return np.array([Ax,Nx,omx,Vx])

import hka_desing as _D    # verified symbolic desingularized flow (for the sonic eigendirection only)

def eigdir(V0, eps=1e-6):
    """Real leading-eigenvalue direction of the desingularized-flow Jacobian at the sonic point,
    from hka_desing (verified). Called once per V0 (cheap: ~5 evals)."""
    Y0, J, w, Vc = _D.sonic_jacobian(V0)
    order=np.argsort(-w.real); vec=Vc[:,order[0]]; vr=np.real(vec); vr/=np.linalg.norm(vr)
    return Y0, vr, w[order[0]]

def rk4(Y,h):
    k1=rhs(Y);k2=rhs(Y+0.5*h*k1);k3=rhs(Y+0.5*h*k2);k4=rhs(Y+h*k3)
    return Y+(h/6)*(k1+2*k2+2*k3+k4)

def inward(V0, es, h=1e-4, nmax=200000):
    """Launch inward, integrate resolved flow in -x. Return trajectory summary."""
    Y0,vr,lam=eigdir(V0)
    Y=Y0+es*abs(1e-5)*vr; x=0.0
    Amin=Y[0]; Amax=Y[0]; Vprev=Y[3]; vz=0
    hit_second_sonic=False; Nmax=Y[1]
    Dprev=(1+Y[1]*Y[3])**2-GM1*(Y[1]+Y[3])**2
    for i in range(nmax):
        Y=rk4(Y,-abs(h)); x-=abs(h)
        if not np.all(np.isfinite(Y)): return ('nan',x,Y,Amin,Amax,vz,Nmax)
        A,N,om,V=Y
        Amin=min(Amin,A);Amax=max(Amax,A);Nmax=max(Nmax,N)
        if V*Vprev<0: vz+=1
        Vprev=V
        D=(1+N*V)**2-GM1*(N+V)**2
        if D*Dprev<0: hit_second_sonic=True
        Dprev=D
        if N>1e5 or abs(A)>1e3 or abs(V)>1e3 or om<-1e-3 or om>1e3:
            break
        if N>100 and abs(V)<1e-3 and abs(A-1)<1e-2:
            return ('center',x,Y,Amin,Amax,vz,Nmax)
    A,N,om,V=Y
    st='center' if (N>1e3 and abs(A-1)<0.3 and abs(V)<0.2) else ('2sonic' if hit_second_sonic else 'diverge')
    return (st,x,Y,Amin,Amax,vz,Nmax)

def classify(V0):
    """Try both eps signs; return the better (closer to center) with a signed 'miss' functional.
    miss>0 vs miss<0 should bracket the EC solution."""
    best=None
    for es in (+1,-1):
        st,x,Y,Amin,Amax,vz,Nmax=inward(V0,es)
        A,N,om,V=Y
        # functional: A at the deepest reliable point relative to 1 (regular => A->1).
        # Use sign of (A-1) at breakpoint combined with V. For a clean bracket use log(A) if N grew.
        score=(Nmax, -abs(A-1))
        if best is None or score>best[0]:
            best=(score,(st,x,Y,Amin,Amax,vz,Nmax,es))
    return best[1]

if __name__=="__main__":
    import sys, time
    t0=time.time()
    print(f"{'V0':>8} {'status':>8} {'x':>8} {'A':>10} {'N':>10} {'om':>9} {'V':>9} {'Amin':>8} {'Vz':>3} {'es':>3}")
    for V0 in np.arange(-0.54,-0.05,0.02):
        st,x,Y,Amin,Amax,vz,Nmax,es=classify(round(V0,4))
        A,N,om,V=Y
        print(f"{V0:8.3f} {st:>8} {x:8.3f} {A:10.4g} {N:10.4g} {om:9.4f} {V:9.4f} {Amin:8.4f} {vz:3d} {es:+d}")
    print(f"# {time.time()-t0:.1f}s")
