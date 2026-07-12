# arbiter_divT.py — GROUND TRUTH: which fluid-slope system (css-rederived vs KHA-printed) actually
# satisfies nabla_a T^{ab} = 0 ?  Build the full (t,r) fields from the self-similar ansatz with the
# candidate x-slopes, compute nabla_a T^{ab} by FINITE DIFFERENCES (independent of any symbolic
# derivation), and report the residual for each candidate. The correct system -> residual ~ 0.
import numpy as np, sympy as sp, pickle, math

# self-similar ansatz: fields depend on x=ln(-r/t). At a chosen (t0,r0) we know x0 and the field
# values (N,A,om,V). To build f(t,r) locally we use f(t,r)=F(x) with x=ln(-r/t); we get F(x) near x0
# from a local quadratic using the candidate slopes F'(x0) (and F''(x0) from d/dx of the slope RHS).
# metric: alpha = N*sqrt(A)*(r/(-t))  [E=r/(-t)? sign]; a=sqrt(A). Use t<0.
# Actually KHA E-def: N = alpha a^{-1} e^{-x}, x=ln(-r/t) => e^{-x} = -t/r => alpha = N*a*(-t/r)*... wait
#   N = alpha/(a) * e^{-x} => alpha = N a e^{x} = N a (-r/t). With t<0, -r/t>0. So alpha=N*a*(-r/t).
# rho = om/(4 pi r^2 A). u^t=1/(alpha sqrt(1-V^2)), u^r=V/(a sqrt(1-V^2)). (a=sqrt A)

def build_field_funcs(x0, vals, slopes, slopes2):
    """Return F_i(x) local quadratic for i in N,A,om,V."""
    def mk(i):
        return lambda x: vals[i]+slopes[i]*(x-x0)+0.5*slopes2[i]*(x-x0)**2
    return [mk(i) for i in range(4)]

def Tup(t,r, Ffun):
    x=math.log(-r/t)
    N=Ffun[0](x); A=Ffun[1](x); om=Ffun[2](x); V=Ffun[3](x)
    a=math.sqrt(A); alpha=N*a*(-r/t)
    rho=om/(4*math.pi*r*r*A); p=rho/3.0
    W=1.0/math.sqrt(1-V*V)
    ut=W/alpha; ur=W*V/a
    u=[ut,ur,0.0,0.0]
    # metric inverse (diagonal): g^tt=-1/alpha^2, g^rr=1/a^2, g^thth=1/r^2, g^phph=1/(r^2 sin^2)
    gi=[-1/alpha**2, 1/a**2, 1/r**2, 0.0]  # phi handled separately; use theta=pi/2
    T=np.zeros((4,4))
    for i in range(4):
        for j in range(4):
            gij= gi[i] if i==j else 0.0
            T[i,j]=(rho+p)*u[i]*u[j]+p*gij
    return T

def metric_diag(t,r,Ffun,theta=math.pi/2):
    x=math.log(-r/t); N=Ffun[0](x); A=Ffun[1](x)
    a=math.sqrt(A); alpha=N*a*(-r/t)
    return np.array([-alpha**2, a**2, r**2, r**2*math.sin(theta)**2])

def christoffel(t,r,Ffun,theta=math.pi/2,ht=1e-6,hr=1e-6):
    # numeric metric and derivatives (diagonal metric). coords: 0=t,1=r,2=theta,3=phi
    def gd(tt,rr,th): return metric_diag(tt,rr,Ffun,th)
    g0=gd(t,r,theta)
    dg=[np.zeros(4) for _ in range(4)]  # dg[c] = d g_ii / d x^c
    dg[0]=(gd(t+ht,r,theta)-gd(t-ht,r,theta))/(2*ht)
    dg[1]=(gd(t,r+hr,theta)-gd(t,r-hr,theta))/(2*hr)
    dth=1e-6; dg[2]=(gd(t,r,theta+dth)-gd(t,r,theta-dth))/(2*dth)
    dg[3]=np.zeros(4)
    gi=1.0/g0
    # Gamma^l_{ij} = 1/2 g^{ll}( d_i g_{lj}+ d_j g_{li} - d_l g_{ij} ), diagonal metric
    G=np.zeros((4,4,4))
    for l in range(4):
        for i in range(4):
            for j in range(4):
                term=0.0
                # d_i g_{lj}: nonzero only if l==j
                if l==j: term+=dg[i][l]
                if l==i: term+=dg[j][l]
                if i==j: term-=dg[l][i]
                G[l,i,j]=0.5*gi[l]*term
    return G

def divT(t,r,Ffun,ht=1e-6,hr=1e-6):
    # nabla_a T^{a b} = d_a T^{ab} + Gamma^a_{ac} T^{cb} + Gamma^b_{ac} T^{ac}
    T0=Tup(t,r,Ffun)
    dT=[np.zeros((4,4)) for _ in range(4)]
    dT[0]=(Tup(t+ht,r,Ffun)-Tup(t-ht,r,Ffun))/(2*ht)
    dT[1]=(Tup(t,r+hr,Ffun)-Tup(t,r-hr,Ffun))/(2*hr)
    # theta,phi derivatives of T^{ab}: T depends on theta only via g^{thth},g^{phph} (p g^{ij}); small
    dth=1e-6
    # build T at theta+-: need metric_diag theta-dep in Tup -> Tup uses theta=pi/2 fixed; approximate d/dtheta=0 for T^{tt},T^{rr}
    dT[2]=np.zeros((4,4)); dT[3]=np.zeros((4,4))
    G=christoffel(t,r,Ffun)
    out=np.zeros(4)
    for b in range(4):
        s=0.0
        for a in range(4):
            s+=dT[a][a,b]
        for a in range(4):
            for c in range(4):
                s+=G[a,a,c]*T0[c,b]
                s+=G[b,a,c]*T0[a,c]
        out[b]=s
    return out

# --- candidate slope systems ---
r=sp.Symbol('r',positive=True); A0,N0,o0,V0=sp.symbols('A0 N0 om0 V0')
d=pickle.load(open("css_syms.pkl","rb"))
css_Ap=sp.lambdify((A0,N0,o0,V0),sp.sympify(d['Axv']),'math')
css_Np=sp.lambdify((A0,N0,o0,V0),sp.sympify(d['Nxv']),'math')
css_wp=sp.lambdify((A0,N0,o0,V0),sp.sympify(d['omx']).subs(r,1),'math')
css_Vp=sp.lambdify((A0,N0,o0,V0),sp.sympify(d['Vxx']).subs(r,1),'math')
# KHA:
N,A,w,V=sp.symbols('N A w V')
Apk=A*(1-A+(2*w/(1-V**2))*(1+V**2/3)); Npk=N*(-2+A-sp.Rational(2,3)*w)
g,u=sp.symbols('g u')
f1=(1+N*V)*g+4*(N+V)*u/(3*(1-V**2))-N*V*Apk/(3*A)+4*V*Npk/3+2*N*(1+4*w/(9*(1-V**2)))
f2=(4*V+N+3*N*V**2)*g+4*(1+V**2+2*N*V)*u/(1-V**2)+N*(1-V**2)*Apk/A+4*(1+V**2)*Npk+2*N*(1+3*V**2)
sol=sp.solve([f1,f2],[g,u],dict=True)[0]
kha_wp=sp.lambdify((N,A,w,V),sp.simplify(w*sol[g]),'math'); kha_Vp=sp.lambdify((N,A,w,V),sp.simplify(sol[u]),'math')
kha_Ap=sp.lambdify((N,A,w,V),Apk,'math'); kha_Np=sp.lambdify((N,A,w,V),Npk,'math')

# Simpler: just test with slopes only (quadratic term small if we use tiny dx window). Use slopes2=0 and
# small FD step so the quadratic doesn't matter much. Evaluate residual scaled by |T|.
def test(vals, sys_name, Ap,Np,wp,Vp, order='ANV'):
    N_,A_,om_,V_=vals
    if order=='css':  # css funcs take (A,N,om,V)
        sl=[Np(A_,N_,om_,V_), Ap(A_,N_,om_,V_), wp(A_,N_,om_,V_), Vp(A_,N_,om_,V_)]
    else:             # kha funcs take (N,A,w,V)
        sl=[Np(N_,A_,om_,V_), Ap(N_,A_,om_,V_), wp(N_,A_,om_,V_), Vp(N_,A_,om_,V_)]
    x0=math.log(2.0)  # arbitrary; pick t,r with -r/t=2 => e.g. t=-1,r=2
    t0=-1.0; r0=2.0
    Ff=build_field_funcs(x0,[N_,A_,om_,V_],sl,[0,0,0,0])
    T0=Tup(t0,r0,Ff); scale=np.max(np.abs(T0))
    res=divT(t0,r0,Ff)
    print(f"  {sys_name}: |divT|/|T| per comp = {np.array2string(np.abs(res)/scale, precision=3)}  max={np.max(np.abs(res))/scale:.3e}")

vals=(1.3,1.2,0.4,0.3)  # N,A,om,V generic (subsonic-ish)
print(f"Ground-truth nabla_a T^ab residual at (N,A,om,V)={vals}, (t,r)=(-1,2):")
test(vals,"css-rederived",css_Ap,css_Np,css_wp,css_Vp,'css')
test(vals,"KHA-printed  ",kha_Ap,kha_Np,kha_wp,kha_Vp,'kha')
print("\n(The system with residual ~0 satisfies conservation and is the correct fluid system.)")
