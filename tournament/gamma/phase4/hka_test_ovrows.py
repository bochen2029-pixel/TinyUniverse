# Isolate the (ombar,V) rows: do they satisfy the gauge mode? Test row2 (ombar') and row3 (V')
# separately, and also test whether a -kappa diagonal shift is needed on these rows.
import numpy as np, math, sympy as sp
import hka_beta4 as B
B.bg(); B.bg_path(); G=4.0/3.0
g=sp.Rational(4,3)
A,N,om,V=sp.symbols('A N om V',real=True); omx,Vx=sp.symbols('om_x V_x',real=True); k=sp.symbols('kappa')
oV2=1-V**2
Ax=1+N*V; Bx=g*(N+V)/oV2; Cx=(g-1)*(N+V); Dx=g*(1+N*V)/oV2
As=sp.Integer(1); Bs=g*V/oV2; Cs=(g-1)*V; Ds=g/oV2
E1=-((g+2)/2)*A*N*V
E2=((6-3*g)/2)*N*V-((2+g)/2)*A*N*V+(2-g)*om*N*V-N*V*omx-g*N*Vx/oV2
E3=(2-g)*om*N*V
E4=((6-3*g)/2)*N-((2+g)/2)*A*N+(2-g)*om*N-N*omx-g*(1+2*N*V+V**2)*Vx/oV2**2
F1=((2-3*g)/2)*A*N
F2=(2-g)*(g-1)*om*N+((7*g-6)/2)*N+((2-3*g)/2)*A*N-(g-1)*N*omx-g*N*V*Vx/oV2
F3=(2-g)*(g-1)*om*N
F4=-(g-1)*omx-g*(N+2*V+N*V**2)*Vx/oV2**2
Mx=sp.Matrix([[Ax,Bx],[Cx,Dx]]); Minv=Mx.inv(); Ms=sp.Matrix([[As,Bs],[Cs,Ds]])
Evec=sp.Matrix([[E1,E2,E3,E4],[F1,F2,F3,F4]])
args=(A,N,om,V,omx,Vx,k)
fMinv=sp.lambdify(args,Minv,'numpy'); fEvec=sp.lambdify(args,Evec,'numpy'); fMs=sp.lambdify(args,Ms,'numpy')

def bgf(x): return B.bg_state(x)
def gauge_vec(x,kbar):
    Ab,Nb,ob,Vb,ox,Vxb=bgf(x); oV2n=1-Vb*Vb
    Axn=Ab*(1-Ab+2*ob*(1+(1/3.0)*Vb*Vb)/oV2n); Nxn=Nb*(-2+Ab-(2-G)*ob)
    return np.array([Axn/Ab, Nxn/Nb+kbar, ox, Vxb],complex)

print("Test (ombar',V') rows on gauge mode. LHS = actual d/dx of (ombar_p,V_p). RHS variants:")
for kbar in [1.0, 2.81]:
    print(f"\n kbar={kbar}:")
    for x in [-4.0,-1.0,-0.3]:
        Psi=gauge_vec(x,kbar); dP=(gauge_vec(x+1e-5,kbar)-gauge_vec(x-1e-5,kbar))/2e-5
        f=bgf(x)
        Minv_=np.array(fMinv(*f,complex(kbar)),complex); Ev=np.array(fEvec(*f,complex(kbar)),complex); Ms_=np.array(fMs(*f,complex(kbar)),complex)
        ovp=Psi[[2,3]]                     # (ombar_p, V_p)
        # variant A (md 5.13): (ombar',V') = Minv[ E.Psi + kappa Ms ovp ]
        rhsA=Minv_.dot(Ev.dot(Psi) + kbar*Ms_.dot(ovp))
        # variant B: with an extra -kappa on the LHS moved to RHS: (ombar',V')+kappa ovp = Minv[...]
        rhsB=rhsA - kbar*ovp
        lhs=dP[[2,3]]
        print(f"  x={x:+.2f}: LHS(ombar',V')={lhs.round(4)}")
        print(f"        variantA Minv[E+kMs] ={rhsA.round(4)}  diff={np.linalg.norm(lhs-rhsA):.2e}")
        print(f"        variantB -kappa shift ={rhsB.round(4)}  diff={np.linalg.norm(lhs-rhsB):.2e}")
