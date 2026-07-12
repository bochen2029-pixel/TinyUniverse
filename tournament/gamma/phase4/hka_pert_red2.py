# Determine the reduced kappa-coupling Ks by requiring the reduced gauge mode exact for ALL kappa.
# L_red = J3 + kappa * Ks. Reduced gauge mode psi_g = (lnN,lnom,V components of the 4D gauge mode)
#        = ( (lnN)'+kappa , (lnom)' , V' ).  [the lnN component carries the +kappa gauge shift]
# psi_g' = L_red(kappa) psi_g must hold for all kappa. Expand: psi_g = D3 + kappa e1, D3=((lnN)',(lnom)',V'),
# e1=(1,0,0). psi_g'=D3'. RHS=J3 D3 + kappa(J3 e1 + Ks D3)+kappa^2 Ks e1.
#  kappa^0: D3'=J3 D3 (holds). kappa^2: Ks e1=0 -> Ks col0=0. kappa^1: J3 e1 + Ks D3=0 -> Ks D3=-J3 e1.
# So Ks must satisfy: Ks[:,0]=0 and Ks.D3=-J3[:,0]. With Ks nonzero only in fluid rows (1,2) x cols(1,2)
# (the MinvMs block), plus possibly a lnN-row coupling. We SOLVE for Ks numerically at each x from these
# constraints + the structural ansatz, and check consistency.
import numpy as np, sympy as sp, math, pickle
import hka_beta4 as B
B.bg(); B.bg_path(); G=4.0/3.0
d=pickle.load(open("hka_pert_reduced.pkl","rb")); lN,lo,V=sp.symbols('lN lo V',real=True)
J=sp.sympify(d['J']); fJ=[[sp.lambdify((lN,lo,V),J[i,j],'numpy') for j in range(3)] for i in range(3)]
# MinvMs block
g=sp.Rational(4,3); oV2=1-V**2; N=sp.exp(lN); om=sp.exp(lo)
Ax_=1+N*V; Bx_=g*(N+V)/oV2; Cx_=(g-1)*(N+V); Dx_=g*(1+N*V)/oV2
As=sp.Integer(1); Bs=g*V/oV2; Cs=(g-1)*V; Ds=g/oV2
MM=sp.Matrix([[Ax_,Bx_],[Cx_,Dx_]]).inv()*sp.Matrix([[As,Bs],[Cs,Ds]])
fMM=[[sp.lambdify((lN,lo,V),MM[i,j],'numpy') for j in range(2)] for i in range(2)]

def y3(x):
    A,N,om,V,ox,Vx=B.bg_state(x); return (math.log(N),math.log(om),V)
def D3(x):
    A,N,om,V,ox,Vx=B.bg_state(x); return np.array([N*(-2+A-(2-G)*om)/N, ox, Vx])
def Jn(x): return np.array([[fJ[i][j](*y3(x)) for j in range(3)] for i in range(3)])
def MMn(x): return np.array([[fMM[i][j](*y3(x)) for j in range(2)] for i in range(2)])

# Constraint check: is J3[:,0] (the response to lnN perturbation) in the span such that a fluid-only
# Ks can satisfy Ks D3 = -J3[:,0] with Ks[:,0]=0? Ks acts: (Ks psi)_row for rows 1,2 = MM.(psi_1,psi_2)*s.
# We want the FLUID s-coupling. Test candidate L_red = J3 + kappa Ks with Ks = [[0,0,0],[0,MM],[.. ]] in
# rows (lnom,V), cols (lnom,V): i.e. Ks[1:,1:]=s*MM, else 0. Then check gauge-mode exactness for all kappa.
def Lred(x,kap,s=+1.0):
    Jm=Jn(x); K=np.zeros((3,3)); mm=MMn(x)
    K[1,1]=s*mm[0,0]; K[1,2]=s*mm[0,1]; K[2,1]=s*mm[1,0]; K[2,2]=s*mm[1,1]
    return Jm+kap*K
def gv(x,kap):
    D=D3(x); return np.array([D[0]+kap, D[1], D[2]],complex)
def gvp(x,kap): return (gv(x+1e-6,kap)-gv(x-1e-6,kap))/2e-6

print("gauge-mode exactness for L_red=J3+kappa*Ks(fluid MM), signs +/-:")
for s in (+1,-1):
    for kap in [0.5,1.0,2.0,2.81]:
        errs=[]
        for x in [-3.0,-1.0,-0.3]:
            res=gvp(x,kap)-Lred(x,complex(kap),s).dot(gv(x,kap))
            errs.append(np.linalg.norm(res)/max(np.linalg.norm(gvp(x,kap)),1e-9))
        print(f"  s={s:+d} kappa={kap}: max rel resid={max(errs):.2e}")
# Also directly reverse-engineer Ks from Ks.D3=-J3[:,0], Ks[:,0]=0 (least-norm) and see its structure.
print("\nreverse-engineered Ks (from Ks.D3=-J3[:,0], Ks[:,0]=0), at x=-1:")
x=-1.0; D=D3(x); J0=Jn(x)[:,0]
# Ks 3x3 with col0=0 => 6 unknowns (cols1,2). Ks.D3 = Ks[:,1]D[1]+Ks[:,2]D[2] = -J0. Underdetermined;
# impose fluid-only (row0 of Ks =0): then Ks[0,:]=0 => -J0[0] must=0? check.
print(f"  -J3[:,0] = {(-J0).round(4)}  (row0 must be 0 if lnN row has no s-coupling)")
