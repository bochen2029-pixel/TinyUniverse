# Find the reduced perturbation operator L_red = J3 + kappa*Ks with FLUID-ONLY Ks (rows lnom,V get
# +/- MinvMs) by testing DIFFERENT gauge-mode forms. The gauge shift structure (which component carries
# +kappa) is the unknown. Test psi_g = D3 + kappa*e_g for e_g in {e_lnN, -e_lnN, e_V, ...}.
import numpy as np, sympy as sp, math, pickle
import hka_beta4 as B
B.bg(); B.bg_path(); G=4.0/3.0
d=pickle.load(open("hka_pert_reduced.pkl","rb")); lN,lo,V=sp.symbols('lN lo V',real=True)
J=sp.sympify(d['J']); fJ=[[sp.lambdify((lN,lo,V),J[i,j],'numpy') for j in range(3)] for i in range(3)]
g=sp.Rational(4,3); oV2=1-V**2; N=sp.exp(lN); om=sp.exp(lo)
MM=sp.Matrix([[1+N*V,g*(N+V)/oV2],[(g-1)*(N+V),g*(1+N*V)/oV2]]).inv()*sp.Matrix([[1,g*V/oV2],[(g-1)*V,g/oV2]])
fMM=[[sp.lambdify((lN,lo,V),MM[i,j],'numpy') for j in range(2)] for i in range(2)]
def y3(x):
    A,N,om,V,ox,Vx=B.bg_state(x); return (math.log(N),math.log(om),V)
def D3(x):
    A,N,om,V,ox,Vx=B.bg_state(x); return np.array([N*(-2+A-(2-G)*om)/N,ox,Vx])
def Jn(x): return np.array([[fJ[i][j](*y3(x)) for j in range(3)] for i in range(3)])
def MMn(x): return np.array([[fMM[i][j](*y3(x)) for j in range(2)] for i in range(2)])
def Lred(x,kap,s=+1.0):
    K=np.zeros((3,3)); mm=MMn(x); K[1,1]=s*mm[0,0];K[1,2]=s*mm[0,1];K[2,1]=s*mm[1,0];K[2,2]=s*mm[1,1]
    return Jn(x)+kap*K

gauges={'e_lnN':np.array([1,0,0.]),'-e_lnN':np.array([-1,0,0.]),'e_lnom':np.array([0,1,0.]),
        'e_V':np.array([0,0,1.]),'lnN+lnom':np.array([1,1,0.])}
print("gauge-mode exactness for L_red=J3 +/- kappa*MM(fluid), various gauge-shift directions e_g:")
for s in (+1,-1):
    for gname,eg in gauges.items():
        # find best over a couple kappa; require exact for all kappa
        ok=True; worst=0
        for kap in [0.5,1.0,2.81]:
            for x in [-2.0,-0.5]:
                def gv(xx): D=D3(xx); return (D+kap*eg).astype(complex)
                dP=(gv(x+1e-6)-gv(x-1e-6))/2e-6
                res=dP-Lred(x,complex(kap),s).dot(gv(x))
                r=np.linalg.norm(res)/max(np.linalg.norm(dP),1e-9); worst=max(worst,r)
        print(f"  s={s:+d} e_g={gname:9s}: worst rel resid={worst:.2e}")
