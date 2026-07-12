# Construct the reduced perturbation kappa-coupling K from first principles:
#   K[:,0]=0 (gauge kappa^2), K D3 = -J3[:,0] (gauge kappa^1), fluid block (rows 1,2 cols 1,2) = s*MinvMs.
# Check consistency of the fluid rows, and solve the lnN-row entries K[0,1],K[0,2] (1-param free ->
# fix by demanding the gate pass at a 2nd kappa). Then VALIDATE the resulting operator on the gauge
# mode at many kappa; if it passes, run the eigenvalue shoot.
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

print("Consistency: does a FLUID block s*MinvMs satisfy rows1,2 of  K.D3 = -J3[:,0]?")
for x in [-3.0,-1.0,-0.3]:
    D=D3(x); J0=Jn(x)[:,0]; mm=MMn(x)
    for s in (+1,-1):
        # fluid rows: K[1,:]=(0,s mm00, s mm01); K[1].D3 = s(mm00 D[1]+mm01 D[2])
        r1=s*(mm[0,0]*D[1]+mm[0,1]*D[2]); r2=s*(mm[1,0]*D[1]+mm[1,1]*D[2])
        print(f"  x={x:+.2f} s={s:+d}: K.D3 rows1,2 = ({r1:+.4f},{r2:+.4f})  target -J3[1:,0]=({-J0[1]:+.4f},{-J0[2]:+.4f})  diff=({r1+J0[1]:+.2e},{r2+J0[2]:+.2e})")
