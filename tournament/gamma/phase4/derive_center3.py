# Center leading-order relations via .coeff (avoids Poly on 1/z expressions).
# Ansatz near center z=e^x->0:  N=m0/z + m1 + ..., A=1+a2 z^2+..., V=v1 z+v2 z^2+...,
# omega=w2 z^2 + w3 z^3 + ...  Determine relations among {m0,m1,a2,v1,v2,w2,w3,...}.
import sympy as sp
z=sp.symbols('z', positive=True)
m0,m1,m2=sp.symbols('m0 m1 m2')
a2,a3,a4=sp.symbols('a2 a3 a4')
v1,v2,v3=sp.symbols('v1 v2 v3')
w2,w3,w4=sp.symbols('w2 w3 w4')
Nser=m0/z+m1+m2*z
Aser=1+a2*z**2+a3*z**3+a4*z**4
Vser=v1*z+v2*z**2+v3*z**3
Wser=w2*z**2+w3*z**3+w4*z**4
def dx(F): return z*sp.diff(F,z)

Ap_ode=Aser*(1-Aser+(2*Wser/(1-Vser**2))*(1+Vser**2/3))
Np_ode=Nser*(-2+Aser-sp.Rational(2,3)*Wser)
N,A,om,V=Nser,Aser,Wser,Vser
omV2=1-V**2
m11=(1+N*V);m12=4*(N+V)/(3*omV2);m21=(4*V+N+3*N*V**2);m22=4*(1+V**2+2*N*V)/omV2
rhs1=-(-N*V*Ap_ode/(3*A)+sp.Rational(4,3)*V*Np_ode+2*N*(1+4*om/(9*omV2)))
rhs2=-( N*omV2*Ap_ode/A+4*(1+V**2)*Np_ode+2*N*(1+3*V**2))
det=m11*m22-m12*m21
u=(rhs1*m22-m12*rhs2)/det
Vp=(m11*rhs2-rhs1*m21)/det
omp=u*om

resA=sp.expand(sp.series(dx(Aser)-Ap_ode,z,0,5).removeO())
resN=sp.expand(sp.series(dx(Nser)-Np_ode,z,0,4).removeO())
resV=sp.expand(sp.series(dx(Vser)-Vp,z,0,5).removeO())
resW=sp.expand(sp.series(dx(Wser)-omp,z,0,6).removeO())

eqs=[]
for R,lo,hi in [(resA,0,5),(resN,-1,4),(resV,0,5),(resW,0,6)]:
    for n in range(lo,hi):
        c=sp.simplify(R.coeff(z,n))
        if c!=0: eqs.append(c)
eqs=[e for e in eqs if e!=0]
print(f"{len(eqs)} equations. Solving...")
unk=[m0,m1,m2,a2,a3,a4,v2,v3,w2,w3,w4]  # treat v1 as free (shoot param)
sol=sp.solve(eqs,unk,dict=True)
print("branches:",len(sol))
for s in sol:
    d={str(k):sp.simplify(v) for k,v in s.items()}
    print("  ---")
    for k in ['m0','m1','m2','a2','a3','a4','v2','v3','w2','w3','w4']:
        if k in d: print(f"    {k} = {sp.nsimplify(d[k],[sp.sqrt(3)])}")
    free=[str(x) for x in unk if x not in s]
    print("    free (besides v1):",free)
