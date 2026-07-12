# Convert nabla_a T^{ab}=0 to KHA self-similar variables and compare to KHA's fluid eqs.
# Strategy: substitute the SELF-SIMILAR ansatz directly. Since CSS: alpha,a,rho,V depend on x
# only (with the correct t,r scaling). KHA vars: N=alpha a^{-1} e^{-x}, A=a^2, omega=4pi r^2 a^2 rho.
# From x=ln(-r/t): e^x=-r/t, so r=-t e^x, e^{-x}=-t/r.
# Scaling of fields (homothety): a=a(x) (dimensionless), alpha=alpha(x) (dimensionless),
#   rho scales as rho= R(x)/r^2 (since omega=4pi r^2 a^2 rho finite => rho ~ 1/r^2). Actually
#   omega(x) finite and =4pi r^2 a^2 rho => rho = omega/(4pi r^2 a^2). So rho depends on r explicitly.
#   V=V(x).
# So the fields as functions of (t,r): a=a(x(t,r)), alpha=alpha(x), V=V(x),
#   rho=omega(x)/(4 pi r^2 a(x)^2).
# Substitute into Dt, Dr; everything reduces to functions of x times powers of r,t. Collect.
import sympy as sp

t,r=sp.symbols('t r', positive=True)   # use t>0 formally; sign handled by e^x=-r/t. We'll keep -t.
x=sp.symbols('x', real=True)
# represent x=ln(-r/t). Work with tt=-t>0 so e^x=r/tt.
tt=sp.symbols('tt', positive=True)     # tt=-t
A=sp.Function('A')(x); N=sp.Function('N')(x); om=sp.Function('om')(x); V=sp.Function('V')(x)
# a=sqrt(A); alpha = N*a*e^x = N sqrt(A) (r/tt). rho=om/(4 pi r^2 A).
a_=sp.sqrt(A)
ex=r/tt
alpha_=N*a_*ex
rho_=om/(4*sp.pi*r**2*A)
V_=V
# x as function of (r,tt): x=ln(r/tt). dx/dr=1/r, dx/dtt=-1/tt (since d ln(r/tt)/dtt=-1/tt),
#   note t=-tt so d/dt=-d/dtt.
def Dx_dr(): return sp.Rational(1,1)/r
def Dx_dt(): return sp.Rational(1,1)/tt   # d/dt of ln(-r/t): t=-tt, x=ln(r/tt), dx/dt=dx/dtt*dtt/dt
# dtt/dt=-1 => dx/dt = (-1/tt)*(-1)=1/tt.
# Build the conservation eqs from scratch in these variables is cleaner than converting the
# giant expression. Recompute Dt,Dr with these explicit substitutions.
th,ph=sp.symbols('theta phi')
coords=[t,r,th,ph]
# metric with explicit field(x(t,r))
# express alpha,a,rho,V as functions of t,r via x=ln(-r/t). Use chain rule through X=x.
# Simpler: define X=ln(-r/t) as sympy expression and compose.
X=sp.log(-r/t)
Af=A.subs(x,X); Nf=N.subs(x,X); omf=om.subs(x,X); Vf=V.subs(x,X)
a_f=sp.sqrt(Af); alpha_f=Nf*a_f*(-r/t); rho_f=omf/(4*sp.pi*r**2*Af)
g=sp.diag(-alpha_f**2, a_f**2, r**2, r**2*sp.sin(th)**2)
gi=g.inv()
n=4
def christ(g,gi,coords):
    n=len(coords); G=[[[0]*n for _ in range(n)] for _ in range(n)]
    for l in range(n):
        for m in range(n):
            for k in range(m,n):
                ss=0
                for d in range(n):
                    if gi[l,d]==0: continue
                    ss+=gi[l,d]*(sp.diff(g[d,m],coords[k])+sp.diff(g[d,k],coords[m])-sp.diff(g[m,k],coords[d]))
                val=sp.simplify(ss/2); G[l][m][k]=val; G[l][k][m]=val
    return G
print("Christoffels...")
G=christ(g,gi,coords)
W=1/sp.sqrt(1-Vf**2)
u=[W/alpha_f, W*Vf/a_f,0,0]
pf=rho_f/3
T=[[ (rho_f+pf)*u[i]*u[j]+pf*gi[i,j] for j in range(n)] for i in range(n)]
def divT(b):
    ss=0
    for ai in range(n):
        ss+=sp.diff(T[ai][b],coords[ai])
        for c in range(n):
            if G[ai][ai][c]!=0: ss+=G[ai][ai][c]*T[c][b]
            if G[b][ai][c]!=0: ss+=G[b][ai][c]*T[ai][c]
    return ss
print("divergence t...")
Dt=divT(0)
print("divergence r...")
Dr=divT(1)
# Now substitute back X->x and derivatives. Replace Derivative(A(log(-r/t)),t) etc.
# Let Ap=A'(x) etc. Use: d/dt f(X)=f'(X)*dX/dt = f'(X)*(1/t)? dX/dt = d ln(-r/t)/dt = -1/t?
#  ln(-r/t)=ln(-r)-ln(t)? for t<0 messy. Just let sympy differentiate log(-r/t):
dXdt=sp.diff(X,t); dXdr=sp.diff(X,r)
print("dX/dt=",dXdt," dX/dr=",dXdr)
# substitute A(X)->A(x) as independent, and its derivatives:
Ax,Nx,ox,Vx=sp.symbols("A_x N_x om_x V_x")
A0,N0,o0,V0=sp.symbols("A0 N0 om0 V0")
def reduce_expr(E):
    E=sp.simplify(E)
    # replace derivatives
    E=E.subs({sp.Derivative(Af,t):Ax*dXdt, sp.Derivative(Nf,t):Nx*dXdt, sp.Derivative(omf,t):ox*dXdt, sp.Derivative(Vf,t):Vx*dXdt,
              sp.Derivative(Af,r):Ax*dXdr, sp.Derivative(Nf,r):Nx*dXdr, sp.Derivative(omf,r):ox*dXdr, sp.Derivative(Vf,r):Vx*dXdr})
    E=E.subs({Af:A0,Nf:N0,omf:o0,Vf:V0})
    return sp.simplify(E)
Dt_r=reduce_expr(Dt); Dr_r=reduce_expr(Dr)
print("\n--- Dt reduced ---"); sp.pprint(sp.simplify(Dt_r))
print("\n--- Dr reduced ---"); sp.pprint(sp.simplify(Dr_r))
# save for comparison
import pickle
with open("fluid_reduced.pkl","wb") as f: pickle.dump((str(Dt_r),str(Dr_r)),f)
