# kha_verbatim.py — implement the KHA-PRINTED autonomous EOM VERBATIM (gr-qc/9503007 eq. EOM),
# with the [VERIFY-BRACES] interpretation, and test: (1) center regularity, (2) sonic point.
# Compare against css_symbols' re-derived full-4D system to see if the prior 'fix' changed the physics.
#
# KHA autonomous (set _,s=0):
#  metric:  A'/A = 1 - A + (2w/(1-V^2))(1+V^2/3)
#           N'/N = -2 + A - 2w/3
#  fluid1 (energy):  (1+NV) w'/w + 4(N+V)V'/(3(1-V^2)) - NV A'/(3A) + 4V N'/3 + 2N(1+4w/(9(1-V^2))) = 0
#  fluid2 (mom):     (4V+N+3NV^2) w'/w + 4(1+V^2+2NV)V'/(1-V^2) + N(1-V^2)A'/A + 4(1+V^2)N' + 2N(1+3V^2) = 0
# Solve fluid1,fluid2 (2x2) for (w'/w, V') given A',N' from metric. Then w' = w*(w'/w).
import sympy as sp, numpy as np, math
from scipy.optimize import brentq

N,A,w,V = sp.symbols('N A w V', real=True)
Ap = A*(1 - A + (2*w/(1-V**2))*(1+V**2/3))       # A'
Np = N*(-2 + A - sp.Rational(2,3)*w)              # N'
# unknowns g=w'/w, u=V'
g,u = sp.symbols('g u', real=True)
f1 = (1+N*V)*g + 4*(N+V)*u/(3*(1-V**2)) - N*V*Ap/(3*A) + 4*V*Np/3 + 2*N*(1+4*w/(9*(1-V**2)))
f2 = (4*V+N+3*N*V**2)*g + 4*(1+V**2+2*N*V)*u/(1-V**2) + N*(1-V**2)*Ap/A + 4*(1+V**2)*Np + 2*N*(1+3*V**2)
sol = sp.solve([f1,f2],[g,u],dict=True)[0]
gexpr = sp.simplify(sol[g]); uexpr = sp.simplify(sol[u])
wp = sp.simplify(w*gexpr)                          # w'
# determinant of the 2x2 (denominator) => sonic locus
M = sp.Matrix([[ (1+N*V), 4*(N+V)/(3*(1-V**2)) ],
               [ (4*V+N+3*N*V**2), 4*(1+V**2+2*N*V)/(1-V**2) ]])
detM = sp.simplify(M.det())
print("KHA fluid 2x2 determinant (num factored):")
print("  det =", sp.factor(sp.numer(sp.together(detM))))
# evaluate at exact sonic state N=2/sqrt3, A=3/2, w=3/4, V=1/sqrt3
S3=sp.sqrt(3)
pt={N:2/S3, A:sp.Rational(3,2), w:sp.Rational(3,4), V:1/S3}
print("  det at exact sonic (V=+1/sqrt3):", sp.simplify(detM.subs(pt)))
ptm={N:2/S3, A:sp.Rational(3,2), w:sp.Rational(3,4), V:-1/S3}
print("  det at V=-1/sqrt3 (ingoing):     ", sp.simplify(detM.subs(ptm)))

f_Ap=sp.lambdify((N,A,w,V),Ap,'math'); f_Np=sp.lambdify((N,A,w,V),Np,'math')
f_wp=sp.lambdify((N,A,w,V),wp,'math'); f_up=sp.lambdify((N,A,w,V),uexpr,'math')
f_det=sp.lambdify((N,A,w,V),detM,'math')

def slopes(Y):
    n,a,om,v=Y
    return (f_Np(n,a,om,v), f_Ap(n,a,om,v), f_wp(n,a,om,v), f_up(n,a,om,v))

# center seed: try BOTH outgoing (V=+z/(2nc)) and ingoing (V=-z/(2nc)) and see which reaches sonic.
def seed(nc,a2,z0,sgn):
    return (nc/z0, 1.0+a2*z0*z0, 1.5*a2*z0*z0, sgn*z0/(2*nc))

def rk4(Y,h):
    k1=slopes(Y); k2=slopes(tuple(Y[i]+0.5*h*k1[i] for i in range(4)))
    k3=slopes(tuple(Y[i]+0.5*h*k2[i] for i in range(4))); k4=slopes(tuple(Y[i]+h*k3[i] for i in range(4)))
    return tuple(Y[i]+(h/6)*(k1[i]+2*k2[i]+2*k3[i]+k4[i]) for i in range(4))

def run(nc,a2,z0,h,sgn):
    Y=seed(nc,a2,z0,sgn); x=math.log(z0); prevD=f_det(*Y); prevY=Y; prevx=x
    n=int((4.0-x)/h)
    for i in range(n):
        try: Y2=rk4(Y,h)
        except Exception: return ("err",x,Y)
        if not all(math.isfinite(v) for v in Y2): return ("blow",x,Y)
        D2=f_det(*Y2)
        if prevD*D2<0 and abs(Y2[3])>1e-3:
            fr=prevD/(prevD-D2); Yh=tuple(prevY[k]+fr*(Y2[k]-prevY[k]) for k in range(4))
            return ("sonic",prevx+fr*h,Yh)
        prevD=D2;prevY=Y2;prevx=x+h;Y=Y2;x+=h
    return ("nocross",x,Y)

print("\nCenter dV/dx behavior (KHA-printed system) as N->inf (A=1,V=0,w->0):")
for Nc in [1e4,1e5,1e6]:
    v=f_up(Nc,1.0,1e-12,0.0)
    print(f"  N={Nc:.0e}: dV/dx={v:.4e}  (regular center needs ->0)")

print("\nShoot from center (both signs), a2 scan, does it reach the exact sonic state?")
for sgn in [+1,-1]:
    print(f" sign(V_center)={sgn:+d}:")
    for a2 in [0.20,0.25,0.2500679,0.30,0.35]:
        st,xh,Y=run(1.0,a2,math.exp(-12),1e-4,sgn)
        print(f"   a2={a2:.7f}: {st:>8} x={xh:8.4f} N={Y[0]:.4f} A={Y[1]:.4f} w={Y[2]:.4f} V={Y[3]:.4f}")
