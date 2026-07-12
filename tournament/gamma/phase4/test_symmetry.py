# test_symmetry.py — is css-rederived = KHA-printed under a CONSISTENT relabeling?
# Test candidate symmetries of the KHA-printed fluid slopes (om',V'):
#   (a) V->-V only:          om'(V) ?= om'_KHA(-V),  V'(V) ?= -V'_KHA(-V)
#   (b) V->-V and x->-x:     under x->-x, d/dx->-d/dx, so om'->-om', V'->-V'; combined with V->-V:
#                            css V'(V) ?= +V'_KHA(-V)*(-1)*(-1)=V'_KHA(-V)? carefully below.
#   (c) N->-N, V->-V:        (radial reflection often flips N too)
# If css == KHA under a consistent relabeling, the spectrum is preserved and 2.81 must exist.
# If css is NOT any relabeling of KHA, then css solves a DIFFERENT system (candidate root cause).
import sympy as sp, pickle
N,A,w,V=sp.symbols('N A w V')
# KHA resolved slopes
Ap=A*(1-A+(2*w/(1-V**2))*(1+V**2/3)); Np=N*(-2+A-sp.Rational(2,3)*w)
g,u=sp.symbols('g u')
f1=(1+N*V)*g+4*(N+V)*u/(3*(1-V**2))-N*V*Ap/(3*A)+4*V*Np/3+2*N*(1+4*w/(9*(1-V**2)))
f2=(4*V+N+3*N*V**2)*g+4*(1+V**2+2*N*V)*u/(1-V**2)+N*(1-V**2)*Ap/A+4*(1+V**2)*Np+2*N*(1+3*V**2)
sol=sp.solve([f1,f2],[g,u],dict=True)[0]
KHA_wp=sp.simplify(w*sol[g]); KHA_Vp=sp.simplify(sol[u]); KHA_Ap=sp.simplify(Ap); KHA_Np=sp.simplify(Np)

d=pickle.load(open("css_syms.pkl","rb")); r=sp.Symbol('r',positive=True)
A0,N0,o0,V0=sp.symbols('A0 N0 om0 V0')
sub={A0:A,N0:N,o0:w,V0:V}
CSS_wp=sp.sympify(d['omx']).subs(r,1).subs(sub); CSS_Vp=sp.sympify(d['Vxx']).subs(r,1).subs(sub)
CSS_Ap=sp.sympify(d['Axv']).subs(sub); CSS_Np=sp.sympify(d['Nxv']).subs(sub)

def chk(name, cond):
    print(f"  {name}: {'YES (0)' if cond else 'no'}")

print("Metric slopes: are css and KHA metric eqs identical? (they should be, both verified vs EC)")
chk("A': css==KHA", sp.simplify(CSS_Ap-KHA_Ap)==0)
chk("N': css==KHA", sp.simplify(CSS_Np-KHA_Np)==0)

print("\n(a) V->-V only  [css_Vp(V) == -KHA_Vp(-V), css_wp(V)==KHA_wp(-V)]:")
chk("V'", sp.simplify(CSS_Vp - (-KHA_Vp.subs(V,-V)))==0)
chk("w'", sp.simplify(CSS_wp - KHA_wp.subs(V,-V))==0)

print("\n(b) x->-x only  [css slopes == -KHA slopes]:")
chk("V'", sp.simplify(CSS_Vp + KHA_Vp)==0)
chk("w'", sp.simplify(CSS_wp + KHA_wp)==0)

print("\n(c) V->-V AND x->-x  [css_Vp(V)==+KHA_Vp(-V), css_wp(V)==-KHA_wp(-V)]:")
chk("V'", sp.simplify(CSS_Vp - KHA_Vp.subs(V,-V))==0)
chk("w'", sp.simplify(CSS_wp + KHA_wp.subs(V,-V))==0)

print("\n(d) N->-N, V->-V  [css_Vp(N,V)==-KHA_Vp(-N,-V), css_wp==KHA_wp(-N,-V)]:")
chk("V'", sp.simplify(CSS_Vp - (-KHA_Vp.subs({N:-N,V:-V})))==0)
chk("w'", sp.simplify(CSS_wp - KHA_wp.subs({N:-N,V:-V}))==0)

print("\n(e) N->-N only  [css_Vp==KHA_Vp(-N)? , A',N' behavior]:")
chk("V'", sp.simplify(CSS_Vp - KHA_Vp.subs(N,-N))==0)
chk("w'", sp.simplify(CSS_wp - KHA_wp.subs(N,-N))==0)

# Direct numeric spot check at a generic point to see how different they are
import numpy as np
fCV=sp.lambdify((N,A,w,V),CSS_Vp,'numpy'); fKV=sp.lambdify((N,A,w,V),KHA_Vp,'numpy')
print("\nNumeric V' at (N=1.3,A=1.2,w=0.4,V=0.3):")
print("  css:",fCV(1.3,1.2,0.4,0.3),"  KHA:",fKV(1.3,1.2,0.4,0.3),"  KHA(-V):",fKV(1.3,1.2,0.4,-0.3))
