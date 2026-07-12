# Generate clean, r-free numerical RHS + sonic-regularity numerators from the derived Euler
# ODEs, and pickle fast callables. This is the authoritative background dynamics module.
import sympy as sp, pickle
en_s,mo_s=pickle.load(open("fluid_manual.pkl","rb"))
en=sp.sympify(en_s); mo=sp.sympify(mo_s)
A0,N0,o0,V0,Ax,Nx,ox,Vx,r=sp.symbols('A0 N0 om0 V0 A_x N_x om_x V_x r')
Axv=A0*(1-A0+(2*o0/(1-V0**2))*(1+V0**2/3)); Nxv=N0*(-2+A0-sp.Rational(2,3)*o0)
e1=sp.expand(en.subs({Ax:Axv,Nx:Nxv})); e2=sp.expand(mo.subs({Ax:Axv,Nx:Nxv}))
a11=e1.coeff(ox);a12=e1.coeff(Vx);a21=e2.coeff(ox);a22=e2.coeff(Vx)
b1=-(e1.subs({ox:0,Vx:0}));b2=-(e2.subs({ox:0,Vx:0}))
det=sp.simplify(a11*a22-a12*a21)
omx=sp.simplify((b1*a22-a12*b2)/det)   # = dom/dx
Vxx=sp.simplify((a11*b2-b1*a21)/det)   # = dV/dx
# regularity numerators (r-free): numerator of numV and numU
numV=sp.together(a11*b2-b1*a21); numU=sp.together(b1*a22-a12*b2)
nV=sp.simplify(sp.numer(numV)); nU=sp.simplify(sp.numer(numU))
# strip r (overall) - substitute r=1 (it's an overall factor in num & denom already cancelled in ratios)
# For omx,Vxx r should already be gone; verify:
print("omx has r?",omx.has(r)," Vxx has r?",Vxx.has(r))
# build callables (set r=1 harmlessly)
f_om=sp.lambdify((A0,N0,o0,V0),omx.subs(r,1),'numpy',cse=True)
f_V =sp.lambdify((A0,N0,o0,V0),Vxx.subs(r,1),'numpy',cse=True)
# regularity: sonic-point analytic passage needs BOTH numerators of the 0/0 to vanish. But at
# det=0, Vxx=nV/det. For finite Vxx we need nV->0 as det->0. The regularity residual we shoot to
# zero is nV evaluated at the sonic locus. Provide nV,nU as r-free polys:
f_nV=sp.lambdify((A0,N0,o0,V0),sp.simplify(nV.subs(r,1)),'numpy',cse=True)
f_nU=sp.lambdify((A0,N0,o0,V0),sp.simplify(nU.subs(r,1)),'numpy',cse=True)
f_det=sp.lambdify((A0,N0,o0,V0),sp.simplify(det.subs(r,1)),'numpy',cse=True)
# also the metric slopes
# save source strings for the C++ port later + pickle callables via cloudpickle-free: save srepr
data={
 'omx':sp.srepr(omx.subs(r,1)),'Vxx':sp.srepr(Vxx.subs(r,1)),
 'nV':sp.srepr(sp.simplify(nV.subs(r,1))),'nU':sp.srepr(sp.simplify(nU.subs(r,1))),
 'det':sp.srepr(sp.simplify(det.subs(r,1))),
}
pickle.dump(data,open("rhs_exprs.pkl","wb"))
# quick center check
import numpy as np
print("center dV/dx (A=1,V=0,om->0,N=1e5):", float(f_V(1.0,1e5,1e-9,0.0)))
print("sonic locus check: nV,nU at a regular sonic guess:")
# Also print omx,Vxx in plain form for C++ (expanded, common-denominator)
print("\n--- omx (dom/dx) ---"); print(sp.ccode(sp.simplify(omx.subs(r,1))))
print("\n--- Vxx (dV/dx) ---"); print(sp.ccode(sp.simplify(Vxx.subs(r,1))))
