# Investigate the HKA constraint (4.2) as a first integral, and the sonic eigenstructure on it.
# (4.2): C(A,N,om,V) := 1 - A + 2 om(1+(g-1)V^2)/(1-V^2) + 2 g N V om/(1-V^2) = 0.
# If dC/dx = 0 along the flow (up to a multiple of C), the physical dynamics is 3D. The complex
# sonic eigenvalues may be an artifact of the off-constraint 4th direction.
import numpy as np, sympy as sp
import hka_desing as D

g=sp.Rational(4,3)
A,N,om,V=sp.symbols('A N om V',real=True)
oV2=1-V**2
C = 1 - A + 2*om*(1+(g-1)*V**2)/oV2 + 2*g*N*V*om/oV2   # (4.2) moved to one side (=0)

# resolved slopes (symbolic, from hka_desing data)
Ax=sp.sympify(D._DATA['Ax']); Nx=sp.sympify(D._DATA['Nx'])
numOm=sp.sympify(D._DATA['numOm']); det=sp.sympify(D._DATA['det'])
numV=sp.sympify(D._DATA['numV'])
omx=numOm/det; Vx=numV/det

dC = sp.diff(C,A)*Ax + sp.diff(C,N)*Nx + sp.diff(C,om)*omx + sp.diff(C,V)*Vx
dC = sp.simplify(dC)
print("dC/dx simplified =", dC)
# is dC/dx proportional to C? compute dC/C
ratio = sp.simplify(dC/C)
print("dC/dx / C =", ratio)

# numeric check along a trajectory: constraint value should stay ~0 if we start on it.
print("\nnumeric: constraint along a launched trajectory")
V0=-0.25
Y0=D.sonic_point(V0)
Cnum=sp.lambdify((A,N,om,V), C, 'numpy')
print(f" C at sonic point (V0={V0}): {float(Cnum(*Y0)):.3e}  (should be 0 by 4.7-4.9)")
