# Diagnose the fluid coefficient matrix at each series order near the sonic point.
import numpy as np, importlib
import hka_background as H

V0 = -0.25
A0,N0,om0,_ = H.sonic_data(V0)
print(f"V0={V0} A0={A0:.5f} N0={N0:.5f} om0={om0:.5f}")

# The fluid 2x2 M(om_x,V_x)=b. At the sonic point it's singular. The order-m coefficient equations
# for (o_{m+1},v_{m+1}) have coefficient matrix = the m-derivative structure. Let's print the sonic
# matrix M and the "R_d,R_e vs (oN,vN)" jacobian at m=0 and m=1 to see the rank.

# reproduce internal resid: monkeypatch by calling sonic_series with order=2 and printing.
# Instead, directly compute M_sonic:
g = H.G; oV2 = 1-V0*V0
M = np.array([[(1+N0*V0)/om0, g*(N0+V0)/oV2],[(g-1)*(N0+V0)/om0, g*(1+N0*V0)/oV2]])
print("M_sonic =", M, " det=", np.linalg.det(M))
# The om_x/V_x jacobian at [x^m] scales the derivative: O'->(m+1)o_{m+1}, so the coeff matrix of
# (o_{m+1},v_{m+1}) in (R_d[x^m],R_e[x^m]) is (m+1)*[[ (1+NV)(1-V^2), g(N+V)om ],[ (g-1)(N+V)(1-V^2), g(1+NV)om ]]
# evaluated at sonic (all higher coeffs 0). That's (m+1)*om*(1-V^2)*M_sonic-ish => STILL singular
# for all m. So the singular direction is NEVER resolved by the local [x^m] fluid eqs alone.
J = np.array([[(1+N0*V0)*oV2, g*(N0+V0)*om0],[(g-1)*(N0+V0)*oV2, g*(1+N0*V0)*om0]])
print("J (coeff of (o_{m+1},v_{m+1}) up to (m+1)) =", J, " det=", np.linalg.det(J))
print("=> singular. The solvability (Fredholm) condition at order m fixes the COMBINATION; the")
print("   free component along ker(J) is set by the NEXT order. This is why direct per-order")
print("   inversion fails. Correct approach: parametrize the singular direction and close via")
print("   the higher-order solvability, i.e. a recursive linear system with a shift.")

# ker and coker of J
U,S,Vt = np.linalg.svd(J)
print("singular values of J:", S)
print("right kernel (o,v dir):", Vt[-1])
print("left kernel (solvability dir):", U[:,-1])
