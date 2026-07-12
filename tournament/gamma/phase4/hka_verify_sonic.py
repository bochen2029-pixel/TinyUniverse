# hka_verify_sonic.py — symbolic verification of the HKA Eq. (4.1) background system.
#
# Purpose (Stage A groundwork, honest): verify FROM the transcribed HKA equations that
#   (a) the sonic condition (4.5) det=0 is the sound cone  (1+NV)^2 = (gamma-1)(N+V)^2,
#   (b) the sonic data (4.7)-(4.9) parameterized by V0 satisfy the L'Hopital condition (4.6)
#       a*f - e*c = 0 AND the algebraic constraint (4.2),
#   (c) the regular-center fixed point (4.11) A=1, M=NV=-2/(3g), om=V=0 is a genuine fixed point.
#
# ALL equation numbers refer to HKA_beta_equations.md (transcribed from gr-qc/9607010).
# gamma = 4/3, sqrt(gamma-1) = 1/sqrt(3).

import sympy as sp

A, N, om, V = sp.symbols('A N om V', real=True)
Ax, Nx, omx, Vx = sp.symbols('A_x N_x om_x V_x', real=True)
g = sp.Rational(4, 3)          # gamma = 4/3 (radiation)
cs = sp.sqrt(g - 1)            # sound speed = 1/sqrt(3)

# ---- Background ODEs, HKA Eq. (4.1) verbatim ----
# (4.1a):  A_x/A = 1 - A + 2 om [1+(g-1)V^2]/(1-V^2)
rhs_Aa = 1 - A + 2*om*(1 + (g-1)*V**2)/(1 - V**2)
# (4.1b):  N_x/N = -2 + A - (2-g) om
rhs_Nb = -2 + A - (2-g)*om
# (4.1d):  (1+NV) om_x/om + g(N+V)V_x/(1-V^2) = (3(2-g)/2)NV - ((2+g)/2)ANV + (2-g)NV om
# (4.1e):  (g-1)(N+V) om_x/om + g(1+NV)V_x/(1-V^2) = (2-g)(g-1)N om + ((7g-6)/2)N + ((2-3g)/2)AN
# Write as [[a,b],[c,d]] . [omx/om, Vx/(1-V^2)]^T ? — careful with the exact variable groupings.
# HKA (4.4): the 2x2 is in (om_x, V_x). Group:
#   coefficient of om_x/om and of V_x explicitly.
# Eq (4.1d): (1+NV)*(om_x/om) + g(N+V)/(1-V^2)*V_x = RHS_d
# Eq (4.1e): (g-1)(N+V)*(om_x/om) + g(1+NV)/(1-V^2)*V_x = RHS_e
a_coef = (1 + N*V)                       # coeff of om_x/om in (4.1d)
b_coef = g*(N+V)/(1-V**2)                # coeff of V_x in (4.1d)
c_coef = (g-1)*(N+V)                     # coeff of om_x/om in (4.1e)
d_coef = g*(1+N*V)/(1-V**2)              # coeff of V_x in (4.1e)
RHS_d = (3*(2-g)/2)*N*V - ((2+g)/2)*A*N*V + (2-g)*N*V*om
RHS_e = (2-g)*(g-1)*N*om + ((7*g-6)/2)*N + ((2-3*g)/2)*A*N

# ---- (4.5) sonic condition: determinant of the 2x2 in (om_x/om, V_x) ----
detM = sp.simplify(a_coef*d_coef - b_coef*c_coef)
print("det M (4.4 in vars om_x/om, V_x) =", detM)
# HKA (4.5): det ∝ g{(1+NV)^2 - (g-1)(N+V)^2}/(1-V^2)
sound_cone = sp.simplify((1+N*V)**2 - (g-1)*(N+V)**2)
print("detM / [g*sound_cone/(1-V^2)] =", sp.simplify(detM/(g*sound_cone/(1-V**2))), " (==1 => matches 4.5)")
print("sound cone factor (1+NV)^2-(g-1)(N+V)^2 =", sp.expand(sound_cone))

# ---- (4.7)-(4.9): sonic data as functions of V0 ----
V0 = sp.symbols('V0', real=True)
N0 = (1 - cs*V0)/(cs - V0)                                                        # (4.7)
A0 = (g**2 + 4*g - 4 + 8*(g-1)**sp.Rational(3,2)*V0 - (3*g-2)*(2-g)*V0**2)/(g**2*(1-V0**2))  # (4.8)
om0 = 2*cs*(cs - V0)*(1 + cs*V0)/(g**2*(1-V0**2))                                 # (4.9)

# CHECK 1: does (N0,V0) sit on the sound cone (4.5)?  (1+N0 V0)^2 = (g-1)(N0+V0)^2
cone_at_sonic = sp.simplify(((1+N0*V0)**2 - (g-1)*(N0+V0)**2))
print("\n[CHECK] sound cone at (N0(V0),V0):", sp.simplify(cone_at_sonic), " (expect 0)")

# CHECK 2: L'Hopital / row-proportionality (4.6): a*f - e*c = 0
#   In HKA (4.4) the RHS vector is (e,f)^T = (RHS_d, RHS_e). Row-prop with singular M:
#   a_coef * RHS_e - c_coef * RHS_d = 0  (so the linear system is consistent at det=0).
lhop = (a_coef*RHS_e - c_coef*RHS_d)
lhop_sonic = sp.simplify(lhop.subs({A: A0, N: N0, om: om0, V: V0}))
print("[CHECK] L'Hopital a*RHS_e - c*RHS_d at sonic:", sp.simplify(lhop_sonic), " (expect 0)")

# CHECK 3: constraint (4.2): 1 - A + 2 om[1+(g-1)V^2]/(1-V^2) = -2 g N V om/(1-V^2)
constraint = sp.simplify((1 - A + 2*om*(1+(g-1)*V**2)/(1-V**2)) - (-2*g*N*V*om/(1-V**2)))
constraint_sonic = sp.simplify(constraint.subs({A: A0, N: N0, om: om0, V: V0}))
print("[CHECK] constraint (4.2) at sonic:", sp.simplify(constraint_sonic), " (expect 0)")

# CHECK 4: regular center fixed point (4.11): A=1, M=NV=-2/(3g), om=V=0.
#   At V->0, om->0, A->1: A_x/A ->0? N_x/N -> -2+1 = -1 (so N~e^{-x}, i.e. N->inf as x->-inf). Good.
#   Need M=NV=-2/(3g). With V->0, N->inf, product finite. Check A_x/A and om_x/om -> consistent.
print("\n[CHECK] center fixed point (4.11):")
print("  A_x/A at (A=1,om=0,V=0):", sp.simplify(rhs_Aa.subs({A:1, om:0, V:0})), " (expect 0)")
print("  N_x/N at (A=1,om=0):", sp.simplify(rhs_Nb.subs({A:1, om:0})), " (expect -1 => N~e^{-x})")
print("  M target NV = -2/(3g) =", sp.nsimplify(-2/(3*g)), "=", float(-2/(3*g)))

# ---- numeric sonic data at a trial V0 (ingoing => V0<0 on the ingoing sound cone) ----
print("\n[NUMERIC] sonic data at trial V0:")
for v0try in [sp.Rational(-1,5), sp.Rational(-3,10), sp.Rational(-1,4), sp.Rational(-2,5)]:
    n0 = float(N0.subs(V0, v0try)); a0 = float(A0.subs(V0, v0try)); o0 = float(om0.subs(V0, v0try))
    ms = 1 - 1/a0   # Misner-Sharp 2m/r
    print(f"  V0={float(v0try):+.3f}: N0={n0:+.5f} A0={a0:.5f} om0={o0:.5f}  2m/r=1-1/A0={ms:.5f}")

# sound speed identity check: at sonic, flow speed should equal cs. Flow speed in these vars:
# The characteristic condition (1+NV)^2=(g-1)(N+V)^2 => (1+NV)/|N+V| = cs. Report.
print("\ncs = sqrt(gamma-1) =", float(cs), " (=1/sqrt(3))")
