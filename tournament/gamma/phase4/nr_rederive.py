# nr_rederive.py — re-derive the fluid perturbation coefficients (M_s: As,Bs,Cs,Ds; Gmat fluid rows
# E1..E4,F1..F4) FROM FIRST PRINCIPLES by symbolically linearizing the primary self-similar EOM
# (rflanl.tex 4.1, 4.2) about the background, keeping d_s -> kappa. Compare to hka_pert_hka99 (=§V).
# A mismatch localizes a §V transcription typo (the operator matches §V verbatim but lacks kappa=2.81).
import sympy as sp

g = sp.Rational(4, 3)
e = sp.symbols('e')                                   # perturbation bookkeeping epsilon
kap = sp.symbols('kappa')
# background (ss) values + their x log-slopes (Abar_x = dlnA/dx etc.), V_x = dV/dx
As, Ns, oms, Vs = sp.symbols('As Ns oms Vs', positive=True)
Abx, Nbx, obx, Vx = sp.symbols('Abx Nbx obx Vx', real=True)
# perturbation p-components (Abar_p,Nbar_p,obar_p,Vp) and their x-derivatives
Ap, Np, op, Vp = sp.symbols('Ap Np op Vp', real=True)
Apx, Npx, opx, Vpx = sp.symbols('Apx Npx opx Vpx', real=True)

# fields with perturbation: Abar=lnA so A=As*(1+e*Ap); V=Vs+e*Vp
A = As*(1 + e*Ap); N = Ns*(1 + e*Np); om = oms*(1 + e*op); V = Vs + e*Vp
oV2 = 1 - V**2
# log-slopes d_x(bar field): dlnA/dx = Abx + e*Apx ; V_x = Vx + e*Vpx ; d_s(bar)= e*kappa*(p)
obar_s = e*kap*op;      V_s = e*kap*Vp
obar_x = obx + e*opx;   V_xf  = Vx + e*Vpx

RHS_d = sp.Rational(3,1)*(2-g)/2*N*V - (2+g)/2*A*N*V + (2-g)*N*V*om
RHS_e = -(g-2)*(g-1)*N*om + (7*g-6)/2*N + (2-3*g)/2*A*N

# eq 4.1 and 4.2 (log form: om_s/om=obar_s, om_x/om=obar_x), written as (LHS - RHS) = 0
eq1 = obar_s + g*V*V_s/oV2 + (1+N*V)*obar_x + g*(N+V)*V_xf/oV2 - RHS_d
eq2 = (g-1)*V*obar_s + g*V_s/oV2 + (g-1)*(N+V)*obar_x + g*(1+N*V)*V_xf/oV2 - RHS_e

# O(e^1) linearized equations
lin1 = sp.expand(sp.series(eq1, e, 0, 2).removeO().coeff(e, 1))
lin2 = sp.expand(sp.series(eq2, e, 0, 2).removeO().coeff(e, 1))

def extract(lin, name):
    # collect: M_s * kappa * (op,Vp) + M_x*(opx,Vpx) + (Gmat fluid row).(Ap,Np,op,Vp) = 0
    # => rearranged to  M_x (opx,Vpx) = (Gmat - kappa M_s)(...)  matching §V's  M_x Psi' = (Gmat-k Ms)Psi
    As_ = lin.coeff(kap).coeff(op); Bs_ = lin.coeff(kap).coeff(Vp)  # kappa-coupling (=M_s row, sign +)
    Ax_ = lin.coeff(opx); Bx_ = lin.coeff(Vpx)                      # d/dx coupling (=M_x row)
    lin0 = sp.expand(lin - kap*(As_*op+Bs_*Vp) - Ax_*opx - Bx_*Vpx) # algebraic remainder
    # Gmat row = -(algebraic remainder coeffs), since eq = M_x Psi' - Gmat Psi + kappa M_s Psi = 0
    G_A = sp.simplify(-lin0.coeff(Ap)); G_N = sp.simplify(-lin0.coeff(Np))
    G_o = sp.simplify(-lin0.coeff(op)); G_V = sp.simplify(-lin0.coeff(Vp))
    return dict(Ms_o=sp.simplify(As_), Ms_V=sp.simplify(Bs_), Mx_o=sp.simplify(Ax_), Mx_V=sp.simplify(Bx_),
                G_A=G_A, G_N=G_N, G_o=G_o, G_V=G_V)

r1 = extract(lin1, 'eq1'); r2 = extract(lin2, 'eq2')

# §V values (hka_pert_hka99), as sympy in the same symbols
def sv():
    oV2b = 1 - Vs**2
    return dict(
        As=sp.Integer(1), Bs=g*Vs/oV2b, Cs=(g-1)*Vs, Ds=g/oV2b,
        Ax=1+Ns*Vs, Bx=g*(Ns+Vs)/oV2b, Cx=(g-1)*(Ns+Vs), Dx=g*(1+Ns*Vs)/oV2b,
        E1=-(g+2)/2*As*Ns*Vs,
        E2=(6-3*g)/2*Ns*Vs-(2+g)/2*As*Ns*Vs+(2-g)*oms*Ns*Vs-Ns*Vs*obx-g*Ns*Vx/oV2b,
        E3=(2-g)*oms*Ns*Vs,
        E4=(6-3*g)/2*Ns-(2+g)/2*As*Ns+(2-g)*oms*Ns-Ns*obx-g*(1+2*Ns*Vs+Vs**2)*Vx/oV2b**2,
        F1=(2-3*g)/2*As*Ns,
        F2=(2-g)*(g-1)*oms*Ns+(7*g-6)/2*Ns+(2-3*g)/2*As*Ns-(g-1)*Ns*obx-g*Ns*Vs*Vx/oV2b,
        F3=(2-g)*(g-1)*oms*Ns,
        F4=-(g-1)*obx-g*(Ns+2*Vs+Ns*Vs**2)*Vx/oV2b**2)
S = sv()

# NOTE: §V uses As (background A) in E1..F4; here As symbol == background A. Map re-derived to §V names.
pairs = [
    ('As (Ms row1,om)', r1['Ms_o'], S['As']), ('Bs (Ms row1,V)', r1['Ms_V'], S['Bs']),
    ('Cs (Ms row2,om)', r2['Ms_o'], S['Cs']), ('Ds (Ms row2,V)', r2['Ms_V'], S['Ds']),
    ('Ax', r1['Mx_o'], S['Ax']), ('Bx', r1['Mx_V'], S['Bx']),
    ('Cx', r2['Mx_o'], S['Cx']), ('Dx', r2['Mx_V'], S['Dx']),
    ('E1', r1['G_A'], S['E1']), ('E2', r1['G_N'], S['E2']), ('E3', r1['G_o'], S['E3']), ('E4', r1['G_V'], S['E4']),
    ('F1', r2['G_A'], S['F1']), ('F2', r2['G_N'], S['F2']), ('F3', r2['G_o'], S['F3']), ('F4', r2['G_V'], S['F4']),
]
print("re-derived (from EOM 4.1/4.2)  vs  §V transcription (hka_pert_hka99):")
for name, red, sv_ in pairs:
    diff = sp.simplify(red - sv_)
    tag = 'OK' if diff == 0 else '  <<<<< MISMATCH'
    print(f"  {name:20s}: diff = {diff}{tag}")
    if diff != 0:
        print(f"       re-derived = {sp.simplify(red)}")
        print(f"       §V value   = {sp.simplify(sv_)}")
