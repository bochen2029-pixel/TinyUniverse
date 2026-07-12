# find_typo.py — compare css_symbols' derived fluid pair (en,mo) TERM BY TERM against the KHA-PRINTED
# fluid pair, both in the SAME form (multiplied out, with A',N' kept symbolic). Reveals exactly which
# term differs (the typo or the derivation error).
import sympy as sp
tau,r,th=sp.symbols('tau r theta',positive=True)
A0,N0,o0,V0=sp.symbols('A0 N0 om0 V0'); Ax,Nx,ox,Vx=sp.symbols('A_x N_x om_x V_x')

# ---- css_symbols derivation (verbatim from css_symbols.py, full-4D divT) ----
E=r/tau; a=sp.sqrt(A0); alpha=N0*a*E; W=1/sp.sqrt(1-V0**2)
rho=o0/(4*sp.pi*r**2*A0); p=rho/3
g=sp.diag(-alpha**2,a**2,r**2,r**2*sp.sin(th)**2); gi=g.inv()
coords=[tau,r,th,sp.symbols('phi')]
def Dk(expr,k):
    e=sp.diff(expr,coords[k]); fac=(-1/tau) if k==0 else ((1/r) if k==1 else 0)
    e+=(sp.diff(expr,A0)*Ax+sp.diff(expr,N0)*Nx+sp.diff(expr,o0)*ox+sp.diff(expr,V0)*Vx)*fac
    return e
def christ(l,i,j):
    ss=sp.S.Zero
    for dd in range(4):
        if gi[l,dd]==0: continue
        ss+=gi[l,dd]*(Dk(g[dd,i],j)+Dk(g[dd,j],i)-Dk(g[i,j],dd))
    return ss/2
uu=[W/alpha,W*V0/a,0,0]
T=[[(rho+p)*uu[i]*uu[j]+p*gi[i,j] for j in range(4)] for i in range(4)]
def divT(b):
    ss=sp.S.Zero
    for ai in range(4):
        ss+=Dk(T[ai][b],ai)
        for c in range(4):
            g1=christ(ai,ai,c); ss+= g1*T[c][b] if g1!=0 else 0
            g2=christ(b,ai,c);  ss+= g2*T[ai][c] if g2!=0 else 0
    return ss
Es=sp.symbols('E',positive=True)
def fin(D):
    D=sp.simplify(D*4*sp.pi*r**2); D=D.subs(tau,r/Es); return sp.expand(sp.simplify(D.subs(Es,1)))
print("deriving css en, mo (full-4D)...", flush=True)
css_en=fin(divT(0)); css_mo=fin(divT(1))

# ---- KHA-printed fluid pair (multiply through to match css's *4 pi r^2 normalization? css multiplied
# the covariant divergence by 4 pi r^2; KHA's are already in the (/omega) autonomous form). To compare,
# put both into "expression = 0" with the same overall normalization. KHA eqs are dimensionless combos.
# KHA energy eq (autonomous, _,s=0):
kha_en=( (1+N0*V0)*ox/o0 + 4*(N0+V0)*Vx/(3*(1-V0**2)) - N0*V0*Ax/(3*A0) + 4*V0*Nx/3 + 2*N0*(1+4*o0/(9*(1-V0**2))) )
kha_mo=( (4*V0+N0+3*N0*V0**2)*ox/o0 + 4*(1+V0**2+2*N0*V0)*Vx/(1-V0**2) + N0*(1-V0**2)*Ax/A0 + 4*(1+V0**2)*Nx + 2*N0*(1+3*V0**2) )

# css_en is the covariant divergence * 4 pi r^2; it should be proportional to a combination of the
# KHA energy&momentum eqs (both are consequences of div T=0). To compare cleanly, note both css_en and
# css_mo are the (t) and (r) projections; KHA's two eqs are specific linear combos. Instead compare the
# RESOLVED slopes (already done in compare_systems). Here: check if css_en is in span{kha_en,kha_mo}.
# Solve: css_en = c1*kha_en + c2*kha_mo for constants c1,c2 (rational in fields)? Try matching the
# coefficients of ox and Vx.
print("\nCoefficient comparison (of om_x, V_x, A_x, N_x, and the source) :")
def coeffs(e):
    return dict(ox=sp.simplify(e.coeff(ox)), Vx=sp.simplify(e.coeff(Vx)),
                Ax=sp.simplify(e.coeff(Ax)), Nx=sp.simplify(e.coeff(Nx)),
                src=sp.simplify(e.subs({ox:0,Vx:0,Ax:0,Nx:0})))
ce=coeffs(css_en); cm=coeffs(css_mo); ke=coeffs(kha_en); km=coeffs(kha_mo)
print("css_en coeffs:", {k:str(v) for k,v in ce.items()})
print("kha_en coeffs:", {k:str(v) for k,v in ke.items()})
# ratio css_en/kha_en per coefficient (should be a single common factor if identical up to scale)
print("\ncss_en / kha_en  per-coefficient ratio (constant across coeffs => same eq up to factor):")
for key in ['ox','Vx','Ax','Nx','src']:
    if ke[key]!=0:
        print(f"  {key}: {sp.simplify(ce[key]/ke[key])}")
    else:
        print(f"  {key}: kha=0, css={ce[key]}")
print("\ncss_mo / kha_mo  per-coefficient ratio:")
for key in ['ox','Vx','Ax','Nx','src']:
    if km[key]!=0:
        print(f"  {key}: {sp.simplify(cm[key]/km[key])}")
    else:
        print(f"  {key}: kha=0, css={cm[key]}")
