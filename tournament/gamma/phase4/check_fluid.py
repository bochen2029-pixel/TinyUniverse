# Definitive check: derive the self-similar reduction of the fluid conservation laws and
# compare to KHA's two fluid equations. Rather than reproduce EC's exact eqs, use the
# covariant conservation nabla_a T^{ab}=0 for a perfect fluid in the KHA metric, projected
# parallel (u_a) and orthogonal to u^a, expressed in KHA's (s,x,N,A,omega,V) variables.
#
# This is the primary derivation. If KHA's fluid eqs match, my transcription+the paper are
# right and the center issue is my analysis. If not, we find the typo.
#
# Setup: metric ds^2=-alpha^2 dt^2 + a^2 dr^2 + r^2 dOmega^2. 4-velocity u^t=a^{-1}? Wait KHA:
#   u_t = a(1-V^2)^{-1/2}?? They wrote u_t = a(1-V^2)^{-1/2}. Hmm that's u_t (lower). And
#   u^r = a V (1-V^2)^{-1/2}. Let's use the standard: u^a=(u^t,u^r,0,0),
#   u^t = 1/(alpha sqrt(1-V^2)), u^r = V/(alpha)*... let's define via V=3-velocity:
#   physical velocity v = (a/alpha) dr/dt = V. Then u^t = 1/(alpha sqrt(1-V^2)),
#   u^r = (V/a) * ... Actually proper: u^mu = gamma_L (1/alpha, V/a, 0,0)?? Let's get u_a u^a=-1.
#   Take u^t = W/alpha, u^r = W V/a with W=1/sqrt(1-V^2). Check: g_tt (u^t)^2+g_rr(u^r)^2
#   = -alpha^2 W^2/alpha^2 + a^2 W^2 V^2/a^2 = -W^2 + W^2 V^2 = -W^2(1-V^2) = -1. GOOD.
#
# T^{ab} = (rho+p) u^a u^b + p g^{ab}, p=rho/3.
# Conservation: nabla_a T^{ab}=0. Two independent components (t and r).
# We assume self-similar: everything function of x=ln(-r/t) (and s=-ln(-t)); reduce with
# _,s=0 for the CSS solution -> ODEs in x. Compare to KHA.
#
# This is algebraically heavy; use sympy with explicit metric and Christoffels.
import sympy as sp

t,r,th,ph=sp.symbols('t r theta phi', real=True)
# fields as functions of (t,r)
alpha=sp.Function('alpha')(t,r)
a=sp.Function('a')(t,r)
rho=sp.Function('rho')(t,r)
Vf=sp.Function('V')(t,r)
p=rho/3
coords=[t,r,th,ph]
g=sp.diag(-alpha**2, a**2, r**2, r**2*sp.sin(th)**2)
gi=g.inv()
n=4
# Christoffel
def christ(g,gi,coords):
    n=len(coords)
    Gamma=[[[0]*n for _ in range(n)] for _ in range(n)]
    for l in range(n):
        for m in range(n):
            for k in range(n):
                s=0
                for d in range(n):
                    s+=gi[l,d]*(sp.diff(g[d,m],coords[k])+sp.diff(g[d,k],coords[m])-sp.diff(g[m,k],coords[d]))
                Gamma[l][m][k]=sp.simplify(s/2)
    return Gamma
Gamma=christ(g,gi,coords)
W=1/sp.sqrt(1-Vf**2)
u=[W/alpha, W*Vf/a, 0, 0]           # u^a (contravariant)
# T^{ab}
T=[[ (rho+p)*u[i]*u[j] + p*gi[i,j] for j in range(n)] for i in range(n)]
# nabla_a T^{ab} = partial_a T^{ab} + Gamma^a_{ac} T^{cb} + Gamma^b_{ac} T^{ac}
def divT(b):
    s=0
    for aidx in range(n):
        s+=sp.diff(T[aidx][b], coords[aidx])
        for c in range(n):
            s+=Gamma[aidx][aidx][c]*T[c][b]
            s+=Gamma[b][aidx][c]*T[aidx][c]
    return sp.simplify(s)

print("computing conservation components (may take ~30s)...")
Dt=divT(0)   # b=t
Dr=divT(1)   # b=r
# Project: parallel u_b D^b (energy) and orthogonal. But easier: Dt and Dr are the two eqs.
# Now impose self-similarity. x=ln(-r/t), s=-ln(-t). For CSS: alpha,a,rho,V functions of x
# only in the reduced sense. Use KHA reduced vars: N=alpha/a * e^{-x}, A=a^2, omega=4 pi r^2 a^2 rho.
# We convert Dt,Dr into these. This is where it gets long. Save Dt,Dr to check structure first.
print("Dt has", sp.count_ops(Dt),"ops; Dr has", sp.count_ops(Dr),"ops")
sp.pprint(sp.nsimplify(Dt))
