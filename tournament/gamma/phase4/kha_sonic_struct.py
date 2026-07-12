# kha_sonic_struct.py — KHA-VERBATIM system: derive the exact sonic-point background series at the
# INGOING sonic point V=-1/sqrt3, and the perturbation operator's Fuchsian residue (indicial exps).
# This tells us whether the KHA operator supports the physical mode kappa~2.81.
import sympy as sp, numpy as np, math, pickle
S3 = sp.sqrt(3)

# ---- KHA background slopes (resolved 2x2) ----
N, A, w, V = sp.symbols('N A w V', real=True)
Ap = A*(1 - A + (2*w/(1-V**2))*(1+V**2/3))
Np = N*(-2 + A - sp.Rational(2,3)*w)
g, u = sp.symbols('g u')
f1 = (1+N*V)*g + 4*(N+V)*u/(3*(1-V**2)) - N*V*Ap/(3*A) + 4*V*Np/3 + 2*N*(1+4*w/(9*(1-V**2)))
f2 = (4*V+N+3*N*V**2)*g + 4*(1+V**2+2*N*V)*u/(1-V**2) + N*(1-V**2)*Ap/A + 4*(1+V**2)*Np + 2*N*(1+3*V**2)
sol = sp.solve([f1,f2],[g,u],dict=True)[0]
wp_expr = sp.simplify(w*sol[g]); Vp_expr = sp.simplify(sol[u])   # w', V' away from sonic
detM = -4*(3*N**2*V**2 - N**2 + 4*N*V - V**2 + 3)                # 2x2 determinant (num)

# sonic point (ingoing): N0=2/sqrt3, A0=3/2, w0=3/4, V0=-1/sqrt3
pt = {N: 2/S3, A: sp.Rational(3,2), w: sp.Rational(3,4), V: -1/S3}
Ap0 = sp.simplify(Ap.subs(pt)); Np0 = sp.simplify(Np.subs(pt))
print("KHA sonic (V=-1/sqrt3) metric slopes: A'=%s  N'=%s" % (Ap0, Np0))

# L'Hopital separatrix slopes (w',V') at the sonic point: numerators of Cramer over det, both ->0/0.
# numerators: solve 2x2 [ [c11,c12],[c21,c22] ] (g,u) = (b1,b2). w'=w*g.
c11=(1+N*V); c12=4*(N+V)/(3*(1-V**2)); c21=(4*V+N+3*N*V**2); c22=4*(1+V**2+2*N*V)/(1-V**2)
b1=-(-N*V*Ap/(3*A)+4*V*Np/3+2*N*(1+4*w/(9*(1-V**2))))
b2=-(N*(1-V**2)*Ap/A+4*(1+V**2)*Np+2*N*(1+3*V**2))
numG = sp.simplify(b1*c22 - c12*b2)   # g*det = numG
numU = sp.simplify(c11*b2 - b1*c21)   # u*det = numU
det_full = sp.simplify(c11*c22 - c12*c21)
# gradients along trajectory Y'=(N',A',w',V'); at sonic det=0, numG=numU=0. L'Hopital:
op, vp = sp.symbols('op vp')          # w'/w? use W1=w', U1=V'
W1, U1 = sp.symbols('W1 U1')
Yp = [Np0, Ap0, W1, U1]               # (N',A',w',V') at sonic (W1=w', U1=V')
def dalong(expr):
    gr = [sp.diff(expr, s) for s in (N, A, w, V)]
    return sum(sp.simplify(gr[i].subs(pt))*Yp[i] for i in range(4))
# g = numG/det with g=W1/w; so W1 = w*numG/det. At sonic use L'Hopital: W1 = w0 * d(numG)/d(det).
DnumG = dalong(numG); DnumU = dalong(numU); Ddet = dalong(det_full)
w0 = sp.Rational(3,4)
eqW = sp.expand(W1*Ddet - w0*DnumG)   # W1*det' = w0 * numG'   (since W1=w0*numG/det -> W1 det = w0 numG)
eqU = sp.expand(U1*Ddet - DnumU)      # U1*det' = numU'
solslope = sp.solve([eqW, eqU], [W1, U1], dict=True)
print("\nKHA sonic separatrix slopes (w', V'):")
roots=[]; allroots=[]
for s in solslope:
    wv=sp.simplify(s[W1]); vv=sp.simplify(s[U1])
    wc=complex(wv); vc=complex(vv)
    real = abs(wc.imag)<1e-9 and abs(vc.imag)<1e-9
    print("   w'=%s , V'=%s   (%s)" % (wc, vc, 'REAL' if real else 'complex'))
    allroots.append((wv,vv))
    if real:
        roots.append((sp.nsimplify(sp.re(wv),rational=False), sp.nsimplify(sp.re(vv),rational=False)))

# 2nd-order series for each branch (mirror sonic_exact)
xs=sp.symbols('x')
def branch_series(w1, v1):
    N2,A2,w2,V2 = sp.symbols('N2 A2 w2 V2')
    Nser = pt[N] + Np0*xs + N2*xs**2/2
    Aser = pt[A] + Ap0*xs + A2*xs**2/2
    Wser = pt[w] + w1*xs + w2*xs**2/2
    Vser = pt[V] + v1*xs + V2*xs**2/2
    sub={N:Nser,A:Aser,w:Wser,V:Vser}
    eN = sp.series(sp.diff(Nser,xs)-Np.subs(sub),xs,0,2).removeO()
    eA = sp.series(sp.diff(Aser,xs)-Ap.subs(sub),xs,0,2).removeO()
    # fluid: numG - det*(w'/w0? ) ... use w' - w*g and V'-u with g=numG/det,u=numU/det:
    # eqs: (w')*det - w*numG =0 ; (V')*det - numU=0 along traj to O(x^2)
    nG=numG.subs(sub); nU=numU.subs(sub); dd=det_full.subs(sub); ww=Wser
    eW = sp.series(sp.diff(Wser,xs)*dd - ww*nG,xs,0,3).removeO()
    eV = sp.series(sp.diff(Vser,xs)*dd - nU,xs,0,3).removeO()
    eqs=[]
    for e in [eN,eA,eW,eV]:
        p=sp.Poly(sp.expand(e),xs)
        for c in p.all_coeffs():
            cc=sp.simplify(c)
            if cc!=0: eqs.append(cc)
    s2=sp.solve(eqs,[N2,A2,w2,V2],dict=True)
    return s2
print("\n2nd-order coeffs per branch:")
for (w1,v1) in roots:
    s2=branch_series(w1,v1)
    if s2:
        ss=s2[0]
        print("  branch (w'=%.4f,V'=%.4f):"%(float(w1),float(v1)), {str(k2):round(float(v2),5) for k2,v2 in ss.items()})
    else:
        print("  branch (w'=%.4f,V'=%.4f): (empty)"%(float(w1),float(v1)))

# ---- KHA perturbation operator residue at the ingoing sonic point ----
print("\n=== KHA perturbation operator Fuchsian residue ===")
d=pickle.load(open('Lmat_kha.pkl','rb')); L=sp.sympify(d['Lmat'])
N0s,A0s,w0s,V0s=sp.symbols('N0 A0 w0 V0'); Nxs,Axs,wxs,Vxs=sp.symbols('N_x A_x w_x V_x'); ksym=sp.symbols('k')
# choose the physical branch (the one whose w',V' match the center-side approach). We'll build the
# residue for BOTH branches and print indicial exps; pick later.
for bi,(w1,v1) in enumerate(roots):
    s2=branch_series(w1,v1)
    if not s2: continue
    ss=s2[0]; N2=ss[list(ss.keys())[0]] if False else None
    N2v=sp.simplify([v for k2,v in ss.items() if str(k2)=='N2'][0])
    A2v=sp.simplify([v for k2,v in ss.items() if str(k2)=='A2'][0])
    w2v=sp.simplify([v for k2,v in ss.items() if str(k2)=='w2'][0])
    V2v=sp.simplify([v for k2,v in ss.items() if str(k2)=='V2'][0])
    t=sp.symbols('t')
    Ns=pt[N]+Np0*t+N2v*t**2/2; As=pt[A]+Ap0*t+A2v*t**2/2
    Ws=pt[w]+w1*t+w2v*t**2/2;   Vs=pt[V]+v1*t+V2v*t**2/2
    sub={N0s:Ns,A0s:As,w0s:Ws,V0s:Vs,Nxs:sp.diff(Ns,t),Axs:sp.diff(As,t),wxs:sp.diff(Ws,t),Vxs:sp.diff(Vs,t)}
    Lt=L.subs(sub)
    R=sp.zeros(4,4)
    for i in range(4):
        for j in range(4):
            ser=sp.series(Lt[i,j]*t,t,0,1).removeO()
            R[i,j]=sp.nsimplify(sp.simplify(ser.subs(t,0))) if ser!=0 else sp.Integer(0)
    lam=sp.symbols('lam')
    cp=sp.factor(sp.simplify(R.charpoly(lam).as_expr()))
    print(f"  branch {bi} (V'={float(v1):.4f}): charpoly(lam) = {cp}")
    # numeric eigs at k=2.81
    Rn=np.array(R.subs(ksym,2.81055255)).astype(np.complex128)
    print(f"    eig(R) at k=2.81: {np.round(np.linalg.eigvals(Rn),5)}")
