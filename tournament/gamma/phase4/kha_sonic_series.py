# kha_sonic_series.py — 2nd-order background series at the KHA-printed ingoing sonic point
# (N=2/sqrt3, A=3/2, w=3/4, V=-1/sqrt3), then the perturbation-operator residue matrix + charpoly.
import sympy as sp
S3=sp.sqrt(3)
N,A,w,V=sp.symbols('N A w V')
Ap_rhs=A*(1-A+(2*w/(1-V**2))*(1+V**2/3)); Np_rhs=N*(-2+A-sp.Rational(2,3)*w)

# resolved fluid slopes away from sonic (KHA), solve 2x2:
g,u=sp.symbols('g u')
en=((0+(1+N*V)*w*g)/w+4*(0+(N+V)*u)/(3*(1-V**2))-N*V*Ap_rhs/(3*A)+4*V*Np_rhs/3+2*N*(1+4*w/(9*(1-V**2))))
mo=((0+(4*V+N+3*N*V**2)*w*g)/w+4*(0+(1+V**2+2*N*V)*u)/(1-V**2)+N*(1-V**2)*Ap_rhs/A+4*(1+V**2)*Np_rhs+2*N*(1+3*V**2))
sol=sp.solve([en,mo],[g,u],dict=True)[0]
wx_rhs=sp.simplify(w*sol[g]); Vx_rhs=sp.simplify(sol[u])

pt={N:2/S3,A:sp.Rational(3,2),w:sp.Rational(3,4),V:-1/S3}
Ap0=sp.simplify(Ap_rhs.subs(pt)); Np0=sp.simplify(Np_rhs.subs(pt))
print("KHA sonic metric slopes: A'=",Ap0," N'=",Np0)

# L'Hopital separatrix slopes (w',V') at sonic: grad-ratio
detf=sp.together(Vx_rhs); den=sp.denom(detf)
# use the standard: numV/det with det=denominator. Build numerators:
numV=sp.numer(sp.together(Vx_rhs)); numW=sp.numer(sp.together(wx_rhs)); detD=sp.denom(sp.together(Vx_rhs))
gV=[sp.diff(numV,f) for f in (N,A,w,V)]; gW=[sp.diff(numW,f) for f in (N,A,w,V)]; gD=[sp.diff(detD,f) for f in (N,A,w,V)]
gV=[sp.simplify(e.subs(pt)) for e in gV]; gW=[sp.simplify(e.subs(pt)) for e in gW]; gD=[sp.simplify(e.subs(pt)) for e in gD]
op,vp=sp.symbols('op vp'); Yp=[Np0,Ap0,op,vp]
DnV=sum(gV[i]*Yp[i] for i in range(4)); DnW=sum(gW[i]*Yp[i] for i in range(4)); Dd=sum(gD[i]*Yp[i] for i in range(4))
sols=sp.solve([sp.expand(vp*Dd-DnV),sp.expand(op*Dd-DnW)],[op,vp],dict=True)
print("KHA separatrix slopes (w',V'):")
roots=[]
for s in sols:
    ov=complex(s[op]); vv=complex(s[vp]); print(f"  w'={sp.nsimplify(s[op])} ({ov.real:+.4f}), V'={sp.nsimplify(s[vp])} ({vv.real:+.4f})")
    if abs(ov.imag)<1e-9 and abs(vv.imag)<1e-9: roots.append((float(ov.real),float(vv.real)))

# pick the branch closest to a mirror of branch1 (w'~4.5): whichever has smaller |w'|
roots.sort(key=lambda t:abs(t[0])); w1,V1=roots[0]
print(f"  -> using branch w'={w1:.4f}, V'={V1:.4f}")

# 2nd-order series (mirror sonic_exact)
xs=sp.symbols('xs'); N2,A2,w2,V2=sp.symbols('N2 A2 w2 V2')
N1=Np0;A1=Ap0
Nser=pt[N]+N1*xs+N2*xs**2/2; Aser=pt[A]+A1*xs+A2*xs**2/2
Wser=pt[w]+sp.nsimplify(w1)*xs+w2*xs**2/2; Vser=pt[V]+sp.nsimplify(V1)*xs+V2*xs**2/2
subT={N:Nser,A:Aser,w:Wser,V:Vser}
eN=sp.series(sp.diff(Nser,xs)-Np_rhs.subs(subT),xs,0,2).removeO()
eA=sp.series(sp.diff(Aser,xs)-Ap_rhs.subs(subT),xs,0,2).removeO()
eV=sp.series(sp.numer(sp.together(Vx_rhs)).subs(subT)-sp.denom(sp.together(Vx_rhs)).subs(subT)*sp.diff(Vser,xs),xs,0,3).removeO()
eW=sp.series(sp.numer(sp.together(wx_rhs)).subs(subT)-sp.denom(sp.together(wx_rhs)).subs(subT)*sp.diff(Wser,xs),xs,0,3).removeO()
eqs=[]
for e in [eN,eA,eV,eW]:
    p=sp.Poly(sp.expand(e),xs)
    for c in p.all_coeffs():
        cc=sp.simplify(c)
        if cc!=0: eqs.append(cc)
s2=sp.solve(eqs,[N2,A2,w2,V2],dict=True)[0]
print("KHA 2nd-order:",{str(kk):float(vv) for kk,vv in s2.items()})

# Now build the KHA perturbation operator's residue matrix using this series.
import pickle
Lm=sp.sympify(pickle.load(open("Lmat_kha.pkl","rb"))['Lmat'])
N0,A0,w0,V0=sp.symbols('N0 A0 w0 V0'); Nx,Ax,wx,Vx=sp.symbols('N_x A_x w_x V_x')
Ns=pt[N]+N1*xs+s2[N2]*xs**2/2; As=pt[A]+A1*xs+s2[A2]*xs**2/2
Ws=pt[w]+sp.nsimplify(w1)*xs+s2[w2]*xs**2/2; Vs=pt[V]+sp.nsimplify(V1)*xs+s2[V2]*xs**2/2
subL={N0:Ns,A0:As,w0:Ws,V0:Vs,Nx:sp.diff(Ns,xs),Ax:sp.diff(As,xs),wx:sp.diff(Ws,xs),Vx:sp.diff(Vs,xs)}
Lx=Lm.subs(subL)
R=sp.zeros(4,4)
for i in range(4):
    for j in range(4):
        ser=sp.series(Lx[i,j]*xs,xs,0,1).removeO()
        R[i,j]=sp.simplify(ser.subs(xs,0)) if ser!=0 else sp.Integer(0)
lam=sp.symbols('lam')
print("\n=== KHA operator residue matrix R ===")
sp.pprint(sp.simplify(R))
print("charpoly:",sp.factor(sp.simplify(R.charpoly(lam).as_expr())))
print("eig(R):",{sp.simplify(e):m for e,m in R.eigenvals().items()})
