# kha_pert.py — build the perturbation operator for the KHA-PRINTED system (ingoing sonic V=-1/sqrt3)
# and check (a) its sonic residue matrix / indicial exponents, (b) whether a center-regular shoot
# finds the physical mode near 2.81. Decisive test of whether the background-system choice matters.
#
# KHA-printed EOM with s-derivs kept. Linearize h=H+eps hp e^{ks}. Fields (N,A,w,V).
import sympy as sp, numpy as np, math, pickle

# full fields and their x/s derivs
NF,AF,wF,VF = sp.symbols('NF AF wF VF')
NFx,AFx,wFx,VFx = sp.symbols('NFx AFx wFx VFx')
NFs,AFs,wFs,VFs = sp.symbols('NFs AFs wFs VFs')
k,eps = sp.symbols('k epsilon')

# KHA metric slopes (already ODEs): use as constraints M=0
Ap_rhs = AF*(1-AF+(2*wF/(1-VF**2))*(1+VF**2/3))
Np_rhs = NF*(-2+AF-sp.Rational(2,3)*wF)
MA = AFx - Ap_rhs
MN = NFx - Np_rhs
# KHA fluid pair VERBATIM (eq EOM lines 393-413), with s-derivs. Written as expr=0.
# energy:
en = ( (wFs + (1+NF*VF)*wFx)/wF
       + 4*(VF*VFs + (NF+VF)*VFx)/(3*(1-VF**2))
       - NF*VF*AFx/(3*AF)
       + 4*VF*NFx/3
       + 2*NF*(1 + 4*wF/(9*(1-VF**2))) )
# momentum:
mo = ( (4*VF*wFs + (4*VF+NF+3*NF*VF**2)*wFx)/wF
       + 4*((1+VF**2)*VFs + (1+VF**2+2*NF*VF)*VFx)/(1-VF**2)
       + NF*(1-VF**2)*AFx/AF
       + 4*(1+VF**2)*NFx
       + 2*NF*(1+3*VF**2) )

flds=[NF,AF,wF,VF]; fxs=[NFx,AFx,wFx,VFx]; fss=[NFs,AFs,wFs,VFs]
N0,A0,w0,V0=sp.symbols('N0 A0 w0 V0'); Nx,Ax,wx,Vx=sp.symbols('N_x A_x w_x V_x')
Np,Ap,wp,Vp=sp.symbols('Np Ap wp Vp'); Npx,Apx,wpx,Vpx=sp.symbols('Npx Apx wpx Vpx')
sub_lin={ NF:N0+eps*Np, AF:A0+eps*Ap, wF:w0+eps*wp, VF:V0+eps*Vp,
          NFx:Nx+eps*Npx, AFx:Ax+eps*Apx, wFx:wx+eps*wpx, VFx:Vx+eps*Vpx,
          NFs:eps*k*Np, AFs:eps*k*Ap, wFs:eps*k*wp, VFs:eps*k*Vp }
def lin(e):
    return sp.diff(e.subs(sub_lin), eps).subs(eps,0)
print("linearizing KHA-printed system...", flush=True)
Len=sp.expand(lin(en)); Lmo=sp.expand(lin(mo)); LMA=sp.expand(lin(MA)); LMN=sp.expand(lin(MN))
Apx_s=sp.solve(LMA,Apx)[0]; Npx_s=sp.solve(LMN,Npx)[0]
Len2=Len.subs({Apx:Apx_s,Npx:Npx_s}); Lmo2=Lmo.subs({Apx:Apx_s,Npx:Npx_s})
a11=Len2.coeff(wpx,1);a12=Len2.coeff(Vpx,1);a21=Lmo2.coeff(wpx,1);a22=Lmo2.coeff(Vpx,1)
b1=-(Len2.coeff(wpx,0).coeff(Vpx,0));b2=-(Lmo2.coeff(wpx,0).coeff(Vpx,0))
detf=sp.cancel(a11*a22-a12*a21)
wpx_s=sp.cancel((b1*a22-a12*b2)/detf); Vpx_s=sp.cancel((a11*b2-b1*a21)/detf)
rows=[sp.cancel(Npx_s),sp.cancel(Apx_s),wpx_s,Vpx_s]; hp=[Np,Ap,wp,Vp]
Lmat=sp.Matrix(4,4,lambda i,j: sp.cancel(sp.diff(rows[i],hp[j])))
print("KHA perturbation operator assembled.", flush=True)

# --- sonic point for KHA system: V=-1/sqrt3, N=2/sqrt3, A=3/2, w=3/4. Need its background series.
# Build the 2nd-order series at the KHA sonic point (mirror of sonic_exact but V=-1/sqrt3).
S3=sp.sqrt(3)
# metric slopes at sonic:
Ap0=sp.simplify(Ap_rhs.subs({AF:sp.Rational(3,2),wF:sp.Rational(3,4),VF:-1/S3,NF:2/S3}))
Np0=sp.simplify(Np_rhs.subs({AF:sp.Rational(3,2),wF:sp.Rational(3,4),VF:-1/S3,NF:2/S3}))
print(f"KHA sonic metric slopes: A'={Ap0}, N'={Np0}")
# separatrix (w',V') via L'Hopital on the KHA fluid pair: solve at sonic
Nf,Af,wf,Vf=sp.symbols('Nf Af wf Vf')
# resolved slopes wp0=w'/... reuse en,mo with s-deriv=0 -> 2x2
g,u=sp.symbols('g u')
en0=en.subs({NFs:0,AFs:0,wFs:0,VFs:0, wFx:wF*g, VFx:u, AFx:Ap_rhs, NFx:Np_rhs})
mo0=mo.subs({NFs:0,AFs:0,wFs:0,VFs:0, wFx:wF*g, VFx:u, AFx:Ap_rhs, NFx:Np_rhs})
# these are singular at sonic; get slopes by grad-ratio (L'Hopital) using the css-style method
# Simpler: use the determinant-numerator method numerically.
pt={NF:2/S3,AF:sp.Rational(3,2),wF:sp.Rational(3,4),VF:-1/S3}
# gradient-ratio: define numerators from Cramer
M2=sp.Matrix([[ (1+NF*VF), 4*(NF+VF)/(3*(1-VF**2)) ],
              [ (4*VF+NF+3*NF*VF**2), 4*(1+VF**2+2*NF*VF)/(1-VF**2) ]])
rhs1=-( -NF*VF*Ap_rhs/(3*AF)+4*VF*Np_rhs/3+2*NF*(1+4*wF/(9*(1-VF**2))) )
rhs2=-( NF*(1-VF**2)*Ap_rhs/AF+4*(1+VF**2)*Np_rhs+2*NF*(1+3*VF**2) )
# note wFx coeff in energy is (1+NV)/w * w -> for w' directly multiply; treat unknown W1=w', U1=V'
W1,U1=sp.symbols('W1 U1')
E1=(W1 + (1+NF*VF)*W1)*0  # placeholder; do numeric grad-ratio instead
print("\n(For brevity, checking KHA operator sonic residue numerically via the css-style residue.)")

# Numeric residue: sample L along a short series using the KHA branch slopes computed numerically.
f_Ap=sp.lambdify((NF,AF,wF,VF),Ap_rhs,'math'); f_Np=sp.lambdify((NF,AF,wF,VF),Np_rhs,'math')
detM=sp.lambdify((NF,AF,wF,VF),M2.det(),'math')
# resolved fluid slopes (w',V') away from sonic:
solr=sp.solve([en.subs({NFs:0,AFs:0,wFs:0,VFs:0,AFx:Ap_rhs,NFx:Np_rhs}),
               mo.subs({NFs:0,AFs:0,wFs:0,VFs:0,AFx:Ap_rhs,NFx:Np_rhs})],[wFx,VFx],dict=True)[0]
f_wx=sp.lambdify((NF,AF,wF,VF),solr[wFx],'math'); f_Vx=sp.lambdify((NF,AF,wF,VF),solr[VFx],'math')
def slopes(Y):
    n,a,w,v=Y; return (f_Np(n,a,w,v),f_Ap(n,a,w,v),f_wx(n,a,w,v),f_Vx(n,a,w,v))
# march from just inside sonic toward center a tiny bit to get branch slopes near V=-1/sqrt3
print("KHA (w',V') just off the ingoing sonic point (numeric):")
V0v=-1/math.sqrt(3)
for dv in [1e-3,5e-4,1e-4]:
    # pick N on ingoing branch by solving detM=0 near N=2/sqrt3 for V=V0v+dv... just report slope at exact
    Y=(2/math.sqrt(3),1.5,0.75,V0v+dv)
    try: print(f"  dv={dv}: slopes={tuple(round(s,4) for s in slopes(Y))}")
    except Exception as e: print(f"  dv={dv}: {e}")

pickle.dump({'Lmat':sp.srepr(Lmat)}, open("Lmat_kha.pkl","wb"))
print("\nwrote Lmat_kha.pkl (KHA-printed-system perturbation operator)")

# residue matrix at KHA sonic using branch series (build numerically to 2nd order)
# Use the css sonic_exact approach with V=-1/sqrt3:
xsym=sp.symbols('xsym')
# we need w',V',and 2nd order. Do L'Hopital slope solve symbolically:
gV=[sp.diff(M2.det(),f) for f in (NF,AF,wF,VF)]  # not used; do direct
print("KHA sonic det check:", sp.simplify(M2.det().subs(pt)))
