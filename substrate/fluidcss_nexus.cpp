// ============================================================================
//  fluidcss_nexus — TINY UNIVERSE crown: radiation-fluid CSS critical exponent
//  Contract: contracts/fluidcss_nexus.contract.md  ·  v1.0.0  (Stage A + Stage B)
//
//  Single-file CPU fp64 (no GPU). gamma = 4/3 (radiation), units G = c = 1.
//
//  STAGE A — the TRUE Evans-Coleman critical CSS background, by HKA's own
//  construction (gr-qc/9607010 sec IV): the sonic point is a ONE-PARAMETER family
//  in V0 (eqs 4.7-4.9 give N0, A0, om0 as functions of V0); shoot sonic -> center;
//  the EC solution sits at the transition between the two failure modes
//  (case 1: A<1; case 2: a second sonic point), with exactly ONE zero of V
//  (Bogoyavlensky index 1).  V0* = 0.1124394014;  Nbar'(sonic) = -0.355699
//  (the paper's sonic-gauge gauge-mode value; its printed "0.35699" is a
//  dropped-digit typo).  D-032.
//
//  HISTORY (D-032 corrects D-027): v0.9.x "Stage A" constructed the collapsing
//  flat FRIEDMANN solution (V0 = -1/sqrt3 = -c_s IS the Friedmann point of the
//  sonic line; its "exact" invariants A=1+2om/3, oi*=3/8, 2m/r=1/3 are FRW
//  radiation identities; V = -sqrt(1-1/A) to 1.9e-10).  That build is retained
//  here as the --friedmann CONTROL face (closed-form exact solution = selftest);
//  the golden `b4f4e463` it froze is superseded (goldens/fluidcss_stageA/NOTE.md).
//
//  STAGE B — the eigenvalue (HKA sec V shoot, on the Stage-A background):
//  sonic Frobenius modes (a0 in ker R; recursion (nI-R)a_n = sum L_j a_{n-1-j})
//  + identity eq:alg-PP + gauge Nbar_p(sonic)=0, norm Abar_p(sonic)=1 ->
//  the unique analytic-at-sonic solution; integrate sonic -> center; the
//  discriminant is V_p(x_c) (the (0,0,0,1)e^{-2x} expanding-mode content).
//  Deterministic bisections:  relevant kappa0 in [2.70,2.90]  -> beta = 1/kappa0;
//  CONTROLS: sonic-gauge mode in [0.30,0.40] (= -Nbar'_ss(sonic), analytic) and
//  origin-gauge mode in [0.98,1.02] (= 1 exactly).
//  Gates: G-ANCHOR |beta - 0.35580192| < 4e-3 · G-CONVERGE · G-UNIQUE.
//  Python reference (this repo): kappa0 = 2.810552374, beta = 0.355801945
//  (lit. 2.8105525488 / 0.35580192).
//
//  Build: cl /std:c++17 /EHsc /O2 /W4 substrate\fluidcss_nexus.cpp /Fe:build\fluidcss_nexus.exe
//  Faces: --stageA | --stageB | --friedmann  x  --json|--golden|--selftest. Exit 0/1/2.
// ============================================================================

#define _CRT_SECURE_NO_WARNINGS
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <array>

static const double GG      = 4.0/3.0;            // adiabatic index (radiation)
static const double GM1     = GG - 1.0;           // 1/3
static const double SG      = 0.5773502691896257; // sqrt(gamma-1) = 1/sqrt3 = c_s
static const double MCENTER = -2.0/(3.0*GG);      // -1/2  (regular-center M = NV)

// ============================================================================
//  BLAKE2b-256 (RFC 7693, host) — lifted VERBATIM from the _nexus family.
// ============================================================================
namespace blake2b {
static const uint64_t IV[8] = {
    0x6a09e667f3bcc908ull, 0xbb67ae8584caa73bull, 0x3c6ef372fe94f82bull,
    0xa54ff53a5f1d36f1ull, 0x510e527fade682d1ull, 0x9b05688c2b3e6c1full,
    0x1f83d9abfb41bd6bull, 0x5be0cd19137e2179ull };
static const uint8_t SIGMA[10][16] = {
    { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15},
    {14,10, 4, 8, 9,15,13, 6, 1,12, 0, 2,11, 7, 5, 3},
    {11, 8,12, 0, 5, 2,15,13,10,14, 3, 6, 7, 1, 9, 4},
    { 7, 9, 3, 1,13,12,11,14, 2, 6, 5,10, 4, 0,15, 8},
    { 9, 0, 5, 7, 2, 4,10,15,14, 1,11,12, 6, 8, 3,13},
    { 2,12, 6,10, 0,11, 8, 3, 4,13, 7, 5,15,14, 1, 9},
    {12, 5, 1,15,14,13, 4,10, 0, 7, 6, 3, 9, 2, 8,11},
    {13,11, 7,14,12, 1, 3, 9, 5, 0,15, 4, 8, 6, 2,10},
    { 6,15,14, 9,11, 3, 0, 8,12, 2,13, 7, 1, 4,10, 5},
    {10, 2, 8, 4, 7, 6, 1, 5,15,11, 9,14, 3,12,13, 0} };
static inline uint64_t rotr(uint64_t x, int n){ return (x >> n) | (x << (64 - n)); }
static inline void Gm(uint64_t* v, int a, int b, int c, int d, uint64_t x, uint64_t y){
    v[a] += v[b] + x; v[d] = rotr(v[d] ^ v[a], 32);
    v[c] += v[d];     v[b] = rotr(v[b] ^ v[c], 24);
    v[a] += v[b] + y; v[d] = rotr(v[d] ^ v[a], 16);
    v[c] += v[d];     v[b] = rotr(v[b] ^ v[c], 63);
}
static void compress(uint64_t h[8], const uint8_t block[128], uint64_t t, bool last){
    uint64_t v[16], m[16];
    for (int i=0;i<8;i++){ v[i]=h[i]; v[i+8]=IV[i]; }
    v[12]^=t; if (last) v[14]=~v[14];
    for (int i=0;i<16;i++){ uint64_t w=0; for (int j=7;j>=0;j--) w=(w<<8)|block[i*8+j]; m[i]=w; }
    for (int r=0;r<12;r++){ const uint8_t* s=SIGMA[r%10];
        Gm(v,0,4, 8,12,m[s[ 0]],m[s[ 1]]); Gm(v,1,5, 9,13,m[s[ 2]],m[s[ 3]]);
        Gm(v,2,6,10,14,m[s[ 4]],m[s[ 5]]); Gm(v,3,7,11,15,m[s[ 6]],m[s[ 7]]);
        Gm(v,0,5,10,15,m[s[ 8]],m[s[ 9]]); Gm(v,1,6,11,12,m[s[10]],m[s[11]]);
        Gm(v,2,7, 8,13,m[s[12]],m[s[13]]); Gm(v,3,4, 9,14,m[s[14]],m[s[15]]); }
    for (int i=0;i<8;i++) h[i]^=v[i]^v[i+8];
}
static std::string hash256_hex(const std::string& in){
    uint64_t h[8]; for (int i=0;i<8;i++) h[i]=IV[i]; h[0]^=0x01010000ull^32ull;
    size_t n=in.size(), off=0; uint8_t block[128];
    while (n-off>128){ memcpy(block,in.data()+off,128); off+=128; compress(h,block,(uint64_t)off,false); }
    size_t rem=n-off; memset(block,0,128); if (rem) memcpy(block,in.data()+off,rem);
    compress(h,block,(uint64_t)n,true);
    char hex[65]; for (int i=0;i<32;i++){ unsigned b=(unsigned)((h[i/8]>>(8*(i%8)))&0xFF); snprintf(hex+2*i,3,"%02x",b); }
    return std::string(hex,64);
}
} // namespace blake2b

static std::string fmt6(double x){ if(!std::isfinite(x)) return "9999999.999999"; if(x==0.0)x=0.0; char b[48]; snprintf(b,sizeof(b),"%.6f",x); return std::string(b); }
static std::string fmt9(double x){ if(!std::isfinite(x)) return "9.999999999e+99"; if(x==0.0)x=0.0; char b[48]; snprintf(b,sizeof(b),"%.9e",x); return std::string(b); }

// ============================================================================
//  Background ODE system (rflanl.tex 4.1-4.4; A eliminated via constraint 4.2)
// ============================================================================
static inline double A_of(double N,double om,double V){
    double oV2=1.0-V*V;
    return 1.0 + 2.0*om*(1.0+GM1*V*V)/oV2 + 2.0*GG*N*V*om/oV2;
}
static inline double Dson(double N,double V){
    double a=1.0+N*V, b=N+V;
    return a*a - GM1*b*b;
}
static inline void fluid_slopes(double A,double N,double om,double V,double& omx,double& Vx){
    double g=GG, oV2=1.0-V*V;
    double RHS_d=(3.0*(2.0-g)/2.0)*N*V - ((2.0+g)/2.0)*A*N*V + (2.0-g)*N*V*om;
    double RHS_e=(2.0-g)*(g-1.0)*N*om + ((7.0*g-6.0)/2.0)*N + ((2.0-3.0*g)/2.0)*A*N;
    double a=(1.0+N*V)/om, b=g*(N+V)/oV2, c=(g-1.0)*(N+V)/om, d=g*(1.0+N*V)/oV2;
    double det=a*d-b*c;
    omx=(d*RHS_d-b*RHS_e)/det;
    Vx =(a*RHS_e-c*RHS_d)/det;
}
static inline void rhs3(double N,double om,double V,double& Nx,double& omx,double& Vx){
    double A=A_of(N,om,V);
    Nx=N*(-2.0+A-(2.0-GG)*om);
    fluid_slopes(A,N,om,V,omx,Vx);
}
static inline double lhop(double A,double N,double om,double V){
    double g=GG;
    double RHS_d=(3.0*(2.0-g)/2.0)*N*V - ((2.0+g)/2.0)*A*N*V + (2.0-g)*N*V*om;
    double RHS_e=(2.0-g)*(g-1.0)*N*om + ((7.0*g-6.0)/2.0)*N + ((2.0-3.0*g)/2.0)*A*N;
    return (1.0+N*V)*RHS_e - (g-1.0)*(N+V)*RHS_d;
}

// Sonic-line values parametrized by V0 (rflanl.tex 4.7-4.9).
static inline void sonic_values(double V0,double& A0,double& N0,double& om0){
    double g=GG, sg=SG;
    N0 = (1.0 - sg*V0)/(sg - V0);
    A0 = (g*g + 4.0*g - 4.0 + 8.0*sg*sg*sg*V0 - (3.0*g-2.0)*(2.0-g)*V0*V0)/(g*g*(1.0-V0*V0));
    om0= 2.0*sg*(sg - V0)*(1.0 + sg*V0)/(g*g*(1.0-V0*V0));
}

// ============================================================================
//  Truncated power-series arithmetic (index k = t^k, t = x - x_s), real fp64.
// ============================================================================
typedef std::vector<double> Ser;
static Ser sconst(int K,double v){ Ser a(K,0.0); a[0]=v; return a; }
static Ser smul(const Ser& a,const Ser& b){
    int K=(int)a.size(); Ser c(K,0.0);
    for (int i=0;i<K;i++){ double ai=a[i]; if (ai==0.0) continue;
        for (int j=0;j<K-i;j++) c[i+j]+=ai*b[j]; }
    return c;
}
static Ser sadd(const Ser& a,const Ser& b){ Ser c=a; for (size_t i=0;i<c.size();i++) c[i]+=b[i]; return c; }
static Ser ssub(const Ser& a,const Ser& b){ Ser c=a; for (size_t i=0;i<c.size();i++) c[i]-=b[i]; return c; }
static Ser sscale(const Ser& a,double s){ Ser c=a; for (size_t i=0;i<c.size();i++) c[i]*=s; return c; }
static Ser sinv(const Ser& a){
    int K=(int)a.size(); Ser b(K,0.0); b[0]=1.0/a[0];
    for (int n=1;n<K;n++){ double s=0.0; for (int k=1;k<=n;k++) s+=a[k]*b[n-k]; b[n]=-s/a[0]; }
    return b;
}
static Ser sderiv(const Ser& a){                       // (a')_k = (k+1) a_{k+1}
    int K=(int)a.size(); Ser d(K,0.0);
    for (int k=0;k<K-1;k++) d[k]=(k+1)*a[k+1];
    return d;
}
static double seval(const Ser& a,double t,int nmax=-1){
    int K=(nmax<0)?(int)a.size():((nmax<(int)a.size())?nmax:(int)a.size());
    double s=0.0, tp=1.0; for (int k=0;k<K;k++){ s+=a[k]*tp; tp*=t; } return s;
}

// A(t) from the constraint 4.2 on series.
static Ser A_series(const Ser& N,const Ser& om,const Ser& V){
    int K=(int)N.size();
    Ser V2=smul(V,V), NV=smul(N,V);
    Ser S=sinv(ssub(sconst(K,1.0),V2));
    Ser t1=sscale(smul(om,smul(sadd(sconst(K,1.0),sscale(V2,GM1)),S)),2.0);
    Ser t2=sscale(smul(NV,smul(om,S)),2.0*GG);
    return sadd(sconst(K,1.0),sadd(t1,t2));
}
// Fluid pair (om-cleared): M(t) (2x2 of series) and rhs b(t) (2 of series).
static void fluid_MB(const Ser& A,const Ser& N,const Ser& om,const Ser& V,
                     Ser M[2][2], Ser b[2]){
    int K=(int)N.size(); double g=GG;
    Ser NV=smul(N,V), V2=smul(V,V);
    Ser S=sinv(ssub(sconst(K,1.0),V2));
    Ser NpV=sadd(N,V), omS=smul(om,S), onePNV=sadd(sconst(K,1.0),NV);
    M[0][0]=onePNV;                    M[0][1]=sscale(smul(NpV,omS),g);
    M[1][0]=sscale(NpV,g-1.0);         M[1][1]=sscale(smul(onePNV,omS),g);
    Ser RHS_d=sadd(sadd(sscale(NV,3.0*(2.0-g)/2.0),sscale(smul(A,NV),-(2.0+g)/2.0)),
                   sscale(smul(NV,om),2.0-g));
    Ser RHS_e=sadd(sadd(sscale(smul(N,om),(2.0-g)*(g-1.0)),sscale(N,(7.0*g-6.0)/2.0)),
                   sscale(smul(A,N),(2.0-3.0*g)/2.0));
    b[0]=smul(om,RHS_d); b[1]=smul(om,RHS_e);
}

// ============================================================================
//  Sonic Taylor series at general V0 (order-1 solvability is QUADRATIC: two
//  branches; branch selected by the root value passed in).  Port of nr_ec2.
// ============================================================================
struct BgSer { Ser A,N,om,V; };

// rank-1 2x2 helpers: right/left null vectors + Moore-Penrose pinv application.
static void rank1_null(const double M0[2][2], double nR[2], double nL[2]){
    // right null: M0 nR = 0 -> nR ~ (M01, -M00) or (M11, -M10) (larger row)
    double r0=std::hypot(M0[0][0],M0[0][1]), r1=std::hypot(M0[1][0],M0[1][1]);
    if (r0>=r1){ nR[0]= M0[0][1]; nR[1]=-M0[0][0]; }
    else       { nR[0]= M0[1][1]; nR[1]=-M0[1][0]; }
    double nn=std::hypot(nR[0],nR[1]); nR[0]/=nn; nR[1]/=nn;
    // left null: nL^T M0 = 0 -> nL ~ (M10, -M00)^T? use columns: nL ~ (c1?) larger col
    double c0=std::hypot(M0[0][0],M0[1][0]), c1=std::hypot(M0[0][1],M0[1][1]);
    if (c0>=c1){ nL[0]= M0[1][0]; nL[1]=-M0[0][0]; }
    else       { nL[0]= M0[1][1]; nL[1]=-M0[0][1]; }
    nn=std::hypot(nL[0],nL[1]); nL[0]/=nn; nL[1]/=nn;
}
static void rank1_pinv_apply(const double M0[2][2], const double r[2], double w[2]){
    // M0 = sigma u v^T; pinv(M0) r = (u.r/sigma) v.  v = unit(bigger row), u = unit(bigger col).
    double v[2], u[2];
    double r0=std::hypot(M0[0][0],M0[0][1]), r1=std::hypot(M0[1][0],M0[1][1]);
    if (r0>=r1){ v[0]=M0[0][0]; v[1]=M0[0][1]; } else { v[0]=M0[1][0]; v[1]=M0[1][1]; }
    double nv=std::hypot(v[0],v[1]); v[0]/=nv; v[1]/=nv;
    double c0=std::hypot(M0[0][0],M0[1][0]), c1=std::hypot(M0[0][1],M0[1][1]);
    if (c0>=c1){ u[0]=M0[0][0]; u[1]=M0[1][0]; } else { u[0]=M0[0][1]; u[1]=M0[1][1]; }
    double nu=std::hypot(u[0],u[1]); u[0]/=nu; u[1]/=nu;
    double sigma = u[0]*(M0[0][0]*v[0]+M0[0][1]*v[1]) + u[1]*(M0[1][0]*v[0]+M0[1][1]*v[1]);
    double ur = u[0]*r[0]+u[1]*r[1];
    w[0]=(ur/sigma)*v[0]; w[1]=(ur/sigma)*v[1];
}

// residual series Phi = M.(om',V') - b at the current truncation
static void phi_of(const Ser& N,const Ser& om,const Ser& V, Ser P[2]){
    Ser A=A_series(N,om,V); Ser M[2][2], b[2];
    fluid_MB(A,N,om,V,M,b);
    Ser wpo=sderiv(om), wpv=sderiv(V);
    P[0]=ssub(sadd(smul(M[0][0],wpo),smul(M[0][1],wpv)),b[0]);
    P[1]=ssub(sadd(smul(M[1][0],wpo),smul(M[1][1],wpv)),b[1]);
}

// Solve the order-1 quadratic solvability: returns number of real roots (0/1/2).
static int order1_roots(double V0,int K,const double nR[2],const double nL[2],
                        const double w1[2], double roots[2]){
    double A0,N0,om0; sonic_values(V0,A0,N0,om0);
    double samp[3]; double alphas[3]={-2.0,0.0,2.0};
    for (int si=0;si<3;si++){
        Ser N=sconst(K,N0), om=sconst(K,om0), V=sconst(K,V0);
        Ser A=A_series(N,om,V);
        // N1 from the metric row (order-0 only):
        Ser rhsN=smul(N,sadd(sadd(sconst(K,-2.0),A),sscale(om,-(2.0-GG))));
        N[1]=rhsN[0];
        om[1]=w1[0]+alphas[si]*nR[0]; V[1]=w1[1]+alphas[si]*nR[1];
        Ser P[2]; phi_of(N,om,V,P);
        samp[si] = -(nL[0]*P[0][1]+nL[1]*P[1][1]);
    }
    double C=samp[1], Aq=(samp[0]+samp[2]-2.0*C)/8.0, Bq=(samp[2]-samp[0])/4.0;
    if (std::fabs(Aq)<1e-12){ roots[0]=-C/Bq; return 1; }
    double disc=Bq*Bq-4.0*Aq*C;
    if (disc<0.0) return 0;
    double sq=std::sqrt(disc);
    double r1=(-Bq-sq)/(2.0*Aq), r2=(-Bq+sq)/(2.0*Aq);
    roots[0]=(r1<r2)?r1:r2; roots[1]=(r1<r2)?r2:r1;
    return 2;
}

// Build the sonic Taylor series for a given branch root (order-1 alpha).
static BgSer bg_series_at(double V0,int order,double branch_root){
    int K=order+2;
    double A0,N0,om0; sonic_values(V0,A0,N0,om0);
    Ser N=sconst(K,N0), om=sconst(K,om0), V=sconst(K,V0);
    // order-0 structures:
    Ser A=A_series(N,om,V); Ser M[2][2], b[2]; fluid_MB(A,N,om,V,M,b);
    double M0m[2][2]={{M[0][0][0],M[0][1][0]},{M[1][0][0],M[1][1][0]}};
    double nR[2],nL[2]; rank1_null(M0m,nR,nL);
    double pnext[2]={0,0};
    for (int m=0;m<=order;m++){
        if (m==1){
            om[1]=pnext[0]+branch_root*nR[0]; V[1]=pnext[1]+branch_root*nR[1];
        } else if (m>=2){
            om[m]=pnext[0]; V[m]=pnext[1];
            Ser P[2]; phi_of(N,om,V,P);
            double q0=-(nL[0]*P[0][m]+nL[1]*P[1][m]);
            om[m]=pnext[0]+nR[0]; V[m]=pnext[1]+nR[1];
            phi_of(N,om,V,P);
            double q1=-(nL[0]*P[0][m]+nL[1]*P[1][m]);
            double den=q1-q0;
            double alpha=(std::fabs(den)>1e-11)? -q0/den : 0.0;
            om[m]=pnext[0]+alpha*nR[0]; V[m]=pnext[1]+alpha*nR[1];
        }
        if (m<order){
            Ser Am=A_series(N,om,V);
            Ser rhsN=smul(N,sadd(sadd(sconst(K,-2.0),Am),sscale(om,-(2.0-GG))));
            N[m+1]=rhsN[m]/(m+1);
            Ser P[2]; phi_of(N,om,V,P);
            double r[2]={-P[0][m],-P[1][m]}, w[2];
            rank1_pinv_apply(M0m,r,w);
            pnext[0]=w[0]/(m+1); pnext[1]=w[1]/(m+1);
            om[m+1]=0.0; V[m+1]=0.0;
        }
    }
    BgSer bs; bs.A=A_series(N,om,V); bs.N=N; bs.om=om; bs.V=V;
    return bs;
}

// ============================================================================
//  STAGE A: the EC root-find (sonic -> center; paper's case1/case2 criterion)
// ============================================================================
struct ShootRes { int side; double x_stop, F, A_end, om_end, nzero_partial; bool ok; };
// side: +1 = case2 (another sonic point / om->0 peel), -1 = case1 (A<1), 0 = unusable

static ShootRes shoot_center_branch(double V0,double branch_root,double delta,double x_end,double dx,int order){
    BgSer bs=bg_series_at(V0,order,branch_root);
    double N=seval(bs.N,-delta), om=seval(bs.om,-delta), V=seval(bs.V,-delta);
    double x=-delta, Dprev=Dson(N,V);
    ShootRes rr; rr.side=0; rr.x_stop=x; rr.F=0; rr.A_end=0; rr.om_end=0; rr.nzero_partial=0; rr.ok=false;
    int nz=0; double Vprev=V;
    long nst=(long)((x - x_end)/dx);
    for (long i=0;i<nst;i++){
        double k1N,k1o,k1V; rhs3(N,om,V,k1N,k1o,k1V);
        double k2N,k2o,k2V; rhs3(N-0.5*dx*k1N,om-0.5*dx*k1o,V-0.5*dx*k1V,k2N,k2o,k2V);
        double k3N,k3o,k3V; rhs3(N-0.5*dx*k2N,om-0.5*dx*k2o,V-0.5*dx*k2V,k3N,k3o,k3V);
        double k4N,k4o,k4V; rhs3(N-dx*k3N,om-dx*k3o,V-dx*k3V,k4N,k4o,k4V);
        double Nn=N-dx/6.0*(k1N+2*k2N+2*k3N+k4N);
        double on=om-dx/6.0*(k1o+2*k2o+2*k3o+k4o);
        double Vn=V-dx/6.0*(k1V+2*k2V+2*k3V+k4V);
        x-=dx;
        if (!std::isfinite(Nn)||!std::isfinite(on)||!std::isfinite(Vn)||std::fabs(Vn)>=1.0||on>1e8){
            // stall/blowup: om->0 = case2-side peel; else classify by F sign
            rr.x_stop=x; rr.A_end=A_of(N,om,V); rr.om_end=om; rr.F=N*V+0.5; rr.nzero_partial=nz;
            rr.side = (om<1e-3||rr.F>0.0)? +1 : -1; rr.ok=true; return rr;
        }
        double Acur=A_of(Nn,on,Vn), Dcur=Dson(Nn,Vn);
        if (Vprev*Vn<0.0) nz++;
        Vprev=Vn;
        if (Acur<1.0){ rr.side=-1; rr.x_stop=x; rr.F=Nn*Vn+0.5; rr.A_end=Acur; rr.om_end=on; rr.nzero_partial=nz; rr.ok=true; return rr; }
        if (i>10 && Dprev<0.0 && Dcur>=0.0){ rr.side=+1; rr.x_stop=x; rr.F=Nn*Vn+0.5; rr.A_end=Acur; rr.om_end=on; rr.nzero_partial=nz; rr.ok=true; return rr; }
        N=Nn; om=on; V=Vn; Dprev=Dcur;
    }
    // reached x_end ("deep"): near-root or the spurious N->0 attractor (A>>1, F=+0.5)
    rr.x_stop=x; rr.A_end=A_of(N,om,V); rr.om_end=om; rr.F=N*V+0.5; rr.nzero_partial=nz;
    if (rr.A_end<2.0 && std::fabs(rr.F-0.5)>1e-3){ rr.side=(rr.F>0.0)?+1:-1; rr.ok=true; }
    else rr.side=0;
    return rr;
}

// Physical branch = the real root whose run goes deepest with a case1/case2 verdict.
static ShootRes shoot_phys(double V0,double delta,double x_end,double dx,int order,double* used_root=nullptr){
    int K=order+2;
    double A0,N0,om0; sonic_values(V0,A0,N0,om0);
    Ser N=sconst(K,N0), om=sconst(K,om0), V=sconst(K,V0);
    Ser A=A_series(N,om,V); Ser M[2][2], b[2]; fluid_MB(A,N,om,V,M,b);
    double M0m[2][2]={{M[0][0][0],M[0][1][0]},{M[1][0][0],M[1][1][0]}};
    double nR[2],nL[2]; rank1_null(M0m,nR,nL);
    Ser P[2]; phi_of(N,om,V,P);
    double r1v[2]={-P[0][0],-P[1][0]}, w1[2];
    rank1_pinv_apply(M0m,r1v,w1);
    double roots[2]; int nr=order1_roots(V0,K,nR,nL,w1,roots);
    ShootRes best; best.side=0; best.ok=false; best.x_stop=0.0;
    for (int i=0;i<nr;i++){
        ShootRes rr=shoot_center_branch(V0,roots[i],delta,x_end,dx,order);
        if (!rr.ok) continue;
        if (!best.ok || rr.x_stop<best.x_stop){ best=rr; if (used_root) *used_root=roots[i]; }
    }
    return best;
}

static double find_V0(double lo,double hi,double delta,double x_end,double dx,int order,int iters){
    // bisection on side: +1 (case2) below the root, -1 (case1) above
    for (int i=0;i<iters;i++){
        double mid=0.5*(lo+hi);
        ShootRes r=shoot_phys(mid,delta,x_end,dx,order);
        int s=r.ok? r.side : 0;
        if (s==0){ // nudge
            r=shoot_phys(mid+1e-9,delta,x_end,dx,order); s=r.ok? r.side:0;
            if (s==0){ r=shoot_phys(mid-1e-9,delta,x_end,dx,order); s=r.ok? r.side:0; }
            if (s==0){ lo=lo+0.25*(mid-lo); continue; }
        }
        if (s>0) lo=mid; else hi=mid;
    }
    return 0.5*(lo+hi);
}

// ============================================================================
//  Background storage (center relaunch) + Hermite interpolation
// ============================================================================
struct Bg {
    double x0, h, xs;                       // grid start, step, sonic x
    std::vector<double> N, om, V;           // node values
    std::vector<double> Nd, omd, Vd;        // node derivatives (rhs3)
    int n;
    void state(double x,double& Nv,double& omv,double& Vv) const {
        double u=(x-x0)/h; int i=(int)u; if (i<0)i=0; if (i>n-2)i=n-2;
        double t=u-i, t2=t*t, t3=t2*t;
        double h00=2*t3-3*t2+1, h10=t3-2*t2+t, h01=-2*t3+3*t2, h11=t3-t2;
        Nv = h00*N[i] +h10*h*Nd[i] +h01*N[i+1] +h11*h*Nd[i+1];
        omv= h00*om[i]+h10*h*omd[i]+h01*om[i+1]+h11*h*omd[i+1];
        Vv = h00*V[i] +h10*h*Vd[i] +h01*V[i+1] +h11*h*Vd[i+1];
    }
};

static Bg relaunch_center(double N_inf,double om_inf,double x0,double h,double xmax){
    Bg bg; bg.x0=x0; bg.h=h;
    double z0=std::exp(x0);
    double N=N_inf/z0, om=om_inf*z0*z0, V=(MCENTER/N_inf)*z0, x=x0;
    double Dprev=Dson(N,V);
    long nmax=(long)((xmax-x0)/h)+2;
    bg.N.reserve(nmax); bg.om.reserve(nmax); bg.V.reserve(nmax);
    bg.Nd.reserve(nmax); bg.omd.reserve(nmax); bg.Vd.reserve(nmax);
    bg.xs=x0;
    for (long i=0;i<nmax;i++){
        double Nx,ox,Vx; rhs3(N,om,V,Nx,ox,Vx);
        bg.N.push_back(N); bg.om.push_back(om); bg.V.push_back(V);
        bg.Nd.push_back(Nx); bg.omd.push_back(ox); bg.Vd.push_back(Vx);
        double k1N=Nx,k1o=ox,k1V=Vx;
        double k2N,k2o,k2V; rhs3(N+0.5*h*k1N,om+0.5*h*k1o,V+0.5*h*k1V,k2N,k2o,k2V);
        double k3N,k3o,k3V; rhs3(N+0.5*h*k2N,om+0.5*h*k2o,V+0.5*h*k2V,k3N,k3o,k3V);
        double k4N,k4o,k4V; rhs3(N+h*k3N,om+h*k3o,V+h*k3V,k4N,k4o,k4V);
        double Nn=N+h/6.0*(k1N+2*k2N+2*k3N+k4N);
        double on=om+h/6.0*(k1o+2*k2o+2*k3o+k4o);
        double Vn=V+h/6.0*(k1V+2*k2V+2*k3V+k4V);
        double Dcur=Dson(Nn,Vn);
        if (Dprev<0.0 && Dcur>=0.0){
            double frac=-Dprev/(Dcur-Dprev);
            bg.xs=x+frac*h;
            break;
        }
        N=Nn; om=on; V=Vn; x+=h; Dprev=Dcur;
    }
    bg.n=(int)bg.N.size();
    return bg;
}

// ============================================================================
//  STAGE B: the perturbation operator (hka_pert_hka99 / nr_rederive-proven
//  coefficients), Frobenius sonic modes, and the shoot.
// ============================================================================
struct Coeffs { double As,Bs,Cs,Ds,Ax,Bx,Cx,Dx,E1,E2,E3,E4,F1,F2,F3,F4,G1,G3,G4,H1,H3; };
static Coeffs coeffs_at(double A,double N,double om,double V,double obx,double Vx){
    double g=GG, oV2=1.0-V*V; Coeffs c;
    c.As=1.0; c.Bs=g*V/oV2; c.Cs=(g-1.0)*V; c.Ds=g/oV2;
    c.Ax=1.0+N*V; c.Bx=g*(N+V)/oV2; c.Cx=(g-1.0)*(N+V); c.Dx=g*(1.0+N*V)/oV2;
    c.E1=-((g+2.0)/2.0)*A*N*V;
    c.E2=((6.0-3.0*g)/2.0)*N*V-((2.0+g)/2.0)*A*N*V+(2.0-g)*om*N*V-N*V*obx-g*N*Vx/oV2;
    c.E3=(2.0-g)*om*N*V;
    c.E4=((6.0-3.0*g)/2.0)*N-((2.0+g)/2.0)*A*N+(2.0-g)*om*N-N*obx-g*(1.0+2.0*N*V+V*V)*Vx/(oV2*oV2);
    c.F1=((2.0-3.0*g)/2.0)*A*N;
    c.F2=(2.0-g)*(g-1.0)*om*N+((7.0*g-6.0)/2.0)*N+((2.0-3.0*g)/2.0)*A*N-(g-1.0)*N*obx-g*N*V*Vx/oV2;
    c.F3=(2.0-g)*(g-1.0)*om*N;
    c.F4=-(g-1.0)*obx-g*(N+2.0*V+N*V*V)*Vx/(oV2*oV2);
    c.G1=-A; c.G3=2.0*(1.0+(g-1.0)*V*V)*om/oV2; c.G4=4.0*g*om*V/(oV2*oV2);
    c.H1=A; c.H3=(g-2.0)*om;
    return c;
}
// L(x;kappa) = Mx^{-1}(Gmat - kappa Ms) as a 4x4 (pointwise).
static void Lnum(double A,double N,double om,double V,double obx,double Vx,double kappa,double L[4][4]){
    Coeffs c=coeffs_at(A,N,om,V,obx,Vx);
    double Gm[4][4]={{c.G1,0,c.G3,c.G4},{c.H1,0,c.H3,0},
                     {c.E1,c.E2,c.E3,c.E4},{c.F1,c.F2,c.F3,c.F4}};
    Gm[2][2]-=kappa*c.As; Gm[2][3]-=kappa*c.Bs;
    Gm[3][2]-=kappa*c.Cs; Gm[3][3]-=kappa*c.Ds;
    // rows 0,1: Mx = I
    for (int j=0;j<4;j++){ L[0][j]=Gm[0][j]; L[1][j]=Gm[1][j]; }
    // fluid block: invert [[Ax,Bx],[Cx,Dx]]
    double det=c.Ax*c.Dx-c.Bx*c.Cx;
    for (int j=0;j<4;j++){
        L[2][j]=( c.Dx*Gm[2][j]-c.Bx*Gm[3][j])/det;
        L[3][j]=(-c.Cx*Gm[2][j]+c.Ax*Gm[3][j])/det;
    }
}

// gauss solve 4x4 (partial pivot), in-place on copies.
static bool gauss4(double A[4][4],double b[4],double x[4]){
    double M[4][5];
    for (int i=0;i<4;i++){ for (int j=0;j<4;j++) M[i][j]=A[i][j]; M[i][4]=b[i]; }
    for (int col=0;col<4;col++){
        int piv=col; for (int r=col+1;r<4;r++) if (std::fabs(M[r][col])>std::fabs(M[piv][col])) piv=r;
        if (std::fabs(M[piv][col])<1e-300) return false;
        if (piv!=col) for (int j=col;j<5;j++){ double t=M[col][j]; M[col][j]=M[piv][j]; M[piv][j]=t; }
        for (int r=col+1;r<4;r++){ double f=M[r][col]/M[col][col];
            for (int j=col;j<5;j++) M[r][j]-=f*M[col][j]; }
    }
    for (int i=3;i>=0;i--){ double s=M[i][4];
        for (int j=i+1;j<4;j++) s-=M[i][j]*x[j];
        x[i]=s/M[i][i]; }
    return true;
}

// Operator Laurent on the branch series: R (residue) and L_0..L_{order}.
struct Laurent { double R[4][4]; std::vector<std::array<double,16>> L; };
static Laurent laurent_exact(const BgSer& bs,double kappa,int order){
    int K=order+4;
    auto pad=[&](const Ser& a){ Ser b(K,0.0); for (size_t i=0;i<a.size()&&(int)i<K;i++) b[i]=a[i]; return b; };
    Ser A=pad(bs.A), N=pad(bs.N), om=pad(bs.om), V=pad(bs.V);
    Ser obx=smul(sderiv(om),sinv(om)), Vx=sderiv(V);
    Ser V2=smul(V,V), NV=smul(N,V), NpV=sadd(N,V);
    Ser oV2=ssub(sconst(K,1.0),V2), S=sinv(oV2), S2=smul(S,S);
    Ser one=sconst(K,1.0);
    Ser cAs=one, cBs=sscale(smul(V,S),GG), cCs=sscale(V,GM1), cDs=sscale(S,GG);
    Ser cAx=sadd(one,NV), cBx=sscale(smul(NpV,S),GG), cCx=sscale(NpV,GM1), cDx=sscale(smul(sadd(one,NV),S),GG);
    Ser AN=smul(A,N), ANV=smul(A,NV), omNV=smul(om,NV), omN=smul(om,N);
    double g=GG;
    Ser E1=sscale(ANV,-(g+2.0)/2.0);
    Ser E2=ssub(ssub(sadd(sadd(sscale(NV,(6.0-3.0*g)/2.0),sscale(ANV,-(2.0+g)/2.0)),sscale(omNV,2.0-g)),
                     smul(NV,obx)), sscale(smul(smul(N,Vx),S),g));
    Ser E3=sscale(omNV,2.0-g);
    Ser E4=ssub(ssub(sadd(sadd(sscale(N,(6.0-3.0*g)/2.0),sscale(AN,-(2.0+g)/2.0)),sscale(omN,2.0-g)),
                     smul(N,obx)),
                sscale(smul(smul(sadd(sadd(one,sscale(NV,2.0)),V2),Vx),S2),g));
    Ser F1=sscale(AN,(2.0-3.0*g)/2.0);
    Ser F2=ssub(ssub(sadd(sadd(sscale(omN,(2.0-g)*(g-1.0)),sscale(N,(7.0*g-6.0)/2.0)),sscale(AN,(2.0-3.0*g)/2.0)),
                     sscale(smul(N,obx),g-1.0)), sscale(smul(smul(NV,Vx),S),g));
    Ser F3=sscale(omN,(2.0-g)*(g-1.0));
    Ser F4=ssub(sscale(obx,-(g-1.0)),
                sscale(smul(smul(sadd(sadd(N,sscale(V,2.0)),smul(N,V2)),Vx),S2),g));
    Ser G1=sscale(A,-1.0), G3=sscale(smul(sadd(one,sscale(V2,g-1.0)),smul(om,S)),2.0), G4=sscale(smul(smul(om,V),S2),4.0*g);
    Ser H1=A, H3=sscale(om,g-2.0);
    Ser Z=sconst(K,0.0);
    // (Gmat - kappa Ms)(t):
    Ser GmM[4][4]={{G1,Z,G3,G4},{H1,Z,H3,Z},
                   {E1,E2,ssub(E3,sscale(cAs,kappa)),ssub(E4,sscale(cBs,kappa))},
                   {F1,F2,ssub(F3,sscale(cCs,kappa)),ssub(F4,sscale(cDs,kappa))}};
    // fluid block inverse as Laurent: Mfl^{-1} = (1/t) P(t)
    Ser det=ssub(smul(cAx,cDx),smul(cBx,cCx));    // ~ t (d1 + ...)
    Ser Dsh(K,0.0); for (int i=0;i<K-1;i++) Dsh[i]=det[i+1];
    Ser invD=sinv(Dsh);
    Ser Pm[2][2]={{smul(cDx,invD),sscale(smul(cBx,invD),-1.0)},
                  {sscale(smul(cCx,invD),-1.0),smul(cAx,invD)}};
    // Q(t) = P(t) . GmM[fluid rows] (2x4)
    Ser Q[2][4];
    for (int i=0;i<2;i++) for (int j=0;j<4;j++)
        Q[i][j]=sadd(smul(Pm[i][0],GmM[2][j]),smul(Pm[i][1],GmM[3][j]));
    Laurent la;
    for (int i=0;i<4;i++) for (int j=0;j<4;j++) la.R[i][j]=0.0;
    for (int j=0;j<4;j++){ la.R[2][j]=Q[0][j][0]; la.R[3][j]=Q[1][j][0]; }
    la.L.resize(order+1);
    for (int m=0;m<=order;m++){
        std::array<double,16> Lm{};
        for (int j=0;j<4;j++){
            Lm[0*4+j]=GmM[0][j][m]; Lm[1*4+j]=GmM[1][j][m];
            Lm[2*4+j]=Q[0][j][m+1]; Lm[3*4+j]=Q[1][j][m+1];
        }
        la.L[m]=Lm;
    }
    return la;
}

// ker(R) basis (R is rank 1: rows 2,3 ~ u_i w^T). 3 orthonormal vectors ⟂ w.
static void kerR_basis(const double R[4][4], double B[3][4]){
    double w[4]; double n2=std::hypot(std::hypot(R[2][0],R[2][1]),std::hypot(R[2][2],R[2][3]));
    double n3=std::hypot(std::hypot(R[3][0],R[3][1]),std::hypot(R[3][2],R[3][3]));
    if (n2>=n3){ for (int j=0;j<4;j++) w[j]=R[2][j]; } else { for (int j=0;j<4;j++) w[j]=R[3][j]; }
    double nw=std::sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]+w[3]*w[3]);
    for (int j=0;j<4;j++) w[j]/=nw;
    int nb=0;
    for (int e=0;e<4 && nb<3;e++){
        double v[4]={0,0,0,0}; v[e]=1.0;
        double d=v[0]*w[0]+v[1]*w[1]+v[2]*w[2]+v[3]*w[3];
        for (int j=0;j<4;j++) v[j]-=d*w[j];
        for (int k=0;k<nb;k++){
            double dd=v[0]*B[k][0]+v[1]*B[k][1]+v[2]*B[k][2]+v[3]*B[k][3];
            for (int j=0;j<4;j++) v[j]-=dd*B[k][j];
        }
        double nv=std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]+v[3]*v[3]);
        if (nv<1e-8) continue;
        for (int j=0;j<4;j++) B[nb][j]=v[j]/nv;
        nb++;
    }
}

// The unique analytic-at-sonic solution at parameter kappa, evaluated at t0.
// (a0 in ker R + identity eq:alg-PP + gauge Nbar_p=0, norm Abar_p=1.)
static bool psi_sonic(const BgSer& bs,double V0,double kappa,double t0,int order,double psi[4]){
    Laurent la=laurent_exact(bs,kappa,order);
    double B[3][4]; kerR_basis(la.R,B);
    double A0,N0,om0; sonic_values(V0,A0,N0,om0);
    double oV2=1.0-V0*V0;
    double idc[4]={kappa-A0, 2.0*GG*N0*V0*om0/oV2,
                   2.0*om0*(1.0+GM1*V0*V0+GG*N0*V0)/oV2,
                   2.0*GG*om0*(N0*(1.0+V0*V0)+2.0*V0)/(oV2*oV2)};
    // solve 3x3: [identity; Nbar=0; Abar=1] on the ker-R combo c
    double M3[3][3], rhs3v[3]={0.0,0.0,1.0};
    for (int i=0;i<3;i++){
        M3[0][i]=idc[0]*B[i][0]+idc[1]*B[i][1]+idc[2]*B[i][2]+idc[3]*B[i][3];
        M3[1][i]=B[i][1];
        M3[2][i]=B[i][0];
    }
    // Cramer 3x3
    auto det3=[&](double m[3][3]){ return m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])
        -m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0]) +m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]); };
    double D=det3(M3); if (std::fabs(D)<1e-300) return false;
    double c[3];
    for (int k=0;k<3;k++){
        double Mk[3][3]; for (int i=0;i<3;i++) for (int j=0;j<3;j++) Mk[i][j]=M3[i][j];
        for (int i=0;i<3;i++) Mk[i][k]=rhs3v[i];
        c[k]=det3(Mk)/D;
    }
    double a0[4]={0,0,0,0};
    for (int k=0;k<3;k++) for (int j=0;j<4;j++) a0[j]+=c[k]*B[k][j];
    // Frobenius recursion a_n, sum at t0
    std::vector<std::array<double,4>> a(order+1);
    for (int j=0;j<4;j++) a[0][j]=a0[j];
    for (int n=1;n<=order;n++){
        double rhs[4]={0,0,0,0};
        for (int j=0;j<n;j++){
            if (j>= (int)la.L.size()) break;
            const std::array<double,16>& Lj=la.L[j];
            for (int r=0;r<4;r++) for (int col=0;col<4;col++) rhs[r]+=Lj[r*4+col]*a[n-1-j][col];
        }
        double Mm[4][4];
        for (int r=0;r<4;r++) for (int col=0;col<4;col++) Mm[r][col]=((r==col)?(double)n:0.0)-la.R[r][col];
        double an[4];
        if (!gauss4(Mm,rhs,an)) return false;
        for (int j=0;j<4;j++) a[n][j]=an[j];
    }
    for (int j=0;j<4;j++){ double s=0.0, tp=1.0; for (int n=0;n<=order;n++){ s+=a[n][j]*tp; tp*=t0; } psi[j]=s; }
    return true;
}

// integrate the eigen-ODE sonic -> center on the stored background; return V_p(xc).
static double shootB(const Bg& bg,const BgSer& bs,double V0,double kappa,
                     double t0,int order,double xc,double dx){
    double y[4];
    if (!psi_sonic(bs,V0,kappa,t0,order,y)) return NAN;
    double x=bg.xs+t0;
    auto deriv=[&](double xx,const double yy[4],double dy[4]){
        double N,om,V; bg.state(xx,N,om,V);
        double A=A_of(N,om,V);
        double omx,Vx; fluid_slopes(A,N,om,V,omx,Vx);
        double L[4][4]; Lnum(A,N,om,V,omx/om,Vx,kappa,L);
        for (int r=0;r<4;r++){ dy[r]=0.0; for (int c2=0;c2<4;c2++) dy[r]+=L[r][c2]*yy[c2]; }
    };
    long nst=(long)((x-xc)/dx);
    for (long i=0;i<nst;i++){
        double k1[4],k2[4],k3[4],k4[4],yt[4];
        deriv(x,y,k1);
        for (int j=0;j<4;j++) yt[j]=y[j]-0.5*dx*k1[j];
        deriv(x-0.5*dx,yt,k2);
        for (int j=0;j<4;j++) yt[j]=y[j]-0.5*dx*k2[j];
        deriv(x-0.5*dx,yt,k3);
        for (int j=0;j<4;j++) yt[j]=y[j]-dx*k3[j];
        deriv(x-dx,yt,k4);
        for (int j=0;j<4;j++) y[j]-=dx/6.0*(k1[j]+2*k2[j]+2*k3[j]+k4[j]);
        x-=dx;
        if (!std::isfinite(y[0])) return NAN;
    }
    return y[3];
}

static double bisect_kappa(const Bg& bg,const BgSer& bs,double V0,double lo,double hi,
                           double t0,int order,double xc,double dx,int iters){
    double flo=shootB(bg,bs,V0,lo,t0,order,xc,dx);
    for (int i=0;i<iters;i++){
        double mid=0.5*(lo+hi);
        double fm=shootB(bg,bs,V0,mid,t0,order,xc,dx);
        if (flo*fm<=0.0){ hi=mid; } else { lo=mid; flo=fm; }
    }
    return 0.5*(lo+hi);
}

// count sign changes of the discriminant over a coarse kappa scan (G-UNIQUE).
static int count_roots(const Bg& bg,const BgSer& bs,double V0,double klo,double khi,double dk,
                       double t0,int order,double xc,double dx){
    int n=0; double prev=NAN;
    for (double k=klo;k<=khi+1e-12;k+=dk){
        double f=shootB(bg,bs,V0,k,t0,order,xc,dx);
        if (std::isfinite(prev)&&std::isfinite(f)&&prev*f<0.0) n++;
        prev=f;
    }
    return n;
}

// ============================================================================
//  Results + JSON (v1.0.0)
// ============================================================================
struct ResA {
    double V0=0, A0=0, N0=0, om0=0, Nbp=0, xs=0, N_inf=0, om_inf=0;
    double close_res=0, conv=0; int nzero=1; bool nan_free=true;
    bool g_v0=false, g_fp=false, g_zero=false, g_close=false, g_conv=false, verdict=false;
};
struct ResB {
    double kappa0=0, beta=0, kg_sonic=0, kg_origin=0, conv=0; int nroots=0;
    bool g_anchor=false, g_conv=false, g_unique=false, g_ctrl=false, verdict=false;
};

static const double BETA_REF = 0.35580192;     // KHA95/HKA99 literature value
static const double NBP_REF  = -0.355699;      // -Nbar'_ss(sonic) fingerprint (analytic from V0*)

struct Pipeline {
    double V0=0; BgSer bs; Bg bg; double N_inf=0, om_inf=0;
    double delta=1e-3, dxA=1e-4, x_endA=-14.0; int orderA=6, itersA=100;
};

static Pipeline buildEC(){
    Pipeline P;
    P.V0 = find_V0(0.110,0.116,P.delta,P.x_endA,P.dxA,P.orderA,P.itersA);
    double used_root=0.0;
    ShootRes rr=shoot_phys(P.V0,P.delta,P.x_endA,P.dxA,P.orderA,&used_root);
    P.bs = bg_series_at(P.V0,18,used_root);
    // tail-extract center data a bit above the stall, then relaunch center->sonic
    // re-integrate the physical branch to the stall to sample the tail:
    {
        BgSer b6=bg_series_at(P.V0,P.orderA,used_root);
        double N=seval(b6.N,-P.delta), om=seval(b6.om,-P.delta), V=seval(b6.V,-P.delta);
        double x=-P.delta, x_fit=rr.x_stop+0.7;
        while (x>x_fit){
            double k1N,k1o,k1V; rhs3(N,om,V,k1N,k1o,k1V);
            double k2N,k2o,k2V; rhs3(N-0.5*P.dxA*k1N,om-0.5*P.dxA*k1o,V-0.5*P.dxA*k1V,k2N,k2o,k2V);
            double k3N,k3o,k3V; rhs3(N-0.5*P.dxA*k2N,om-0.5*P.dxA*k2o,V-0.5*P.dxA*k2V,k3N,k3o,k3V);
            double k4N,k4o,k4V; rhs3(N-P.dxA*k3N,om-P.dxA*k3o,V-P.dxA*k3V,k4N,k4o,k4V);
            N-=P.dxA/6.0*(k1N+2*k2N+2*k3N+k4N);
            om-=P.dxA/6.0*(k1o+2*k2o+2*k3o+k4o);
            V-=P.dxA/6.0*(k1V+2*k2V+2*k3V+k4V);
            x-=P.dxA;
        }
        P.om_inf=om*std::exp(-2.0*x); P.N_inf=N*std::exp(x);
    }
    P.bg = relaunch_center(P.N_inf,P.om_inf,-16.0,2.0e-4,1.0);
    return P;
}

static ResA runStageA(const Pipeline& P){
    ResA R;
    R.V0=P.V0; sonic_values(P.V0,R.A0,R.N0,R.om0);
    R.Nbp=-2.0+R.A0-(2.0-GG)*R.om0;
    R.xs=P.bg.xs; R.N_inf=P.N_inf; R.om_inf=P.om_inf;
    // closure: relaunched sonic state vs sonic_values(V0)
    double Ns,oms,Vs; P.bg.state(P.bg.xs-1e-6,Ns,oms,Vs);
    double As=A_of(Ns,oms,Vs);
    R.close_res=std::fabs(As-R.A0)+std::fabs(Ns-R.N0)+std::fabs(oms-R.om0)+std::fabs(Vs-P.V0);
    // one V zero along the relaunched profile
    int nz=0;
    for (int i=1;i<P.bg.n;i++) if (P.bg.V[i-1]*P.bg.V[i]<0.0) nz++;
    R.nzero=nz;
    // convergence: V0 with dx/2
    double V0f=find_V0(0.110,0.116,P.delta,P.x_endA,P.dxA*0.5,P.orderA,P.itersA);
    R.conv=std::fabs(P.V0-V0f);
    R.nan_free=std::isfinite(R.V0)&&std::isfinite(R.close_res);
    R.g_v0   = std::fabs(R.V0-0.1124394014)<5e-4;        // vs the Python-validated root
    R.g_fp   = std::fabs(R.Nbp-NBP_REF)<5e-4;            // the paper's sonic-gauge fingerprint
    R.g_zero = (R.nzero==1);                             // Evans-Coleman = Bogoyavlensky index 1
    R.g_close= (R.close_res<5e-3);                       // two-sided closure of the orbit
    R.g_conv = (R.conv<5e-5);
    R.verdict=R.nan_free&&R.g_v0&&R.g_fp&&R.g_zero&&R.g_close&&R.g_conv;
    return R;
}

static ResB runStageB(const Pipeline& P){
    ResB R;
    double t0=-0.02; int order=18; double xc=-12.0, dx=2.0e-4; int iters=60;
    R.kappa0  = bisect_kappa(P.bg,P.bs,P.V0,2.70,2.90,t0,order,xc,dx,iters);
    R.beta    = 1.0/R.kappa0;
    R.kg_sonic= bisect_kappa(P.bg,P.bs,P.V0,0.30,0.40,t0,order,xc,dx,iters);
    R.kg_origin=bisect_kappa(P.bg,P.bs,P.V0,0.98,1.02,t0,order,xc,dx,iters);
    // G-CONVERGE: kappa0 under dx/2 and t0=-0.03
    double k_f = bisect_kappa(P.bg,P.bs,P.V0,2.70,2.90,t0,order,xc,dx*0.5,iters);
    double k_t = bisect_kappa(P.bg,P.bs,P.V0,2.70,2.90,-0.03,order,xc,dx,iters);
    R.conv = std::fabs(R.kappa0-k_f)+std::fabs(R.kappa0-k_t);
    // G-UNIQUE: coarse scan sign changes in (1.1, 3.6) — must be exactly 1 (the relevant)
    R.nroots = count_roots(P.bg,P.bs,P.V0,1.10,3.60,0.05,t0,order,xc,dx);
    R.g_anchor = std::fabs(R.beta-BETA_REF)<4e-3;
    R.g_conv   = (R.conv<1e-4);
    R.g_unique = (R.nroots==1);
    double A0s,N0s,om0s; sonic_values(P.V0,A0s,N0s,om0s);
    double kg_ref = -(-2.0+A0s-(2.0-GG)*om0s);           // = -Nbar'_ss(sonic), analytic
    R.g_ctrl   = std::fabs(R.kg_sonic-kg_ref)<5e-4 && std::fabs(R.kg_origin-1.0)<5e-3;
    R.verdict  = R.g_anchor&&R.g_conv&&R.g_unique&&R.g_ctrl;
    return R;
}

static std::string declaredJsonA(const ResA& R){
    std::string s; s.reserve(1024);
    s += "{\"tool\":\"fluidcss_nexus\",\"version\":\"1.0.0\",\"units\":\"G=c=1\",\"stage\":\"A\"";
    s += ",\"params\":{\"gamma\":" + fmt6(GG) + ",\"delta\":1.0e-03,\"dx\":1.0e-04,\"x_end\":-14.0,\"order\":6,\"iters\":100}";
    s += ",\"background\":{";
    s += "\"V0_star\":" + fmt9(R.V0);
    s += ",\"sonic_A0\":" + fmt9(R.A0) + ",\"sonic_N0\":" + fmt9(R.N0) + ",\"sonic_om0\":" + fmt9(R.om0);
    s += ",\"Nbar_prime_sonic\":" + fmt9(R.Nbp) + ",\"fingerprint_ref\":-0.355699000";
    s += ",\"x_sonic\":" + fmt6(R.xs) + ",\"N_inf\":" + fmt9(R.N_inf) + ",\"om_inf\":" + fmt9(R.om_inf);
    s += ",\"closure_res\":" + fmt9(R.close_res) + ",\"V_zeros\":" + std::to_string(R.nzero);
    s += ",\"V0_converge\":" + fmt9(R.conv);
    s += ",\"nan_free\":" + std::string(R.nan_free?"1":"0");
    s += "},\"gates\":{";
    s += "\"G_V0\":" + std::string(R.g_v0?"true":"false");
    s += ",\"G_FINGERPRINT\":" + std::string(R.g_fp?"true":"false");
    s += ",\"G_ONE_ZERO\":" + std::string(R.g_zero?"true":"false");
    s += ",\"G_CLOSURE\":" + std::string(R.g_close?"true":"false");
    s += ",\"G_CONVERGE\":" + std::string(R.g_conv?"true":"false");
    s += "},\"note\":\"true Evans-Coleman background (D-032; supersedes the v0.9.x Friedmann-as-EC build)\"";
    s += ",\"verdict\":\"" + std::string(R.verdict?"pass":"fail") + "\"";
    return s;
}

static std::string declaredJsonB(const ResA& RA,const ResB& R){
    std::string s; s.reserve(1024);
    s += "{\"tool\":\"fluidcss_nexus\",\"version\":\"1.0.0\",\"units\":\"G=c=1\",\"stage\":\"B\"";
    s += ",\"params\":{\"t0\":-0.02,\"order\":18,\"x_c\":-12.0,\"dx\":2.0e-04,\"iters\":60}";
    s += ",\"background_V0_star\":" + fmt9(RA.V0);
    s += ",\"eigen\":{";
    s += "\"kappa0\":" + fmt9(R.kappa0) + ",\"beta\":" + fmt9(R.beta);
    s += ",\"beta_ref\":" + fmt9(BETA_REF);
    s += ",\"gauge_sonic\":" + fmt9(R.kg_sonic) + ",\"gauge_origin\":" + fmt9(R.kg_origin);
    s += ",\"converge_spread\":" + fmt9(R.conv) + ",\"roots_in_1p1_3p6\":" + std::to_string(R.nroots);
    s += "},\"gates\":{";
    s += "\"G_ANCHOR\":" + std::string(R.g_anchor?"true":"false");
    s += ",\"G_CONVERGE\":" + std::string(R.g_conv?"true":"false");
    s += ",\"G_UNIQUE\":" + std::string(R.g_unique?"true":"false");
    s += ",\"G_CONTROLS\":" + std::string(R.g_ctrl?"true":"false");
    s += "},\"verdict\":\"" + std::string(R.verdict?"pass":"fail") + "\"";
    return s;
}

// ---------- legacy Friedmann control (the v0.9.x construction, now correctly labeled) ----
static void runFriedmann(){
    // center->sonic shoot targeting V=-c_s: converges to the FRIEDMANN point (3/2,2/sqrt3,3/4)
    // -- closed-form exact; a control face, not a golden.
    double lo=0.35, hi=0.40, N_inf=1.0, x0=-12.0, dx=1e-3, xmax=0.5;
    auto sh=[&](double oi,double& Vs,double& As,double& Ns,double& oms)->bool{
        double z0=std::exp(x0);
        double N=N_inf/z0, om=oi*z0*z0, V=(MCENTER/N_inf)*z0, Dprev=Dson(N,V);
        long nst=(long)((xmax-x0)/dx);
        for (long i=0;i<nst;i++){
            double k1N,k1o,k1V; rhs3(N,om,V,k1N,k1o,k1V);
            double k2N,k2o,k2V; rhs3(N+0.5*dx*k1N,om+0.5*dx*k1o,V+0.5*dx*k1V,k2N,k2o,k2V);
            double k3N,k3o,k3V; rhs3(N+0.5*dx*k2N,om+0.5*dx*k2o,V+0.5*dx*k2V,k3N,k3o,k3V);
            double k4N,k4o,k4V; rhs3(N+dx*k3N,om+dx*k3o,V+dx*k3V,k4N,k4o,k4V);
            double Nn=N+dx/6.0*(k1N+2*k2N+2*k3N+k4N);
            double on=om+dx/6.0*(k1o+2*k2o+2*k3o+k4o);
            double Vn=V+dx/6.0*(k1V+2*k2V+2*k3V+k4V);
            double Dcur=Dson(Nn,Vn);
            if (Dprev<0.0&&Dcur>=0.0){
                double f=-Dprev/(Dcur-Dprev);
                Vs=V+f*(Vn-V); Ns=N+f*(Nn-N); oms=om+f*(on-om); As=A_of(Ns,oms,Vs);
                return true;
            }
            N=Nn; om=on; V=Vn; Dprev=Dcur;
        }
        return false;
    };
    double Vs=0,As=0,Ns=0,oms=0, flo;
    sh(lo,Vs,As,Ns,oms); flo=Vs+SG;
    for (int i=0;i<60;i++){
        double mid=0.5*(lo+hi), fm;
        sh(mid,Vs,As,Ns,oms); fm=Vs+SG;
        if (flo*fm<=0.0) hi=mid; else { lo=mid; flo=fm; }
    }
    double oi=0.5*(lo+hi);
    sh(oi,Vs,As,Ns,oms);
    printf("[friedmann control] the v0.9.x construction (sonic criterion V=-c_s):\n");
    printf("  oi* = %.9f   (Friedmann exact 3/8 = 0.375; (3/2)H^2t^2)\n", oi);
    printf("  sonic (A,N,om,V) = (%.6f, %.6f, %.6f, %.6f)   (Friedmann point of the sonic line)\n", As,Ns,oms,Vs);
    printf("  This solution is the collapsing flat FRIEDMANN spacetime, NOT Evans-Coleman (D-032).\n");
    bool ok = std::fabs(oi-0.375)<1e-3 && std::fabs(As-1.5)<2e-3;
    printf("VERDICT: %s\n", ok?"PASS":"FAIL");
    exit(ok?0:1);
}

static void printHumanA(const ResA& R){
    printf("fluidcss_nexus v1.0.0 - STAGE A: the TRUE Evans-Coleman CSS background (D-032)\n");
    printf("-------------------------------------------------------\n");
    printf("  sonic parameter V0*        %.10f   [gate |.-0.1124394|<5e-4]  %s\n", R.V0, R.g_v0?"PASS":"FAIL");
    printf("  sonic (A0,N0,om0)          (%.6f, %.6f, %.6f)\n", R.A0,R.N0,R.om0);
    printf("  Nbar'(sonic)               %.6f   (paper fingerprint -0.355699)  %s\n", R.Nbp, R.g_fp?"PASS":"FAIL");
    printf("  zeros of V (EC index)      %d   (exactly 1)  %s\n", R.nzero, R.g_zero?"PASS":"FAIL");
    printf("  center-relaunch closure    %.2e   [<5e-3]  %s\n", R.close_res, R.g_close?"PASS":"FAIL");
    printf("  V0 grid-convergence        %.2e   [<5e-5]  %s\n", R.conv, R.g_conv?"PASS":"FAIL");
    printf("  center data: N_inf=%.6f om_inf=%.6f   sonic at x_s=%.6f\n", R.N_inf, R.om_inf, R.xs);
    printf("-------------------------------------------------------\n");
    printf("VERDICT: %s\n", R.verdict?"PASS":"FAIL");
}
static void printHumanB(const ResA& RA,const ResB& R){
    printf("fluidcss_nexus v1.0.0 - STAGE B: the relevant eigenvalue -> beta (D-032)\n");
    printf("-------------------------------------------------------\n");
    printf("  kappa0 (relevant)   %.9f   (lit 2.8105525488)\n", R.kappa0);
    printf("  beta = 1/kappa0     %.9f   (lit %.8f)  [G-ANCHOR |dbeta|<4e-3]  %s\n",
           R.beta, BETA_REF, R.g_anchor?"PASS":"FAIL");
    printf("  CONTROL sonic-gauge %.9f   (analytic -Nbar'_ss = %.6f)  \n", R.kg_sonic, -RA.Nbp);
    printf("  CONTROL origin-gauge %.9f  (exact 1)   [both]  %s\n", R.kg_origin, R.g_ctrl?"PASS":"FAIL");
    printf("  G-CONVERGE spread   %.2e   [<1e-4]  %s\n", R.conv, R.g_conv?"PASS":"FAIL");
    printf("  G-UNIQUE roots in (1.1,3.6)  %d   (exactly 1)  %s\n", R.nroots, R.g_unique?"PASS":"FAIL");
    printf("-------------------------------------------------------\n");
    printf("VERDICT: %s\n", R.verdict?"PASS":"FAIL");
}

static int goldenFace(const std::string& declared,const char* path){
    std::string hash=blake2b::hash256_hex(declared);
    FILE* f=fopen(path,"rb");
    if (!f){ fprintf(stderr,"GOLDEN NOT FROZEN %s\n",hash.c_str()); printf("%s\n",hash.c_str()); return 2; }
    char want[128]={0};
    if (fscanf(f,"%127s",want)!=1){ fclose(f); fprintf(stderr,"GOLDEN FILE UNREADABLE\n"); return 2; }
    fclose(f);
    if (hash==std::string(want)){ fprintf(stderr,"GOLDEN OK %.8s\n",hash.c_str()); printf("%s\n",hash.c_str()); return 0; }
    fprintf(stderr,"GOLDEN MISMATCH have=%.8s want=%.8s\n",hash.c_str(),want);
    printf("%s\n",hash.c_str()); return 1;
}

int main(int argc, char** argv){
    bool json=false, selftest=false, golden=false, stageB=false, fried=false;
    for (int i=1;i<argc;i++){
        std::string a=argv[i];
        if      (a=="--json") json=true;
        else if (a=="--selftest") selftest=true;
        else if (a=="--golden") golden=true;
        else if (a=="--stageA") {/*default*/}
        else if (a=="--stageB") stageB=true;
        else if (a=="--friedmann") fried=true;
        else { fprintf(stderr,"usage: fluidcss_nexus [--stageA|--stageB|--friedmann] [--json] [--golden] [--selftest]\n"); return 2; }
    }

    if (selftest){
        // closed-form checks: (a) sonic_values at the Friedmann point V0=-1/sqrt3 must give
        // (3/2, 2/sqrt3, 3/4); (b) Dson=0 there; (c) A_of consistency.
        double A0,N0,om0; sonic_values(-SG,A0,N0,om0);
        double N0e=2.0/std::sqrt(3.0);
        bool ok1=std::fabs(A0-1.5)<1e-12&&std::fabs(N0-N0e)<1e-12&&std::fabs(om0-0.75)<1e-12;
        bool ok2=std::fabs(Dson(N0e,-SG))<1e-12;
        bool ok3=std::fabs(A_of(N0e,0.75,-SG)-1.5)<1e-9;
        printf("[selftest] sonic_values(Friedmann) [%s]  Dson=0 [%s]  A_of match [%s]\n",
               ok1?"PASS":"FAIL", ok2?"PASS":"FAIL", ok3?"PASS":"FAIL");
        printf("VERDICT: %s\n",(ok1&&ok2&&ok3)?"PASS":"FAIL");
        return (ok1&&ok2&&ok3)?0:1;
    }
    if (fried){ runFriedmann(); }

    Pipeline P=buildEC();
    ResA RA=runStageA(P);
    if (!stageB){
        std::string declared=declaredJsonA(RA);
        if (golden) return goldenFace(declared,"goldens/fluidcss_stageA/golden.hash");
        if (json){ printf("%s,\"notes\":\"hash=%.8s\"}\n",declared.c_str(),blake2b::hash256_hex(declared).c_str()); return RA.verdict?0:1; }
        printHumanA(RA);
        printf("declared hash: %s\n", blake2b::hash256_hex(declared).c_str());
        return RA.verdict?0:1;
    }
    ResB RB=runStageB(P);
    std::string declared=declaredJsonB(RA,RB);
    if (golden) return goldenFace(declared,"goldens/fluidcss_stageB/golden.hash");
    if (json){ printf("%s,\"notes\":\"hash=%.8s\"}\n",declared.c_str(),blake2b::hash256_hex(declared).c_str()); return RB.verdict?0:1; }
    printHumanB(RA,RB);
    printf("declared hash: %s\n", blake2b::hash256_hex(declared).c_str());
    return RB.verdict?0:1;
}
