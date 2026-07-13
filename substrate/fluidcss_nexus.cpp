// ============================================================================
//  fluidcss_nexus — TINY UNIVERSE crown (Stage A): radiation-fluid CSS background
//  Contract: contracts/fluidcss_nexus.contract.md (Stage-A object; Stage-B/beta BLOCKED)
//
//  Single-file CPU fp64 (no GPU). Constructs the Evans-Coleman critical continuously-
//  self-similar (CSS) radiation-fluid collapse BACKGROUND (Hara-Koike-Adachi
//  gr-qc/9607010; Evans-Coleman gr-qc/9402041). Deterministic C++ port of the
//  overnight-run's verified `hka_ec.py`/`hka_background.py` (which landed at machine
//  precision). gamma=4/3, c_s=1/sqrt3.  (Supersedes the walled v0.9.0 4D reduction:
//  that ODE lacked a regular center on the ingoing sound cone; this one has it.)
//
//  STAGE A (LANDED): shoot the central density oi from a REGULAR CENTER (ingoing sound
//  cone, V<0) outward to the SONIC point (Dson=0), tuning oi so the sonic velocity
//  reaches V=-c_s (the analytic Evans-Coleman crossing). The result EMERGES (nothing
//  tuned to a target):
//     oi* = 3/8,  sonic (A0,N0,om0,V0) = (3/2, 2/sqrt3, 3/4, -1/sqrt3) EXACT,
//     2m/r = 1 - 1/A0 = 1/3,  exact invariants  N = N_inf e^{-x},  A = 1 + (2/3) om.
//
//  Stage B (the perturbation eigenvalue -> beta = 1/Re(kappa0) = 0.35580192) remains an
//  HONEST WALL (the perturbation dS-coupling; RESULTS_hka_beta.md). beta is NOT reported
//  here (D-016/D-021). This tool banks the verified Stage-A background as a golden.
//
//  Build: cl /std:c++17 /EHsc /O2 /W4 substrate\fluidcss_nexus.cpp /Fe:build\fluidcss_nexus.exe
//  Faces: --stageA (default) --json|--golden|--selftest. Exit 0/1/2. Units G=c=1.
// ============================================================================

#define _CRT_SECURE_NO_WARNINGS
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <string>

static const double GG      = 4.0/3.0;            // adiabatic index (radiation)
static const double GM1     = GG - 1.0;           // 1/3
static const double CS      = 0.5773502691896257; // 1/sqrt(3) = sqrt(gamma-1)
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
//  HKA/Evans-Coleman reduced ODE system (verbatim from hka_ec.py)
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

// Integrate center -> sonic (Dson: - -> +) with fixed-step RK4. Returns the crossing
// state (interpolated) + max invariant residuals along the way.
struct Cross { bool ok; double x,N,om,V,A,lh; double invA,invN; };
static Cross shoot(double N_inf,double oi,double x0,double dx,double xmax){
    double z0=std::exp(x0);
    double N=N_inf/z0, om=oi*z0*z0, V=(MCENTER/N_inf)*z0, x=x0;
    double Dprev=Dson(N,V), invA=0.0, invN=0.0;
    Cross cr; cr.ok=false; cr.x=cr.N=cr.om=cr.V=cr.A=cr.lh=0.0;
    long nsteps=(long)((xmax-x0)/dx);
    for (long i=0;i<nsteps;i++){
        double A=A_of(N,om,V);
        double ra=std::fabs(A-(1.0+(2.0/3.0)*om));
        double rn=std::fabs(N*std::exp(x)/N_inf - 1.0);
        if (ra>invA) invA=ra; if (rn>invN) invN=rn;
        double k1N,k1o,k1V; rhs3(N,om,V,k1N,k1o,k1V);
        double k2N,k2o,k2V; rhs3(N+0.5*dx*k1N,om+0.5*dx*k1o,V+0.5*dx*k1V,k2N,k2o,k2V);
        double k3N,k3o,k3V; rhs3(N+0.5*dx*k2N,om+0.5*dx*k2o,V+0.5*dx*k2V,k3N,k3o,k3V);
        double k4N,k4o,k4V; rhs3(N+dx*k3N,om+dx*k3o,V+dx*k3V,k4N,k4o,k4V);
        double Nn=N+dx/6.0*(k1N+2*k2N+2*k3N+k4N);
        double on=om+dx/6.0*(k1o+2*k2o+2*k3o+k4o);
        double Vn=V+dx/6.0*(k1V+2*k2V+2*k3V+k4V);
        double Dcurr=Dson(Nn,Vn);
        if (Dprev<0.0 && Dcurr>=0.0){
            double frac=-Dprev/(Dcurr-Dprev);
            cr.x=x+frac*dx; cr.N=N+frac*(Nn-N); cr.om=om+frac*(on-om); cr.V=V+frac*(Vn-V);
            cr.A=A_of(cr.N,cr.om,cr.V); cr.lh=lhop(cr.A,cr.N,cr.om,cr.V);
            cr.invA=invA; cr.invN=invN; cr.ok=true; return cr;
        }
        N=Nn; om=on; V=Vn; x+=dx; Dprev=Dcurr;
        if (!std::isfinite(N)||!std::isfinite(om)||!std::isfinite(V)) break;
    }
    cr.invA=invA; cr.invN=invN; return cr;
}

// Shoot oi so the sonic-crossing velocity reaches V = -c_s (the analytic EC crossing).
// Deterministic bisection (fixed iters) on obj(oi) = V_sonic(oi) + c_s.
static double find_oi(double N_inf,double x0,double dx,double xmax,double lo,double hi,int iters, Cross& out){
    Cross rlo=shoot(N_inf,lo,x0,dx,xmax);
    double flo=rlo.ok? (rlo.V+CS) : 1.0;
    for (int i=0;i<iters;i++){
        double mid=0.5*(lo+hi);
        Cross rm=shoot(N_inf,mid,x0,dx,xmax);
        double fm=rm.ok? (rm.V+CS) : 1.0;
        if (flo*fm <= 0.0){ hi=mid; } else { lo=mid; flo=fm; }
    }
    double oistar=0.5*(lo+hi);
    out=shoot(N_inf,oistar,x0,dx,xmax);
    return oistar;
}

// ============================================================================
//  Result + JSON
// ============================================================================
struct Result {
    double N_inf=1.0, x0=0, dx=0, xmax=0; int biters=0;
    double oi=0, A0=0,N0=0,om0=0,V0=0, mr=0, xs=0, lh=0, invA=0, invN=0, oi_fine=0, conv=0;
    bool nan_free=true;
    bool g_oi=false, g_sonic=false, g_mr=false, g_inv=false, g_conv=false, verdict=false;
};

static std::string declaredJson(const Result& R){
    double N0e=2.0/std::sqrt(3.0), V0e=-1.0/std::sqrt(3.0);
    std::string s; s.reserve(1536);
    s += "{\"tool\":\"fluidcss_nexus\",\"version\":\"0.9.1\",\"units\":\"G=c=1\",\"stage\":\"A\"";
    s += ",\"params\":{\"gamma\":" + fmt6(GG) + ",\"N_inf\":" + fmt6(R.N_inf);
    s += ",\"x0\":" + fmt6(R.x0) + ",\"dx\":" + fmt9(R.dx) + ",\"xmax\":" + fmt6(R.xmax);
    s += ",\"bisect_iters\":" + std::to_string(R.biters) + "}";
    s += ",\"background\":{";
    s += "\"oi_star\":" + fmt9(R.oi) + ",\"oi_exact_3_8\":0.375000000";
    s += ",\"sonic_A0\":" + fmt9(R.A0) + ",\"sonic_N0\":" + fmt9(R.N0);
    s += ",\"sonic_om0\":" + fmt9(R.om0) + ",\"sonic_V0\":" + fmt9(R.V0);
    s += ",\"sonic_exact\":[1.500000000," + fmt9(N0e) + ",0.750000000," + fmt9(V0e) + "]";
    s += ",\"mr_2m_over_r\":" + fmt9(R.mr) + ",\"mr_exact_1_3\":0.333333333";
    s += ",\"x_sonic\":" + fmt6(R.xs) + ",\"lhop_sonic\":" + fmt9(R.lh);
    s += ",\"inv_A_resid\":" + fmt9(R.invA) + ",\"inv_N_resid\":" + fmt9(R.invN);
    s += ",\"oi_star_fine\":" + fmt9(R.oi_fine) + ",\"converge_spread\":" + fmt9(R.conv);
    s += ",\"nan_free\":" + std::string(R.nan_free?"1":"0");
    s += "},\"gates\":{";
    s += "\"G_OI\":"    + std::string(R.g_oi?"true":"false");
    s += ",\"G_SONIC\":" + std::string(R.g_sonic?"true":"false");
    s += ",\"G_MR\":"    + std::string(R.g_mr?"true":"false");
    s += ",\"G_INV\":"   + std::string(R.g_inv?"true":"false");
    s += ",\"G_CONVERGE\":" + std::string(R.g_conv?"true":"false");
    s += "},\"beta\":\"BLOCKED (Stage B honest wall; D-016/D-021)\"";
    s += ",\"verdict\":\"" + std::string(R.verdict?"pass":"fail") + "\"";
    return s;
}

static Result runStageA(){
    Result R; R.N_inf=1.0; R.x0=-12.0; R.dx=1.0e-3; R.xmax=0.5; R.biters=60;
    Cross c; R.oi = find_oi(R.N_inf,R.x0,R.dx,R.xmax, 0.35,0.40, R.biters, c);
    R.A0=c.A; R.N0=c.N; R.om0=c.om; R.V0=c.V; R.mr=1.0-1.0/c.A; R.xs=c.x; R.lh=c.lh;
    R.invA=c.invA; R.invN=c.invN;
    Cross cf; R.oi_fine = find_oi(R.N_inf,R.x0,R.dx*0.5,R.xmax, 0.35,0.40, R.biters, cf);
    R.conv = std::fabs(R.oi - R.oi_fine);
    double N0e=2.0/std::sqrt(3.0), V0e=-1.0/std::sqrt(3.0);
    R.nan_free = c.ok && std::isfinite(R.oi) && std::isfinite(R.A0);
    R.g_oi    = std::fabs(R.oi - 0.375) < 1e-3;
    R.g_sonic = std::fabs(R.A0-1.5)<2e-3 && std::fabs(R.N0-N0e)<2e-3
             && std::fabs(R.om0-0.75)<2e-3 && std::fabs(R.V0-V0e)<2e-3;
    R.g_mr    = std::fabs(R.mr - 1.0/3.0) < 1e-3;
    R.g_inv   = (R.invA < 1e-4) && (R.invN < 1e-4);
    R.g_conv  = (R.conv < 5e-4);
    R.verdict = R.nan_free && R.g_oi && R.g_sonic && R.g_mr && R.g_inv && R.g_conv;
    return R;
}

static void printHuman(const Result& R){
    double N0e=2.0/std::sqrt(3.0), V0e=-1.0/std::sqrt(3.0);
    printf("fluidcss_nexus v0.9.1 - radiation-fluid CSS critical BACKGROUND (Stage A)\n");
    printf("  HKA/Evans-Coleman, gamma=4/3, c_s=1/sqrt3, G=c=1.  (x0=%.1f dx=%.0e biters=%d)\n", R.x0, R.dx, R.biters);
    printf("-------------------------------------------------------\n");
    printf("  central density oi*     %.10f   (exact 3/8 = 0.3750000000)  [gate <1e-3]  %s\n", R.oi, R.g_oi?"PASS":"FAIL");
    printf("  sonic point (EMERGES, nothing tuned):\n");
    printf("    A0  = %.8f   exact 3/2      = 1.50000000\n", R.A0);
    printf("    N0  = %.8f   exact 2/sqrt3  = %.8f\n", R.N0, N0e);
    printf("    om0 = %.8f   exact 3/4      = 0.75000000\n", R.om0);
    printf("    V0  = %.8f   exact -1/sqrt3 = %.8f  (= -c_s: sonic on the sound cone)\n", R.V0, V0e);
    printf("    [gate |.-exact|<2e-3]  %s\n", R.g_sonic?"PASS":"FAIL");
    printf("  Misner-Sharp 2m/r = 1-1/A0 = %.10f   (exact 1/3)  [gate <1e-3]  %s\n", R.mr, R.g_mr?"PASS":"FAIL");
    printf("  exact invariants:  max|A-(1+2om/3)|=%.2e  max|N e^x/N_inf-1|=%.2e  [gate<1e-4]  %s\n",
           R.invA, R.invN, R.g_inv?"PASS":"FAIL");
    printf("  grid-convergence: oi*(dx)-oi*(dx/2) = %.2e   [gate<5e-4]  %s\n", R.conv, R.g_conv?"PASS":"FAIL");
    printf("  sonic at x_s=%.5f   lhop(sonic)=%.2e\n", R.xs, R.lh);
    printf("  beta = 1/Re(kappa0): BLOCKED (Stage B honest wall, D-016/D-021 -- not reported)\n");
    printf("-------------------------------------------------------\n");
    printf("VERDICT: %s\n", R.verdict?"PASS":"FAIL");
}

int main(int argc, char** argv){
    bool json=false, selftest=false, golden=false;
    for (int i=1;i<argc;i++){
        std::string a=argv[i];
        if      (a=="--json") json=true;
        else if (a=="--selftest") selftest=true;
        else if (a=="--golden") golden=true;
        else if (a=="--stageA") {/*default*/}
        else { fprintf(stderr,"usage: fluidcss_nexus [--stageA] [--json] [--golden] [--selftest]\n"); return 2; }
    }

    if (selftest){
        double N0e=2.0/std::sqrt(3.0), V0e=-1.0/std::sqrt(3.0), A0e=1.5, om0e=0.75;
        double d=Dson(N0e,V0e);
        double ainv=std::fabs(A0e-(1.0+(2.0/3.0)*om0e));
        double amatch=std::fabs(A_of(N0e,om0e,V0e)-A0e);
        bool ok1=std::fabs(d)<1e-12, ok2=ainv<1e-12, ok3=amatch<1e-9;
        printf("[selftest] exact sonic Dson=%.2e [%s]  A=1+2om/3 resid=%.2e [%s]  A_of match=%.2e [%s]\n",
               d,ok1?"PASS":"FAIL", ainv,ok2?"PASS":"FAIL", amatch,ok3?"PASS":"FAIL");
        printf("VERDICT: %s\n", (ok1&&ok2&&ok3)?"PASS":"FAIL");
        return (ok1&&ok2&&ok3)?0:1;
    }

    Result R = runStageA();
    std::string declared = declaredJson(R);
    std::string hash = blake2b::hash256_hex(declared);

    if (golden){
        const char* path="goldens/fluidcss_stageA/golden.hash";
        FILE* f=fopen(path,"rb");
        if (!f){ fprintf(stderr,"GOLDEN NOT FROZEN %s\n",hash.c_str()); printf("%s\n",hash.c_str()); return 2; }
        char want[128]={0};
        if (fscanf(f,"%127s",want)!=1){ fclose(f); fprintf(stderr,"GOLDEN FILE UNREADABLE\n"); return 2; }
        fclose(f);
        if (hash==std::string(want)){ fprintf(stderr,"GOLDEN OK %.8s\n",hash.c_str()); printf("%s\n",hash.c_str()); return 0; }
        fprintf(stderr,"GOLDEN MISMATCH have=%.8s want=%.8s\n",hash.c_str(),want);
        printf("%s\n",hash.c_str()); return 1;
    }
    if (json){
        std::string out = declared + ",\"notes\":\"non-declared; hash=" + hash.substr(0,8) + "\"}";
        printf("%s\n", out.c_str());
        return R.verdict?0:1;
    }
    printHuman(R);
    printf("declared hash: %s\n", hash.c_str());
    return R.verdict?0:1;
}
