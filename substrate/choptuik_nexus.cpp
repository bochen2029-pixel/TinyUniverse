// choptuik_nexus.cpp — direct Choptuik mass-scaling oracle (CPU fp64, stdlib only).
// Contract: contracts/choptuik.contract.md v1.0.0. Faithful port of the research
// evolver tournament/gamma/dss/nr_evolve.py::Ev (session-6 run-4 observables):
// polar-areal massless scalar, staggered grid, trapezoid log-metric (a^2 lags one
// node), Phi_t = d_r(w Pi) plain gradient (even ghost), Pi_t = d_r(w Phi) + 2 w Phi/r
// (odd ghost; outgoing-ish outer row), KO eps=0.5, RK4, CFL 0.4, phi(0,t) EVOLVED.
// Masses: freeze-peak (raw argmax at the freeze step) + first-crossing M70/M65 at
// 2m/r = 0.70/0.65 with SUB-CELL parabolic peak interpolation (the dr staircase fix).
// Oracle: the committed Python table gamma_scaling.npy (see contract) + the
// literature gamma band. gamma to +-0.001 is AMR-gated (D-021 triple) — NOT claimed.
// Build (BUILD.md): cl /std:c++17 /EHsc /O2 /W4 substrate\choptuik_nexus.cpp
//                      /Fe:build\choptuik_nexus.exe
#include <cstdio>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

// ----------------------------------------------------------------- blake2b-256
// (house idiom, verbatim class of substrate/fluidcss_nexus.cpp)
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

static std::string fmt9(double x){ if(!std::isfinite(x)) return "9.999999999e+99"; if(x==0.0)x=0.0; char b[48]; snprintf(b,sizeof(b),"%.9e",x); return std::string(b); }

// ---------------------------------------------------------------- constants
static const double PSTAR_800  = 0.03751655962597;
static const double PSTAR_1600 = 0.03732817692976;
static const double PSTAR_3200 = 0.03739102496155509;
static const double PI_ = 3.14159265358979323846;

// oracle anchors: tournament/gamma/dss/gamma_scaling.npy (run 4, committed) — see contract
struct Anchor { double dp, M70, M65, Mfrz; };
static const Anchor ANCH[4] = {
    {0.0010964781961431851, 0.16446089316636592, 0.15673331402621063, 0.1758678147376308 },
    {0.0022908676527677745, 0.20794151004571293, 0.19959739389584485, 0.23191592433047195},
    {0.0047863009232263802, 0.31810886500533275, 0.34128180133591657, 0.33960666294128899},
    {0.01,                  0.55267430178984989, 0.55294965651816241, 0.54268339619480543},
};

// ---------------------------------------------------------------- evolver
struct RunOut {
    int fate = 2;                 // 0 disperse, 1 bh, 2 tmax, 3 nan
    double t=0, tG=0;
    double m2r_frz=0, rH_frz=0, M_frz=0;
    double r70=0, M70=0, r65=0, M65=0;
    bool has70=false, has65=false;
};

struct Ev {
    int N; double rmax, dr, cfl, eko;
    std::vector<double> r;
    // scratch (rhs outputs/staging)
    std::vector<double> a, al, wPi, wPhi, dwPhi, koTmp, aS, alS;
    Ev(int N_, double rmax_=60.0): N(N_), rmax(rmax_), cfl(0.4), eko(0.5){
        dr = rmax/N;
        r.resize(N); for (int j=0;j<N;j++) r[j]=(j+0.5)*dr;
        a.resize(N); al.resize(N); wPi.resize(N); wPhi.resize(N);
        dwPhi.resize(N); koTmp.resize(N); aS.resize(N); alS.resize(N);
    }
    void slice_metric(const double* Phi, const double* Pi){
        double a2=1.0, fpa=0.0, fpl=0.0, acca=0.0, accl=0.0;
        for (int j=0;j<N;j++){
            double S = 2.0*PI_*r[j]*(Pi[j]*Pi[j] + Phi[j]*Phi[j]);
            double fa = (1.0-a2)/(2.0*r[j]) + S;
            double fl = (a2-1.0)/(2.0*r[j]) + S;
            double h = (j>0)? dr : r[0];
            acca += 0.5*(fpa+fa)*h;
            accl += 0.5*(fpl+fl)*h;
            a[j]=std::exp(acca); al[j]=std::exp(accl);
            a2 = std::exp(2.0*acca);
            fpa=fa; fpl=fl;
        }
        double norm = 1.0/(al[N-1]*a[N-1]);
        for (int j=0;j<N;j++) al[j]*=norm;
    }
    void rhs(const double* Phi, const double* Pi, double* dPhi, double* dPi){
        slice_metric(Phi,Pi);
        for (int j=0;j<N;j++){ double w=al[j]/a[j]; wPi[j]=w*Pi[j]; wPhi[j]=w*Phi[j]; }
        for (int j=1;j<N-1;j++) dPhi[j]=(wPi[j+1]-wPi[j-1])/(2.0*dr);
        dPhi[0]   =(wPi[1]-wPi[0])/(2.0*dr);            // even ghost
        dPhi[N-1] =(wPi[N-1]-wPi[N-2])/dr;
        for (int j=1;j<N-1;j++) dwPhi[j]=(wPhi[j+1]-wPhi[j-1])/(2.0*dr);
        dwPhi[0]   =(wPhi[1]+wPhi[0])/(2.0*dr);         // odd ghost
        dwPhi[N-1] =(wPhi[N-1]-wPhi[N-2])/dr;
        for (int j=0;j<N;j++) dPi[j]=dwPhi[j]+2.0*wPhi[j]/r[j];
        dPi[N-1] = -(wPi[N-1])/r[N-1] - (wPi[N-1]-wPi[N-2])/dr;   // outgoing-ish
    }
    void ko(double* f, double sgn){
        auto GH=[&](int i)->double{
            if (i==-2) return sgn*f[1];
            if (i==-1) return sgn*f[0];
            if (i>=N)  return f[N-1];
            return f[i];
        };
        for (int j=0;j<N;j++)
            koTmp[j]=f[j]-(eko/16.0)*(GH(j-2)-4.0*GH(j-1)+6.0*GH(j)-4.0*GH(j+1)+GH(j+2));
        memcpy(f,koTmp.data(),(size_t)N*sizeof(double));
    }
    // sub-cell parabolic peak (the dr-staircase fix; identical to the Python)
    void peak(const double* aArr, int jmax, double& rp, double& mp) const {
        double y1 = 1.0-1.0/(aArr[jmax]*aArr[jmax]);
        if (jmax>0 && jmax<N-1){
            double y0 = 1.0-1.0/(aArr[jmax-1]*aArr[jmax-1]);
            double y2 = 1.0-1.0/(aArr[jmax+1]*aArr[jmax+1]);
            double d = y0-2.0*y1+y2;
            double s = (std::fabs(d)>1e-30)? 0.5*(y0-y2)/d : 0.0;
            s = std::min(0.5, std::max(-0.5, s));
            rp = r[jmax]+s*dr;
            mp = y1-0.25*(y0-y2)*s;
        } else { rp=r[jmax]; mp=y1; }
    }
    RunOut run(double p, double tmax=200.0){
        std::vector<double> phi0(N), Phi(N), Pi(N,0.0);
        for (int j=0;j<N;j++){ double u=(r[j]-12.0)/2.0; phi0[j]=p*std::exp(-u*u); }
        Phi[0]=(phi0[1]-phi0[0])/dr;                    // np.gradient edges
        for (int j=1;j<N-1;j++) Phi[j]=(phi0[j+1]-phi0[j-1])/(2.0*dr);
        Phi[N-1]=(phi0[N-1]-phi0[N-2])/dr;
        std::vector<double> k1P(N),k1p(N),k2P(N),k2p(N),k3P(N),k3p(N),k4P(N),k4p(N),tP(N),tp(N);
        RunOut o;
        double t=0.0, tG=0.0;
        while (t<tmax){
            rhs(Phi.data(),Pi.data(),k1P.data(),k1p.data());
            aS=a; alS=al;                               // step-start metric (Python semantics)
            double wmax=0.0; for (int j=0;j<N;j++) wmax=std::max(wmax, al[j]/a[j]);
            double dt = cfl*dr/std::max(wmax,1e-9);
            for (int j=0;j<N;j++){ tP[j]=Phi[j]+0.5*dt*k1P[j]; tp[j]=Pi[j]+0.5*dt*k1p[j]; }
            rhs(tP.data(),tp.data(),k2P.data(),k2p.data());
            for (int j=0;j<N;j++){ tP[j]=Phi[j]+0.5*dt*k2P[j]; tp[j]=Pi[j]+0.5*dt*k2p[j]; }
            rhs(tP.data(),tp.data(),k3P.data(),k3p.data());
            for (int j=0;j<N;j++){ tP[j]=Phi[j]+dt*k3P[j]; tp[j]=Pi[j]+dt*k3p[j]; }
            rhs(tP.data(),tp.data(),k4P.data(),k4p.data());
            for (int j=0;j<N;j++){
                Phi[j]+=dt/6.0*(k1P[j]+2.0*k2P[j]+2.0*k3P[j]+k4P[j]);
                Pi[j] +=dt/6.0*(k1p[j]+2.0*k2p[j]+2.0*k3p[j]+k4p[j]);
            }
            ko(Phi.data(),-1.0); ko(Pi.data(),+1.0);
            t+=dt; tG+=alS[0]*dt;
            bool finite=true;
            for (int j=0;j<N;j++) if(!std::isfinite(Phi[j])){ finite=false; break; }
            if (!finite){ o.fate=3; break; }
            double mx=-1e30; int jmax=0;
            for (int j=0;j<N;j++){ double m2=1.0-1.0/(aS[j]*aS[j]); if (m2>mx){ mx=m2; jmax=j; } }
            if (!o.has65 && mx>=0.65){ o.has65=true; double rp,mp; peak(aS.data(),jmax,rp,mp); o.r65=rp; o.M65=0.5*rp*mp; }
            if (!o.has70 && mx>=0.70){ o.has70=true; double rp,mp; peak(aS.data(),jmax,rp,mp); o.r70=rp; o.M70=0.5*rp*mp; }
            if (mx>0.90 || alS[0]<0.02){
                o.fate=1; o.m2r_frz=mx; o.rH_frz=r[jmax]; o.M_frz=0.5*r[jmax]*mx;   // raw argmax (Python)
                break;
            }
            if (t>40.0){
                double imx=-1e30;
                for (int j=0;j<N/3;j++) imx=std::max(imx, 1.0-1.0/(aS[j]*aS[j]));
                if (imx<1e-3){ o.fate=0; break; }
            }
        }
        o.t=t; o.tG=tG;
        return o;
    }
};

// ---------------------------------------------------------------- fits + faces
struct LRow { double dp, M70, M65, Mfrz, r70; int fate; };

static double slope_lnln(const std::vector<double>& x, const std::vector<double>& y){
    // least-squares slope of ln y vs ln x (2x2 normal equations)
    double n=(double)x.size(), sx=0, sy=0, sxx=0, sxy=0;
    for (size_t i=0;i<x.size();i++){
        double lx=std::log(x[i]), ly=std::log(y[i]);
        sx+=lx; sy+=ly; sxx+=lx*lx; sxy+=lx*ly;
    }
    return (n*sxy-sx*sy)/(n*sxx-sx*sx);
}

static int goldenFace(const std::string& declared, const char* path){
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

static int faceScaling(bool golden, bool json){
    // frozen golden config = the campaign run-4 config EXACTLY (N=1600, p*(1600),
    // 26 pts log-spaced dp in [1e-4, 1e-2]) — the full window partially averages the
    // fine-structure wiggle (~1.3 periods); any sub-period window is wiggle-phase
    // biased by up to +-0.2 (measured: [1e-3,1e-2] alone gives a 0.55 plain slope).
    // With the bit-exact port, the golden table IS the committed oracle table.
    Ev ev(1600);
    std::vector<LRow> rows;
    for (int i=0;i<26;i++){
        double ldp = -4.0 + (double)i*(2.0/25.0);
        double dp = std::pow(10.0, ldp);
        RunOut o = ev.run(PSTAR_1600 + dp);
        rows.push_back({dp, o.M70, o.M65, o.M_frz, o.r70, o.fate});
    }
    bool gFate=true; for (auto& rw: rows) gFate = gFate && (rw.fate==1);
    std::vector<double> xs, ys; int nkeep=0;
    for (auto& rw: rows) if (rw.fate==1 && rw.r70>=8.0*ev.dr){ xs.push_back(rw.dp); ys.push_back(rw.M70); nkeep++; }
    double gam = (nkeep>=3)? slope_lnln(xs,ys) : 0.0;
    bool gGamma = (gam>0.30 && gam<0.45);
    bool gCut = (nkeep==24);      // idx 0,1 sit below the 8dr resolvability cut (measured)
    bool verdict = gFate && gGamma && gCut;
    std::string d = "{\"module\":\"choptuik_nexus\",\"ver\":\"1.0.0\",\"face\":\"scaling\",\"N\":1600,\"npts\":26";
    d += ",\"pstar\":"+fmt9(PSTAR_1600)+",\"rows\":[";
    for (size_t i=0;i<rows.size();i++){
        d += (i?",":"");
        d += "["+fmt9(rows[i].dp)+","+fmt9(rows[i].M70)+","+fmt9(rows[i].M65)+","+fmt9(rows[i].Mfrz)+"]";
    }
    d += "],\"gamma_M70\":"+fmt9(gam);
    d += std::string(",\"gates\":{\"fate\":")+(gFate?"1":"0")+",\"gamma\":"+(gGamma?"1":"0")+",\"cut\":"+(gCut?"1":"0")+"}";
    d += std::string(",\"verdict\":")+(verdict?"1":"0")+"}";
    if (golden) return goldenFace(d, "goldens/choptuik_scaling/golden.hash");
    if (json){ printf("%s\n", d.c_str()); return verdict?0:1; }
    printf("choptuik_nexus v1.0.0 - SCALING (direct Choptuik mass scaling, N=1600)\n");
    printf("-------------------------------------------------------\n");
    for (auto& rw: rows)
        printf("  dp=%.4e  M70=%.6f  M65=%.6f  Mfrz=%.6f  %s\n",
               rw.dp, rw.M70, rw.M65, rw.Mfrz, rw.fate==1?"bh":"NOT-BH");
    printf("  gamma[M70, plain slope, rH>=8dr] = %.6f   [G-GAMMA 0.30..0.45]  %s\n", gam, gGamma?"PASS":"FAIL");
    printf("  (band is wide BY DESIGN: the 2-decade window holds ~1.3 fine-structure periods =>\n");
    printf("   residual wiggle-phase bias ~+-0.05; the exact value is FROZEN in the golden hash.\n");
    printf("   Campaign quote: gamma = 0.37 +- 0.02 at N=3200, RESULTS_dss.md; +-0.001 = AMR-gated.)\n");
    printf("  G-FATE all-bh %s   G-CUT 24/26 kept (%d) %s\n", gFate?"PASS":"FAIL", nkeep, gCut?"PASS":"FAIL");
    printf("-------------------------------------------------------\n");
    printf("VERDICT: %s\n", verdict?"PASS":"FAIL");
    printf("declared hash: %s\n", blake2b::hash256_hex(d).c_str());
    return verdict?0:1;
}

static int faceCross(bool golden, bool json){
    // oracle gate: C++ vs the committed Python research table at 4 anchor dps
    Ev ev(1600);
    double worst=0.0; bool gFate=true;
    std::string d = "{\"module\":\"choptuik_nexus\",\"ver\":\"1.0.0\",\"face\":\"cross\",\"N\":1600,\"rows\":[";
    double vals[4][3];
    for (int i=0;i<4;i++){
        RunOut o = ev.run(PSTAR_1600 + ANCH[i].dp);
        gFate = gFate && (o.fate==1);
        vals[i][0]=o.M70; vals[i][1]=o.M65; vals[i][2]=o.M_frz;
        double d0=std::fabs(o.M70 /ANCH[i].M70 -1.0);
        double d1=std::fabs(o.M65 /ANCH[i].M65 -1.0);
        double d2=std::fabs(o.M_frz/ANCH[i].Mfrz-1.0);
        worst=std::max(worst,std::max(d0,std::max(d1,d2)));
        d += (i?",":"");
        d += "["+fmt9(ANCH[i].dp)+","+fmt9(o.M70)+","+fmt9(o.M65)+","+fmt9(o.M_frz)+"]";
    }
    bool gOracle = worst < 1e-6;
    bool verdict = gFate && gOracle;
    d += "],\"worst_reldev\":"+fmt9(worst);
    d += std::string(",\"gates\":{\"fate\":")+(gFate?"1":"0")+",\"oracle\":"+(gOracle?"1":"0")+"}";
    d += std::string(",\"verdict\":")+(verdict?"1":"0")+"}";
    if (golden) return goldenFace(d, "goldens/choptuik_cross/golden.hash");
    if (json){ printf("%s\n", d.c_str()); return verdict?0:1; }
    printf("choptuik_nexus v1.0.0 - CROSS (C++ vs committed Python research table)\n");
    printf("-------------------------------------------------------\n");
    for (int i=0;i<4;i++)
        printf("  dp=%.10e  M70=%.12f (py %.12f)  M65=%.12f  Mfrz=%.12f\n",
               ANCH[i].dp, vals[i][0], ANCH[i].M70, vals[i][1], vals[i][2]);
    printf("  worst relative deviation = %.3e   [G-ORACLE < 1e-6]  %s\n", worst, gOracle?"PASS":"FAIL");
    printf("  G-FATE all-bh %s\n", gFate?"PASS":"FAIL");
    printf("-------------------------------------------------------\n");
    printf("VERDICT: %s\n", verdict?"PASS":"FAIL");
    printf("declared hash: %s\n", blake2b::hash256_hex(d).c_str());
    return verdict?0:1;
}

static int faceSelftest(){
    // (a) flat-space metric at machine zero
    Ev ev(800);
    std::vector<double> z(800,0.0);
    ev.slice_metric(z.data(),z.data());
    double da=0, dl=0;
    for (int j=0;j<800;j++){ da=std::max(da,std::fabs(ev.a[j]-1.0)); dl=std::max(dl,std::fabs(ev.al[j]-1.0)); }
    bool ok1 = (da<1e-14 && dl<1e-14);
    // (b) parabolic peak interpolator exact on a synthetic parabola in 2m/r
    // choose a(r) so that m2(r) = 0.8 - 0.3*(r-rc)^2 near jmax; peak must recover (rc, 0.8)
    {
        int jc=400; double rc=ev.r[jc]+0.31*ev.dr;
        std::vector<double> aa(800,1.0);
        for (int j=jc-2;j<=jc+2;j++){
            double m2 = 0.8-0.3*(ev.r[j]-rc)*(ev.r[j]-rc);
            aa[j]=1.0/std::sqrt(1.0-m2);
        }
        double rp,mp; ev.peak(aa.data(),jc,rp,mp);
        bool okr = std::fabs(rp-rc)<1e-9, okm = std::fabs(mp-0.8)<1e-9;
        ok1 = ok1 && okr && okm;
        printf("[selftest] flat metric max|a-1|=%.1e |al-1|=%.1e  peak dr=%.1e dm=%.1e\n", da, dl, std::fabs(rp-rc), std::fabs(mp-0.8));
    }
    // (c) vacuum run disperses
    RunOut o = ev.run(0.0, 60.0);
    bool ok2 = (o.fate==0);
    printf("[selftest] vacuum run fate=%s (t=%.1f)\n", o.fate==0?"disperse":"OTHER", o.t);
    bool ok = ok1 && ok2;
    printf("VERDICT: %s\n", ok?"PASS":"FAIL");
    return ok?0:1;
}

int main(int argc, char** argv){
    bool json=false, golden=false, selftest=false, scaling=false, cross=false, ladder=false, bisect=false;
    int N=1600, npts=26, iters=40;
    double lo=-4.0, hi=-2.0, pstar=PSTAR_1600;
    std::vector<std::string> pos;
    for (int i=1;i<argc;i++){
        std::string s=argv[i];
        if      (s=="--json") json=true;
        else if (s=="--golden") golden=true;
        else if (s=="--selftest") selftest=true;
        else if (s=="--scaling") scaling=true;
        else if (s=="--cross") cross=true;
        else if (s=="--ladder") ladder=true;
        else if (s=="--bisect") bisect=true;
        else pos.push_back(s);
    }
    if (selftest) return faceSelftest();
    if (scaling)  return faceScaling(golden,json);
    if (cross)    return faceCross(golden,json);
    if (ladder){
        // --ladder N lo hi npts pstar   (research face, ungated)
        if (pos.size()>0) N=atoi(pos[0].c_str());
        if (pos.size()>1) lo=atof(pos[1].c_str());
        if (pos.size()>2) hi=atof(pos[2].c_str());
        if (pos.size()>3) npts=atoi(pos[3].c_str());
        if (pos.size()>4) pstar=atof(pos[4].c_str());
        Ev ev(N);
        std::vector<double> xs, ys;
        for (int i=0;i<npts;i++){
            double ldp = lo + (npts>1? (double)i*(hi-lo)/(npts-1) : 0.0);
            double dp = std::pow(10.0, ldp);
            RunOut o = ev.run(pstar+dp);
            printf("[%2d] dp=%.4e  M70=%.5f r70=%.4f | M65=%.5f | Mfrz=%.5f m2r=%.3f  fate=%d  t=%.1f\n",
                   i, dp, o.M70, o.r70, o.M65, o.M_frz, o.m2r_frz, o.fate, o.t);
            if (o.fate==1 && o.r70>=8.0*ev.dr){ xs.push_back(dp); ys.push_back(o.M70); }
        }
        if (xs.size()>=3) printf("gamma[M70, rH>=8dr] = %.4f  (%zu pts)\n", slope_lnln(xs,ys), xs.size());
        return 0;
    }
    if (bisect){
        // --bisect N [iters]   (reproduce p*, ungated)
        if (pos.size()>0) N=atoi(pos[0].c_str());
        if (pos.size()>1) iters=atoi(pos[1].c_str());
        Ev ev(N);
        double blo=0.003, bhi=0.02;
        RunOut rl=ev.run(blo), rh=ev.run(bhi);
        if (rl.fate!=0 || rh.fate!=1){ fprintf(stderr,"bracket invalid (fates %d,%d)\n",rl.fate,rh.fate); return 2; }
        for (int i=0;i<iters;i++){
            double mid=0.5*(blo+bhi);
            RunOut o=ev.run(mid);
            if (o.fate==1) bhi=mid; else blo=mid;
            if ((i+1)%8==0) printf("  bisect %d/%d: p* in [%.12f, %.12f]\n", i+1, iters, blo, bhi);
        }
        printf("p*(N=%d) = [%.17g, %.17g]\n", N, blo, bhi);
        return 0;
    }
    fprintf(stderr,"usage: choptuik_nexus --selftest | --scaling [--golden|--json] | --cross [--golden|--json] |\n"
                   "                      --ladder N lo hi npts pstar | --bisect N [iters]\n");
    return 2;
}
