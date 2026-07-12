// ============================================================================
//  curve_nexus — TINY UNIVERSE v2 N3 "curve" geometry-curves oracle tool
//  Contract: contracts/curve.contract.md v1.0.0
//
//  Single-file CPU fp64 (no GPU — like N0 substrate_nexus; runs under any card
//  contention). Geodesics through the substrate's weak-field metric
//    ds^2 = -(1 + 2 Phi/c^2) c^2 dt^2 + (1 - 2 Phi/c^2)(dx^2+dy^2+dz^2),  Phi=-GM/r
//  reproduce the exact-GR curvature observables:
//    deflect : light deflection 4GM/(b c^2) -- DECOMPOSED into the N2-lapse half
//              (time, 2GM/bc^2) + the N3-space half (the 1919 factor of 2).
//    precess : perihelion precession 6 pi GM/(c^2 a(1-e^2)) per orbit.
//
//  Envelope face (--scenario X --json|--golden|--selftest). Exit 0/1/2.
//
//  Build (BUILD.md CPU path):
//    cl /std:c++17 /EHsc /O2 /W4 substrate\curve_nexus.cpp /Fe:build\curve_nexus.exe
//
//  Dials (nexus v0, frozen): c=20, G=2e-3, M=1e5 (r_s=2GM/c^2=1 su). hbar/m/dt inert.
//  Determinism: fp64 RK4, fixed steps -> byte-identical declared JSON + blake2b.
// ============================================================================

#define _CRT_SECURE_NO_WARNINGS
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <string>

static constexpr double PI = 3.14159265358979323846;

struct Dials { double c = 20.0; double G = 2.0e-3; double hbar = 0.5; double m = 1.0; };

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
    for (int i = 0; i < 8; i++){ v[i] = h[i]; v[i + 8] = IV[i]; }
    v[12] ^= t;
    if (last) v[14] = ~v[14];
    for (int i = 0; i < 16; i++){
        uint64_t w = 0;
        for (int j = 7; j >= 0; j--) w = (w << 8) | block[i*8 + j];
        m[i] = w;
    }
    for (int r = 0; r < 12; r++){
        const uint8_t* s = SIGMA[r % 10];
        Gm(v, 0, 4,  8, 12, m[s[ 0]], m[s[ 1]]);
        Gm(v, 1, 5,  9, 13, m[s[ 2]], m[s[ 3]]);
        Gm(v, 2, 6, 10, 14, m[s[ 4]], m[s[ 5]]);
        Gm(v, 3, 7, 11, 15, m[s[ 6]], m[s[ 7]]);
        Gm(v, 0, 5, 10, 15, m[s[ 8]], m[s[ 9]]);
        Gm(v, 1, 6, 11, 12, m[s[10]], m[s[11]]);
        Gm(v, 2, 7,  8, 13, m[s[12]], m[s[13]]);
        Gm(v, 3, 4,  9, 14, m[s[14]], m[s[15]]);
    }
    for (int i = 0; i < 8; i++) h[i] ^= v[i] ^ v[i + 8];
}
static std::string hash256_hex(const std::string& in){
    uint64_t h[8];
    for (int i = 0; i < 8; i++) h[i] = IV[i];
    h[0] ^= 0x01010000ull ^ 32ull;
    size_t n = in.size(), off = 0;
    uint8_t block[128];
    while (n - off > 128){ memcpy(block, in.data() + off, 128); off += 128; compress(h, block, (uint64_t)off, false); }
    size_t rem = n - off;
    memset(block, 0, 128);
    if (rem) memcpy(block, in.data() + off, rem);
    compress(h, block, (uint64_t)n, true);
    char hex[65];
    for (int i = 0; i < 32; i++){ unsigned byte = (unsigned)((h[i / 8] >> (8 * (i % 8))) & 0xFF); snprintf(hex + 2*i, 3, "%02x", byte); }
    return std::string(hex, 64);
}
} // namespace blake2b

static std::string fmt6(double x){ if(!std::isfinite(x)) return "9999999.999999"; if(x==0.0)x=0.0; char b[48]; snprintf(b,sizeof(b),"%.6f",x); return std::string(b); }
static std::string fmt9(double x){ if(!std::isfinite(x)) return "9.999999999e+99"; if(x==0.0)x=0.0; char b[48]; snprintf(b,sizeof(b),"%.9e",x); return std::string(b); }

// ============================================================================
//  GEODESIC INTEGRATORS (fp64 RK4)
// ============================================================================

// Light deflection of a photon at impact parameter b through the effective index
// of refraction n(r) = 1 + K*GM/(r c^2)  [K=2 full metric (space+time), K=1 lapse-
// only (time = N2's contribution)]. Eikonal ray ODE, arc-length s:
//   dx/ds = t_hat ,  dt_hat/ds = (grad n - (t_hat . grad n) t_hat)/n .
// Ray launched at (-X0, b) with t_hat=(1,0), integrated to x=+X0 (X0 >> b). The
// deflection is the angle t_hat has turned. Returns the (positive) bend angle.
static double rayDeflection(double GM, double c, double b, double K, double X0, long steps){
    double x=-X0, y=b, tx=1.0, ty=0.0;
    double ds = (2.0*X0)/(double)steps;
    double c2 = c*c;
    auto acc = [&](double x,double y,double tx,double ty,double& dtx,double& dty){
        double r2=x*x+y*y; double r=std::sqrt(r2); double r3=r2*r;
        double n = 1.0 + K*GM/(r*c2);
        double gnx = -K*GM/(r3*c2)*x;         // grad n
        double gny = -K*GM/(r3*c2)*y;
        double tg = tx*gnx + ty*gny;
        dtx = (gnx - tg*tx)/n;
        dty = (gny - tg*ty)/n;
    };
    for (long i=0;i<steps;i++){
        // RK4 on (x,y,tx,ty); dx/ds=tx, dy/ds=ty, dt/ds=acc
        double k1x=tx,k1y=ty,k1tx,k1ty; acc(x,y,tx,ty,k1tx,k1ty);
        double k2x=tx+0.5*ds*k1tx, k2y=ty+0.5*ds*k1ty, k2tx,k2ty;
        acc(x+0.5*ds*k1x, y+0.5*ds*k1y, tx+0.5*ds*k1tx, ty+0.5*ds*k1ty, k2tx,k2ty);
        double k3x=tx+0.5*ds*k2tx, k3y=ty+0.5*ds*k2ty, k3tx,k3ty;
        acc(x+0.5*ds*k2x, y+0.5*ds*k2y, tx+0.5*ds*k2tx, ty+0.5*ds*k2ty, k3tx,k3ty);
        double k4x=tx+ds*k3tx, k4y=ty+ds*k3ty, k4tx,k4ty;
        acc(x+ds*k3x, y+ds*k3y, tx+ds*k3tx, ty+ds*k3ty, k4tx,k4ty);
        x  += ds/6.0*(k1x + 2*k2x + 2*k3x + k4x);
        y  += ds/6.0*(k1y + 2*k2y + 2*k3y + k4y);
        tx += ds/6.0*(k1tx+ 2*k2tx+ 2*k3tx+ k4tx);
        ty += ds/6.0*(k1ty+ 2*k2ty+ 2*k3ty+ k4ty);
        double tn = std::sqrt(tx*tx+ty*ty); tx/=tn; ty/=tn;   // keep t_hat unit
    }
    return std::fabs(std::atan2(ty, tx));   // turn angle from the +x axis
}

// Perihelion precession per orbit of a bound timelike geodesic, weak-field metric.
// Orbit equation (1PN): u'' + u = GM/L^2 + (3GM/c^2) u^2, u=1/r, ' = d/dphi.
// L^2 = GM a(1-e^2). Integrate phi from perihelion; the advance to the next
// perihelion beyond 2*pi is the precession Delta-omega = 6 pi GM/(c^2 a(1-e^2)).
static double orbitPrecession(double GM, double c, double a, double e, long steps){
    double p  = a*(1.0-e*e);
    double C  = 1.0/p;                       // GM/L^2 with L^2=GM p
    double eps= 3.0*GM/(c*c);
    double u = 1.0/(a*(1.0-e));               // perihelion: u max
    double w = 0.0;                            // du/dphi = 0 at perihelion
    double dphi = (2.0*PI + 0.6)/(double)steps;
    auto dw = [&](double u){ return C + eps*u*u - u; };
    double phi=0.0, prevw=w, prevphi=0.0;
    for (long i=0;i<steps;i++){
        // RK4 on (u,w): du/dphi=w, dw/dphi=C+eps u^2 - u
        double k1u=w,            k1w=dw(u);
        double k2u=w+0.5*dphi*k1w, k2w=dw(u+0.5*dphi*k1u);
        double k3u=w+0.5*dphi*k2w, k3w=dw(u+0.5*dphi*k2u);
        double k4u=w+dphi*k3w,     k4w=dw(u+dphi*k3u);
        double un = u + dphi/6.0*(k1u+2*k2u+2*k3u+k4u);
        double wn = w + dphi/6.0*(k1w+2*k2w+2*k3w+k4w);
        double phin = phi + dphi;
        // next perihelion = w crosses + -> - with phi > pi
        if (phin > PI && prevw > 0.0 && wn <= 0.0){
            double frac = prevw/(prevw - wn);          // linear interp of the w zero
            double phi_cross = prevphi + frac*(phin - prevphi);
            return phi_cross - 2.0*PI;
        }
        prevw = wn; prevphi = phin; u=un; w=wn; phi=phin;
    }
    return std::nan("");   // no perihelion found (should not happen)
}

// ============================================================================
//  Result + JSON
// ============================================================================
struct Result {
    std::string scenario;
    double c=0, GM=0, r_s=0;
    // deflect
    double b1=0, b2=0, X0=0; long dsteps=0;
    double df1_full=0, df1_lapse=0, df2_full=0, df2_lapse=0;   // measured deflections
    double rf1=0, rl1=0, rf2=0, rl2=0;                         // ratios vs GR (full/4, lapse/2)
    double doubling1=0, doubling2=0;                           // full/lapse (~2)
    // precess
    double a=0, e=0; long psteps=0;
    double pr_meas=0, pr_exact=0, pr_rel=0;
    bool nan_free=true, gate_primary=false, gate_secondary=false, gate_nan=false, verdict=false;
};

static std::string declaredJson(const Dials& D, uint64_t seed, const Result& R){
    std::string s; s.reserve(1536);
    s += "{\"tool\":\"curve_nexus\",\"version\":\"1.0.0\",\"seed\":" + std::to_string(seed);
    s += ",\"params\":{\"scenario\":\"" + R.scenario + "\"";
    s += ",\"c\":" + fmt6(D.c) + ",\"G\":" + fmt9(D.G) + ",\"GM\":" + fmt6(R.GM) + ",\"r_s\":" + fmt6(R.r_s);
    if (R.scenario == "deflect"){
        s += ",\"b1\":" + fmt6(R.b1) + ",\"b2\":" + fmt6(R.b2) + ",\"X0\":" + fmt6(R.X0) + ",\"steps\":" + std::to_string(R.dsteps);
        s += "},\"result\":{";
        s += "\"defl_full_b1\":"  + fmt9(R.df1_full)  + ",\"defl_lapse_b1\":" + fmt9(R.df1_lapse);
        s += ",\"defl_full_b2\":" + fmt9(R.df2_full)  + ",\"defl_lapse_b2\":" + fmt9(R.df2_lapse);
        s += ",\"ratio_full_b1\":"  + fmt6(R.rf1) + ",\"ratio_lapse_b1\":" + fmt6(R.rl1);
        s += ",\"ratio_full_b2\":"  + fmt6(R.rf2) + ",\"ratio_lapse_b2\":" + fmt6(R.rl2);
        s += ",\"doubling_b1\":" + fmt6(R.doubling1) + ",\"doubling_b2\":" + fmt6(R.doubling2);
    } else {
        s += ",\"a\":" + fmt6(R.a) + ",\"e\":" + fmt6(R.e) + ",\"steps\":" + std::to_string(R.psteps);
        s += "},\"result\":{";
        s += "\"precess_meas\":" + fmt9(R.pr_meas) + ",\"precess_exact\":" + fmt9(R.pr_exact);
        s += ",\"precess_rel\":" + fmt9(R.pr_rel);
    }
    s += ",\"nan_free\":" + std::string(R.nan_free?"1":"0");
    s += "},\"gates\":{\"primary\":" + std::string(R.gate_primary?"true":"false");
    s += ",\"secondary\":" + std::string(R.gate_secondary?"true":"false");
    s += ",\"nan\":" + std::string(R.gate_nan?"true":"false");
    s += "},\"verdict\":\"" + std::string(R.verdict?"pass":"fail") + "\"";
    return s;
}

// ============================================================================
//  SCENARIOS
// ============================================================================
static Result runDeflect(const Dials& D){
    Result R; R.scenario = "deflect";
    const double M = 1.0e5; R.GM = D.G*M; R.r_s = 2.0*R.GM/(D.c*D.c);
    R.b1 = 100.0; R.b2 = 200.0; R.X0 = 5000.0; R.dsteps = 2000000;
    double GM=R.GM, c=D.c;
    R.df1_full  = rayDeflection(GM, c, R.b1, 2.0, R.X0, R.dsteps);
    R.df1_lapse = rayDeflection(GM, c, R.b1, 1.0, R.X0, R.dsteps);
    R.df2_full  = rayDeflection(GM, c, R.b2, 2.0, R.X0, R.dsteps);
    R.df2_lapse = rayDeflection(GM, c, R.b2, 1.0, R.X0, R.dsteps);
    double gr1 = 4.0*GM/(R.b1*c*c), gr2 = 4.0*GM/(R.b2*c*c);   // full GR deflection
    R.rf1 = R.df1_full/gr1;  R.rl1 = R.df1_lapse/(0.5*gr1);     // lapse-only vs 2GM/bc^2
    R.rf2 = R.df2_full/gr2;  R.rl2 = R.df2_lapse/(0.5*gr2);
    R.doubling1 = R.df1_full/R.df1_lapse;  R.doubling2 = R.df2_full/R.df2_lapse;
    R.nan_free = std::isfinite(R.df1_full)&&std::isfinite(R.df2_full)&&std::isfinite(R.df1_lapse)&&std::isfinite(R.df2_lapse);
    double tol = 0.02;
    R.gate_primary = std::fabs(R.rf1-1)<tol && std::fabs(R.rf2-1)<tol;   // full = exact GR 4GM/bc^2
    R.gate_secondary = std::fabs(R.rl1-1)<tol && std::fabs(R.rl2-1)<tol; // lapse = exactly half (space doubles it)
    R.gate_nan = R.nan_free;
    R.verdict = R.gate_primary && R.gate_secondary && R.gate_nan;
    return R;
}

static Result runPrecess(const Dials& D){
    Result R; R.scenario = "precess";
    const double M = 1.0e5; R.GM = D.G*M; R.r_s = 2.0*R.GM/(D.c*D.c);
    R.a = 1000.0; R.e = 0.5; R.psteps = 2000000;
    double GM=R.GM, c=D.c, p=R.a*(1.0-R.e*R.e);
    R.pr_meas  = orbitPrecession(GM, c, R.a, R.e, R.psteps);
    R.pr_exact = 6.0*PI*GM/(c*c*p);
    R.pr_rel   = std::fabs(R.pr_meas/R.pr_exact - 1.0);
    R.nan_free = std::isfinite(R.pr_meas);
    R.gate_primary = (R.pr_rel < 0.02);
    R.gate_secondary = true;
    R.gate_nan = R.nan_free;
    R.verdict = R.gate_primary && R.gate_nan;
    return R;
}

static void printHuman(const Result& R){
    printf("curve_nexus v1.0.0 - TINY UNIVERSE v2 N3 curve (geometry curves)\n");
    printf("scenario: %s   GM=%.3f  r_s=2GM/c^2=%.4f su\n", R.scenario.c_str(), R.GM, R.r_s);
    printf("-------------------------------------------------------\n");
    if (R.scenario == "deflect"){
        printf("  light deflection: full metric vs exact GR 4GM/bc^2; lapse-only vs 2GM/bc^2\n");
        printf("    b=%.0f:  full  %.6e (x%.4f GR)   lapse %.6e (x%.4f half)\n", R.b1, R.df1_full, R.rf1, R.df1_lapse, R.rl1);
        printf("    b=%.0f:  full  %.6e (x%.4f GR)   lapse %.6e (x%.4f half)\n", R.b2, R.df2_full, R.rf2, R.df2_lapse, R.rl2);
        printf("    space DOUBLES the bending:  full/lapse = %.4f, %.4f  (GR: 2.0 exactly)\n", R.doubling1, R.doubling2);
        printf("    full = exact GR 4GM/bc^2      [gate |ratio-1|<0.02]  %s\n", R.gate_primary?"PASS":"FAIL");
        printf("    lapse = half (the N2 time part) [gate |ratio-1|<0.02]  %s\n", R.gate_secondary?"PASS":"FAIL");
    } else {
        printf("  perihelion precession (a=%.0f, e=%.2f) vs exact GR 6*pi*GM/(c^2 a(1-e^2)):\n", R.a, R.e);
        printf("    measured  %.6e rad/orbit\n", R.pr_meas);
        printf("    exact GR  %.6e rad/orbit\n", R.pr_exact);
        printf("    rel error %.3e   [gate < 0.02]  %s\n", R.pr_rel, R.gate_primary?"PASS":"FAIL");
    }
    printf("  nan_free      %s\n", R.nan_free?"yes":"NO");
    printf("-------------------------------------------------------\n");
    printf("VERDICT: %s\n", R.verdict?"PASS":"FAIL");
}

int main(int argc, char** argv){
    Dials D; uint64_t seed = 20260711ull;
    std::string scenario = "deflect";
    bool json=false, selftest=false, golden=false;
    for (int i=1;i<argc;i++){
        std::string a=argv[i];
        if      (a=="--json") json=true;
        else if (a=="--selftest") selftest=true;
        else if (a=="--golden") golden=true;
        else if (a=="--scenario" && i+1<argc) scenario=argv[++i];
        else if (a=="--seed" && i+1<argc) seed=strtoull(argv[++i],nullptr,10);
        else { fprintf(stderr,"usage: curve_nexus --scenario deflect|precess [--seed N] [--json] [--golden] [--selftest]\n"); return 2; }
    }

    if (selftest){
        // flat (GM=0): the ray goes straight (zero deflection) and the orbit closes.
        double defl = rayDeflection(0.0, D.c, 100.0, 2.0, 5000.0, 200000);
        double prec = orbitPrecession(0.0, D.c, 100.0, 0.5, 200000);
        bool ok1 = (std::fabs(defl) < 1e-9);
        bool ok2 = (std::isfinite(prec) && std::fabs(prec) < 1e-6);
        printf("[selftest] flat-space ray straight   deflection=%.3e  [%s]\n", defl, ok1?"PASS":"FAIL");
        printf("[selftest] flat-space orbit closed   precession=%.3e  [%s]\n", prec, ok2?"PASS":"FAIL");
        printf("VERDICT: %s\n", (ok1&&ok2)?"PASS":"FAIL");
        return (ok1&&ok2)?0:1;
    }

    if (golden){ D = Dials(); seed = 20260711ull; }

    Result R = (scenario=="deflect") ? runDeflect(D) : (scenario=="precess") ? runPrecess(D)
             : (fprintf(stderr,"error: unknown scenario '%s' (deflect|precess)\n", scenario.c_str()), std::exit(2), Result());
    std::string declared = declaredJson(D, seed, R);
    std::string hash = blake2b::hash256_hex(declared);

    if (golden){
        std::string path = "goldens/curve_" + R.scenario + "/golden.hash";
        FILE* f = fopen(path.c_str(), "rb");
        if (!f){ fprintf(stderr, "GOLDEN NOT FROZEN %s\n", hash.c_str()); printf("%s\n", hash.c_str()); return 2; }
        char want[128]={0};
        if (fscanf(f, "%127s", want)!=1){ fclose(f); fprintf(stderr, "GOLDEN FILE UNREADABLE\n"); return 2; }
        fclose(f);
        if (hash == std::string(want)){ fprintf(stderr, "GOLDEN OK %.8s\n", hash.c_str()); printf("%s\n", hash.c_str()); return 0; }
        fprintf(stderr, "GOLDEN MISMATCH have=%.8s want=%.8s\n", hash.c_str(), want);
        printf("%s\n", hash.c_str()); return 1;
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
