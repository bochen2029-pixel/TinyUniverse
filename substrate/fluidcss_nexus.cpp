// ============================================================================
//  fluidcss_nexus — radiation-fluid CSS critical exponent (eigenvalue route)
//  Contract: contracts/fluidcss_nexus.contract.md v0.9.0
//
//  Single file, C++17, fp64, stdlib only, NO GPU. Continuously-self-similar
//  (Evans-Coleman / Koike-Hara-Adachi) critical collapse of a p=rho/3 radiation
//  fluid. Target (KHA95, gr-qc/9503007): Re kappa0 = 2.81055255, beta = 0.35580192.
//
//  Units: geometric G = c = 1 (beta is a dimensionless universal number).
//
//  STATUS (honesty, D-016/D-021): the sonic-point regularity condition (the §1.6
//  derivation) is DERIVED and VERIFIED (see RESULTS_fluidcss.md). Stage A builds
//  the background critical CSS solution by shooting from the regular center through
//  the sonic point. Stage B (the perturbation eigenvalue that yields beta) is NOT
//  yet landing beta; this tool freezes the VERIFIED Stage-A object and reports the
//  measured background quantities. beta is emitted as NOT-MEASURED (nan) until the
//  eigenvalue shoot converges with G-ANCHOR/G-CONVERGE/G-UNIQUE firing. No beta is
//  faked (a faked eigenvalue poisons the oracle farm).
//
//  The fluid ODE RHS is the FULL 4D covariant reduction of nabla_a T^{ab}=0 in the
//  KHA metric (the angular-pressure terms Gamma^r_{th th}T^{th th}+... included);
//  this is verified to (a) share the KHA sonic locus 3N^2V^2-N^2+4NV-V^2+3=0 and
//  (b) give the correct regular center dV/dx->0 as N->inf (the transcribed KHA eq.18
//  fluid pair does NOT — see RESULTS). Metric slopes match Evans-Coleman eqs (4),(5).
//
//  Build (MSVC, golden platform):
//    cl /std:c++17 /EHsc /O2 /W4 substrate\fluidcss_nexus.cpp /Fe:build\fluidcss_nexus.exe
//
//  Exit: 0 all gates pass / 1 a declared gate fired (incl. verdict "blocked") / 2 error.
//  Determinism: (params) -> byte-identical declared JSON (notes excluded).
// ============================================================================
#define _CRT_SECURE_NO_WARNINGS
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

static constexpr double PI = 3.14159265358979323846;

// ----------------------------------------------------------------------------
// BLAKE2b-256 (RFC 7693) — lifted verbatim from substrate_nexus.cpp (D-020).
// ----------------------------------------------------------------------------
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
static inline void G(uint64_t* v, int a, int b, int c, int d, uint64_t x, uint64_t y){
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
        G(v, 0, 4,  8, 12, m[s[ 0]], m[s[ 1]]);
        G(v, 1, 5,  9, 13, m[s[ 2]], m[s[ 3]]);
        G(v, 2, 6, 10, 14, m[s[ 4]], m[s[ 5]]);
        G(v, 3, 7, 11, 15, m[s[ 6]], m[s[ 7]]);
        G(v, 0, 5, 10, 15, m[s[ 8]], m[s[ 9]]);
        G(v, 1, 6, 11, 12, m[s[10]], m[s[11]]);
        G(v, 2, 7,  8, 13, m[s[12]], m[s[13]]);
        G(v, 3, 4,  9, 14, m[s[14]], m[s[15]]);
    }
    for (int i = 0; i < 8; i++) h[i] ^= v[i] ^ v[i + 8];
}
static std::string hash256_hex(const std::string& in){
    uint64_t h[8];
    for (int i = 0; i < 8; i++) h[i] = IV[i];
    h[0] ^= 0x01010000ull ^ 32ull;
    size_t n = in.size(), off = 0;
    uint8_t block[128];
    while (n - off > 128){
        memcpy(block, in.data() + off, 128);
        off += 128;
        compress(h, block, (uint64_t)off, false);
    }
    size_t rem = n - off;
    memset(block, 0, 128);
    if (rem) memcpy(block, in.data() + off, rem);
    compress(h, block, (uint64_t)n, true);
    char hex[65];
    for (int i = 0; i < 32; i++){
        unsigned byte = (unsigned)((h[i / 8] >> (8 * (i % 8))) & 0xFF);
        snprintf(hex + 2*i, 3, "%02x", byte);
    }
    return std::string(hex, 64);
}
} // namespace blake2b

static std::string fmtnum(double x){
    if (!std::isfinite(x)) return "9.999999999e+99";
    if (x == 0.0) x = 0.0;
    char buf[40];
    snprintf(buf, sizeof(buf), "%.9e", x);
    return std::string(buf);
}

// ============================================================================
// The autonomous CSS ODE system (KHA reduced variables N,A,omega,V; independent x).
//   State Y = {N, A, omega, V}.  Center x->-inf: A=1, V=0, omega->0, N~e^{-x}.
//   Sonic point: D(N,V) = 3 N^2 V^2 - N^2 + 4 N V - V^2 + 3 = 0  (gauge x=0).
// ----------------------------------------------------------------------------
// Metric slopes (VERIFIED vs Evans-Coleman eqs 4,5):
static inline double dAdx(double /*N*/, double A, double om, double V){
    return A*(1.0 - A + (2.0*om/(1.0 - V*V))*(1.0 + V*V/3.0));
}
static inline double dNdx(double N, double A, double om, double /*V*/){
    return N*(-2.0 + A - (2.0/3.0)*om);
}
// Sonic locus:
static inline double Dloc(double N, double V){
    return 3.0*N*N*V*V - N*N + 4.0*N*V - V*V + 3.0;
}
// Fluid slopes (FULL 4D covariant reduction of nabla_a T^{ab}=0; CSE-generated).
// Returns d(omega)/dx and d(V)/dx via the out-params. Singular where D->0.
static void dFluid(double N0, double A0, double om0, double V0, double& dom, double& dV){
    const double x0 = pow(V0, 4);
    const double x1 = pow(V0, 2);
    const double x2 = 4*x1;
    const double x3 = pow(V0, 3);
    const double x4 = 4*N0;
    const double x5 = pow(N0, 2);
    const double x6 = V0*x4 + x5 - 3;
    const double x7 = A0*V0;
    const double x8 = 2*om0;
    const double x9 = x1*x5;
    const double x10 = 24*x9;
    const double x11 = A0*N0;
    const double x12 = 2*x11;
    const double x13 = A0*x3;
    const double x14 = N0*x1;
    const double x15 = 3*x5;
    const double x16 = (2.0/3.0)*om0;
    dom = (1.0/3.0)*om0*(12*A0*N0*x3 - 15*A0*x0*x5 - 3*A0*x0 + 24*A0*x1*x5 - 9*A0*x5 + 3*A0
          + 24*N0*V0*om0 + 8*N0*om0*x3 - 12*N0*x7 + 6*om0*x0*x5 - 8*om0*x1 - om0*x10 + 2*om0*x5
          - 6*om0 + 9*x0*x5 - x0*x8 + 3*x0 - x10 + 15*x5 - 3)
          /(3*x0*x5 - x0 - x2*x5 + x2 - x3*x4 + x6);
    dV = ((3.0/2.0)*N0*x0 + (7.0/2.0)*N0 - V0*x15 + V0*x16*x5 + V0*x8 + V0 - 16.0/3.0*om0*x14
          - x0*x12 + x11*x2 - x12 - x13*x5 + x13 - 5*x14 + x15*x3 + x16*x3 + x3*x5*x8 - x3
          + x5*x7 - x7)
          /(x1 + x6 - 3*x9);
}
// Full RHS dY/dx:
static void rhs(const double Y[4], double dY[4]){
    double N = Y[0], A = Y[1], om = Y[2], V = Y[3];
    dY[0] = dNdx(N, A, om, V);
    dY[1] = dAdx(N, A, om, V);
    dFluid(N, A, om, V, dY[2], dY[3]);
}
static void rk4step(double Y[4], double h){
    double k1[4], k2[4], k3[4], k4[4], t[4];
    rhs(Y, k1);
    for (int i=0;i<4;i++) t[i]=Y[i]+0.5*h*k1[i]; rhs(t,k2);
    for (int i=0;i<4;i++) t[i]=Y[i]+0.5*h*k2[i]; rhs(t,k3);
    for (int i=0;i<4;i++) t[i]=Y[i]+h*k3[i];     rhs(t,k4);
    for (int i=0;i<4;i++) Y[i]+=(h/6.0)*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
}

// Center seed (regular center, INGOING per Evans-Coleman): at z=e^{x0}->0,
//   N = nc/z,  A = 1 + a2 z^2,  omega = 1.5 a2 z^2,  V = -(1.5/nc) z.
// (Relations VERIFIED against the RHS to machine precision: dA/dx=(2w-a2)z^2 with w=1.5a2,
//  dV/dx = -v1 z, dN/dx = -N.) nc is an x-translation gauge; a2 (central density) is the
//  physical shooting parameter.
struct Seed { double nc, a2, z0; };
static void seed_state(const Seed& s, double Y[4]){
    double z = s.z0;
    Y[0] = s.nc / z;
    Y[1] = 1.0 + s.a2*z*z;
    Y[2] = 1.5*s.a2*z*z;
    Y[3] = -(1.5/s.nc)*z;
}

// Integrate from the center seed; record the state at the first sonic-locus crossing
// (D sign change on the outgoing V>thr branch) and diagnostics. Returns:
//   code 0 = reached sonic; 1 = blew up (collapse) before sonic; 2 = ran to xmax w/o sonic.
struct BgResult {
    int code = 2;
    double x_sonic=0, N_s=0, A_s=0, om_s=0, V_s=0;   // state at sonic crossing
    double Vp_near=0;                                  // dV/dx just before sonic (analyticity probe)
    double x_end=0, N_e=0, A_e=0, om_e=0, V_e=0;       // final state
    double Vmin=0;                                     // most-ingoing V (should be <=0)
};
static BgResult integrate_bg(const Seed& s, double xmax, double h, double Dthr){
    BgResult r;
    double Y[4]; seed_state(s, Y);
    double x = std::log(s.z0);
    r.Vmin = Y[3];
    double Dprev = Dloc(Y[0], Y[3]);
    long nmax = (long)((xmax - x)/h);
    for (long i=0;i<nmax;i++){
        double dY[4]; rhs(Y,dY);
        // analyticity probe: dV/dx as we approach the sonic locus on the outgoing branch
        double Dcur = Dloc(Y[0], Y[3]);
        if (Dcur > 0 && Dcur < Dthr && Y[3] > 0.05){
            r.code = 0; r.x_sonic = x;
            r.N_s=Y[0]; r.A_s=Y[1]; r.om_s=Y[2]; r.V_s=Y[3];
            r.Vp_near = dY[3];
            r.x_end=x; r.N_e=Y[0]; r.A_e=Y[1]; r.om_e=Y[2]; r.V_e=Y[3];
            return r;
        }
        rk4step(Y, h); x += h;
        if (Y[3] < r.Vmin) r.Vmin = Y[3];
        if (!std::isfinite(Y[0])||!std::isfinite(Y[3])||std::fabs(Y[3])>=1.0||Y[1]<=0.0){
            r.code = 1; r.x_end=x; r.N_e=Y[0]; r.A_e=Y[1]; r.om_e=Y[2]; r.V_e=Y[3];
            return r;
        }
        (void)Dprev; Dprev = Dcur;
    }
    r.code = 2; r.x_end=x; r.N_e=Y[0]; r.A_e=Y[1]; r.om_e=Y[2]; r.V_e=Y[3];
    return r;
}

// ============================================================================
// Result plumbing (tiny_nexus idiom)
// ============================================================================
struct Metric { std::string k; double v; };

struct Report {
    // sonic-locus coefficients (declared, exact): 3 N^2 V^2 - N^2 + 4 N V - V^2 + 3
    // background critical solution (Stage A)
    bool   bg_reached_sonic = false;
    double a2_star = 0, nc = 0;
    double x_sonic=0, N_sonic=0, A_sonic=0, om_sonic=0, V_sonic=0, Dresid=0;
    double Vmin=0;
    // convergence probe (Stage A): sonic V under grid refinement
    double conv_spread = 0;
    // Stage B (eigenvalue) — NOT yet measured
    double beta = std::nan("");
    double re_kappa0 = std::nan("");
    double im_kappa0 = std::nan("");
    int    n_relevant_modes = -1;   // -1 = not computed
    double gauge_mode_kappa = 0.35699;  // KHA-declared spurious mode (reference)
    // gates
    bool G_ANCHOR=false, G_CONVERGE=false, G_UNIQUE=false;
    std::string verdict = "blocked";
};

// verify the sonic-point regularity locus + metric slopes vs Evans-Coleman (selftest anchors)
static bool selftest(std::string& why){
    // 1. metric slope A_,x/A at a probe vs Evans-Coleman (4): A_,x/A = 1-A+2 om[4/(3(1-V^2))-1/3]
    {
        double N=2.0,A=1.3,om=0.2,V=0.25;
        double lhs = dAdx(N,A,om,V)/A;
        double rhs_ec = 1.0 - A + 2.0*om*(4.0/(3.0*(1.0-V*V)) - 1.0/3.0);
        if (std::fabs(lhs - rhs_ec) > 1e-12){ why="metric A-slope != Evans-Coleman eq(4)"; return false; }
    }
    // 2. N slope vs Evans-Coleman (5) reduced form: N_,x/N = -2 + A - 2 om/3
    {
        double N=2.0,A=1.3,om=0.2,V=0.25;
        double lhs = dNdx(N,A,om,V)/N;
        double rhs_ec = -2.0 + A - (2.0/3.0)*om;
        if (std::fabs(lhs - rhs_ec) > 1e-12){ why="metric N-slope != Evans-Coleman eq(5)"; return false; }
    }
    // 3. sonic locus vanishes at a known point on it (V=0 => N=sqrt3)
    {
        if (std::fabs(Dloc(std::sqrt(3.0), 0.0)) > 1e-12){ why="sonic locus D(sqrt3,0)!=0"; return false; }
    }
    // 4. regular center: dV/dx -> 0 as N->inf at (A=1,V=0,om->0)
    {
        double dom,dV; dFluid(1e6, 1.0, 1e-10, 0.0, dom, dV);
        if (std::fabs(dV) > 1e-4){ why="center dV/dx not ->0 (regularity broken)"; return false; }
    }
    // 5. blake2b known-answer (empty string)
    {
        std::string h = blake2b::hash256_hex("");
        if (h.substr(0,16) != "0e5751c026e543b2"){ why="blake2b KAT failed"; return false; }
    }
    why = "ok";
    return true;
}

// Bracket + bisect the background shoot on a2 for analytic sonic passage.
// Returns a2_star and fills the sonic-point state; sets bg_reached_sonic.
static void run_stageA(Report& R, double nc, double z0, double a2lo, double a2hi,
                       int nscan, double xmax, double h, double Dthr){
    R.nc = nc;
    // scan for a min in |dV/dx| at sonic (the critical solution has the smoothest passage),
    // or a sign change; here we locate where dV/dx at approach is minimized (tangency), which
    // for this system is the physically-selected analytic separatrix candidate.
    double best_a2 = a2lo, best_abs = 1e300; bool any=false;
    double Vs_lo=0, Vs_hi=0; bool have_lo=false, have_hi=false;
    for (int i=0;i<nscan;i++){
        double a2 = a2lo + (a2hi-a2lo)*i/(double)(nscan-1);
        Seed s{nc, a2, z0};
        BgResult r = integrate_bg(s, xmax, h, Dthr);
        if (r.code == 0){
            any = true;
            double av = std::fabs(r.Vp_near);
            if (av < best_abs){ best_abs = av; best_a2 = a2;
                R.x_sonic=r.x_sonic; R.N_sonic=r.N_s; R.A_sonic=r.A_s; R.om_sonic=r.om_s;
                R.V_sonic=r.V_s; R.Dresid=Dloc(r.N_s,r.V_s); R.Vmin=r.Vmin; }
            if (i==0){ Vs_lo=r.V_s; have_lo=true; }
            Vs_hi=r.V_s; have_hi=true;
        }
    }
    R.a2_star = best_a2;
    R.bg_reached_sonic = any;
    // convergence probe: recompute sonic V at half the step; spread = |dV|
    if (any){
        Seed s{nc, best_a2, z0};
        BgResult rc = integrate_bg(s, xmax, h*0.5, Dthr);
        if (rc.code==0) R.conv_spread = std::fabs(rc.V_s - R.V_sonic);
    }
    (void)Vs_lo; (void)Vs_hi; (void)have_lo; (void)have_hi;
}

// declared JSON (hash domain: everything below; notes appended outside)
static std::string declared_json(const Report& R, double nc, double z0, double a2lo,
                                 double a2hi, int nscan){
    std::string s; s.reserve(2048);
    s += "{\"tool\":\"fluidcss_nexus\",\"version\":\"0.9.0\",\"units\":\"G=c=1\",\"params\":{";
    s += "\"nc\":" + fmtnum(nc) + ",\"z0\":" + fmtnum(z0);
    s += ",\"a2lo\":" + fmtnum(a2lo) + ",\"a2hi\":" + fmtnum(a2hi);
    s += ",\"nscan\":" + std::to_string(nscan) + "},";
    s += "\"sonic_locus_coeffs\":[3,-1,4,-1,3],";
    s += "\"background\":{";
    s += "\"reached_sonic\":"; s += R.bg_reached_sonic ? "true":"false";
    s += ",\"a2_star\":" + fmtnum(R.a2_star);
    s += ",\"x_sonic\":" + fmtnum(R.x_sonic);
    s += ",\"N_sonic\":" + fmtnum(R.N_sonic);
    s += ",\"A_sonic\":" + fmtnum(R.A_sonic);
    s += ",\"omega_sonic\":" + fmtnum(R.om_sonic);
    s += ",\"V_sonic\":" + fmtnum(R.V_sonic);
    s += ",\"D_residual\":" + fmtnum(R.Dresid);
    s += ",\"Vmin_ingoing\":" + fmtnum(R.Vmin);
    s += ",\"conv_spread\":" + fmtnum(R.conv_spread);
    s += "},\"eigenvalue\":{";
    s += "\"beta\":" + fmtnum(R.beta);
    s += ",\"re_kappa0\":" + fmtnum(R.re_kappa0);
    s += ",\"im_kappa0\":" + fmtnum(R.im_kappa0);
    s += ",\"n_relevant_modes\":" + std::to_string(R.n_relevant_modes);
    s += ",\"gauge_mode_kappa\":" + fmtnum(R.gauge_mode_kappa);
    s += "},\"gates\":{";
    s += "\"G_ANCHOR\":"; s += R.G_ANCHOR?"true":"false";
    s += ",\"G_CONVERGE\":"; s += R.G_CONVERGE?"true":"false";
    s += ",\"G_UNIQUE\":"; s += R.G_UNIQUE?"true":"false";
    s += "},\"verdict\":\"" + R.verdict + "\"}";
    return s;
}

int main(int argc, char** argv){
    bool json=false, dostest=false, golden=false, stageA=false;
    // golden params (declared)
    double nc = 1.5, z0 = std::exp(-10.0);
    double a2lo = 0.02, a2hi = 1.0;
    int    nscan = 40;
    double xmax = 3.0, h = 2e-4, Dthr = 0.05;

    for (int i=1;i<argc;i++){
        std::string a = argv[i];
        if (a=="--json") json=true;
        else if (a=="--selftest") dostest=true;
        else if (a=="--golden") golden=true;
        else if (a=="--stageA") stageA=true;
        else if (a=="--nc" && i+1<argc) nc=strtod(argv[++i],nullptr);
        else if (a=="--nscan" && i+1<argc) nscan=atoi(argv[++i]);
        else if (a=="--h" && i+1<argc) h=strtod(argv[++i],nullptr);
        else {
            fprintf(stderr,"usage: fluidcss_nexus [--json] [--selftest] [--golden] [--stageA] "
                           "[--nc V] [--nscan N] [--h H]\n");
            return 2;
        }
    }
    (void)stageA;

    if (dostest){
        std::string why;
        bool ok = selftest(why);
        printf("selftest: %s (%s)\n", ok?"PASS":"FAIL", why.c_str());
        return ok?0:1;
    }

    // sanity: selftest must pass before we trust anything
    {
        std::string why;
        if (!selftest(why)){ fprintf(stderr,"SELFTEST FAILED: %s\n", why.c_str()); return 2; }
    }

    Report R;
    run_stageA(R, nc, z0, a2lo, a2hi, nscan, xmax, h, Dthr);

    // Gates. Stage B (beta) not computed -> beta is nan -> G_ANCHOR cannot pass.
    // G_CONVERGE here checks the STAGE-A sonic V is stable under step refinement (a real,
    // meaningful metamorphic check on the piece we DID compute).
    R.G_CONVERGE = R.bg_reached_sonic && (R.conv_spread < 5e-3);
    R.G_UNIQUE   = false;                 // eigenvalue box not scanned yet
    R.G_ANCHOR   = std::isfinite(R.beta) && (std::fabs(R.beta - 0.35580192) < 4e-3);
    // verdict: "blocked" while beta not measured; never "pass" without G-ANCHOR.
    R.verdict = R.G_ANCHOR ? "pass" : "blocked";

    std::string declared = declared_json(R, nc, z0, a2lo, a2hi, nscan);
    std::string hash = blake2b::hash256_hex(declared);

    if (golden){
        FILE* f = fopen("goldens/fluidcss/golden.hash","rb");
        if (!f){ fprintf(stderr,"GOLDEN NOT FROZEN %s\n", hash.c_str());
                 printf("%s\n", hash.c_str()); return 2; }
        char want[128]={0};
        if (fscanf(f,"%127s",want)!=1){ fclose(f); fprintf(stderr,"GOLDEN FILE UNREADABLE\n"); return 2; }
        fclose(f);
        if (hash==want){ fprintf(stderr,"GOLDEN OK %.8s\n", hash.c_str()); printf("%s\n", hash.c_str()); return 0; }
        fprintf(stderr,"GOLDEN MISMATCH have=%.8s want=%.8s\n", hash.c_str(), want);
        printf("%s\n", hash.c_str());
        return 1;
    }

    if (json){
        std::string out = declared;
        out.pop_back();  // drop closing brace
        out += ",\"notes\":\"beta NOT measured (Stage B eigenvalue open); Stage-A background+sonic "
               "verified. hash=" + hash.substr(0,8) + "\"}";
        printf("%s\n", out.c_str());
        // exit 1 for the declared 'blocked' verdict (a real negative result, not an error)
        return R.verdict=="pass" ? 0 : 1;
    }

    // human
    printf("fluidcss_nexus v0.9.0 - radiation-fluid CSS critical exponent (eigenvalue route)\n");
    printf("--------------------------------------------------------------------------------\n");
    printf("  units G=c=1 ; target Re kappa0=2.81055255 beta=0.35580192 (KHA95)\n");
    printf("  sonic locus: 3 N^2 V^2 - N^2 + 4 N V - V^2 + 3 = 0  [VERIFIED]\n");
    printf("  STAGE A (background critical CSS solution):\n");
    printf("    reached sonic point : %s\n", R.bg_reached_sonic?"yes":"no");
    printf("    a2* (density param) : %.6f  (nc gauge=%.3f)\n", R.a2_star, R.nc);
    printf("    sonic point         : x=%.4f N=%.5f A=%.5f omega=%.5f V=%.5f\n",
           R.x_sonic, R.N_sonic, R.A_sonic, R.om_sonic, R.V_sonic);
    printf("    D residual @ sonic  : %.3e   Vmin(ingoing): %.4f\n", R.Dresid, R.Vmin);
    printf("    conv spread (h/2)   : %.3e\n", R.conv_spread);
    printf("  STAGE B (perturbation eigenvalue -> beta):\n");
    printf("    beta                : %s\n", std::isfinite(R.beta)?fmtnum(R.beta).c_str():"NOT MEASURED (open)");
    printf("    Re kappa0           : %s\n", std::isfinite(R.re_kappa0)?fmtnum(R.re_kappa0).c_str():"NOT MEASURED (open)");
    printf("  GATES: G-ANCHOR=%s  G-CONVERGE=%s  G-UNIQUE=%s\n",
           R.G_ANCHOR?"PASS":"fail", R.G_CONVERGE?"PASS":"fail", R.G_UNIQUE?"PASS":"fail");
    printf("--------------------------------------------------------------------------------\n");
    printf("VERDICT: %s\n", R.verdict.c_str());
    printf("declared hash: %s\n", hash.c_str());
    return R.verdict=="pass" ? 0 : 1;
}
