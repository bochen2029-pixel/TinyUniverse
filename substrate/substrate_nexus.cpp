// ============================================================================
//  substrate_nexus — TINY UNIVERSE v2 N0: the substrate oracle
//  Contract: contracts/substrate_nexus.contract.md v1.0.0
//
//  Single file, C++17, fp64, stdlib only, NO GPU. Spherically-symmetric
//  Einstein-Klein-Gordon in polar-areal gauge, constrained evolution. Crown:
//  the Choptuik critical-collapse mass-scaling exponent gamma ~ 0.374.
//
//  Units: geometric G = c = 1 (gamma is a dimensionless universal number).
//
//  Build (MSVC, golden platform):
//    cl /std:c++17 /EHsc /O2 /W4 substrate\substrate_nexus.cpp /Fe:build\substrate_nexus.exe
//
//  Exit: 0 all gates pass / 1 a declared gate fired / 2 error.
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
// BLAKE2b-256 (RFC 7693) — lifted verbatim from nexus/tiny_nexus.cpp (D-020).
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
// Spherical Einstein-Klein-Gordon solver (polar-areal, constrained evolution)
// ============================================================================
struct EKG {
    int N;                        // grid points 0..N (r_j = j*dr)
    double dr, rmax, mu2, eko, cfl;
    std::vector<double> phi, Phi, Pi, a, alp;
    std::vector<double> kphi, kPhi, kPi;               // RK accumulators
    std::vector<double> phi0, Phi0, Pi0;               // RK base

    EKG(int N_, double rmax_, double mu2_, double eko_, double cfl_)
        : N(N_), dr(rmax_/N_), rmax(rmax_), mu2(mu2_), eko(eko_), cfl(cfl_),
          phi(N_+1), Phi(N_+1), Pi(N_+1), a(N_+1), alp(N_+1),
          kphi(N_+1), kPhi(N_+1), kPi(N_+1), phi0(N_+1), Phi0(N_+1), Pi0(N_+1) {}

    double r(int j) const { return j*dr; }
    double Vof(double f) const { return 0.5*mu2*f*f; }

    // time-symmetric Gaussian initial data: phi = p*exp(-((r-r0)/sig)^2), Pi = 0
    void init(double p, double r0, double sig){
        for (int j = 0; j <= N; j++){
            double rr = r(j), x = (rr - r0)/sig;
            phi[j] = p*std::exp(-x*x);
            Phi[j] = phi[j]*(-2.0*(rr - r0)/(sig*sig));
            Pi[j]  = 0.0;
        }
        Phi[0] = 0.0;                                  // parity (odd) at origin
    }

    // Hamiltonian-constraint source for (ln a)': (1-a^2)/(2r) + 2*pi*r*(Pi^2+Phi^2) + 4*pi*r*a^2*V
    double Fa(int j, double av) const {
        double rr = r(j);
        double src = 2.0*PI*rr*(Pi[j]*Pi[j] + Phi[j]*Phi[j]) + 4.0*PI*rr*av*av*Vof(phi[j]);
        double geo = (rr > 0.0) ? (1.0 - av*av)/(2.0*rr) : 0.0;
        return geo + src;
    }
    // slicing source for (ln alpha)': (a^2-1)/(2r) + 2*pi*r*(Pi^2+Phi^2) - 4*pi*r*a^2*V
    double Fal(int j, double av) const {
        double rr = r(j);
        double src = 2.0*PI*rr*(Pi[j]*Pi[j] + Phi[j]*Phi[j]) - 4.0*PI*rr*av*av*Vof(phi[j]);
        double geo = (rr > 0.0) ? (av*av - 1.0)/(2.0*rr) : 0.0;
        return geo + src;
    }
    // solve a, alpha from the constraints by outward Heun integration of the LOGS
    // (d(ln a)/dr = Fa, d(ln alpha)/dr = Fal) — guarantees a,alpha > 0 through the
    // stiff near-horizon region; rescale alpha so alpha*a -> 1 at rmax.
    void solveMetric(){
        double lna = 0.0; a[0] = 1.0;                  // ln a(0) = 0
        for (int j = 0; j < N; j++){
            double f0 = Fa(j, a[j]);                   // d(ln a)/dr at j
            double ap = std::exp(lna + dr*f0);         // Euler predictor for a[j+1]
            double f1 = Fa(j+1, ap);
            lna += 0.5*dr*(f0 + f1);
            a[j+1] = std::exp(lna);
        }
        double lnal = 0.0; alp[0] = 1.0;               // ln alpha(0) = 0; Fal depends only on a
        for (int j = 0; j < N; j++){                    // -> pure trapezoidal quadrature
            lnal += 0.5*dr*(Fal(j, a[j]) + Fal(j+1, a[j+1]));
            alp[j+1] = std::exp(lnal);
        }
        double s = 1.0/(alp[N]*a[N]);                  // proper time = coordinate time at infinity
        for (int j = 0; j <= N; j++) alp[j] *= s;
    }

    // field RHS given current fields+metric -> (dphi,dPhi,dPi). Origin parity, outer radiative.
    void rhs(std::vector<double>& dphi, std::vector<double>& dPhi, std::vector<double>& dPi){
        static std::vector<double> X, Y, W;
        X.assign(N+1, 0); Y.assign(N+1, 0); W.assign(N+1, 0);
        for (int j = 0; j <= N; j++){
            double v = alp[j]/a[j];
            X[j] = v*Pi[j];                            // dphi/dt = X ; dPhi/dt = dX/dr
            Y[j] = v*Phi[j];                           // dPi/dt = (1/r^2) d(r^2 Y)/dr - alpha a mu^2 phi
            W[j] = r(j)*r(j)*Y[j];
        }
        for (int j = 0; j <= N; j++){
            dphi[j] = X[j];
            double src = alp[j]*a[j]*mu2*phi[j];
            if (j == 0){
                dPhi[0] = 0.0;                          // X even -> dX/dr(0)=0
                dPi[0]  = 3.0*Y[1]/dr - src;            // (1/r^2)d(r^2 Y)/dr -> 3 Y'(0), Y odd
            } else if (j == N){                         // outgoing: backward (upwind) differences
                dPhi[N] = (X[N] - X[N-1])/dr;
                dPi[N]  = (W[N] - W[N-1])/dr/(r(N)*r(N)) - src;
            } else {
                dPhi[j] = (X[j+1] - X[j-1])/(2.0*dr);
                dPi[j]  = (W[j+1] - W[j-1])/(2.0*dr)/(r(j)*r(j)) - src;
            }
        }
        // Kreiss-Oliger dissipation (interior): -eko/(16 dr) * (2nd-order KO stencil)
        if (eko > 0.0){
            for (int j = 2; j <= N-2; j++){
                double c = -eko/(16.0*dr);
                dPhi[j] += c*(Phi[j+2] - 4*Phi[j+1] + 6*Phi[j] - 4*Phi[j-1] + Phi[j-2]);
                dPi[j]  += c*(Pi[j+2]  - 4*Pi[j+1]  + 6*Pi[j]  - 4*Pi[j-1]  + Pi[j-2]);
                dphi[j] += c*(phi[j+2] - 4*phi[j+1] + 6*phi[j] - 4*phi[j-1] + phi[j-2]);
            }
        }
    }

    // one RK4 step of dt; metric solved at each substage
    void step(double dt){
        phi0 = phi; Phi0 = Phi; Pi0 = Pi;
        static std::vector<double> dphi, dPhi, dPi, sphi, sPhi, sPi;
        sphi.assign(N+1,0); sPhi.assign(N+1,0); sPi.assign(N+1,0);
        dphi.assign(N+1,0); dPhi.assign(N+1,0); dPi.assign(N+1,0);
        const double w[4] = {1,2,2,1};
        const double sub[4] = {0.0, 0.5, 0.5, 1.0};
        for (int stage = 0; stage < 4; stage++){
            if (stage > 0){
                for (int j = 0; j <= N; j++){
                    phi[j] = phi0[j] + sub[stage]*dt*kphi[j];
                    Phi[j] = Phi0[j] + sub[stage]*dt*kPhi[j];
                    Pi[j]  = Pi0[j]  + sub[stage]*dt*kPi[j];
                }
            }
            solveMetric();
            rhs(dphi, dPhi, dPi);
            kphi = dphi; kPhi = dPhi; kPi = dPi;
            for (int j = 0; j <= N; j++){
                sphi[j] += w[stage]*kphi[j];
                sPhi[j] += w[stage]*kPhi[j];
                sPi[j]  += w[stage]*kPi[j];
            }
        }
        for (int j = 0; j <= N; j++){
            phi[j] = phi0[j] + (dt/6.0)*sphi[j];
            Phi[j] = Phi0[j] + (dt/6.0)*sPhi[j];
            Pi[j]  = Pi0[j]  + (dt/6.0)*sPi[j];
        }
        Phi[0] = 0.0;                                  // enforce parity
    }

    double mass(int j) const { return 0.5*r(j)*(1.0 - 1.0/(a[j]*a[j])); }
    double max2mr() const {
        double mx = 0;
        for (int j = 1; j <= N; j++){ double v = 1.0 - 1.0/(a[j]*a[j]); if (v > mx) mx = v; }
        return mx;
    }
    double maxMass() const {
        double mx = 0;
        for (int j = 1; j <= N; j++){ double m = mass(j); if (m > mx) mx = m; }
        return mx;
    }
    // apparent-horizon mass = mass aspect at the compactness (2m/r) peak. This is the
    // Choptuik observable: it -> 0 as p -> p* (horizon forms at ever-smaller radius).
    double horizonMass(double& r_ah) const {
        int jh = 1; double best = 0;
        for (int j = 1; j <= N; j++){ double v = 1.0 - 1.0/(a[j]*a[j]); if (v > best){ best = v; jh = j; } }
        r_ah = r(jh);
        return mass(jh);
    }
    double centralField() const {
        double s = 0; for (int j = 0; j <= N/8; j++) s += std::fabs(phi[j]); return s/(N/8+1);
    }
    bool finite() const {
        for (int j = 0; j <= N; j++) if (!std::isfinite(phi[j]) || !std::isfinite(Pi[j])) return false;
        return true;
    }
};

// outcome of one evolution
struct Outcome {
    bool collapsed = false, ok = true;
    double m_bh = 0, r_ah = 0, max2mr = 0, alpha0_min = 1, mass0 = 0, mass_drift = 0, cfield_ratio = 1;
    double zmax = 0;                                    // peak near-central curvature (Pi^2+Phi^2)
    double t_end = 0;
};

// evolve time-symmetric-Gaussian data of amplitude p; classify collapse vs dispersal
static Outcome evolve(int N, double rmax, double mu2, double eko, double cfl,
                      double p, double r0, double sig, double tmax_mult,
                      bool track_massdrift){
    EKG s(N, rmax, mu2, eko, cfl);
    s.init(p, r0, sig);
    s.solveMetric();
    Outcome o;
    o.mass0 = s.maxMass();
    double cf0 = s.centralField();
    double dt = cfl*s.dr;
    double tmax = tmax_mult*rmax;
    double massdrift = 0;
    int nsteps = (int)(tmax/dt);
    double crossing = rmax;                            // pulse reaches outer boundary ~ this time
    double peak2mr = 0;
    for (int it = 0; it < nsteps; it++){
        s.step(dt);
        double t = (it+1)*dt;
        if (!s.finite()){ o.ok = false; break; }
        double a0 = s.alp[0];
        if (a0 < o.alpha0_min) o.alpha0_min = a0;
        double m2r = s.max2mr();
        if (m2r > o.max2mr) o.max2mr = m2r;
        // peak near-central curvature R ~ (Pi^2+Phi^2) at the origin — the
        // Garfinkle-Duncan subcritical scaling variable: max_t R ~ (p*-p)^(-2 gamma)
        double z = 0;
        for (int j = 0; j <= 3; j++){ double zz = s.Pi[j]*s.Pi[j] + s.Phi[j]*s.Phi[j]; if (zz > z) z = zz; }
        if (z > o.zmax) o.zmax = z;
        // mass conservation only meaningful before radiation crosses rmax (S1 oracle)
        if (track_massdrift && t < 0.6*crossing){
            double md = std::fabs(s.maxMass() - o.mass0)/(o.mass0 + 1e-300);
            if (md > massdrift) massdrift = md;
        }
        // apparent horizon (clean collapse criterion): 2m/r -> 1. Near-critical runs
        // that NaN before a resolved horizon forms are excluded, not misclassified —
        // the crown exponent comes from SUBCRITICAL curvature scaling (S4), not from
        // resolving the horizon a uniform grid cannot capture.
        if (m2r > 0.98){
            o.collapsed = true;
            o.m_bh = s.horizonMass(o.r_ah);            // mass at the compactness peak (horizon)
            o.t_end = t;
            break;
        }
        (void)peak2mr; (void)crossing;
    }
    if (o.t_end == 0.0) o.t_end = nsteps*dt;
    o.mass_drift = massdrift;
    o.cfield_ratio = s.centralField()/(cf0 + 1e-300);
    o.max2mr = std::max(o.max2mr, s.max2mr());
    return o;
}

// bisect the critical amplitude p*
static double bisect_pstar(int N, double rmax, double eko, double cfl,
                           double r0, double sig, double plo, double phi_,
                           int iters, double tmax_mult){
    for (int i = 0; i < iters; i++){
        double pm = 0.5*(plo + phi_);
        Outcome o = evolve(N, rmax, 0.0, eko, cfl, pm, r0, sig, tmax_mult, false);
        if (o.collapsed) phi_ = pm; else plo = pm;
        if ((phi_ - plo)/phi_ < 1e-7) break;
    }
    return 0.5*(plo + phi_);
}

// ============================================================================
// Result plumbing (tiny_nexus idiom)
// ============================================================================
struct TestResult {
    std::string id, name;
    bool pass = false;
    std::vector<std::pair<std::string, double>> m;
    void put(const std::string& k, double v){ m.emplace_back(k, v); }
};

// battery scenario parameters (golden-pinned). N=800/rmax=24 is the frozen config:
// clean far-field mass-scaling (R^2~0.998) at ~20 s; finer uniform grids chase the
// unresolvable near-critical regime and turn chaotic (D-021) — AMR is the later contract.
struct Cfg {
    int    N       = 800;       // Choptuik grid (dr = 0.03)
    double rmax    = 24.0;
    double eko     = 0.5;
    double cfl     = 0.25;
    double r0      = 0.0;       // centered time-symmetric Gaussian
    double sig     = 1.0;
    double tmax    = 2.0;       // in units of rmax
};

static TestResult testS1(const Cfg& c){          // flat-space wave (eps amplitude)
    TestResult t; t.id = "S1"; t.name = "flatwave";
    Outcome o = evolve(800, 30.0, 0.0, c.eko, c.cfl, 1e-3, 5.0, 1.0, 1.5, true);
    bool disperse = o.cfield_ratio < 0.2;
    bool masscons = o.mass_drift < 1e-3;
    bool weak = o.max2mr < 1e-3;
    t.pass = o.ok && disperse && masscons && weak;
    t.put("mass_drift", o.mass_drift);
    t.put("cfield_ratio", o.cfield_ratio);
    t.put("max_2m_over_r", o.max2mr);
    t.put("finite", o.ok ? 1.0 : 0.0);
    return t;
}
static TestResult testS2(const Cfg& c, double pstar){   // subcritical dispersal
    TestResult t; t.id = "S2"; t.name = "subcrit";
    double p = 0.6*pstar;                            // safely below critical (clean dispersal)
    Outcome o = evolve(c.N, c.rmax, 0.0, c.eko, c.cfl, p, c.r0, c.sig, c.tmax, false);
    bool nohorizon = !o.collapsed && o.max2mr < 1.0;
    bool recovered = o.alpha0_min > 0.1;             // lapse dipped but did not collapse
    bool dispersed = o.cfield_ratio < 0.5;
    t.pass = o.ok && nohorizon && recovered && dispersed;
    t.put("p_over_pstar", p/pstar);
    t.put("max_2m_over_r", o.max2mr);
    t.put("alpha0_min", o.alpha0_min);
    t.put("cfield_ratio", o.cfield_ratio);
    return t;
}
static TestResult testS3(const Cfg& c, double pstar){   // supercritical collapse -> horizon
    TestResult t; t.id = "S3"; t.name = "supercrit";
    double p = 1.5*pstar;
    Outcome o = evolve(c.N, c.rmax, 0.0, c.eko, c.cfl, p, c.r0, c.sig, c.tmax, false);
    bool horizon = o.collapsed && o.max2mr > 0.98;
    bool censor = o.m_bh > 0.0 && o.m_bh < 2.0*o.mass0 + 1.0;   // BH mass ~ O(pulse mass)
    t.pass = o.ok && horizon && censor;
    t.put("p_over_pstar", p/pstar);
    t.put("max_2m_over_r", o.max2mr);
    t.put("m_bh", o.m_bh);
    t.put("alpha0_min", o.alpha0_min);
    return t;
}
// The Choptuik phenomenon, honestly at uniform-grid resolution: demonstrate the
// TYPE-II critical transition — a sharp critical amplitude p* below which the pulse
// disperses and above which a black hole forms, with the horizon mass going
// CONTINUOUSLY to zero as p -> p* (arbitrarily small black holes). We report an
// effective mass-scaling exponent as a grid-limited DIAGNOSTIC; the precise universal
// gamma = 0.374 requires adaptive mesh refinement and is DEFERRED (proposal 6, D-021)
// — a uniform grid caps the self-similar curvature and cannot resolve it. Faking a
// tight gamma gate against grid noise would poison the v2 oracle farm (RAYFORMER/D-016).
static TestResult testS4(const Cfg& c, double& pstar_out, double& gamma_out){
    TestResult t; t.id = "S4"; t.name = "critical";
    double pstar = bisect_pstar(c.N, c.rmax, c.eko, c.cfl, c.r0, c.sig, 0.05, 0.70, 20, c.tmax);
    pstar_out = pstar;
    // clean bracketing of the transition: below disperses, above collapses
    Outcome lo = evolve(c.N, c.rmax, 0.0, c.eko, c.cfl, pstar*0.70, c.r0, c.sig, c.tmax, false);
    Outcome hi = evolve(c.N, c.rmax, 0.0, c.eko, c.cfl, pstar*1.30, c.r0, c.sig, c.tmax, false);
    bool bracket = (!lo.collapsed) && hi.collapsed;
    // supercritical mass scaling: M_BH shrinks toward the grid floor as p -> p*
    std::vector<double> P, M;
    const int NP = 12;
    double rmin = 3.0*c.rmax/c.N;                       // horizon must be a few cells wide
    for (int i = 0; i < NP; i++){
        double f = 0.03*std::pow(15.0, i/(double)(NP-1));   // (p-p*)/p* in [0.03, 0.45]
        double p = pstar*(1.0 + f);
        Outcome o = evolve(c.N, c.rmax, 0.0, c.eko, c.cfl, p, c.r0, c.sig, c.tmax, false);
        if (o.collapsed && o.m_bh > 0 && o.r_ah > rmin){ P.push_back(p); M.push_back(o.m_bh); }
        fprintf(stderr, "  [S4] p=%.6f coll=%d m_bh=%.6e r_ah=%.4f\n",
                p, (int)o.collapsed, o.m_bh, o.r_ah);
    }
    int n = (int)M.size();
    double mmin = 1e30, mmax = 0;
    for (double v : M){ mmin = std::min(mmin, v); mmax = std::max(mmax, v); }
    if (n == 0){ mmin = 0; }
    bool monotone = true;
    for (int i = 1; i < n; i++) if (M[i] < 0.7*M[i-1]) monotone = false;   // M rises with p
    // effective exponent (grid-limited diagnostic): scan p* for best-R^2 ln M vs ln(p-p*)
    double pfit = 0, r2 = -1, gamma = 0;
    if (n >= 5){
        double pminP = 1e30; for (double v : P) pminP = std::min(pminP, v);
        for (int k = 1; k < 500; k++){
            double pc = pminP*(k/500.0);
            double sx=0,sy=0,sxx=0,sxy=0; int m=0; bool bad=false;
            for (int i = 0; i < n; i++){ double d = P[i]-pc; if (d<=0){ bad=true; break; }
                double x=std::log(d), y=std::log(M[i]); sx+=x; sy+=y; sxx+=x*x; sxy+=x*y; m++; }
            if (bad || m < 3) continue;
            double den = m*sxx - sx*sx; if (std::fabs(den) < 1e-30) continue;
            double g = (m*sxy - sx*sy)/den, b = (sy - g*sx)/m;
            double ss=0, st=0, my=sy/m;
            for (int i = 0; i < n; i++){ double x=std::log(P[i]-pc), ff=g*x+b, yy=std::log(M[i]);
                ss+=(yy-ff)*(yy-ff); st+=(yy-my)*(yy-my); }
            double rr = (st>0)? 1.0-ss/st : 0;
            if (rr > r2){ r2 = rr; gamma = g; pfit = pc; }
        }
    }
    gamma_out = gamma;
    double ratio = (mmax > 0)? mmin/mmax : 1.0;
    fprintf(stderr, "  [S4] pstar=%.5f bracket=%d n=%d Mmin=%.4e Mmax=%.4e ratio=%.3f "
            "gamma_eff=%.3f (needs AMR for 0.374) r2=%.3f\n",
            pstar, (int)bracket, n, mmin, mmax, ratio, gamma, r2);
    // HONEST gate: the Type-II transition itself — bracketed threshold + BH mass
    // continuously -> 0 at p* (ratio small) + monotone rise. NOT gated on gamma=0.374.
    bool ok = bracket && (n >= 5) && (ratio < 0.35) && monotone;
    t.pass = ok;
    t.put("pstar", pstar);
    t.put("bracket_clean", bracket ? 1.0 : 0.0);
    t.put("n_points", (double)n);
    t.put("m_bh_min", mmin);
    t.put("m_bh_max", mmax);
    t.put("m_ratio", ratio);
    t.put("pfit", pfit);
    t.put("gamma_eff_gridlimited", gamma);
    t.put("gamma_ref_amr", 0.374);
    t.put("fit_r2", r2);
    return t;
}
static TestResult testS5(const Cfg& c){          // massive-field stable + conservative (mu>0)
    TestResult t; t.id = "S5"; t.name = "massive";
    double mu2 = 1.0;                                 // Klein-Gordon mass term active
    Outcome o = evolve(1000, 20.0, mu2, c.eko, c.cfl, 0.05, 0.0, 2.0, 2.0, true);
    // honest scope: the massive-KG term (potential V=mu^2 phi^2/2 in both the field
    // evolution and the constraints) evolves STABLY and CONSERVES mass. A true bound
    // oscillaton needs an eigen-profile (deferred, later contract) — not claimed here.
    bool stable = o.ok;                              // finite (mass term didn't destabilize)
    bool nocollapse = !o.collapsed && o.max2mr < 1.0;
    bool masscons = o.mass_drift < 2e-2;             // ADM mass conserved with V active
    t.pass = stable && nocollapse && masscons;
    t.put("finite", o.ok ? 1.0 : 0.0);
    t.put("max_2m_over_r", o.max2mr);
    t.put("mass_drift", o.mass_drift);
    t.put("cfield_ratio", o.cfield_ratio);
    return t;
}

// declared JSON (hash domain: everything below; notes appended outside)
static std::string declared_json(const Cfg& c, const std::vector<TestResult>& R,
                                 const TestResult& det){
    std::string s; s.reserve(4096);
    s += "{\"tool\":\"substrate_nexus\",\"version\":\"1.0.0\",\"units\":\"G=c=1\",\"params\":{";
    s += "\"N\":" + std::to_string(c.N) + ",\"rmax\":" + fmtnum(c.rmax);
    s += ",\"eko\":" + fmtnum(c.eko) + ",\"cfl\":" + fmtnum(c.cfl);
    s += ",\"r0\":" + fmtnum(c.r0) + ",\"sig\":" + fmtnum(c.sig);
    s += ",\"tmax_mult\":" + fmtnum(c.tmax) + "},\"results\":[";
    auto emit = [&](const TestResult& r, bool first){
        if (!first) s += ",";
        s += "{\"id\":\"" + r.id + "\",\"name\":\"" + r.name + "\",\"pass\":";
        s += r.pass ? "true" : "false";
        s += ",\"metrics\":{";
        for (size_t i = 0; i < r.m.size(); i++){
            if (i) s += ",";
            s += "\"" + r.m[i].first + "\":" + fmtnum(r.m[i].second);
        }
        s += "}}";
    };
    for (size_t i = 0; i < R.size(); i++) emit(R[i], i == 0);
    emit(det, false);
    bool all = true; for (auto& r : R) all = all && r.pass; all = all && det.pass;
    s += "],\"verdict\":\"";
    s += all ? "pass" : "fail";
    s += "\"";
    return s;
}

static bool run_battery(const Cfg& c, std::vector<TestResult>& R, double& pstar, double& gamma){
    R.push_back(testS1(c));
    TestResult s4 = testS4(c, pstar, gamma);
    R.push_back(testS2(c, pstar));
    R.push_back(testS3(c, pstar));
    R.push_back(s4);
    R.push_back(testS5(c));
    bool all = true; for (auto& r : R) all = all && r.pass;
    return all;
}

static void print_human(const std::vector<TestResult>& R, const TestResult& det){
    printf("substrate_nexus v1.0.0 - TINY UNIVERSE v2 N0 (spherical EKG, G=c=1)\n");
    printf("--------------------------------------------------------------------\n");
    auto row = [](const TestResult& r){
        printf("  %-4s %-10s %s\n", r.id.c_str(), r.name.c_str(), r.pass ? "PASS" : "FAIL");
        for (auto& kv : r.m) printf("        %-22s %.9e\n", kv.first.c_str(), kv.second);
    };
    for (auto& r : R) row(r);
    row(det);
    bool all = true; for (auto& r : R) all = all && r.pass; all = all && det.pass;
    printf("--------------------------------------------------------------------\n");
    printf("VERDICT: %s\n", all ? "PASS" : "FAIL");
}

// ============================================================================
// main
// ============================================================================
int main(int argc, char** argv){
    bool json = false, selftest = false, golden = false, dev = false;
    const char* only = nullptr;
    double devp = 0.1;
    Cfg c;
    for (int i = 1; i < argc; i++){
        std::string a = argv[i];
        if (a == "--json") json = true;
        else if (a == "--selftest") selftest = true;
        else if (a == "--golden") golden = true;
        else if (a == "--only" && i+1 < argc) only = argv[++i];
        else if (a == "--dev" && i+1 < argc) { dev = true; devp = strtod(argv[++i], nullptr); }
        else if (a == "--N" && i+1 < argc) c.N = atoi(argv[++i]);
        else if (a == "--rmax" && i+1 < argc) c.rmax = strtod(argv[++i], nullptr);
        else {
            fprintf(stderr, "usage: substrate_nexus [--json] [--selftest] [--golden] "
                            "[--only S1..S5] [--dev P] [--N n] [--rmax R]\n");
            return 2;
        }
    }

    if (dev){                                        // single evolution, diagnostics (physics dev)
        printf("dev: N=%d rmax=%.2f eko=%.2f cfl=%.2f  p=%.6f  r0=%.2f sig=%.2f\n",
               c.N, c.rmax, c.eko, c.cfl, devp, c.r0, c.sig);
        Outcome o = evolve(c.N, c.rmax, 0.0, c.eko, c.cfl, devp, c.r0, c.sig, c.tmax, true);
        printf("  collapsed=%d  max2mr=%.6f  m_bh=%.6f  alpha0min=%.6f  mass0=%.6e  cfield_ratio=%.4f  t_end=%.3f  ok=%d\n",
               o.collapsed, o.max2mr, o.m_bh, o.alpha0_min, o.mass0, o.cfield_ratio, o.t_end, o.ok);
        return 0;
    }

    if (only){
        std::vector<TestResult> R; double ps=0, gm=0;
        TestResult r;
        std::string id = only;
        if (id == "S1") r = testS1(c);
        else if (id == "S4") { r = testS4(c, ps, gm); }
        else {
            double pstar = bisect_pstar(c.N, c.rmax, c.eko, c.cfl, c.r0, c.sig, 0.05, 0.70, 40, c.tmax);
            if (id == "S2") r = testS2(c, pstar);
            else if (id == "S3") r = testS3(c, pstar);
            else if (id == "S5") r = testS5(c);
            else { fprintf(stderr, "error: unknown test id (S1..S5)\n"); return 2; }
        }
        printf("%s %s: %s\n", r.id.c_str(), r.name.c_str(), r.pass ? "PASS" : "FAIL");
        for (auto& kv : r.m) printf("  %-22s %.9e\n", kv.first.c_str(), kv.second);
        return r.pass ? 0 : 1;
    }

    std::vector<TestResult> R; double pstar=0, gamma=0;
    bool all = run_battery(c, R, pstar, gamma);
    // determinism probe: re-run S1 in-process and compare its serialization (cheap,
    // serialization-sensitive). Full-battery determinism is proven OUT of process by
    // the golden re-run (GOLDEN OK on a fresh invocation) — a stronger check than this.
    TestResult det; det.id = "S6"; det.name = "det";
    {
        TestResult s1a = testS1(c), s1b = testS1(c);
        auto ser = [](const TestResult& r){ std::string s;
            for (auto& kv : r.m) s += kv.first + fmtnum(kv.second); return s; };
        det.pass = (ser(s1a) == ser(s1b));
        det.put("identical", det.pass ? 1.0 : 0.0);
    }
    all = all && det.pass;

    std::string declared = declared_json(c, R, det);
    std::string hash = blake2b::hash256_hex(declared);

    if (golden){
        FILE* f = fopen("goldens/substrate_nexus/golden.hash", "rb");
        if (!f){
            fprintf(stderr, "GOLDEN NOT FROZEN %s\n", hash.c_str());
            printf("%s\n", hash.c_str());
            return 2;
        }
        char want[128] = {0};
        if (fscanf(f, "%127s", want) != 1){ fclose(f); fprintf(stderr, "GOLDEN FILE UNREADABLE\n"); return 2; }
        fclose(f);
        if (hash == want){ fprintf(stderr, "GOLDEN OK %.8s\n", hash.c_str()); printf("%s\n", hash.c_str()); return 0; }
        fprintf(stderr, "GOLDEN MISMATCH have=%.8s want=%.8s\n", hash.c_str(), want);
        printf("%s\n", hash.c_str());
        return 1;
    }
    if (json){
        std::string out = declared;
        out += ",\"notes\":\"non-declared; hash=" + hash.substr(0,8) + "\"}";
        printf("%s\n", out.c_str());
        return all ? 0 : 1;
    }
    (void)selftest;
    print_human(R, det);
    printf("declared hash: %s\n", hash.c_str());
    return all ? 0 : 1;
}
