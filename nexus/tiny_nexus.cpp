// ============================================================================
//  tiny_nexus — TINY UNIVERSE M0 composition proof & standing oracle
//  Contract: contracts/nexus.contract.md v1.0.0 (+ D-011 errata)
//
//  Single file, C++17, fp64, stdlib only. Proves the four regimes (quantum,
//  classical, relativistic, compact) compose without contradiction under one
//  dial set, and freezes the dial defaults. Battery N1..N11.
//
//  Build (MSVC, golden platform):
//    cl /std:c++17 /EHsc /O2 /W4 nexus\tiny_nexus.cpp /Fe:build\tiny_nexus.exe
//  Parity: g++/clang++ -std=c++17 -O2 -Wall -o tiny_nexus tiny_nexus.cpp
//
//  Exit: 0 all gates pass · 1 a declared gate fired · 2 error.
//  Determinism: (dials, seed) -> byte-identical declared JSON (notes excluded).
// ============================================================================

#define _CRT_SECURE_NO_WARNINGS
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <random>
#include <complex>
#include <algorithm>

static constexpr double PI = 3.14159265358979323846;

// ----------------------------------------------------------------------------
// Dials (contract v0 defaults)
// ----------------------------------------------------------------------------
struct Dials {
    double m_p   = 1.0;
    double c     = 20.0;
    double hbar  = 0.5;
    double G     = 2.0e-3;
    double dt    = 1.0 / 240.0;
    double L_box = 512.0;
    double k_B   = 1.0;
};

// ----------------------------------------------------------------------------
// splitmix64 (purpose-keying, liborrery idiom hand-inlined per contract)
// ----------------------------------------------------------------------------
static uint64_t splitmix64(uint64_t x){
    x += 0x9E3779B97F4A7C15ull;
    uint64_t z = x;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
}

// ----------------------------------------------------------------------------
// BLAKE2b-256 (RFC 7693, self-contained). Envelope byte-compat with liborrery
// is deferred to the M1 lift (MODULE.md); the golden freezes THIS hash.
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
    v[12] ^= t;           // t_lo (inputs < 2^32 bytes here)
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
    h[0] ^= 0x01010000ull ^ 32ull;      // digest_length=32, fanout=1, depth=1
    size_t n = in.size(), off = 0;
    uint8_t block[128];
    // process all but the final block
    while (n - off > 128){
        memcpy(block, in.data() + off, 128);
        off += 128;
        compress(h, block, (uint64_t)off, false);
    }
    size_t rem = n - off;               // 0..128 (0 only when input empty)
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

// ----------------------------------------------------------------------------
// Canonical JSON number formatting: %.9e, -0 normalized. Deterministic per
// platform; the golden is MSVC-pinned (contract determinism clause).
// ----------------------------------------------------------------------------
static std::string fmtnum(double x){
    if (!std::isfinite(x)) return "9.999999999e+99";   // sentinel; gates catch it
    if (x == 0.0) x = 0.0;                              // normalize -0
    char buf[40];
    snprintf(buf, sizeof(buf), "%.9e", x);
    return std::string(buf);
}

// ----------------------------------------------------------------------------
// Result plumbing
// ----------------------------------------------------------------------------
struct TestResult {
    std::string id, name;
    bool pass = false;
    std::vector<std::pair<std::string, double>> m;   // ordered metrics
    void put(const std::string& k, double v){ m.emplace_back(k, v); }
};
struct Battery {
    std::vector<TestResult> results;
    std::vector<std::pair<std::string, double>> derived;
    bool all_pass() const {
        for (auto& r : results) if (!r.pass) return false;
        return true;
    }
};

// ----------------------------------------------------------------------------
// FFT — iterative radix-2, fp64 (N power of two)
// ----------------------------------------------------------------------------
typedef std::complex<double> cplx;
static void fft(std::vector<cplx>& a, bool inverse){
    const size_t n = a.size();
    for (size_t i = 1, j = 0; i < n; i++){
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }
    for (size_t len = 2; len <= n; len <<= 1){
        double ang = 2.0 * PI / (double)len * (inverse ? 1.0 : -1.0);
        cplx wl(std::cos(ang), std::sin(ang));
        for (size_t i = 0; i < n; i += len){
            cplx w(1.0, 0.0);
            for (size_t k = 0; k < len / 2; k++){
                cplx u = a[i + k], v = a[i + k + len/2] * w;
                a[i + k] = u + v;
                a[i + k + len/2] = u - v;
                w *= wl;
            }
        }
    }
    if (inverse){
        double inv = 1.0 / (double)n;
        for (auto& x : a) x *= inv;
    }
}

// 1D split-step workspace ------------------------------------------------------
struct Psi1D {
    int N; double L, dx, m, hbar;
    std::vector<cplx> psi;
    std::vector<double> k;
    Psi1D(int N_, double L_, double m_, double hbar_) : N(N_), L(L_), m(m_), hbar(hbar_){
        dx = L / N;
        psi.assign(N, cplx(0, 0));
        k.assign(N, 0.0);
        for (int j = 0; j < N; j++){
            int jj = (j <= N/2) ? j : j - N;
            k[j] = 2.0 * PI * jj / L;
        }
    }
    double x(int i) const { return -L/2 + i * dx; }
    void normalize(){
        double s = 0;
        for (auto& p : psi) s += std::norm(p);
        s = std::sqrt(s * dx);
        for (auto& p : psi) p /= s;
    }
    void gaussian(double x0, double sigma0, double k0){
        for (int i = 0; i < N; i++){
            double d = x(i) - x0;
            psi[i] = std::exp(cplx(-d*d/(4.0*sigma0*sigma0), k0 * x(i)));
        }
        normalize();
    }
    void step_free_real(double dt, int nsteps){        // V=0: exact in k-space
        fft(psi, false);
        for (int s = 0; s < nsteps; s++)
            for (int j = 0; j < N; j++){
                double ph = -hbar * k[j]*k[j] / (2.0*m) * dt;
                psi[j] *= cplx(std::cos(ph), std::sin(ph));
            }
        fft(psi, true);
    }
    void step_sho_imag(double tau, int iters, double omega){   // Strang, e^{-H tau/hbar}
        std::vector<double> ek(N), ev(N);
        for (int j = 0; j < N; j++) ek[j] = std::exp(-hbar*k[j]*k[j]/(2.0*m) * tau * 0.5);
        for (int i = 0; i < N; i++){
            double V = 0.5 * m * omega*omega * x(i)*x(i);
            ev[i] = std::exp(-V * tau / hbar);
        }
        for (int it = 0; it < iters; it++){
            fft(psi, false);
            for (int j = 0; j < N; j++) psi[j] *= ek[j];
            fft(psi, true);
            for (int i = 0; i < N; i++) psi[i] *= ev[i];
            fft(psi, false);
            for (int j = 0; j < N; j++) psi[j] *= ek[j];
            fft(psi, true);
            normalize();
        }
    }
    double sigma() const {
        double s = 0, sx = 0, sxx = 0;
        for (int i = 0; i < N; i++){
            double p = std::norm(psi[i]);
            s += p; sx += p * x(i); sxx += p * x(i)*x(i);
        }
        sx /= s; sxx /= s;
        return std::sqrt(sxx - sx*sx);
    }
    double energy(double omega) const {                 // <T> (k-space) + <V> (x-space)
        std::vector<cplx> ph = psi;
        fft(ph, false);
        double T = 0, norm2 = 0;
        for (int j = 0; j < N; j++){
            T += hbar*hbar*k[j]*k[j]/(2.0*m) * std::norm(ph[j]);
            norm2 += std::norm(ph[j]);
        }
        T /= norm2;
        double V = 0, s = 0;
        for (int i = 0; i < N; i++){
            double p = std::norm(psi[i]);
            V += p * 0.5 * m * omega*omega * x(i)*x(i);
            s += p;
        }
        return T + V / s;
    }
};

// ----------------------------------------------------------------------------
// Small vector helpers (2D orbital tests)
// ----------------------------------------------------------------------------
struct V2 { double x = 0, y = 0; };
static inline V2 operator+(V2 a, V2 b){ return {a.x + b.x, a.y + b.y}; }
static inline V2 operator-(V2 a, V2 b){ return {a.x - b.x, a.y - b.y}; }
static inline V2 operator*(double s, V2 a){ return {s * a.x, s * a.y}; }
static inline double dot(V2 a, V2 b){ return a.x*b.x + a.y*b.y; }
static inline double len(V2 a){ return std::sqrt(dot(a, a)); }

// LRL eccentricity vector (relative 2D orbit, mu = G*Mtot)
static V2 lrl(V2 r, V2 v, double mu){
    double L = r.x*v.y - r.y*v.x;
    double rn = len(r);
    return { v.y*L/mu - r.x/rn, -v.x*L/mu - r.y/rn };
}

// ============================================================================
// N1 — dials audit
// ============================================================================
static TestResult testN1(const Dials& D, std::vector<std::pair<std::string,double>>& derived){
    TestResult t; t.id = "N1"; t.name = "dials";
    double c2 = D.c * D.c;
    double lambda = D.hbar / (D.m_p * 1.0);                       // v = 1 su/s
    double rs1e5  = 2.0 * D.G * (1e5 * D.m_p) / c2;
    double rs1e2  = 2.0 * D.G * (1e2 * D.m_p) / c2;
    double vcirc  = std::sqrt(D.G * (1e4 * D.m_p) / 10.0) / D.c;
    auto TH = [&](double M){ return D.hbar * c2 * D.c / (8.0 * PI * D.G * M * D.k_B); };
    double th1e5 = TH(1e5 * D.m_p), th1e3 = TH(1e3 * D.m_p);
    double th_ratio = th1e3 / th1e5;
    double n_classical = 30.0, n_star = 1e4;                      // declared nominals
    double n_bh = c2 * 1.0 / (2.0 * D.G * D.m_p);                 // N where r_s = 1 su
    double span = std::log10(n_bh / n_classical);

    bool ok = true;
    ok &= (lambda >= 0.1 && lambda <= 10.0);
    ok &= (rs1e5 >= 0.5 && rs1e5 <= 5.0);
    ok &= (rs1e2 < 0.01);
    ok &= (vcirc >= 0.03 && vcirc <= 0.3);
    ok &= (std::fabs(th_ratio - 100.0) <= 1e-9);
    ok &= (span <= 4.0);

    t.put("lambda_dB_v1", lambda);
    t.put("r_s_1e5", rs1e5);
    t.put("r_s_1e2", rs1e2);
    t.put("v_circ_1e4_r10_over_c", vcirc);
    t.put("T_H_1e5", th1e5);
    t.put("T_H_1e3", th1e3);
    t.put("T_H_ratio", th_ratio);
    t.put("n_classical_nominal", n_classical);
    t.put("n_star_nominal", n_star);
    t.put("n_bh_rs1su", n_bh);
    t.put("span_decades", span);
    t.pass = ok;

    derived.emplace_back("lambda_dB_v1", lambda);
    derived.emplace_back("r_s_1e5", rs1e5);
    derived.emplace_back("r_s_1e2", rs1e2);
    derived.emplace_back("v_circ_1e4_r10_over_c", vcirc);
    derived.emplace_back("T_H_1e5", th1e5);
    derived.emplace_back("T_H_1e3", th1e3);
    derived.emplace_back("span_decades", span);
    return t;
}

// ============================================================================
// N2 — kepler: two-body KDK leapfrog vs analytic, e=0.6, 100 orbits
// Gates: |dE/E| < 1e-9, |dL/L| < 1e-11 (D-011 erratum), peri drift < 1e-4 rad/orbit
// ============================================================================
static TestResult testN2(const Dials& D){
    TestResult t; t.id = "N2"; t.name = "kepler";
    const double m1 = 1e4 * D.m_p, m2 = 1.0 * D.m_p, Mtot = m1 + m2;
    const double mu = D.G * Mtot, a = 150.0, e = 0.6, h = D.dt;
    const double rap = a * (1.0 + e);
    const double vap = std::sqrt(mu * (2.0/rap - 1.0/a));
    // CM frame: relative state -> body states
    V2 r1{ -(m2/Mtot)*rap, 0 }, r2{ (m1/Mtot)*rap, 0 };
    V2 v1{ 0, -(m2/Mtot)*vap }, v2{ 0, (m1/Mtot)*vap };
    const double T = 2.0*PI*std::sqrt(a*a*a/mu);
    const long long nsteps = llround(100.0 * T / h);

    auto accel = [&](V2 R1, V2 R2, V2& a1, V2& a2){
        V2 d = R2 - R1;
        double r = len(d), r3 = r*r*r;
        a1 = ( D.G*m2/r3) * d;
        a2 = (-D.G*m1/r3) * d;
    };
    auto energy = [&](){
        double r = len(r2 - r1);
        return 0.5*m1*dot(v1,v1) + 0.5*m2*dot(v2,v2) - D.G*m1*m2/r;
    };
    auto angmom = [&](){
        return m1*(r1.x*v1.y - r1.y*v1.x) + m2*(r2.x*v2.y - r2.y*v2.x);
    };

    const double E0 = energy(), L0 = angmom();
    V2 e0 = lrl(r2 - r1, v2 - v1, mu);
    const double ecc0 = len(e0), ang0 = std::atan2(e0.y, e0.x);

    V2 a1, a2; accel(r1, r2, a1, a2);
    double dEmax = 0, dLmax = 0;
    for (long long s = 0; s < nsteps; s++){
        v1 = v1 + (0.5*h)*a1;  v2 = v2 + (0.5*h)*a2;
        r1 = r1 + h*v1;        r2 = r2 + h*v2;
        accel(r1, r2, a1, a2);
        v1 = v1 + (0.5*h)*a1;  v2 = v2 + (0.5*h)*a2;
        if ((s & 4095) == 0){
            dEmax = std::max(dEmax, std::fabs((energy() - E0)/E0));
            dLmax = std::max(dLmax, std::fabs((angmom() - L0)/L0));
        }
    }
    dEmax = std::max(dEmax, std::fabs((energy() - E0)/E0));
    dLmax = std::max(dLmax, std::fabs((angmom() - L0)/L0));
    V2 e1 = lrl(r2 - r1, v2 - v1, mu);
    double ang1 = std::atan2(e1.y, e1.x);
    double dang = ang1 - ang0;
    while (dang >  PI) dang -= 2.0*PI;
    while (dang < -PI) dang += 2.0*PI;
    double drift = std::fabs(dang) / 100.0;

    t.put("de_max", dEmax);
    t.put("dl_max", dLmax);
    t.put("peri_drift_per_orbit", drift);
    t.put("ecc0", ecc0);
    t.put("orbits", 100.0);
    t.pass = (dEmax < 1e-9) && (dLmax < 1e-11) && (drift < 1e-4)
          && (std::fabs(ecc0 - e) < 1e-6);
    return t;
}

// ============================================================================
// N3 — pn1: 1PN precession vs 6*pi*GM/(a(1-e^2)c^2), within 5%
// Test particle, harmonic-gauge 1PN accel; RK4; 30 perihelion passages.
// ============================================================================
static TestResult testN3(const Dials& D){
    TestResult t; t.id = "N3"; t.name = "pn1";
    const double GM = D.G * 1e4 * D.m_p, c2 = D.c*D.c;
    const double a = 10.0, e = 0.6, h = D.dt;
    const double prec_exact = 6.0*PI*GM / (a * (1.0 - e*e) * c2);

    auto acc = [&](V2 r, V2 v){
        double rn = len(r), r3 = rn*rn*rn;
        double v2 = dot(v, v), rv = dot(r, v);
        V2 aN = (-GM/r3) * r;
        V2 a1 = (GM/(c2*r3)) * ( ((4.0*GM/rn - v2)) * r + (4.0*rv) * v );
        return aN + a1;
    };
    // RK4 on (r, v)
    auto rk4 = [&](V2& r, V2& v){
        V2 k1r = v,               k1v = acc(r, v);
        V2 k2r = v + (0.5*h)*k1v, k2v = acc(r + (0.5*h)*k1r, v + (0.5*h)*k1v);
        V2 k3r = v + (0.5*h)*k2v, k3v = acc(r + (0.5*h)*k2r, v + (0.5*h)*k2v);
        V2 k4r = v + h*k3v,       k4v = acc(r + h*k3r, v + h*k3v);
        r = r + (h/6.0)*(k1r + 2.0*k2r + 2.0*k3r + k4r);
        v = v + (h/6.0)*(k1v + 2.0*k2v + 2.0*k3v + k4v);
    };

    double rap = a*(1.0 + e);
    V2 r{rap, 0}, v{0, std::sqrt(GM*(2.0/rap - 1.0/a))};
    double rdot_prev = dot(r, v)/len(r);
    std::vector<double> peri_angles;
    long long maxsteps = (long long)(35.0 * 2.0*PI*std::sqrt(a*a*a/GM) / h);
    for (long long s = 0; s < maxsteps && (int)peri_angles.size() < 30; s++){
        rk4(r, v);
        double rdot = dot(r, v)/len(r);
        if (rdot_prev < 0.0 && rdot >= 0.0){          // perihelion passage
            V2 ev = lrl(r, v, GM);
            peri_angles.push_back(std::atan2(ev.y, ev.x));
        }
        rdot_prev = rdot;
    }
    bool ok = false;
    double slope = 0, rel = 1;
    if (peri_angles.size() >= 30){
        // unwrap
        for (size_t i = 1; i < peri_angles.size(); i++){
            while (peri_angles[i] < peri_angles[i-1] - PI) peri_angles[i] += 2.0*PI;
            while (peri_angles[i] > peri_angles[i-1] + PI) peri_angles[i] -= 2.0*PI;
        }
        slope = (peri_angles.back() - peri_angles.front()) / (double)(peri_angles.size() - 1);
        rel = std::fabs(slope/prec_exact - 1.0);
        ok = rel < 0.05;
    }
    t.put("prec_meas_rad_per_orbit", slope);
    t.put("prec_exact_rad_per_orbit", prec_exact);
    t.put("rel_err", rel);
    t.put("passages", (double)peri_angles.size());
    t.pass = ok;
    return t;
}

// ============================================================================
// N4 — pw: Paczynski-Wiita ISCO at 6GM/c^2 within 1%; orbit at 5.5 GM/c^2 captured
// ============================================================================
static TestResult testN4(const Dials& D){
    TestResult t; t.id = "N4"; t.name = "pw";
    const double GM = D.G * 1e5 * D.m_p, c2 = D.c*D.c;
    const double rs = 2.0*GM/c2, h = D.dt;
    const double isco_exact = 6.0*GM/c2;

    // circular-orbit angular momentum L(r) = sqrt(GM r^3)/(r - rs); minimize
    auto Lc = [&](double r){ return std::sqrt(GM*r*r*r) / (r - rs); };
    double lo = 1.05*rs, hi = 30.0*rs;
    for (int it = 0; it < 300; it++){
        double m1 = lo + (hi - lo)/3.0, m2 = hi - (hi - lo)/3.0;
        if (Lc(m1) < Lc(m2)) hi = m2; else lo = m1;
    }
    double isco_meas = 0.5*(lo + hi);
    double rel = std::fabs(isco_meas/isco_exact - 1.0);

    // capture from unstable circular orbit at 5.5 GM/c^2 with tiny inward kick
    auto acc = [&](V2 r){
        double rn = len(r), d = rn - rs;
        double amag = GM / (d*d);
        return (-amag/rn) * r;
    };
    auto rk4 = [&](V2& r, V2& v){
        V2 k1r = v,               k1v = acc(r);
        V2 k2r = v + (0.5*h)*k1v, k2v = acc(r + (0.5*h)*k1r);
        V2 k3r = v + (0.5*h)*k2v, k3v = acc(r + (0.5*h)*k2r);
        V2 k4r = v + h*k3v,       k4v = acc(r + h*k3r);
        r = r + (h/6.0)*(k1r + 2.0*k2r + 2.0*k3r + k4r);
        v = v + (h/6.0)*(k1v + 2.0*k2v + 2.0*k3v + k4v);
    };
    double r0 = 5.5*GM/c2;
    double vphi = std::sqrt(GM*r0)/(r0 - rs);
    V2 r{r0, 0}, v{-1e-3*vphi, vphi};
    bool captured = false, bad = false;
    double t_cap = -1;
    long long maxsteps = llround(60.0 / h);
    for (long long s = 0; s < maxsteps; s++){
        rk4(r, v);
        double rn = len(r);
        if (!std::isfinite(rn)){ bad = true; break; }
        if (rn < 1.02*rs){ captured = true; t_cap = (s+1)*h; break; }
        if (rn > 50.0*rs){ break; }                    // escaped: fail
    }
    t.put("r_isco_meas", isco_meas);
    t.put("r_isco_exact", isco_exact);
    t.put("rel_err", rel);
    t.put("captured", captured ? 1.0 : 0.0);
    t.put("t_capture_s", t_cap);
    t.pass = (rel < 0.01) && captured && !bad;
    return t;
}

// ============================================================================
// N5 — qm: split-step free-packet sigma(t) within 1%; SHO E0 = hbar*w/2 within 1e-3 rel
// ============================================================================
static TestResult testN5(const Dials& D){
    TestResult t; t.id = "N5"; t.name = "qm";
    const double m = D.m_p, hb = D.hbar;

    // (a) free packet: sigma doubles
    Psi1D freeP(1024, 64.0, m, hb);
    const double s0 = 1.0;
    freeP.gaussian(-8.0, s0, 2.0);
    double t_target = 2.0*std::sqrt(3.0)*m*s0*s0/hb;
    long long nq = llround(t_target / D.dt);
    double t_actual = nq * D.dt;
    freeP.step_free_real(D.dt, (int)nq);
    double sig_meas = freeP.sigma();
    double u = hb*t_actual/(2.0*m*s0*s0);
    double sig_exp = s0*std::sqrt(1.0 + u*u);
    double sig_rel = std::fabs(sig_meas/sig_exp - 1.0);

    // (b) SHO ground state by imaginary time
    const double omega = 1.0;
    Psi1D sho(512, 32.0, m, hb);
    sho.gaussian(0.3, 1.5, 0.0);
    sho.step_sho_imag(0.002, 20000, omega);
    double e_meas = sho.energy(omega);
    double e_exp = 0.5*hb*omega;
    double e_rel = std::fabs(e_meas/e_exp - 1.0);

    t.put("sigma_meas", sig_meas);
    t.put("sigma_exp", sig_exp);
    t.put("sigma_rel_err", sig_rel);
    t.put("e0_meas", e_meas);
    t.put("e0_exp", e_exp);
    t.put("e0_rel_err", e_rel);
    t.pass = (sig_rel < 0.01) && (e_rel < 1e-3);
    return t;
}

// ============================================================================
// N6 — ratchet: MC vs min(1, p/((1-p)rho))^R at (0.3, 0.6, R in {1,4,16}), 1e7 trials
// Birth-death walk: down w.p. p/(p+(1-p)rho), up otherwise; absorb 0, cap R+30.
// ============================================================================
static TestResult testN6(const Dials&, uint64_t seed){
    TestResult t; t.id = "N6"; t.name = "ratchet";
    const double p = 0.3, rho = 0.6;
    const double lam = p / ((1.0 - p) * rho);            // 0.714286 (subcritical unwrite)
    const double pdown = p / (p + (1.0 - p)*rho);
    const int Rs[3] = {1, 4, 16};
    const long long trials = 10000000;

    std::mt19937_64 rng(splitmix64(seed ^ 6ull));
    // integer-threshold bernoulli; each 64-bit draw yields two 32-bit draws
    // (same RNG source; threshold quantization bias 2^-32, five orders below
    // the 1% gate) — keeps the full battery inside the 30 s selftest budget
    const uint32_t th32 = (uint32_t)llround(pdown * 4294967296.0);
    bool ok = true;
    for (int ri = 0; ri < 3; ri++){
        int R = Rs[ri], cap = R + 30;
        long long unwrote = 0;
        uint64_t buf = 0; int have = 0;
        for (long long tr = 0; tr < trials; tr++){
            int n = R;
            while (n > 0 && n < cap){
                if (!have){ buf = rng(); have = 2; }
                uint32_t u = (uint32_t)buf; buf >>= 32; have--;
                n += (u < th32) ? -1 : +1;
            }
            if (n == 0) unwrote++;
        }
        double pm = (double)unwrote / (double)trials;
        double pe = std::pow(std::min(1.0, lam), R);
        double rel = std::fabs(pm/pe - 1.0);
        char k1[32], k2[32], k3[32];
        snprintf(k1, sizeof k1, "p_meas_r%d", R);
        snprintf(k2, sizeof k2, "p_exact_r%d", R);
        snprintf(k3, sizeof k3, "rel_err_r%d", R);
        t.put(k1, pm); t.put(k2, pe); t.put(k3, rel);
        ok &= (rel < 0.01);
    }
    t.put("trials_per_r", (double)trials);
    t.pass = ok;
    return t;
}

// ============================================================================
// N7 — rapidity: v<c ladder; cosh path exact; naive path's cancellation demonstrated
// (D-011 erratum: 1e-12 agreement gated at gamma <= 100)
// ============================================================================
static TestResult testN7(const Dials& D){
    TestResult t; t.id = "N7"; t.name = "rapidity";
    const double c = D.c, m = D.m_p;
    bool ok = true;

    // (a) v(p) < c strictly for p = 10^0..10^9
    double worst_gap = 1e300;
    for (int e = 0; e <= 9; e++){
        double pmag = std::pow(10.0, e);
        double E = std::sqrt(pmag*pmag*c*c + m*m*c*c*c*c);
        double v = pmag*c*c/E;
        worst_gap = std::min(worst_gap, c - v);
        ok &= (v < c);
    }

    // (b) cosh path == naive path to 1e-12 for gamma <= 100
    double agree_max = 0;
    for (double g : {2.0, 10.0, 100.0}){
        double w = std::acosh(g);
        double beta = std::tanh(w);
        double g_cosh = std::cosh(w);
        double g_naive = 1.0/std::sqrt(1.0 - beta*beta);
        agree_max = std::max(agree_max, std::fabs(g_cosh/g_naive - 1.0));
    }
    ok &= (agree_max < 1e-12);

    // (c) cosh path correct at gamma = 1e7; naive path fails somewhere on the ladder
    double g_target = 1e7;
    double w7 = std::acosh(g_target);
    double g_cosh7 = std::cosh(w7);
    double cosh_rel = std::fabs(g_cosh7/g_target - 1.0);
    ok &= (cosh_rel < 1e-12);

    // ladder extends to 1e8, where beta rounds to exactly 1.0 in fp64 and the
    // naive path deterministically diverges (D-011: demonstration made
    // deterministic instead of rounding-luck-dependent)
    double naive_max_err = 0, naive_err_1e7 = 0;
    for (double g : {1e5, 3e5, 1e6, 3e6, 1e7, 3e7, 1e8}){
        double beta = std::tanh(std::acosh(g));
        double gn = 1.0/std::sqrt(1.0 - beta*beta);
        double err = std::isfinite(gn) ? std::fabs(gn/g - 1.0) : 1.0;
        if (g == 1e7) naive_err_1e7 = err;
        naive_max_err = std::max(naive_max_err, err);
    }
    double err_ratio = naive_max_err / std::max(cosh_rel, 1e-18);
    ok &= (naive_max_err > 1e-3) && (err_ratio > 1e6);

    // (d) the 8-nines ceiling: beta with 8 digits of precision cannot express 1e7
    double beta8 = 0.99999999;
    double g8 = 1.0/std::sqrt(1.0 - beta8*beta8);       // ~7071, nowhere near 1e7
    ok &= (g8 < 1e4);

    t.put("v_lt_c_min_gap", worst_gap);
    t.put("agree_max_rel_gamma_le_100", agree_max);
    t.put("cosh_rel_err_1e7", cosh_rel);
    t.put("naive_max_err_ladder", naive_max_err);
    t.put("naive_err_1e7", naive_err_1e7);
    t.put("err_ratio_naive_over_cosh", err_ratio);
    t.put("gamma_beta_8nines", g8);
    t.put("omega_1e7", w7);
    t.pass = ok;
    return t;
}

// ============================================================================
// N8 — tau: proper-time rules lock (SR exact; weak-field static; combined circular)
// ============================================================================
static TestResult testN8(const Dials& D){
    TestResult t; t.id = "N8"; t.name = "tau";
    const double c = D.c, c2 = c*c;
    bool ok = true;

    // (a) SR: sqrt(1 - v^2/c^2) vs 1/gamma(rapidity) at beta = 0.6
    double beta = 0.6;
    double tau_rule = std::sqrt(1.0 - beta*beta);
    double tau_rap = 1.0/std::cosh(std::atanh(beta));
    double sr_err = std::fabs(tau_rule/tau_rap - 1.0);
    ok &= (sr_err < 1e-12);

    // (b) weak-field static: sqrt(1 - 2GM/(r c^2)) vs Schwarzschild sqrt(1 - rs/r), r >= 20 rs
    const double GM = D.G * 1e5 * D.m_p;
    double rs = 2.0*GM/c2;
    double static_err = 0;
    for (double rr : {20.0*rs, 50.0*rs, 200.0*rs}){
        double a1 = std::sqrt(1.0 - 2.0*GM/(rr*c2));
        double a2 = std::sqrt(1.0 - rs/rr);
        static_err = std::max(static_err, std::fabs(a1/a2 - 1.0));
    }
    ok &= (static_err < 1e-6);

    // (c) circular orbit: sim rule sqrt(1 - v^2/c^2 - 2GM/(r c^2)) with v^2 = GM/r
    //     vs exact Schwarzschild circular sqrt(1 - 3GM/(r c^2)); r = 20 rs
    double rr = 20.0*rs;
    double v2 = GM/rr;
    double rule = std::sqrt(1.0 - v2/c2 - 2.0*GM/(rr*c2));
    double exact = std::sqrt(1.0 - 3.0*GM/(rr*c2));
    double circ_err = std::fabs(rule/exact - 1.0);
    ok &= (circ_err < 1e-4);

    t.put("sr_err", sr_err);
    t.put("static_err", static_err);
    t.put("circ_err", circ_err);
    t.pass = ok;
    return t;
}

// ============================================================================
// N9 — t_emit: retarded-time solver vs closed form; monotone lag; (1 - v/c) rate
// ============================================================================
static TestResult testN9(const Dials& D){
    TestResult t; t.id = "N9"; t.name = "t_emit";
    const double c = D.c, vobs = 4.0, D0 = 100.0, Om = 1.0, rorb = 2.0;
    bool ok = true;

    // generic retarded-time fixed-point: t_e = t_r - |x_obs(t_r) - x_em(t_e)|/c
    auto solve_te = [&](double tr, auto&& xem){
        double xo = D0 + vobs*tr;
        double te = tr - std::fabs(xo)/c;
        for (int it = 0; it < 80; it++){
            double ex, ey; xem(te, ex, ey);
            double dx = xo - ex, dy = -ey;
            te = tr - std::sqrt(dx*dx + dy*dy)/c;
        }
        return te;
    };

    // (i) static emitter at origin: closed form t_e = t_r - (D0 + v t_r)/c
    double solver_err = 0;
    for (double tr : {0.0, 10.0, 25.0, 50.0}){
        double te_solver = solve_te(tr, [](double, double& x, double& y){ x = 0; y = 0; });
        double te_closed = tr - (D0 + vobs*tr)/c;
        solver_err = std::max(solver_err, std::fabs(te_solver - te_closed));
    }
    ok &= (solver_err < 1e-9);

    // apparent clock rate for static emitter: dt_e/dt_r = 1 - v/c exactly
    double te_a = 0 - (D0)/c, te_b = 50.0 - (D0 + vobs*50.0)/c;
    double rate_static = (te_b - te_a)/50.0;
    ok &= (std::fabs(rate_static - (1.0 - vobs/c)) < 1e-12);

    // (ii) orbiting emitter (radius 2, Omega = 1): monotone lag + mean rate within 1%
    auto orbit = [&](double te, double& x, double& y){
        x = rorb*std::cos(Om*te); y = rorb*std::sin(Om*te);
    };
    double prev_lag = -1e300;
    bool monotone = true;
    double te0 = 0, te1 = 0;
    for (int i = 0; i <= 10; i++){
        double tr = 5.0*i;
        double te = solve_te(tr, orbit);
        if (i == 0) te0 = te;
        if (i == 10) te1 = te;
        double lag = Om*tr - Om*te;
        if (lag <= prev_lag) monotone = false;
        prev_lag = lag;
    }
    double mean_rate = (te1 - te0)/50.0;
    double mean_rel = std::fabs(mean_rate/(1.0 - vobs/c) - 1.0);
    ok &= monotone && (mean_rel < 0.01);

    t.put("solver_vs_closed", solver_err);
    t.put("rate_static", rate_static);
    t.put("rate_exact", 1.0 - vobs/c);
    t.put("mean_rate_orbiting", mean_rate);
    t.put("mean_rate_rel_err", mean_rel);
    t.put("lag_monotone", monotone ? 1.0 : 0.0);
    t.pass = ok;
    return t;
}

// ============================================================================
// N10 — compose: one scene, all regimes, 1e4 ticks. Subsystems are dynamically
// independent by declaration (cross-gravity is the CUDA sim's job); composition
// here = shared clock, shared typing, shared bookkeeping, zero contradictions.
// ============================================================================
static const uint32_t R_QUANTUM  = 0x01;
static const uint32_t R_CLASSICAL= 0x02;
static const uint32_t R_REL      = 0x04;
static const uint32_t R_COMPACT  = 0x08;
static const uint32_t R_BOUND    = 0x10;
static const uint32_t R_MASSLESS = 0x20;

static TestResult testN10(const Dials& D){
    TestResult t; t.id = "N10"; t.name = "compose";
    const double c = D.c, c2 = c*c, h = D.dt;
    const long long nticks = 10000;
    bool ok = true, nan_free = true;

    // --- subsystem A: bound binary (N2 params) ---
    const double m1 = 1e4*D.m_p, m2 = 1.0*D.m_p, Mtot = m1 + m2, mu = D.G*Mtot;
    const double a = 150.0, e = 0.6, rap = a*(1.0 + e);
    const double vap = std::sqrt(mu*(2.0/rap - 1.0/a));
    V2 br1{-(m2/Mtot)*rap, 0}, br2{(m1/Mtot)*rap, 0};
    V2 bv1{0, -(m2/Mtot)*vap}, bv2{0, (m1/Mtot)*vap};
    auto bE = [&](){
        double r = len(br2 - br1);
        return 0.5*m1*dot(bv1,bv1) + 0.5*m2*dot(bv2,bv2) - D.G*m1*m2/r;
    };
    double bE0 = bE();

    // --- subsystem B: P-W BH (M = 1e5) + orbiter at r = 8 (ISCO = 3, rs = 1) ---
    const double GMbh = D.G*1e5*D.m_p, rsbh = 2.0*GMbh/c2;
    auto pwacc = [&](V2 r){
        double rn = len(r), d = rn - rsbh;
        return (-GMbh/(d*d)/rn) * r;
    };
    double orb_r = 8.0;
    V2 orx{orb_r, 0}, orv{0, std::sqrt(GMbh*orb_r)/(orb_r - rsbh)};
    auto oE = [&](){
        double rn = len(orx);
        return 0.5*dot(orv,orv) - GMbh/(rn - rsbh);
    };
    double oE0 = oE();

    // --- subsystem C: photon (straight in nexus; lensing is a sim/renderer tier) ---
    V2 phx{0, -100.0}; V2 phdir{1.0/std::sqrt(2.0), 1.0/std::sqrt(2.0)};
    double pmag = 5.0;                                     // photon momentum
    double phE0 = pmag*c;

    // --- subsystem D: free quantum packet ---
    Psi1D pk(2048, 128.0, D.m_p, D.hbar);
    pk.gaussian(0.0, 1.0, 0.0);
    double H0 = pk.energy(0.0);

    // regime classifier (thresholds declared here; the frame contract freezes them at M1)
    auto classify = [&](double v, double mass, double dist_bh, bool has_psi, bool bound){
        uint32_t m = 0;
        if (has_psi) m |= R_QUANTUM;
        if (mass == 0.0) m |= R_MASSLESS | R_REL;
        else if (v > 0.05*c) m |= R_REL;
        else m |= R_CLASSICAL;
        if (dist_bh < 10.0*rsbh) m |= R_COMPACT;
        if (bound) m |= R_BOUND;
        return m;
    };
    auto masks_now = [&](uint32_t out[4]){
        out[0] = classify(len(bv1), m1, 1e9, false, bE() < 0);            // binary primary
        out[1] = classify(len(bv2), m2, 1e9, false, bE() < 0);            // binary secondary
        out[2] = classify(c, 0.0, 1e9, false, false);                     // photon
        out[3] = classify(len(orv), D.m_p, len(orx), false, oE() < 0);    // BH orbiter
    };
    uint32_t masks_ref[4] = {0,0,0,0}, masks_end[4] = {0,0,0,0};

    V2 ba1, ba2;
    {
        V2 d = br2 - br1; double r = len(d), r3 = r*r*r;
        ba1 = ( D.G*m2/r3)*d; ba2 = (-D.G*m1/r3)*d;
    }
    double bEdrift = 0, oEdrift = 0;
    for (long long s = 0; s < nticks; s++){
        // binary: KDK
        bv1 = bv1 + (0.5*h)*ba1; bv2 = bv2 + (0.5*h)*ba2;
        br1 = br1 + h*bv1;       br2 = br2 + h*bv2;
        {
            V2 d = br2 - br1; double r = len(d), r3 = r*r*r;
            ba1 = ( D.G*m2/r3)*d; ba2 = (-D.G*m1/r3)*d;
        }
        bv1 = bv1 + (0.5*h)*ba1; bv2 = bv2 + (0.5*h)*ba2;
        // orbiter: RK4 in P-W potential
        {
            V2 k1r = orv,               k1v = pwacc(orx);
            V2 k2r = orv + (0.5*h)*k1v, k2v = pwacc(orx + (0.5*h)*k1r);
            V2 k3r = orv + (0.5*h)*k2v, k3v = pwacc(orx + (0.5*h)*k2r);
            V2 k4r = orv + h*k3v,       k4v = pwacc(orx + h*k3r);
            orx = orx + (h/6.0)*(k1r + 2.0*k2r + 2.0*k3r + k4r);
            orv = orv + (h/6.0)*(k1v + 2.0*k2v + 2.0*k3v + k4v);
        }
        // photon: straight at c
        phx = phx + (c*h)*phdir;

        if (s == 100) masks_now(masks_ref);
        if ((s & 511) == 0){
            bEdrift = std::max(bEdrift, std::fabs((bE() - bE0)/bE0));
            oEdrift = std::max(oEdrift, std::fabs((oE() - oE0)/oE0));
            nan_free &= std::isfinite(br1.x) && std::isfinite(br2.y)
                     && std::isfinite(orx.x) && std::isfinite(orv.y)
                     && std::isfinite(phx.x);
        }
    }
    // quantum packet: evolve the same span in one exact free-space pass
    pk.step_free_real(h, (int)nticks);
    double H1 = pk.energy(0.0);
    double Hdrift = std::fabs((H1 - H0)/H0);
    double tq = nticks*h;
    double uu = D.hbar*tq/(2.0*D.m_p*1.0);
    double sig_exp = std::sqrt(1.0 + uu*uu);
    double sig_rel = std::fabs(pk.sigma()/sig_exp - 1.0);

    masks_now(masks_end);
    bool masks_stable = true;
    for (int i = 0; i < 4; i++) masks_stable &= (masks_ref[i] == masks_end[i]);
    // expected typing (frozen)
    masks_stable &= (masks_ref[0] == (R_CLASSICAL | R_BOUND));
    masks_stable &= (masks_ref[1] == (R_CLASSICAL | R_BOUND));
    masks_stable &= (masks_ref[2] == (R_MASSLESS | R_REL));
    masks_stable &= (masks_ref[3] == (R_REL | R_COMPACT | R_BOUND));

    // photon dispersion in lab and boosted frame (beta = 0.5 x-boost)
    double E = pmag*c, px = pmag*phdir.x, py = pmag*phdir.y;
    double beta = 0.5, gam = 1.0/std::sqrt(1.0 - beta*beta);
    double Eb = gam*(E - beta*c*px), pxb = gam*(px - beta*E/c);
    double disp_lab = std::fabs(E*E - c*c*(px*px + py*py))/(E*E);
    double disp_boost = std::fabs(Eb*Eb - c*c*(pxb*pxb + py*py))/(Eb*Eb);
    bool frames_ok = (disp_lab < 1e-12) && (disp_boost < 1e-12);
    double phEdrift = std::fabs((pmag*c - phE0)/phE0);     // photon E const by construction

    ok = nan_free && (bEdrift < 1e-6) && (oEdrift < 1e-6) && (Hdrift < 1e-6)
       && (phEdrift < 1e-12) && (sig_rel < 0.01) && masks_stable && frames_ok;

    t.put("e_binary_drift", bEdrift);
    t.put("e_orbiter_drift", oEdrift);
    t.put("e_photon_drift", phEdrift);
    t.put("h_packet_drift", Hdrift);
    t.put("packet_sigma_rel_err", sig_rel);
    t.put("masks_stable", masks_stable ? 1.0 : 0.0);
    t.put("mask_binary", (double)masks_ref[0]);
    t.put("mask_photon", (double)masks_ref[2]);
    t.put("mask_orbiter", (double)masks_ref[3]);
    t.put("photon_frames_ok", frames_ok ? 1.0 : 0.0);
    t.put("nan_free", nan_free ? 1.0 : 0.0);
    t.put("ticks", (double)nticks);
    t.pass = ok;
    return t;
}

// ============================================================================
// Battery + serialization
// ============================================================================
static Battery run_battery(const Dials& D, uint64_t seed){
    Battery B;
    B.results.push_back(testN1(D, B.derived));
    B.results.push_back(testN2(D));
    B.results.push_back(testN3(D));
    B.results.push_back(testN4(D));
    B.results.push_back(testN5(D));
    B.results.push_back(testN6(D, seed));
    B.results.push_back(testN7(D));
    B.results.push_back(testN8(D));
    B.results.push_back(testN9(D));
    B.results.push_back(testN10(D));
    return B;
}

// Declared JSON (hash domain: everything below; "notes" is appended outside it)
static std::string declared_json(const Dials& D, uint64_t seed, const Battery& B,
                                 const TestResult& n11){
    std::string s;
    s.reserve(8192);
    s += "{\"tool\":\"tiny_nexus\",\"version\":\"1.0.0\",\"seed\":";
    s += std::to_string(seed);
    s += ",\"dials\":{";
    s += "\"m_p\":"  + fmtnum(D.m_p)  + ",\"c\":"     + fmtnum(D.c);
    s += ",\"hbar\":" + fmtnum(D.hbar) + ",\"G\":"     + fmtnum(D.G);
    s += ",\"dt\":"   + fmtnum(D.dt)   + ",\"L_box\":" + fmtnum(D.L_box);
    s += ",\"k_B\":"  + fmtnum(D.k_B)  + "},\"derived\":{";
    for (size_t i = 0; i < B.derived.size(); i++){
        if (i) s += ",";
        s += "\"" + B.derived[i].first + "\":" + fmtnum(B.derived[i].second);
    }
    s += "},\"results\":[";
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
    for (size_t i = 0; i < B.results.size(); i++) emit(B.results[i], i == 0);
    emit(n11, false);
    bool all = B.all_pass() && n11.pass;
    s += "],\"verdict\":\"";
    s += all ? "pass" : "fail";
    s += "\"";
    return s;   // caller closes with notes + '}'
}

static void print_human(const Battery& B, const TestResult& n11){
    printf("tiny_nexus v1.0.0 — TINY UNIVERSE composition proof\n");
    printf("---------------------------------------------------\n");
    auto row = [](const TestResult& r){
        printf("  %-4s %-9s %s\n", r.id.c_str(), r.name.c_str(), r.pass ? "PASS" : "FAIL");
        if (!r.pass)
            for (auto& kv : r.m)
                printf("        %-28s %.9e\n", kv.first.c_str(), kv.second);
    };
    for (auto& r : B.results) row(r);
    row(n11);
    bool all = B.all_pass() && n11.pass;
    printf("---------------------------------------------------\n");
    printf("VERDICT: %s\n", all ? "PASS" : "FAIL");
}

// ----------------------------------------------------------------------------
// Dials file (flat JSON, partial overrides allowed)
// ----------------------------------------------------------------------------
static bool load_dials(const char* path, Dials& D){
    FILE* f = fopen(path, "rb");
    if (!f) return false;
    std::string txt;
    char buf[4096]; size_t n;
    while ((n = fread(buf, 1, sizeof buf, f)) > 0) txt.append(buf, n);
    fclose(f);
    auto grab = [&](const char* key, double& out){
        std::string pat = std::string("\"") + key + "\"";
        size_t p = txt.find(pat);
        if (p == std::string::npos) return;
        p = txt.find(':', p);
        if (p == std::string::npos) return;
        out = strtod(txt.c_str() + p + 1, nullptr);
    };
    grab("m_p", D.m_p); grab("c", D.c); grab("hbar", D.hbar); grab("G", D.G);
    grab("dt", D.dt); grab("L_box", D.L_box); grab("k_B", D.k_B);
    return true;
}

// ----------------------------------------------------------------------------
// main
// ----------------------------------------------------------------------------
int main(int argc, char** argv){
    Dials D;
    uint64_t seed = 20260711ull;
    bool json = false, selftest = false, golden = false;
    const char* only = nullptr;

    for (int i = 1; i < argc; i++){
        std::string a = argv[i];
        if (a == "--json") json = true;
        else if (a == "--selftest") selftest = true;
        else if (a == "--golden") golden = true;
        else if (a == "--seed" && i + 1 < argc) seed = strtoull(argv[++i], nullptr, 10);
        else if (a == "--only" && i + 1 < argc) only = argv[++i];
        else if (a == "--dials" && i + 1 < argc){
            if (!load_dials(argv[++i], D)){
                fprintf(stderr, "error: cannot read dials file\n");
                return 2;
            }
        }
        else {
            fprintf(stderr, "usage: tiny_nexus [--dials PATH] [--seed N] [--json] "
                            "[--selftest] [--golden] [--only ID]\n");
            return 2;
        }
    }

    if (only){
        Dials Dd = D;
        std::string id = only;
        TestResult r;
        std::vector<std::pair<std::string,double>> derived;
        if      (id == "N1") r = testN1(Dd, derived);
        else if (id == "N2") r = testN2(Dd);
        else if (id == "N3") r = testN3(Dd);
        else if (id == "N4") r = testN4(Dd);
        else if (id == "N5") r = testN5(Dd);
        else if (id == "N6") r = testN6(Dd, seed);
        else if (id == "N7") r = testN7(Dd);
        else if (id == "N8") r = testN8(Dd);
        else if (id == "N9") r = testN9(Dd);
        else if (id == "N10") r = testN10(Dd);
        else { fprintf(stderr, "error: unknown test id (N1..N10)\n"); return 2; }
        printf("%s %s: %s\n", r.id.c_str(), r.name.c_str(), r.pass ? "PASS" : "FAIL");
        for (auto& kv : r.m) printf("  %-28s %.9e\n", kv.first.c_str(), kv.second);
        return r.pass ? 0 : 1;
    }

    if (golden){ D = Dials(); seed = 20260711ull; }     // golden = frozen defaults

    // battery, twice (N11 determinism gate compares declared strings in-process)
    Battery B1 = run_battery(D, seed);
    Battery B2 = run_battery(D, seed);
    TestResult probe; probe.id = "N11"; probe.name = "det"; probe.pass = true;   // placeholder
    std::string s1 = declared_json(D, seed, B1, probe);
    std::string s2 = declared_json(D, seed, B2, probe);
    TestResult n11; n11.id = "N11"; n11.name = "det";
    bool identical = (s1 == s2);
    n11.put("identical", identical ? 1.0 : 0.0);
    n11.put("inprocess_repeats", 2.0);
    n11.pass = identical;

    std::string declared = declared_json(D, seed, B1, n11);
    std::string hash = blake2b::hash256_hex(declared);
    bool all = B1.all_pass() && n11.pass;

    if (golden){
        FILE* f = fopen("goldens/nexus/golden.hash", "rb");
        if (!f){
            fprintf(stderr, "GOLDEN NOT FROZEN %s\n", hash.c_str());
            printf("%s\n", hash.c_str());
            return 2;
        }
        char want[128] = {0};
        if (fscanf(f, "%127s", want) != 1) { fclose(f); fprintf(stderr, "GOLDEN FILE UNREADABLE\n"); return 2; }
        fclose(f);
        if (hash == want){ fprintf(stderr, "GOLDEN OK %.8s\n", hash.c_str()); printf("%s\n", hash.c_str()); return 0; }
        fprintf(stderr, "GOLDEN MISMATCH have=%.8s want=%.8s\n", hash.c_str(), want);
        printf("%s\n", hash.c_str());
        return 1;
    }

    if (json){
        std::string out = declared;
        out += ",\"notes\":\"non-declared; hash=" + hash.substr(0, 8) + "\"}";
        printf("%s\n", out.c_str());
        return all ? 0 : 1;
    }

    // default and --selftest: human-readable
    (void)selftest;
    print_human(B1, n11);
    printf("declared hash: %s\n", hash.c_str());
    return all ? 0 : 1;
}
