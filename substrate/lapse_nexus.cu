// ============================================================================
//  lapse_nexus — TINY UNIVERSE v2 N2 "lapse" the-clock substrate oracle tool
//  Contract: contracts/lapse.contract.md v1.0.0
//
//  Single-file CUDA. Turns the substrate's Newtonian potential Phi into a CLOCK:
//    lapse   alpha(x) = sqrt(1 + 2 Phi/c^2)      (static-observer clock rate)
//    proper  tau(x)  = integral alpha dt          (declared, hashed per-cell field)
//  Gravitational redshift z = 1/sqrt(1 - r_s/r) - 1 (r_s = 2GM/c^2). For a point
//  mass Phi=-GM/r the lapse is EXACTLY the Schwarzschild static lapse sqrt(1-r_s/r),
//  so the redshift oracle is exact (temporal metric only; light-bending is N3).
//
//  Reuses N1's PM cuFFT-Poisson (kGreen lineage, VERBATIM) for the redshiftPM weld.
//  Envelope face (--scenario X --json|--golden|--selftest). Exit 0/1/2, preflight 3.
//
//  Scenarios:
//    redshift   : ANALYTIC Phi=-GM/r -> alpha,tau -> z(r) vs exact Schwarzschild.
//                 The lapse+proper-time machinery oracle (like N1 freepacket/sho3d).
//    redshiftPM : PM-Poisson Phi of a resolved mass -> lapse -> fit well depth A=GM.
//                 The substrate weld: substrate gravity -> correct gravitational
//                 redshift (like N1 soliton / --poissontest).
//
//  Build (BUILD.md; no fast-math — Invariant 6):
//    nvcc -O3 -arch=sm_89 -Xcompiler "/O2" -o build\lapse_nexus.exe
//         substrate\lapse_nexus.cu cufft.lib
//
//  Dials (nexus v0, frozen): c=20 (LIVE), G=2e-3, m=1, hbar=0.5(inert), dt=1/240,
//  L_box=512. Exit: 0 pass · 1 gate fired · 2 error · 3 GPU contended.
// ============================================================================

#define _CRT_SECURE_NO_WARNINGS
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <cuda_runtime.h>
#include <cufft.h>

#define PI_F 3.14159265358979f

// ---- dials (nexus v0 defaults, frozen). c is LIVE in N2. ---------------------
struct Dials {
    double c     = 20.0;      // su/s — LIVE: lapse = sqrt(1 + 2 Phi/c^2)
    double hbar  = 0.5;       // inert in N2 (no psi evolution) — declared, unused
    double m     = 1.0;
    double G     = 2.0e-3;
    double dt    = 1.0 / 240.0;
    double L_box = 512.0;
};

// ============================================================================
//  BLAKE2b-256 (RFC 7693) — lifted VERBATIM from field_nexus.cu / tiny_nexus.cpp.
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
static std::string hash256_bytes(const void* data, size_t n){
    std::string s((const char*)data, n);
    return hash256_hex(s);
}
} // namespace blake2b

// ---- canonical number formatting (fmt6: %.6f ; fmt9: %.9e) -------------------
static std::string fmt6(double x){
    if (!std::isfinite(x)) return "9999999.999999";
    if (x == 0.0) x = 0.0;
    char buf[48]; snprintf(buf, sizeof(buf), "%.6f", x); return std::string(buf);
}
static std::string fmt9(double x){
    if (!std::isfinite(x)) return "9.999999999e+99";
    if (x == 0.0) x = 0.0;
    char buf[48]; snprintf(buf, sizeof(buf), "%.9e", x); return std::string(buf);
}

// ---- fixed-point accumulator (Invariant 4; liborrery reduce.cuh idiom) -------
__host__ __device__ inline long long fp_encode(double v){ return (long long)(v * 4294967296.0); }
__host__ __device__ inline double     fp_decode(long long q){ return (double)q / 4294967296.0; }
__device__ inline void fp_atomic_add(unsigned long long* acc, double v){
    atomicAdd(acc, (unsigned long long)fp_encode(v));
}

#define CUDA_CHECK(call) do {                                                  \
    cudaError_t _e = (call);                                                   \
    if (_e != cudaSuccess) {                                                   \
        fprintf(stderr, "CUDA error %s (%d) at %s:%d\n  %s\n",                 \
                cudaGetErrorName(_e), (int)_e, __FILE__, __LINE__, #call);     \
        std::exit(2);                                                          \
    }                                                                          \
} while (0)
#define CUFFT_CHECK(call) do { cufftResult _r = (call); if (_r != CUFFT_SUCCESS){ \
    fprintf(stderr, "cuFFT error %d at %s:%d\n", (int)_r, __FILE__, __LINE__); std::exit(2);} } while (0)

// ============================================================================
//  KERNELS
// ============================================================================

// analytic point-mass potential Phi(x) = -GM / sqrt(r^2 + eps^2), r from box center
__global__ void kFillPointPhi(float* phi, int LN, float gL, float dx, float GM, float eps){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    int N3 = LN*LN*LN; if (idx >= N3) return;
    int ix = idx % LN, iy = (idx/LN) % LN, iz = idx/(LN*LN);
    float x = -gL*0.5f + (ix + 0.5f)*dx;
    float y = -gL*0.5f + (iy + 0.5f)*dx;
    float z = -gL*0.5f + (iz + 0.5f)*dx;
    float r2 = x*x + y*y + z*z;
    phi[idx] = -GM / sqrtf(r2 + eps*eps);
}

// real Gaussian source psi = exp(-r^2/(4 sigma^2)), phase 0 (redshiftPM mass blob)
__global__ void kGauss3(cufftComplex* q, int LN, float gL, float dx, float s0){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    int N3 = LN*LN*LN; if (idx >= N3) return;
    int ix = idx % LN, iy = (idx/LN) % LN, iz = idx/(LN*LN);
    float x = -gL*0.5f + (ix + 0.5f)*dx;
    float y = -gL*0.5f + (iy + 0.5f)*dx;
    float z = -gL*0.5f + (iz + 0.5f)*dx;
    float env = expf(-0.25f*(x*x + y*y + z*z)/(s0*s0));
    q[idx].x = env; q[idx].y = 0.0f;
}

__global__ void kScale(cufftComplex* q, int N3, float s){
    int idx = blockIdx.x*blockDim.x + threadIdx.x; if (idx >= N3) return;
    q[idx].x *= s; q[idx].y *= s;
}

// psi-mass deposit -> fixed-point uint64 grid (1:1 gather, VERBATIM from N1). Inv 4.
__global__ void kPsiDeposit(const cufftComplex* q, unsigned long long* g, int N3, double mcell){
    int idx = blockIdx.x*blockDim.x + threadIdx.x; if (idx >= N3) return;
    double rho = ((double)q[idx].x*q[idx].x + (double)q[idx].y*q[idx].y)*mcell;
    fp_atomic_add(&g[idx], rho);
}
__global__ void kFixToReal(const unsigned long long* g, float* r, int N3){
    int i = blockIdx.x*blockDim.x + threadIdx.x; if (i < N3) r[i] = (float)fp_decode((long long)g[i]);
}
__global__ void kZeroFix(unsigned long long* g, int n){
    int i = blockIdx.x*blockDim.x + threadIdx.x; if (i < n) g[i] = 0ull;
}
// PM Green multiply (VERBATIM from N1 kGreen; k=0 zeroed = mean does not gravitate)
__global__ void kGreen(cufftComplex* s, int LN, int FNZC, float gL, float dx, float G){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    int total = LN*LN*FNZC; if (i >= total) return;
    int z = i % FNZC, y = (i/FNZC) % LN, x = i/(FNZC*LN);
    int fx = (x <= LN/2) ? x : x - LN;
    int fy = (y <= LN/2) ? y : y - LN;
    int fz = z;
    float kf = 2.0f*PI_F/gL;
    float k2 = kf*kf*(float)(fx*fx + fy*fy + fz*fz);
    float f = (k2 > 0.0f) ? -4.0f*PI_F*G/(k2 * (dx*dx*dx) * (float)(LN*LN*LN)) : 0.0f;
    s[i].x *= f; s[i].y *= f;
}

// lapse: alpha = sqrt( max(1 + 2 Phi/c^2, floor^2) ). Phi<0 in wells -> alpha<1.
// floor clamps the horizon interior (1+2Phi/c^2 <= 0) to a tiny positive alpha
// (declared) so no NaN; the gates probe only r > r_s where the argument is > 0.
__global__ void kLapse(float* alpha, const float* phi, int N3, float invc2, float floor2){
    int i = blockIdx.x*blockDim.x + threadIdx.x; if (i >= N3) return;
    float a = 1.0f + 2.0f*phi[i]*invc2;
    alpha[i] = sqrtf(fmaxf(a, floor2));
}
// proper-time accumulate one tick: tau += alpha*dt (per-cell clock integrator)
__global__ void kAccumTau(float* tau, const float* alpha, int N3, float dt){
    int i = blockIdx.x*blockDim.x + threadIdx.x; if (i >= N3) return;
    tau[i] += alpha[i]*dt;
}

// ============================================================================
//  Grid + host helpers
// ============================================================================
struct Grid {
    int LN, N3, FNZC;
    float gL, dx;
    float* dPhi = nullptr;      // potential (analytic fill) OR density->Phi (PM)
    float* dAlpha = nullptr;    // lapse
    float* dTau = nullptr;      // proper time (declared field)
    // PM-only working set
    cufftComplex* dQ = nullptr;
    unsigned long long* dFix = nullptr;
    cufftComplex* dSpec = nullptr;
    cufftHandle planR2C, planC2R;
    bool poisson = false;
};
static dim3 gridBlocks(int n, int b){ return dim3((n + b - 1)/b); }

static void gridAlloc(Grid& g, bool withPoisson){
    g.N3 = g.LN*g.LN*g.LN; g.FNZC = g.LN/2 + 1; g.dx = g.gL/g.LN;
    CUDA_CHECK(cudaMalloc(&g.dPhi,   (size_t)g.N3*sizeof(float)));
    CUDA_CHECK(cudaMalloc(&g.dAlpha, (size_t)g.N3*sizeof(float)));
    CUDA_CHECK(cudaMalloc(&g.dTau,   (size_t)g.N3*sizeof(float)));
    CUDA_CHECK(cudaMemset(g.dTau, 0, (size_t)g.N3*sizeof(float)));
    g.poisson = withPoisson;
    if (withPoisson){
        CUDA_CHECK(cudaMalloc(&g.dQ,   (size_t)g.N3*sizeof(cufftComplex)));
        CUDA_CHECK(cudaMalloc(&g.dFix, (size_t)g.N3*sizeof(unsigned long long)));
        CUDA_CHECK(cudaMalloc(&g.dSpec,(size_t)g.LN*g.LN*g.FNZC*sizeof(cufftComplex)));
        CUFFT_CHECK(cufftPlan3d(&g.planR2C, g.LN, g.LN, g.LN, CUFFT_R2C));
        CUFFT_CHECK(cufftPlan3d(&g.planC2R, g.LN, g.LN, g.LN, CUFFT_C2R));
    }
}
static void gridFree(Grid& g){
    cudaFree(g.dPhi); cudaFree(g.dAlpha); cudaFree(g.dTau);
    if (g.poisson){ cudaFree(g.dQ); cudaFree(g.dFix); cudaFree(g.dSpec);
        cufftDestroy(g.planR2C); cufftDestroy(g.planC2R); }
}

// integrate proper time: tau = 0; N ticks of tau += alpha*dt (per-cell). Static
// lapse -> tau = alpha*N*dt, computed as a genuine per-tick accumulation (the clock).
static void integrateTau(Grid& g, int nsteps, float dt){
    dim3 b(256), gr = gridBlocks(g.N3, 256);
    CUDA_CHECK(cudaMemset(g.dTau, 0, (size_t)g.N3*sizeof(float)));
    for (int s = 0; s < nsteps; s++) kAccumTau<<<gr, b>>>(g.dTau, g.dAlpha, g.N3, dt);
    CUDA_CHECK(cudaDeviceSynchronize());
}
static std::string captureTauHash(Grid& g){
    std::vector<float> h((size_t)g.N3);
    CUDA_CHECK(cudaMemcpy(h.data(), g.dTau, (size_t)g.N3*sizeof(float), cudaMemcpyDeviceToHost));
    return blake2b::hash256_bytes(h.data(), (size_t)g.N3*sizeof(float));
}
static double hostNormMass(Grid& g){
    std::vector<cufftComplex> h((size_t)g.N3);
    CUDA_CHECK(cudaMemcpy(h.data(), g.dQ, (size_t)g.N3*sizeof(cufftComplex), cudaMemcpyDeviceToHost));
    double s = 0.0; for (int i=0;i<g.N3;i++){ double a=h[i].x,b=h[i].y; s += a*a+b*b; }
    return s * (double)g.dx*g.dx*g.dx;
}

// ============================================================================
//  Result + JSON
// ============================================================================
struct Result {
    std::string scenario;
    int grid=0; long long steps=0;
    double gL=0, M=0, r_s=0;
    // redshift (analytic)
    double z_inner=0, z_inner_ana=0, r_inner=0;   // strong-field probe (~2 su)
    double z_mid=0, z_mid_ana=0, r_mid=0;         // mid probe (~8 su)
    double redshift_rel=0;                         // max |alpha_meas/alpha_ana - 1|
    double tau_rel=0;                              // |tau/(alpha*N*dt) - 1| max
    double alpha_min=0;
    // redshiftPM (weld)
    double A_fit=0, A_expect=0, A_rel=0, z_at8=0;
    double mass0=0;
    bool nan_free=true;
    bool gate_primary=false, gate_norm=false, gate_nan=false, verdict=false;
    std::string state_b2b;
};

static std::string declaredJson(const Dials& D, uint64_t seed, const Result& R){
    std::string s; s.reserve(2048);
    s += "{\"tool\":\"lapse_nexus\",\"version\":\"1.0.0\",\"seed\":";
    s += std::to_string(seed);
    s += ",\"params\":{";
    s += "\"scenario\":\"" + R.scenario + "\"";
    s += ",\"grid\":" + std::to_string(R.grid);
    s += ",\"steps\":" + std::to_string(R.steps);
    s += ",\"c\":" + fmt6(D.c) + ",\"G\":" + fmt9(D.G) + ",\"m\":" + fmt6(D.m);
    s += ",\"hbar\":" + fmt6(D.hbar) + ",\"dt\":" + fmt9(D.dt) + ",\"L_box\":" + fmt6(D.L_box);
    s += ",\"gL\":" + fmt6(R.gL) + ",\"M\":" + fmt6(R.M) + ",\"r_s\":" + fmt6(R.r_s);
    s += "},\"result\":{";
    s += "\"state_b2b\":\"" + R.state_b2b + "\"";
    if (R.scenario == "redshift"){
        s += ",\"r_inner\":" + fmt6(R.r_inner) + ",\"z_inner\":" + fmt6(R.z_inner);
        s += ",\"z_inner_ana\":" + fmt6(R.z_inner_ana);
        s += ",\"r_mid\":" + fmt6(R.r_mid) + ",\"z_mid\":" + fmt6(R.z_mid);
        s += ",\"z_mid_ana\":" + fmt6(R.z_mid_ana);
        s += ",\"redshift_rel\":" + fmt9(R.redshift_rel);
        s += ",\"tau_rel\":" + fmt9(R.tau_rel);
        s += ",\"alpha_min\":" + fmt6(R.alpha_min);
    } else if (R.scenario == "redshiftPM"){
        s += ",\"A_fit\":" + fmt6(R.A_fit) + ",\"A_expect\":" + fmt6(R.A_expect);
        s += ",\"A_rel\":" + fmt9(R.A_rel) + ",\"z_at8\":" + fmt6(R.z_at8);
        s += ",\"mass0\":" + fmt6(R.mass0) + ",\"alpha_min\":" + fmt6(R.alpha_min);
    }
    s += ",\"nan_free\":" + std::string(R.nan_free ? "1":"0");
    s += "},\"gates\":{";
    s += "\"primary\":" + std::string(R.gate_primary ? "true":"false");
    s += ",\"norm\":"   + std::string(R.gate_norm ? "true":"false");
    s += ",\"nan\":"    + std::string(R.gate_nan ? "true":"false");
    s += "},\"verdict\":\"";
    s += R.verdict ? "pass" : "fail";
    s += "\"";
    return s;
}

// ============================================================================
//  SCENARIO: redshift — ANALYTIC Phi=-GM/r -> lapse -> tau -> z(r) vs Schwarzschild.
//  Probes the cell nearest each target radius (uses the cell's EXACT r), so the
//  comparison is cell-exact (no shell-averaging bias on the steep alpha(r)).
// ============================================================================
static Result runRedshift(const Dials& D, int gridReq){
    Result R; R.scenario = "redshift";
    const int LN = (gridReq>0? gridReq : 128);
    Grid g; g.LN = LN; g.gL = 64.0f; gridAlloc(g, false);
    R.grid = LN; R.gL = g.gL;
    const double M = 1.0e5;                       // r_s = 2GM/c^2 = 1 su
    const double GM = D.G*M;
    const double r_s = 2.0*GM/(D.c*D.c);
    const double eps = 0.01;                       // core softening (< smallest probe)
    R.M = M; R.r_s = r_s;
    const int Nsteps = 480; R.steps = Nsteps;
    dim3 b(256), gr = gridBlocks(g.N3, 256);
    kFillPointPhi<<<gr,b>>>(g.dPhi, LN, g.gL, g.dx, (float)GM, (float)eps);
    kLapse<<<gr,b>>>(g.dAlpha, g.dPhi, g.N3, (float)(1.0/(D.c*D.c)), 1e-8f);
    integrateTau(g, Nsteps, (float)D.dt);

    std::vector<float> ha((size_t)g.N3), ht((size_t)g.N3);
    CUDA_CHECK(cudaMemcpy(ha.data(), g.dAlpha, (size_t)g.N3*sizeof(float), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(ht.data(), g.dTau,   (size_t)g.N3*sizeof(float), cudaMemcpyDeviceToHost));

    // probe radii: strong -> weak, plus a far reference. Find the cell whose center
    // r is closest to each target; use that cell's EXACT r for the analytic value.
    const double targets[5] = {2.0, 4.0, 8.0, 16.0, 28.0};
    double rcell[5]={0}, acell[5]={0}, tcell[5]={0}, bestd[5]={1e30,1e30,1e30,1e30,1e30};
    double amin = 1e30; double Ttot = Nsteps*D.dt;
    double gL=g.gL, dx=g.dx; int N=LN;
    for (int iz=0; iz<N; iz++){ double z=-gL*0.5+(iz+0.5)*dx;
      for (int iy=0; iy<N; iy++){ double y=-gL*0.5+(iy+0.5)*dx;
        for (int ix=0; ix<N; ix++){ double x=-gL*0.5+(ix+0.5)*dx;
            size_t idx=((size_t)iz*N+iy)*N+ix;
            double a=ha[idx]; if (a<amin) amin=a;
            double r=std::sqrt(x*x+y*y+z*z);
            for (int p=0;p<5;p++){ double dd=std::fabs(r-targets[p]);
                if (dd<bestd[p]){ bestd[p]=dd; rcell[p]=r; acell[p]=a; tcell[p]=ht[idx]; } }
        }}}
    double maxrel=0, maxtaurel=0;
    auto zana=[&](double r){ return 1.0/std::sqrt(1.0 - r_s/r) - 1.0; };
    for (int p=0;p<5;p++){
        double a_ana = std::sqrt(1.0 - r_s/rcell[p]);
        double rel = std::fabs(acell[p]/a_ana - 1.0);
        if (rel>maxrel) maxrel=rel;
        double taurel = std::fabs(tcell[p]/(acell[p]*Ttot) - 1.0);
        if (taurel>maxtaurel) maxtaurel=taurel;
    }
    R.r_inner=rcell[0]; R.z_inner=1.0/acell[0]-1.0; R.z_inner_ana=zana(rcell[0]);
    R.r_mid=rcell[2];   R.z_mid=1.0/acell[2]-1.0;   R.z_mid_ana=zana(rcell[2]);
    R.redshift_rel=maxrel; R.tau_rel=maxtaurel; R.alpha_min=amin;
    R.nan_free = std::isfinite(maxrel) && std::isfinite(amin);
    R.state_b2b = captureTauHash(g);
    R.gate_primary = (maxrel < 1e-3);
    R.gate_norm    = (maxtaurel < 5e-5);       // the proper-time integrator (fp32 accum floor ~7e-6)
    R.gate_nan     = R.nan_free;
    R.verdict = R.gate_primary && R.gate_norm && R.gate_nan;
    gridFree(g);
    return R;
}

// ============================================================================
//  SCENARIO: redshiftPM — the substrate weld. PM-Poisson Phi of a resolved mass
//  blob -> lapse -> invert to Phi_meas=1/2 c^2(alpha^2-1) -> fit -A/r+C over the
//  far field; A should equal G*M (Newtonian well depth) -> the redshift magnitude
//  is right. Same PM solver as N1 (kGreen verbatim); c goes live.
// ============================================================================
static Result runRedshiftPM(const Dials& D, int gridReq){
    Result R; R.scenario = "redshiftPM";
    const int LN = (gridReq>0? gridReq : 128);
    Grid g; g.LN = LN; g.gL = 64.0f; gridAlloc(g, true);
    R.grid = LN; R.gL = g.gL;
    const double M = 1.0e5, sigma = 2.0;           // compact resolved source (poissontest regime)
    const double GM = D.G*M;
    const double r_s = 2.0*GM/(D.c*D.c);
    R.M = M; R.r_s = r_s;
    const int Nsteps = 480; R.steps = Nsteps;
    dim3 b(256), gr = gridBlocks(g.N3, 256);
    int nspec = g.LN*g.LN*g.FNZC; dim3 grS = gridBlocks(nspec,256), grF = gridBlocks(g.N3,256);
    // source psi = Gaussian, normalized so integral m|psi|^2 dx^3 = M
    kGauss3<<<gr,b>>>(g.dQ, LN, g.gL, g.dx, (float)sigma);
    { double mass=hostNormMass(g); double sc=std::sqrt((M/D.m)/mass);
      kScale<<<gr,b>>>(g.dQ, g.N3, (float)sc); }
    R.mass0 = M;
    // PM Poisson (VERBATIM N1 weld path): deposit -> rho -> Phi (true, M2 convention)
    kZeroFix<<<grF,b>>>(g.dFix, g.N3);
    kPsiDeposit<<<grF,b>>>(g.dQ, g.dFix, g.N3, (double)D.m*g.dx*g.dx*g.dx);
    kFixToReal<<<grF,b>>>(g.dFix, g.dPhi, g.N3);           // dPhi = density
    CUFFT_CHECK(cufftExecR2C(g.planR2C, g.dPhi, g.dSpec));
    kGreen<<<grS,b>>>(g.dSpec, g.LN, g.FNZC, g.gL, g.dx, (float)D.G);
    CUFFT_CHECK(cufftExecC2R(g.planC2R, g.dSpec, g.dPhi));  // dPhi = Phi (true)
    // lapse + proper time
    kLapse<<<gr,b>>>(g.dAlpha, g.dPhi, g.N3, (float)(1.0/(D.c*D.c)), 1e-8f);
    integrateTau(g, Nsteps, (float)D.dt);

    std::vector<float> ha((size_t)g.N3);
    CUDA_CHECK(cudaMemcpy(ha.data(), g.dAlpha, (size_t)g.N3*sizeof(float), cudaMemcpyDeviceToHost));
    // radial shell-mean of the lapse-inverted potential Phi_meas = 1/2 c^2 (alpha^2-1)
    double gL=g.gL, dx=g.dx; int N=LN; double c2=D.c*D.c;
    const int NB=256; std::vector<double> psum(NB,0), pcnt(NB,0); double amin=1e30;
    for (int iz=0; iz<N; iz++){ double z=-gL*0.5+(iz+0.5)*dx;
      for (int iy=0; iy<N; iy++){ double y=-gL*0.5+(iy+0.5)*dx;
        for (int ix=0; ix<N; ix++){ double x=-gL*0.5+(ix+0.5)*dx;
            size_t idx=((size_t)iz*N+iy)*N+ix; double a=ha[idx]; if(a<amin)amin=a;
            double phim = 0.5*c2*((double)a*a - 1.0);
            double r=std::sqrt(x*x+y*y+z*z); int bi=(int)(r/(0.5*gL)*NB); if(bi>=NB)bi=NB-1;
            psum[bi]+=phim; pcnt[bi]+=1; }}}
    auto phiAtR=[&](double rq)->double{ int bi=(int)(rq/(0.5*gL)*NB); if(bi>=NB)bi=NB-1;
        return pcnt[bi]>0? psum[bi]/pcnt[bi]:0.0; };
    // fit Phi_meas = -A/r + C at two far radii (r >> sigma, r << box/2): A = GM.
    double r1=8.0, r2=16.0, p1=phiAtR(r1), p2=phiAtR(r2);
    double A = (p2-p1)/(1.0/r1 - 1.0/r2);
    double A_expect = GM;
    R.A_fit=A; R.A_expect=A_expect; R.A_rel=std::fabs(A/A_expect - 1.0);
    // implied clock rate at r=8 from the fitted well: alpha=sqrt(1+2 Phi/c^2), z=1/alpha-1
    double phi8 = -A/8.0;   // use the fitted well (C cancels for the depth)
    double a8 = std::sqrt(std::fmax(1.0 + 2.0*phi8/c2, 1e-8));
    R.z_at8 = 1.0/a8 - 1.0;
    R.alpha_min = amin;
    R.state_b2b = captureTauHash(g);
    R.nan_free = std::isfinite(A) && std::isfinite(amin);
    R.gate_primary = (R.A_rel < 0.05);          // Newtonian well depth (poissontest regime)
    R.gate_norm    = true;                        // static source; mass is the deposited M
    R.gate_nan     = R.nan_free;
    R.verdict = R.gate_primary && R.gate_nan;
    gridFree(g);
    return R;
}

static void printHuman(const Result& R){
    printf("lapse_nexus v1.0.0 — TINY UNIVERSE v2 N2 lapse (the clock)\n");
    printf("scenario: %s   grid: %d^3   steps: %lld\n", R.scenario.c_str(), R.grid, R.steps);
    printf("  M=%.0f  r_s=2GM/c^2=%.4f su\n", R.M, R.r_s);
    printf("-------------------------------------------------------\n");
    if (R.scenario == "redshift"){
        printf("  gravitational time dilation vs exact Schwarzschild z=1/sqrt(1-r_s/r)-1:\n");
        printf("    r=%.3f (strong):  z_meas=%.6f  z_ana=%.6f\n", R.r_inner, R.z_inner, R.z_inner_ana);
        printf("    r=%.3f (mid):     z_meas=%.6f  z_ana=%.6f\n", R.r_mid, R.z_mid, R.z_mid_ana);
        printf("    max |alpha_meas/alpha_ana - 1|  %.3e   [gate < 1e-3]  %s\n",
               R.redshift_rel, R.gate_primary?"PASS":"FAIL");
        printf("    proper-time integrator |tau/(alpha*T)-1|  %.3e   [gate < 5e-5]  %s\n",
               R.tau_rel, R.gate_norm?"PASS":"FAIL");
        printf("    alpha_min (deepest clock)  %.6f\n", R.alpha_min);
    } else if (R.scenario == "redshiftPM"){
        printf("  substrate weld: PM-Poisson well through the lapse (A should = G*M):\n");
        printf("    fitted well depth A   %.6f\n", R.A_fit);
        printf("    expected G*M          %.6f\n", R.A_expect);
        printf("    A/GM rel error        %.3e   [gate < 5e-2]  %s\n", R.A_rel, R.gate_primary?"PASS":"FAIL");
        printf("    implied redshift z(r=8)  %.6f\n", R.z_at8);
        printf("    alpha_min             %.6f\n", R.alpha_min);
    }
    printf("  nan_free      %s\n", R.nan_free?"yes":"NO");
    printf("-------------------------------------------------------\n");
    printf("VERDICT: %s\n", R.verdict ? "PASS" : "FAIL");
}

// ---- GPU preflight (harness exit-3 lineage) ---------------------------------
static int gpuPreflight(size_t need_bytes){
    size_t freeB=0, totalB=0;
    cudaError_t e = cudaMemGetInfo(&freeB, &totalB);
    if (e != cudaSuccess){ fprintf(stderr, "[preflight] cudaMemGetInfo failed: %s\n", cudaGetErrorString(e)); return 2; }
    size_t min_free = need_bytes + (size_t)512*1024*1024;
    if (freeB < min_free){
        fprintf(stderr, "[preflight] GPU busy: %.0f MB free (< %.0f MB needed). Exit 3.\n",
                freeB/1048576.0, min_free/1048576.0);
        return 3;
    }
    return 0;
}
static size_t needBytes(const std::string& sc, int LN){
    size_t N3=(size_t)LN*LN*LN;
    size_t base = N3*sizeof(float)*3;                  // Phi + alpha + tau
    if (sc == "redshiftPM"){
        base += N3*sizeof(cufftComplex) + N3*sizeof(unsigned long long)
              + (size_t)LN*LN*(LN/2+1)*sizeof(cufftComplex)*3;
    }
    return base;
}
static Result dispatch(const std::string& sc, const Dials& D, int grid){
    if (sc == "redshift")   return runRedshift(D, grid);
    if (sc == "redshiftPM") return runRedshiftPM(D, grid);
    fprintf(stderr, "error: unknown scenario '%s' (redshift|redshiftPM)\n", sc.c_str());
    std::exit(2);
}

int main(int argc, char** argv){
    Dials D;
    uint64_t seed = 20260711ull;
    std::string scenario = "redshift";
    int grid = 128;
    bool json=false, selftest=false, golden=false;
    for (int i=1;i<argc;i++){
        std::string a = argv[i];
        if      (a == "--json") json = true;
        else if (a == "--selftest") selftest = true;
        else if (a == "--golden") golden = true;
        else if (a == "--scenario" && i+1<argc) scenario = argv[++i];
        else if (a == "--grid" && i+1<argc) grid = atoi(argv[++i]);
        else if (a == "--seed" && i+1<argc) seed = strtoull(argv[++i], nullptr, 10);
        else { fprintf(stderr, "usage: lapse_nexus --scenario redshift|redshiftPM "
                       "[--grid N=128] [--seed N] [--json] [--golden] [--selftest]\n"); return 2; }
    }

    if (golden){ D = Dials(); seed = 20260711ull; grid = 128; }

    int pf = gpuPreflight(needBytes(scenario, grid));
    if (pf != 0) return pf;

    if (selftest){
        // flatlapse: Phi=0 -> alpha=1 exactly; tau = N*dt exactly.
        Grid g; g.LN = 64; g.gL = 64.0f; gridAlloc(g, false);
        dim3 b(256), gr = gridBlocks(g.N3, 256);
        CUDA_CHECK(cudaMemset(g.dPhi, 0, (size_t)g.N3*sizeof(float)));
        kLapse<<<gr,b>>>(g.dAlpha, g.dPhi, g.N3, (float)(1.0/(D.c*D.c)), 1e-8f);
        const int Ns = 100; integrateTau(g, Ns, (float)D.dt);
        std::vector<float> ha((size_t)g.N3), ht((size_t)g.N3);
        CUDA_CHECK(cudaMemcpy(ha.data(), g.dAlpha, (size_t)g.N3*sizeof(float), cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(ht.data(), g.dTau,   (size_t)g.N3*sizeof(float), cudaMemcpyDeviceToHost));
        double amax=0, tmax=0; double Texp=Ns*D.dt;
        for (int i=0;i<g.N3;i++){ amax=std::max(amax, std::fabs((double)ha[i]-1.0));
            tmax=std::max(tmax, std::fabs((double)ht[i]/Texp - 1.0)); }
        gridFree(g);
        bool ok1=(amax<1e-6), ok2=(tmax<1e-6);
        printf("[selftest] flat-lapse max|alpha-1|=%.3e  [%s]\n", amax, ok1?"PASS":"FAIL");
        printf("[selftest] proper-time  max|tau/(N*dt)-1|=%.3e  [%s]\n", tmax, ok2?"PASS":"FAIL");
        printf("VERDICT: %s\n", (ok1&&ok2)?"PASS":"FAIL");
        return (ok1&&ok2)?0:1;
    }

    Result R = dispatch(scenario, D, grid);
    std::string declared = declaredJson(D, seed, R);
    std::string hash = blake2b::hash256_hex(declared);

    if (golden){
        std::string path = "goldens/lapse_" + R.scenario + "/golden.hash";
        FILE* f = fopen(path.c_str(), "rb");
        if (!f){ fprintf(stderr, "GOLDEN NOT FROZEN %s\n", hash.c_str()); printf("%s\n", hash.c_str()); return 2; }
        char want[128]={0};
        if (fscanf(f, "%127s", want)!=1){ fclose(f); fprintf(stderr, "GOLDEN FILE UNREADABLE\n"); return 2; }
        fclose(f);
        if (hash == std::string(want)){ fprintf(stderr, "GOLDEN OK %.8s\n", hash.c_str()); printf("%s\n", hash.c_str()); return 0; }
        fprintf(stderr, "GOLDEN MISMATCH have=%.8s want=%.8s\n", hash.c_str(), want);
        printf("%s\n", hash.c_str());
        return 1;
    }
    if (json){
        std::string out = declared;
        out += ",\"notes\":\"non-declared; hash=" + hash.substr(0,8) + "\"}";
        printf("%s\n", out.c_str());
        return R.verdict ? 0 : 1;
    }
    printHuman(R);
    printf("declared hash: %s\n", hash.c_str());
    return R.verdict ? 0 : 1;
}
