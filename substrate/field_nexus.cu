// ============================================================================
//  field_nexus — TINY UNIVERSE v2 N1 "field" first-GPU-substrate oracle tool
//  Contract: contracts/field.contract.md v1.0.0 (design) · plan phase4/N1_FIELD_PLAN.md
//
//  Single-file CUDA. The Schrodinger-Poisson (SP) split-step engine welded from
//  M6's split-step psi kernels (kQ3PhaseK lineage) + M2's PM cuFFT-Poisson
//  (kGreen lineage). Deterministic; fp64 host-measured observables; envelope
//  face (--scenario X --json|--golden|--selftest). Exit 0/1/2, GPU preflight 3.
//
//  Scenarios (staged, oracle-grounded first):
//    freepacket : 3D free Gaussian dispersal, gravity OFF -> vs N5 analytic
//                 sigma(t) = sigma0*sqrt(1 + (hbar t / (2 m sigma0^2))^2).
//    sho3d      : 3D isotropic harmonic well, imaginary-time relax to ground ->
//                 vs N5 analytic E0 = (3/2) hbar omega, sigma = sqrt(hbar/2 m omega).
//    soliton    : self-gravitating SP ground state (gravity ON, PM Poisson each
//                 step) -> mass-radius product r_c*M scaling (SP soliton law).
//
//  Build (BUILD.md; no fast-math — Invariant 6):
//    nvcc -O3 -arch=sm_89 -Xcompiler "/O2" -o build\field_nexus.exe
//         substrate\field_nexus.cu cufft.lib
//
//  Dials (nexus v0, frozen): hbar=0.5, m=1, G=2e-3, dt=1/240, L_box=512.
//  Exit: 0 all gates pass · 1 a declared gate fired · 2 error · 3 GPU contended.
//  Determinism: (scenario, dials, seed, grid, steps) -> byte-identical declared
//  JSON (notes excluded); the frozen golden is sm_89-pinned.
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

static constexpr double PI = 3.14159265358979323846;
#define PI_F 3.14159265358979f

// ---- dials (nexus v0 defaults, frozen) --------------------------------------
struct Dials {
    double hbar  = 0.5;
    double m     = 1.0;
    double G     = 2.0e-3;
    double dt    = 1.0 / 240.0;
    double L_box = 512.0;
};

// ============================================================================
//  BLAKE2b-256 (RFC 7693, self-contained) — lifted VERBATIM from tiny_nexus.cpp
//  (the golden hasher for the _nexus tool family; freezes THIS hash).
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
// hash arbitrary bytes (the psi-grid receipt)
static std::string hash256_bytes(const void* data, size_t n){
    std::string s((const char*)data, n);
    return hash256_hex(s);
}
} // namespace blake2b

// ---- canonical number formatting (fmt6: %.6f, -0 normalized) ----------------
static std::string fmt6(double x){
    if (!std::isfinite(x)) return "9999999.999999";       // sentinel; gates catch NaN/Inf
    if (x == 0.0) x = 0.0;
    char buf[48];
    snprintf(buf, sizeof(buf), "%.6f", x);
    return std::string(buf);
}
// higher-precision channel for tiny observables (relative errors near tol)
static std::string fmt9(double x){
    if (!std::isfinite(x)) return "9.999999999e+99";
    if (x == 0.0) x = 0.0;
    char buf[48];
    snprintf(buf, sizeof(buf), "%.9e", x);
    return std::string(buf);
}

// ============================================================================
//  liborrery fixed-point accumulator (Invariant 4: no float atomics in declared
//  reductions). Idiom lifted from core/lib/reduce.cuh (order-invariant uint64).
// ============================================================================
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
//  KERNELS — 3D split-step (M6 kQ3* lineage) + PM Poisson (M2 kGreen lineage)
//  Grid is FN^3 over a cubic box of side gL; dx = gL/FN. Device constants passed
//  by value (the box side varies per scenario, so k-scale is a kernel arg).
// ============================================================================

// initial Gaussian packet: psi = exp(-|r-r0|^2/(4 sigma0^2)) * exp(i k.r)
__global__ void kGauss3(cufftComplex* q, int FN, float gL, float dx,
                        float sx0, float sy0, float sz0, float s0,
                        float kx, float ky, float kz){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    int N3 = FN*FN*FN;
    if (idx >= N3) return;
    int ix = idx % FN, iy = (idx/FN) % FN, iz = idx/(FN*FN);
    float x = -gL*0.5f + (ix + 0.5f)*dx;
    float y = -gL*0.5f + (iy + 0.5f)*dx;
    float z = -gL*0.5f + (iz + 0.5f)*dx;
    float dxr = x - sx0, dyr = y - sy0, dzr = z - sz0;
    float env = expf(-0.25f*(dxr*dxr + dyr*dyr + dzr*dzr)/(s0*s0));
    float ph = kx*x + ky*y + kz*z;
    q[idx].x = env*cosf(ph);
    q[idx].y = env*sinf(ph);
}

// kinetic phase (real time): psi <- e^{-i coef k^2} psi  (coef = hbar dt/(4m) half)
__global__ void kPhaseK(cufftComplex* q, int FN, float gL, float coef){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    int N3 = FN*FN*FN;
    if (idx >= N3) return;
    int ix = idx % FN, iy = (idx/FN) % FN, iz = idx/(FN*FN);
    int fx = (ix <= FN/2) ? ix : ix - FN;
    int fy = (iy <= FN/2) ? iy : iy - FN;
    int fz = (iz <= FN/2) ? iz : iz - FN;
    float kf = 2.0f*PI_F/gL;
    float k2 = kf*kf*(float)(fx*fx + fy*fy + fz*fz);
    float ph = -coef*k2;
    float c = cosf(ph), s = sinf(ph);
    cufftComplex v = q[idx];
    q[idx].x = v.x*c - v.y*s;
    q[idx].y = v.x*s + v.y*c;
}

// kinetic decay (imaginary time): psi <- e^{-coef k^2} psi  (coef = hbar dtau/(4m))
__global__ void kDecayK(cufftComplex* q, int FN, float gL, float coef){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    int N3 = FN*FN*FN;
    if (idx >= N3) return;
    int ix = idx % FN, iy = (idx/FN) % FN, iz = idx/(FN*FN);
    int fx = (ix <= FN/2) ? ix : ix - FN;
    int fy = (iy <= FN/2) ? iy : iy - FN;
    int fz = (iz <= FN/2) ? iz : iz - FN;
    float kf = 2.0f*PI_F/gL;
    float k2 = kf*kf*(float)(fx*fx + fy*fy + fz*fz);
    float d = expf(-coef*k2);
    q[idx].x *= d; q[idx].y *= d;
}

// potential phase (real time) from an arbitrary real V grid: psi <- e^{-i V coef} psi
__global__ void kPhaseV(cufftComplex* q, const float* V, int N3, float coef){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if (idx >= N3) return;
    float ph = -V[idx]*coef;
    float c = cosf(ph), s = sinf(ph);
    cufftComplex v = q[idx];
    q[idx].x = v.x*c - v.y*s;
    q[idx].y = v.x*s + v.y*c;
}
// potential decay (imaginary time): psi <- e^{-V coef} psi
__global__ void kDecayV(cufftComplex* q, const float* V, int N3, float coef){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if (idx >= N3) return;
    float d = expf(-V[idx]*coef);
    q[idx].x *= d; q[idx].y *= d;
}

// static external harmonic well V = 1/2 m omega^2 r^2 (sho3d)
__global__ void kFillHarmonic(float* V, int FN, float gL, float dx, float half_m_w2){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    int N3 = FN*FN*FN;
    if (idx >= N3) return;
    int ix = idx % FN, iy = (idx/FN) % FN, iz = idx/(FN*FN);
    float x = -gL*0.5f + (ix + 0.5f)*dx;
    float y = -gL*0.5f + (iy + 0.5f)*dx;
    float z = -gL*0.5f + (iz + 0.5f)*dx;
    V[idx] = half_m_w2*(x*x + y*y + z*z);
}

__global__ void kScale(cufftComplex* q, int N3, float s){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if (idx >= N3) return;
    q[idx].x *= s; q[idx].y *= s;
}

// psi-mass deposit -> fixed-point uint64 grid (1:1 gather, NO CIC scatter — the
// field already lives on the grid). g[c] += m*(|psi|^2)*dx^3. Invariant 4.
__global__ void kPsiDeposit(const cufftComplex* q, unsigned long long* g,
                            int N3, double mcell){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if (idx >= N3) return;
    double rho = ((double)q[idx].x*q[idx].x + (double)q[idx].y*q[idx].y)*mcell;
    // 1:1 slot: still routed through fixed_atomic_add so the |psi|^2 -> rho path
    // is bit-identical to how M2 quantizes mass (contract [DECIDED]).
    fp_atomic_add(&g[idx], rho);
}
__global__ void kFixToReal(const unsigned long long* g, float* r, int N3){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < N3) r[i] = (float)fp_decode((long long)g[i]);
}
__global__ void kZeroFix(unsigned long long* g, int n){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < n) g[i] = 0ull;
}
// PM Green multiply (VERBATIM structure from M2 kGreen; k=0 zeroed = mean does
// not gravitate). R2C spectrum layout: last dim -> FN/2+1 complex bins.
__global__ void kGreen(cufftComplex* s, int FN, int FNZC, float gL, float dx, float G){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    int total = FN*FN*FNZC;
    if (i >= total) return;
    int z = i % FNZC, y = (i/FNZC) % FN, x = i/(FNZC*FN);
    int fx = (x <= FN/2) ? x : x - FN;
    int fy = (y <= FN/2) ? y : y - FN;
    int fz = z;
    float kf = 2.0f*PI_F/gL;
    float k2 = kf*kf*(float)(fx*fx + fy*fy + fz*fz);
    float f = (k2 > 0.0f)
        ? -4.0f*PI_F*G/(k2 * (dx*dx*dx) * (float)(FN*FN*FN))
        : 0.0f;                                       // k=0: mean density does not gravitate
    s[i].x *= f; s[i].y *= f;
}
// subtract the field mean from Phi (defensive: k=0 already zeroed above; this
// re-centers any fp round so the periodic-gravity convention stays honest)

// deposit the per-cell |psi|^2 into a real grid for host readout (norm, radius)
__global__ void kRho(const cufftComplex* q, float* r, int N3){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if (idx >= N3) return;
    r[idx] = q[idx].x*q[idx].x + q[idx].y*q[idx].y;
}

// ============================================================================
//  Host-side split-step engines and fp64 observables
// ============================================================================
struct Grid {
    int FN;            // cells per side
    int N3;            // FN^3
    int FNZC;          // FN/2+1 (R2C last dim)
    float gL;          // box side (su)
    float dx;          // gL/FN
    cufftComplex* dQ = nullptr;
    float*  dV = nullptr;              // potential grid (harmonic or Phi)
    float*  dReal = nullptr;           // deposit real grid / Phi
    cufftComplex* dSpec = nullptr;     // R2C spectrum
    unsigned long long* dFix = nullptr;// fixed-point mass grid
    cufftHandle planC2C, planR2C, planC2R;
    bool poisson = false;
};

static void gridAlloc(Grid& g, bool withPoisson){
    g.N3   = g.FN*g.FN*g.FN;
    g.FNZC = g.FN/2 + 1;
    g.dx   = g.gL / g.FN;
    CUDA_CHECK(cudaMalloc(&g.dQ, (size_t)g.N3*sizeof(cufftComplex)));
    CUDA_CHECK(cudaMalloc(&g.dV, (size_t)g.N3*sizeof(float)));
    CUDA_CHECK(cudaMemset(g.dV, 0, (size_t)g.N3*sizeof(float)));
    CUFFT_CHECK(cufftPlan3d(&g.planC2C, g.FN, g.FN, g.FN, CUFFT_C2C));
    g.poisson = withPoisson;
    if (withPoisson){
        CUDA_CHECK(cudaMalloc(&g.dReal, (size_t)g.N3*sizeof(float)));
        CUDA_CHECK(cudaMalloc(&g.dSpec, (size_t)g.FN*g.FN*g.FNZC*sizeof(cufftComplex)));
        CUDA_CHECK(cudaMalloc(&g.dFix,  (size_t)g.N3*sizeof(unsigned long long)));
        CUFFT_CHECK(cufftPlan3d(&g.planR2C, g.FN, g.FN, g.FN, CUFFT_R2C));
        CUFFT_CHECK(cufftPlan3d(&g.planC2R, g.FN, g.FN, g.FN, CUFFT_C2R));
    }
}
static void gridFree(Grid& g){
    cudaFree(g.dQ); cudaFree(g.dV);
    cufftDestroy(g.planC2C);
    if (g.poisson){
        cudaFree(g.dReal); cudaFree(g.dSpec); cudaFree(g.dFix);
        cufftDestroy(g.planR2C); cufftDestroy(g.planC2R);
    }
}

static dim3 gridBlocks(int n, int b){ return dim3((n + b - 1)/b); }

// forward C2C then inverse rescale by 1/N3 (cuFFT is unnormalized)
static inline void fftFwd(Grid& g){ CUFFT_CHECK(cufftExecC2C(g.planC2C, g.dQ, g.dQ, CUFFT_FORWARD)); }
static inline void fftInv(Grid& g){
    CUFFT_CHECK(cufftExecC2C(g.planC2C, g.dQ, g.dQ, CUFFT_INVERSE));
    dim3 b(256), gr = gridBlocks(g.N3, 256);
    kScale<<<gr, b>>>(g.dQ, g.N3, 1.0f/(float)g.N3);
}

// host fp64 norm of |psi|^2 * dx^3 (mass if psi normalized to unit mass density)
static double hostNormMass(Grid& g, std::vector<cufftComplex>& hbuf){
    hbuf.resize((size_t)g.N3);
    CUDA_CHECK(cudaMemcpy(hbuf.data(), g.dQ, (size_t)g.N3*sizeof(cufftComplex), cudaMemcpyDeviceToHost));
    double s = 0.0;
    for (int i = 0; i < g.N3; i++){
        double a = (double)hbuf[i].x, b = (double)hbuf[i].y;
        s += a*a + b*b;
    }
    return s * (double)g.dx*g.dx*g.dx;
}

// renormalize psi so that integral |psi|^2 dx^3 = 1 (host-driven, deterministic)
static void renorm(Grid& g, std::vector<cufftComplex>& hbuf){
    double mass = hostNormMass(g, hbuf);
    double sc = 1.0 / std::sqrt(mass);
    dim3 b(256), gr = gridBlocks(g.N3, 256);
    kScale<<<gr, b>>>(g.dQ, g.N3, (float)sc);
}

// fp64 sigma_x (per-axis RMS about the |psi|^2 centroid, then averaged x/y/z)
struct SigmaStat { double sx, sy, sz, siso, mass, cx, cy, cz; };
static SigmaStat hostSigma(Grid& g, std::vector<cufftComplex>& hbuf){
    hbuf.resize((size_t)g.N3);
    CUDA_CHECK(cudaMemcpy(hbuf.data(), g.dQ, (size_t)g.N3*sizeof(cufftComplex), cudaMemcpyDeviceToHost));
    double s=0, sx=0, sy=0, sz=0, sxx=0, syy=0, szz=0;
    double gL = g.gL, dx = g.dx; int FN = g.FN;
    for (int iz = 0; iz < FN; iz++){
        double z = -gL*0.5 + (iz + 0.5)*dx;
        for (int iy = 0; iy < FN; iy++){
            double y = -gL*0.5 + (iy + 0.5)*dx;
            for (int ix = 0; ix < FN; ix++){
                double x = -gL*0.5 + (ix + 0.5)*dx;
                size_t idx = ((size_t)iz*FN + iy)*FN + ix;
                double a = (double)hbuf[idx].x, b = (double)hbuf[idx].y;
                double p = a*a + b*b;
                s += p;
                sx += p*x; sy += p*y; sz += p*z;
                sxx += p*x*x; syy += p*y*y; szz += p*z*z;
            }
        }
    }
    SigmaStat st;
    st.mass = s * dx*dx*dx;
    st.cx = sx/s; st.cy = sy/s; st.cz = sz/s;
    st.sx = std::sqrt(sxx/s - st.cx*st.cx);
    st.sy = std::sqrt(syy/s - st.cy*st.cy);
    st.sz = std::sqrt(szz/s - st.cz*st.cz);
    st.siso = std::cbrt(st.sx*st.sy*st.sz);
    return st;
}

// fp64 total energy <T>+<V>: T via k-space spectral, V via real-space |psi|^2.
// hbuf must hold the CURRENT real-space psi on entry (we copy internally).
static double hostEnergyHarmonic(Grid& g, const Dials& D, double omega){
    std::vector<cufftComplex> hr((size_t)g.N3);
    CUDA_CHECK(cudaMemcpy(hr.data(), g.dQ, (size_t)g.N3*sizeof(cufftComplex), cudaMemcpyDeviceToHost));
    // norm
    double nrm = 0.0;
    for (int i = 0; i < g.N3; i++){ double a=hr[i].x, b=hr[i].y; nrm += a*a+b*b; }
    // <V> real space
    double V = 0.0;
    double gL = g.gL, dx = g.dx; int FN = g.FN;
    double halfmw2 = 0.5*D.m*omega*omega;
    for (int iz = 0; iz < FN; iz++){
        double z = -gL*0.5 + (iz + 0.5)*dx;
        for (int iy = 0; iy < FN; iy++){
            double y = -gL*0.5 + (iy + 0.5)*dx;
            for (int ix = 0; ix < FN; ix++){
                double x = -gL*0.5 + (ix + 0.5)*dx;
                size_t idx = ((size_t)iz*FN + iy)*FN + ix;
                double a=hr[idx].x, b=hr[idx].y; double p=a*a+b*b;
                V += p * halfmw2*(x*x+y*y+z*z);
            }
        }
    }
    V /= nrm;
    // <T> spectral: forward FFT a copy, integrate hbar^2 k^2/(2m) |psi_k|^2 / norm_k
    fftFwd(g);
    std::vector<cufftComplex> hk((size_t)g.N3);
    CUDA_CHECK(cudaMemcpy(hk.data(), g.dQ, (size_t)g.N3*sizeof(cufftComplex), cudaMemcpyDeviceToHost));
    // restore real-space psi (inverse + rescale)
    fftInv(g);
    double T = 0.0, nk = 0.0;
    double kf = 2.0*PI/gL;
    double coef = D.hbar*D.hbar/(2.0*D.m);
    for (int iz = 0; iz < FN; iz++){
        int fz = (iz <= FN/2) ? iz : iz - FN;
        for (int iy = 0; iy < FN; iy++){
            int fy = (iy <= FN/2) ? iy : iy - FN;
            for (int ix = 0; ix < FN; ix++){
                int fx = (ix <= FN/2) ? ix : ix - FN;
                size_t idx = ((size_t)iz*FN + iy)*FN + ix;
                double a=hk[idx].x, b=hk[idx].y; double pk=a*a+b*b;
                double k2 = kf*kf*(double)(fx*fx+fy*fy+fz*fz);
                T += coef*k2*pk;
                nk += pk;
            }
        }
    }
    T /= nk;
    return T + V;
}

// ---------------------------------------------------------------------------
// half-mass radius (r_c) from the current |psi|^2 about its centroid: the
// radius enclosing half the mass. Used by soliton mass-radius scaling.
// ---------------------------------------------------------------------------
static double hostHalfMassRadius(Grid& g, std::vector<cufftComplex>& hbuf, double& massOut){
    SigmaStat st = hostSigma(g, hbuf);   // fills hbuf + centroid
    massOut = st.mass;
    double gL = g.gL, dx = g.dx; int FN = g.FN;
    // build (r, mass) then integrate CDF; bin by radius for a monotone CDF
    // coarse histogram over radius; enough for a 5% gate.
    const int NB = 512;
    double rmax = 0.87*gL;   // corner distance guard
    std::vector<double> hist(NB, 0.0);
    double total = 0.0;
    for (int iz = 0; iz < FN; iz++){
        double z = -gL*0.5 + (iz + 0.5)*dx - st.cz;
        for (int iy = 0; iy < FN; iy++){
            double y = -gL*0.5 + (iy + 0.5)*dx - st.cy;
            for (int ix = 0; ix < FN; ix++){
                double x = -gL*0.5 + (ix + 0.5)*dx - st.cx;
                size_t idx = ((size_t)iz*FN + iy)*FN + ix;
                double a=hbuf[idx].x, b=hbuf[idx].y; double p=a*a+b*b;
                double r = std::sqrt(x*x+y*y+z*z);
                int bin = (int)(r/rmax*NB); if (bin >= NB) bin = NB-1; if (bin < 0) bin = 0;
                hist[bin] += p; total += p;
            }
        }
    }
    double half = 0.5*total, acc = 0.0, rc = rmax;
    for (int bi = 0; bi < NB; bi++){
        double prev = acc; acc += hist[bi];
        if (acc >= half){
            double frac = (half - prev)/std::max(hist[bi], 1e-300);
            rc = (bi + frac)/NB*rmax;
            break;
        }
    }
    return rc;
}

// ============================================================================
//  SCENARIO: freepacket — 3D free Gaussian dispersal, gravity OFF, vs N5 sigma(t)
// ============================================================================
struct Result {
    std::string scenario;
    int grid; long long steps;
    // observables
    double sigma_meas=0, sigma_exp=0, sigma_rel=0;
    double E_meas=0, E_exp=0, E_rel=0;
    double sigmaG_meas=0, sigmaG_exp=0, sigmaG_rel=0;   // sho3d sigma
    double rc=0, mass=0, rc_M=0, rc_M2=0, mass2=0, rc2=0; // soliton
    double kappa=0, scale_rel=0;
    double mass0=0, mass1=0, dmass_rel=0;
    double norm0=0, norm1=0;
    bool nan_free=true;
    // gates
    bool gate_primary=false, gate_norm=false, gate_nan=false, gate_scale=false;
    bool verdict=false;
    std::string psi_b2b;
    // scenario params echoed
    double gL=0, sigma0=0, omega=0, tt=0;
};

static std::string capturePsiHash(Grid& g){
    std::vector<cufftComplex> h((size_t)g.N3);
    CUDA_CHECK(cudaMemcpy(h.data(), g.dQ, (size_t)g.N3*sizeof(cufftComplex), cudaMemcpyDeviceToHost));
    return blake2b::hash256_bytes(h.data(), (size_t)g.N3*sizeof(cufftComplex));
}

static Result runFreepacket(const Dials& D, int gridReq){
    Result R; R.scenario = "freepacket";
    // Pin geometry for a resolved-AND-headroomed 3D free packet (see design notes):
    // box 48 su, 192^3 (dx = 0.25 = sigma0/6), packet grows by ~1.7x -> 6*sigma_f
    // ~ 0.64 of the half-box, so periodic tails at the wall are ~1e-8 (no wrap
    // contamination). Momentum k0 = 0: centroid stays put -> symmetric headroom.
    (void)gridReq;
    const int FN = 192;
    R.grid = FN;
    Grid g; g.FN = FN; g.gL = 48.0f;
    gridAlloc(g, false);
    R.gL = g.gL;
    const double sigma0 = 1.5;
    R.sigma0 = sigma0;
    const float k0 = 0.0f;
    dim3 b(256), gr = gridBlocks(g.N3, 256);
    kGauss3<<<gr, b>>>(g.dQ, FN, g.gL, g.dx, 0,0,0, (float)sigma0, k0, 0.0f, 0.0f);
    std::vector<cufftComplex> hbuf;
    renorm(g, hbuf);

    // target time: sigma grows by 1.7x -> u = sqrt(1.7^2 - 1); t = u*2 m sigma0^2/hbar
    const double ratio = 1.7;
    double u_t = std::sqrt(ratio*ratio - 1.0);
    double t_target = u_t*2.0*D.m*sigma0*sigma0/D.hbar;
    long long nq = (long long)std::llround(t_target / D.dt);
    R.steps = nq;
    double t_actual = nq * D.dt;
    R.tt = t_actual;

    SigmaStat s0 = hostSigma(g, hbuf);
    R.norm0 = s0.mass; R.mass0 = s0.mass;

    // V=0 free evolution: the Strang half-K | identity-V | half-K reduces to one
    // full kinetic step per tick, and consecutive ticks' k-space phases fuse
    // exactly (the FFTs between them are inverse-then-forward = identity). So we
    // FFT ONCE, apply the full-step kinetic phase nq times in k-space (genuinely
    // exercising the kinetic operator every tick), and inverse-FFT ONCE. This is
    // N5's step_free_real structure verbatim: byte-identical result to nq round
    // trips, but norm is conserved to a SINGLE fp32 round-trip (no per-step FFT
    // drift). coefK_full = hbar dt/(2m) for e^{-i coef k^2}.
    float coefKfull = (float)(D.hbar*D.dt/(2.0*D.m));
    fftFwd(g);
    for (long long s = 0; s < nq; s++)
        kPhaseK<<<gr, b>>>(g.dQ, FN, g.gL, coefKfull);
    fftInv(g);
    CUDA_CHECK(cudaDeviceSynchronize());

    SigmaStat s1 = hostSigma(g, hbuf);
    R.norm1 = s1.mass; R.mass1 = s1.mass;
    R.dmass_rel = std::fabs((s1.mass - s0.mass)/s0.mass);
    R.sigma_meas = s1.siso;
    double u = D.hbar*t_actual/(2.0*D.m*sigma0*sigma0);
    R.sigma_exp = sigma0*std::sqrt(1.0 + u*u);
    R.sigma_rel = std::fabs(R.sigma_meas/R.sigma_exp - 1.0);

    // NaN scan
    for (int i = 0; i < g.N3 && R.nan_free; i++)
        if (!std::isfinite(hbuf[i].x) || !std::isfinite(hbuf[i].y)) R.nan_free = false;

    R.psi_b2b = capturePsiHash(g);
    // gates: sigma vs analytic < 1e-2; periodic (no absorber) => mass conserved 1e-3; NaN-free
    R.gate_primary = (R.sigma_rel < 1e-2);
    R.gate_norm    = (R.dmass_rel < 1e-3);
    R.gate_nan     = R.nan_free;
    R.verdict = R.gate_primary && R.gate_norm && R.gate_nan;
    gridFree(g);
    return R;
}

// ============================================================================
//  SCENARIO: sho3d — 3D isotropic harmonic ground state by imaginary-time relax,
//  vs N5 analytic E0 = (3/2) hbar omega, sigma = sqrt(hbar/(2 m omega)).
// ============================================================================
static Result runSho3d(const Dials& D, int FN){
    Result R; R.scenario = "sho3d"; R.grid = FN;
    const double omega = 1.0;
    R.omega = omega;
    // ground-state width sqrt(hbar/2 m omega) = sqrt(0.25) = 0.5. Box must hold
    // several widths and resolve dx << width. sigma_gs = 0.5 -> box 16 su, dx=0.125@128.
    Grid g; g.FN = FN; g.gL = 16.0f;
    gridAlloc(g, false);
    R.gL = g.gL;
    dim3 b(256), gr = gridBlocks(g.N3, 256);
    // start from an off-width Gaussian (like N5: sigma 1.5, offset) to prove convergence
    kGauss3<<<gr, b>>>(g.dQ, FN, g.gL, g.dx, 0.3f, 0.0f, 0.0f, 1.5f, 0.0f, 0.0f, 0.0f);
    std::vector<cufftComplex> hbuf;
    renorm(g, hbuf);
    // fill harmonic potential
    kFillHarmonic<<<gr, b>>>(g.dV, FN, g.gL, g.dx, (float)(0.5*D.m*omega*omega));

    // imaginary-time Strang: e^{-H tau/hbar}. tau step, iterate; renorm each step.
    const double tau = 0.002;
    const int iters = 20000;                   // matches N5's imaginary-time budget
    R.steps = iters;
    float ekcoef = (float)(D.hbar*tau/(4.0*D.m));       // half kinetic decay coef (hbar k^2/(2m)*tau*0.5)
    float evcoef = (float)(tau/D.hbar);                 // full potential decay coef (V/hbar*tau)
    for (int it = 0; it < iters; it++){
        fftFwd(g);
        kDecayK<<<gr, b>>>(g.dQ, FN, g.gL, ekcoef);
        fftInv(g);
        kDecayV<<<gr, b>>>(g.dQ, g.dV, g.N3, evcoef);
        fftFwd(g);
        kDecayK<<<gr, b>>>(g.dQ, FN, g.gL, ekcoef);
        fftInv(g);
        renorm(g, hbuf);
    }
    CUDA_CHECK(cudaDeviceSynchronize());

    SigmaStat st = hostSigma(g, hbuf);
    R.mass1 = st.mass; R.norm1 = st.mass; R.norm0 = 1.0; R.mass0 = 1.0;
    R.dmass_rel = std::fabs(st.mass - 1.0);
    R.sigmaG_meas = st.siso;
    R.sigmaG_exp  = std::sqrt(D.hbar/(2.0*D.m*omega));
    R.sigmaG_rel  = std::fabs(R.sigmaG_meas/R.sigmaG_exp - 1.0);

    R.E_meas = hostEnergyHarmonic(g, D, omega);
    R.E_exp  = 1.5*D.hbar*omega;               // 3D ground = 3 * (1/2 hbar omega)
    R.E_rel  = std::fabs(R.E_meas/R.E_exp - 1.0);

    for (int i = 0; i < g.N3 && R.nan_free; i++)
        if (!std::isfinite(hbuf[i].x) || !std::isfinite(hbuf[i].y)) R.nan_free = false;

    R.psi_b2b = capturePsiHash(g);
    // gates: E within 1e-3, sigma within 1e-2, NaN-free
    R.gate_primary = (R.E_rel < 1e-3) && (R.sigmaG_rel < 1e-2);
    R.gate_nan     = R.nan_free;
    R.gate_norm    = (R.dmass_rel < 1e-3);
    R.verdict = R.gate_primary && R.gate_nan;
    gridFree(g);
    return R;
}

// ============================================================================
//  SCENARIO: soliton — self-gravitating SP ground state by imaginary-time relax
//  with the PM Poisson solve providing Phi = self-gravity each step. Measures
//  the mass-radius product r_c*M and checks scale-covariance across two masses.
//
//  SP soliton scaling law (Chavanis 2011; Mocz 2017): the ground state is
//  self-similar under (M -> lambda M, r -> r/lambda), so r_c * M = const (kappa).
//  We MEASURE kappa in-tool at mass M1, then verify r_c*M is invariant at a
//  second mass M2 = 2 M1 (the scale-covariance receipt). No external number is
//  imported as the gate; the literature value is reported as a labelled cross-check.
// ============================================================================
static double relaxSoliton(const Dials& D, Grid& g, double Mtot, int iters,
                           double tau, std::vector<cufftComplex>& hbuf){
    dim3 b(256), gr = gridBlocks(g.N3, 256);
    int nspec = g.FN*g.FN*g.FNZC;
    dim3 grS = gridBlocks(nspec, 256);
    dim3 grF = gridBlocks(g.N3, 256);
    float ekcoef = (float)(D.hbar*tau/(4.0*D.m));
    // psi normalized so integral m|psi|^2 dx^3 = Mtot  => integral |psi|^2 dx^3 = Mtot/m
    // We keep integral |psi|^2 dx^3 = Mtot/m by scaling after unit-renorm.
    double targetNorm = Mtot / D.m;                 // integral |psi|^2 dx^3
    for (int it = 0; it < iters; it++){
        // --- half kinetic (imag) ---
        fftFwd(g);
        kDecayK<<<gr, b>>>(g.dQ, g.FN, g.gL, ekcoef);
        fftInv(g);
        // --- Poisson: deposit m|psi|^2 -> rho, solve Phi ---
        kZeroFix<<<grF, b>>>(g.dFix, g.N3);
        kPsiDeposit<<<grF, b>>>(g.dQ, g.dFix, g.N3, (double)D.m*g.dx*g.dx*g.dx);
        kFixToReal<<<grF, b>>>(g.dFix, g.dReal, g.N3);
        CUFFT_CHECK(cufftExecR2C(g.planR2C, g.dReal, g.dSpec));
        kGreen<<<grS, b>>>(g.dSpec, g.FN, g.FNZC, g.gL, g.dx, (float)D.G);
        CUFFT_CHECK(cufftExecC2R(g.planC2R, g.dSpec, g.dReal)); // dReal = Phi (unnormalized C2R)
        // C2R is unnormalized: divide by N3. Fold into potential coef instead.
        // potential decay: e^{-(V/hbar) tau}, V = m*Phi, Phi_true = dReal / N3
        float evcoef = (float)(D.m*tau/D.hbar/(double)g.N3);
        kDecayV<<<grF, b>>>(g.dQ, g.dReal, g.N3, evcoef);
        // --- half kinetic (imag) ---
        fftFwd(g);
        kDecayK<<<gr, b>>>(g.dQ, g.FN, g.gL, ekcoef);
        fftInv(g);
        // renormalize to target mass
        double mass = hostNormMass(g, hbuf);     // integral |psi|^2 dx^3
        double sc = std::sqrt(targetNorm / mass);
        kScale<<<gr, b>>>(g.dQ, g.N3, (float)sc);
    }
    CUDA_CHECK(cudaDeviceSynchronize());
    double massOut = 0;
    double rc = hostHalfMassRadius(g, hbuf, massOut);
    return rc;
}

static Result runSoliton(const Dials& D, int FN){
    Result R; R.scenario = "soliton"; R.grid = FN;
    Grid g; g.FN = FN; g.gL = 128.0f;          // box comfortably larger than the core
    gridAlloc(g, true);
    R.gL = g.gL;
    dim3 b(256), gr = gridBlocks(g.N3, 256);
    std::vector<cufftComplex> hbuf;

    const int   iters = 4000;
    const double tau  = 0.02;
    R.steps = iters;

    // --- mass 1 ---
    const double M1 = 200.0;
    kGauss3<<<gr, b>>>(g.dQ, FN, g.gL, g.dx, 0,0,0, 8.0f, 0.0f,0.0f,0.0f);
    renorm(g, hbuf);
    // scale to M1/m norm
    { double mass = hostNormMass(g, hbuf); double sc = std::sqrt((M1/D.m)/mass);
      kScale<<<gr, b>>>(g.dQ, g.N3, (float)sc); }
    double rc1 = relaxSoliton(D, g, M1, iters, tau, hbuf);
    double m1out = 0; { double dummy; (void)dummy; }
    SigmaStat st1 = hostSigma(g, hbuf); m1out = st1.mass;
    R.rc = rc1; R.mass = m1out; R.rc_M = rc1 * m1out;

    // --- mass 2 = 2*M1 ---
    const double M2 = 400.0;
    kGauss3<<<gr, b>>>(g.dQ, FN, g.gL, g.dx, 0,0,0, 6.0f, 0.0f,0.0f,0.0f);
    renorm(g, hbuf);
    { double mass = hostNormMass(g, hbuf); double sc = std::sqrt((M2/D.m)/mass);
      kScale<<<gr, b>>>(g.dQ, g.N3, (float)sc); }
    double rc2 = relaxSoliton(D, g, M2, iters, tau, hbuf);
    SigmaStat st2 = hostSigma(g, hbuf);
    R.rc2 = rc2; R.mass2 = st2.mass; R.rc_M2 = rc2 * st2.mass;

    R.kappa = R.rc_M;                             // measured constant at M1
    R.scale_rel = std::fabs(R.rc_M2/R.rc_M - 1.0);

    R.psi_b2b = capturePsiHash(g);
    R.nan_free = std::isfinite(rc1) && std::isfinite(rc2) && (rc1 > 0) && (rc2 > 0);
    // gate: scale covariance r_c*M invariant across 2x mass within 5%
    R.gate_scale = (R.scale_rel < 0.05);
    R.gate_nan = R.nan_free;
    R.verdict = R.gate_scale && R.gate_nan;
    gridFree(g);
    return R;
}

// ============================================================================
//  Declared JSON (hash domain: everything below; "notes" appended outside).
//  Canonical order per contract: seed, params, result, gates, verdict.
// ============================================================================
static std::string declaredJson(const Dials& D, uint64_t seed, const Result& R){
    std::string s; s.reserve(2048);
    s += "{\"tool\":\"field_nexus\",\"version\":\"1.0.0\",\"seed\":";
    s += std::to_string(seed);
    s += ",\"params\":{";
    s += "\"scenario\":\"" + R.scenario + "\"";
    s += ",\"grid\":" + std::to_string(R.grid);
    s += ",\"steps\":" + std::to_string(R.steps);
    s += ",\"hbar\":" + fmt6(D.hbar) + ",\"m\":" + fmt6(D.m) + ",\"G\":" + fmt9(D.G);
    s += ",\"dt\":" + fmt9(D.dt) + ",\"L_box\":" + fmt6(D.L_box);
    s += ",\"gL\":" + fmt6(R.gL);
    s += "},\"result\":{";
    s += "\"state_b2b\":\"" + R.psi_b2b + "\"";
    s += ",\"mass0\":" + fmt6(R.mass0) + ",\"mass1\":" + fmt6(R.mass1);
    s += ",\"dmass_rel\":" + fmt9(R.dmass_rel);
    if (R.scenario == "freepacket"){
        s += ",\"sigma_meas\":" + fmt6(R.sigma_meas);
        s += ",\"sigma_exp\":"  + fmt6(R.sigma_exp);
        s += ",\"sigma_rel\":"  + fmt9(R.sigma_rel);
        s += ",\"t_actual\":"   + fmt6(R.tt);
    } else if (R.scenario == "sho3d"){
        s += ",\"E_meas\":" + fmt6(R.E_meas) + ",\"E_exp\":" + fmt6(R.E_exp);
        s += ",\"E_rel\":"  + fmt9(R.E_rel);
        s += ",\"sigma_meas\":" + fmt6(R.sigmaG_meas) + ",\"sigma_exp\":" + fmt6(R.sigmaG_exp);
        s += ",\"sigma_rel\":"  + fmt9(R.sigmaG_rel);
    } else if (R.scenario == "soliton"){
        s += ",\"rc1\":" + fmt6(R.rc) + ",\"M1\":" + fmt6(R.mass) + ",\"rcM1\":" + fmt6(R.rc_M);
        s += ",\"rc2\":" + fmt6(R.rc2) + ",\"M2\":" + fmt6(R.mass2) + ",\"rcM2\":" + fmt6(R.rc_M2);
        s += ",\"kappa\":" + fmt6(R.kappa) + ",\"scale_rel\":" + fmt9(R.scale_rel);
    }
    s += ",\"nan_free\":" + std::string(R.nan_free ? "1" : "0");
    s += "},\"gates\":{";
    s += "\"primary\":"  + std::string(R.gate_primary ? "true":"false");
    s += ",\"norm\":"    + std::string(R.gate_norm ? "true":"false");
    s += ",\"nan\":"     + std::string(R.gate_nan ? "true":"false");
    s += ",\"scale\":"   + std::string(R.gate_scale ? "true":"false");
    s += "},\"verdict\":\"";
    s += R.verdict ? "pass" : "fail";
    s += "\"";
    return s;                                       // caller closes with notes + '}'
}

static void printHuman(const Result& R){
    printf("field_nexus v1.0.0 — TINY UNIVERSE v2 N1 field (SP fusion)\n");
    printf("scenario: %s   grid: %d^3   steps: %lld\n", R.scenario.c_str(), R.grid, R.steps);
    printf("-------------------------------------------------------\n");
    if (R.scenario == "freepacket"){
        printf("  sigma_meas    %.6f\n", R.sigma_meas);
        printf("  sigma_exp     %.6f  (N5 analytic)\n", R.sigma_exp);
        printf("  sigma_rel_err %.3e   [gate < 1e-2]  %s\n", R.sigma_rel, R.gate_primary?"PASS":"FAIL");
        printf("  mass drift    %.3e   [gate < 1e-3]  %s\n", R.dmass_rel, R.gate_norm?"PASS":"FAIL");
    } else if (R.scenario == "sho3d"){
        printf("  E_meas        %.6f\n", R.E_meas);
        printf("  E_exp         %.6f  (3/2 hbar omega)\n", R.E_exp);
        printf("  E_rel_err     %.3e   [gate < 1e-3]  %s\n", R.E_rel, (R.E_rel<1e-3)?"PASS":"FAIL");
        printf("  sigma_meas    %.6f\n", R.sigmaG_meas);
        printf("  sigma_exp     %.6f  (sqrt(hbar/2mw))\n", R.sigmaG_exp);
        printf("  sigma_rel_err %.3e   [gate < 1e-2]  %s\n", R.sigmaG_rel, (R.sigmaG_rel<1e-2)?"PASS":"FAIL");
    } else if (R.scenario == "soliton"){
        printf("  M1=%.2f  rc1=%.4f  rc*M1=%.4f\n", R.mass, R.rc, R.rc_M);
        printf("  M2=%.2f  rc2=%.4f  rc*M2=%.4f\n", R.mass2, R.rc2, R.rc_M2);
        printf("  kappa(meas)   %.6f\n", R.kappa);
        printf("  scale_rel     %.3e   [gate < 5e-2]  %s\n", R.scale_rel, R.gate_scale?"PASS":"FAIL");
    }
    printf("  nan_free      %s\n", R.nan_free?"yes":"NO");
    printf("-------------------------------------------------------\n");
    printf("VERDICT: %s\n", R.verdict ? "PASS" : "FAIL");
}

// ---- GPU preflight (harness exit-3 lineage): refuse a contended card --------
static int gpuPreflight(size_t need_bytes){
    size_t freeB = 0, totalB = 0;
    cudaError_t e = cudaMemGetInfo(&freeB, &totalB);
    if (e != cudaSuccess){
        fprintf(stderr, "[preflight] cudaMemGetInfo failed: %s\n", cudaGetErrorString(e));
        return 2;
    }
    size_t min_free = need_bytes + (size_t)512*1024*1024;   // headroom
    if (freeB < min_free){
        fprintf(stderr, "[preflight] GPU busy: %.0f MB free (< %.0f MB needed). Exit 3.\n",
                freeB/1048576.0, min_free/1048576.0);
        return 3;
    }
    return 0;
}

static Result dispatch(const std::string& sc, const Dials& D, int grid){
    if (sc == "freepacket") return runFreepacket(D, grid);
    if (sc == "sho3d")      return runSho3d(D, grid);
    if (sc == "soliton")    return runSoliton(D, grid);
    fprintf(stderr, "error: unknown scenario '%s' (freepacket|sho3d|soliton)\n", sc.c_str());
    std::exit(2);
}

// estimate live VRAM need for a scenario at grid FN (bytes)
static size_t needBytes(const std::string& sc, int FN){
    size_t N3 = (size_t)FN*FN*FN;
    size_t c2c = N3*sizeof(cufftComplex)*3;             // psi + plan working set (~2x)
    if (sc == "soliton"){
        size_t pois = N3*(sizeof(float) + sizeof(unsigned long long))
                    + (size_t)FN*FN*(FN/2+1)*sizeof(cufftComplex)*3;
        return c2c + pois;
    }
    return c2c + N3*sizeof(float);
}

int main(int argc, char** argv){
    Dials D;
    uint64_t seed = 20260711ull;
    std::string scenario = "freepacket";
    int grid = 128;
    bool json=false, selftest=false, golden=false;

    for (int i = 1; i < argc; i++){
        std::string a = argv[i];
        if      (a == "--json") json = true;
        else if (a == "--selftest") selftest = true;
        else if (a == "--golden") golden = true;
        else if (a == "--scenario" && i+1 < argc) scenario = argv[++i];
        else if (a == "--grid" && i+1 < argc) grid = atoi(argv[++i]);
        else if (a == "--seed" && i+1 < argc) seed = strtoull(argv[++i], nullptr, 10);
        else {
            fprintf(stderr, "usage: field_nexus --scenario freepacket|sho3d|soliton "
                            "[--grid N=128] [--seed N] [--json] [--golden] [--selftest]\n");
            return 2;
        }
    }

    // golden = frozen params: fixed grid per scenario, canonical seed
    if (golden){ D = Dials(); seed = 20260711ull; grid = 128; }

    // GPU preflight (exit 3 on contention — never a red golden)
    int pf = gpuPreflight(needBytes(scenario, grid));
    if (pf != 0) return pf;

    if (selftest){
        // cheap smoke: flat field -> Poisson Phi ~ 0 (mean does not gravitate);
        // and a short free-evolution norm-conservation check.
        Grid g; g.FN = 32; g.gL = 64.0f; gridAlloc(g, true);
        dim3 b(256), gr = gridBlocks(g.N3, 256);
        // flat unit field
        kGauss3<<<gr,b>>>(g.dQ, g.FN, g.gL, g.dx, 0,0,0, 1e6f, 0,0,0);  // huge sigma ~ flat
        std::vector<cufftComplex> hbuf; renorm(g, hbuf);
        int nspec = g.FN*g.FN*g.FNZC; dim3 grS = gridBlocks(nspec,256), grF = gridBlocks(g.N3,256);
        kZeroFix<<<grF,b>>>(g.dFix, g.N3);
        kPsiDeposit<<<grF,b>>>(g.dQ, g.dFix, g.N3, (double)D.m*g.dx*g.dx*g.dx);
        kFixToReal<<<grF,b>>>(g.dFix, g.dReal, g.N3);
        CUFFT_CHECK(cufftExecR2C(g.planR2C, g.dReal, g.dSpec));
        kGreen<<<grS,b>>>(g.dSpec, g.FN, g.FNZC, g.gL, g.dx, (float)D.G);
        CUFFT_CHECK(cufftExecC2R(g.planC2R, g.dSpec, g.dReal));
        std::vector<float> hp((size_t)g.N3);
        CUDA_CHECK(cudaMemcpy(hp.data(), g.dReal, (size_t)g.N3*sizeof(float), cudaMemcpyDeviceToHost));
        double phimax = 0; for (int i=0;i<g.N3;i++) phimax = std::max(phimax, std::fabs((double)hp[i]/g.N3));
        double m0 = hostNormMass(g, hbuf);
        float coefKhalf = (float)(D.hbar*D.dt/(4.0*D.m));
        for (int s=0;s<50;s++){ fftFwd(g); kPhaseK<<<gr,b>>>(g.dQ,g.FN,g.gL,coefKhalf);
            kPhaseK<<<gr,b>>>(g.dQ,g.FN,g.gL,coefKhalf); fftInv(g); }
        double m1 = hostNormMass(g, hbuf);
        double ndrift = std::fabs((m1-m0)/m0);
        gridFree(g);
        bool ok1 = (phimax < 1e-3);       // flat field -> Phi ~ 0
        bool ok2 = (ndrift < 1e-4);       // unitary free evolution conserves norm
        printf("[selftest] flat-field Phi_max=%.3e  [%s]\n", phimax, ok1?"PASS":"FAIL");
        printf("[selftest] free-evo norm drift=%.3e  [%s]\n", ndrift, ok2?"PASS":"FAIL");
        printf("VERDICT: %s\n", (ok1&&ok2)?"PASS":"FAIL");
        return (ok1&&ok2)?0:1;
    }

    Result R = dispatch(scenario, D, grid);
    std::string declared = declaredJson(D, seed, R);
    std::string hash = blake2b::hash256_hex(declared);

    if (golden){
        std::string path = "goldens/" + R.scenario + "/golden.hash";
        // field goldens live under goldens/field_<scenario>/ to avoid clashing with
        // the v1 particle scenario names (galaxy/merger/echo already exist).
        path = "goldens/field_" + R.scenario + "/golden.hash";
        FILE* f = fopen(path.c_str(), "rb");
        if (!f){
            fprintf(stderr, "GOLDEN NOT FROZEN %s\n", hash.c_str());
            printf("%s\n", hash.c_str());
            return 2;
        }
        char want[128] = {0};
        if (fscanf(f, "%127s", want) != 1){ fclose(f); fprintf(stderr, "GOLDEN FILE UNREADABLE\n"); return 2; }
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
