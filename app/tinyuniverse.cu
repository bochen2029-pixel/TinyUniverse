// ============================================================================
//  TINY UNIVERSE — app v0.1.0 "first light" (M1 canvas)
//
//  1M inert drifting particles (galaxy scenario), rendered entirely in CUDA:
//  HDR float4 accumulation -> mip-chain threshold-free bloom -> log-luminance
//  auto-exposure -> AgX (default) | ACES-Narkowicz (GARGANTUA parity) -> grade
//  -> triangular dither -> 8-bit. Presentation = thin GL blit (D-012 P0).
//  Two CUDA streams: sim (low prio, fixed dt ticks, ping-pong publish) and
//  present (high prio) — presentation never waits beyond the published tick.
//
//  Contracts: contracts/frame.contract.md v1.0.0 · CINEMATIC.md gate
//  liborrery consumer: core/lib/rng.cuh (counter-keyed init; D-005 lift)
//
//  Build (BUILD.md; no fast-math — Invariant 6):
//    nvcc -O3 -arch=sm_89 -Xcompiler "/O2" -o build\tinyuniverse.exe
//         app\tinyuniverse.cu user32.lib gdi32.lib opengl32.lib
//
//  Run:  tinyuniverse.exe [--w 1920 --h 1080] [--n 1000000] [--seed N]
//                         [--ssaa] [--bench FRAMES] [--shot PATH.bmp [--frames N]]
//  Keys: L-drag orbit · wheel zoom · 1-3 presets · O auto-orbit · P pause
//        C cinematic/physical · T tonemap AgX/ACES · B bloom · A auto-exposure
//        +/- exposure EV · F1 HUD · V vsync · Esc quit
// ============================================================================

#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#include <windowsx.h>
#include <GL/gl.h>
#include <cuda_runtime.h>
#include <cuda_gl_interop.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <vector>
#include <string>
#include <cufft.h>
#include "../core/lib/rng.cuh"
#include "../core/lib/reduce.cuh"
#include "../core/lib/envelope.h"

#define PI_F 3.14159265358979f

// --- M2 newton: dials + PM constants (contracts/newton.contract.md) ---
#define G_DIAL   2.0e-3f
#define PMN      128
#define PMNZC    (PMN/2 + 1)
#define CELL     4.0f
#define BOXL     512.0f
#define REL_V2   1.0f            // (0.05c)^2 = 1 su^2/s^2 — regime threshold
enum Scenario { SC_GALAXY = 0, SC_KEPLER, SC_THREEBODY, SC_CLOUD };
enum Solver   { SOLV_PM = 0, SOLV_DIRECT, SOLV_TINY };
static const char* SC_NAME[4] = {"galaxy", "kepler", "threebody", "cloud"};

#define CUDA_CHECK(call) do {                                                  \
    cudaError_t _e = (call);                                                   \
    if (_e != cudaSuccess) {                                                   \
        fprintf(stderr, "CUDA error %s (%d) at %s:%d\n  %s\n",                 \
                cudaGetErrorName(_e), (int)_e, __FILE__, __LINE__, #call);     \
        exit(2);                                                               \
    }                                                                          \
} while (0)

// ----------------------------------------------------------------------------
// Params (by-value to kernels)
// ----------------------------------------------------------------------------
struct RP {
    int   RW, RH;          // render (internal) resolution
    int   OW, OH;          // output resolution
    int   ss;              // supersample factor (1 or 2)
    float cpos[3], fwd[3], rgt[3], upv[3];
    float tf, aspect;
    float exposure, bloomAmt, satCine;
    int   mode;            // 0 cinematic, 1 physical
    int   tonemap;         // 0 AgX, 1 ACES
    int   bloomOn;
    unsigned frame;
};

// ----------------------------------------------------------------------------
// Color: blackbody (Tanner-Helland fit — lifted from GARGANTUA), AgX, ACES
// ----------------------------------------------------------------------------
__device__ __forceinline__ float clampf(float x, float a, float b){
    return x < a ? a : (x > b ? b : x);
}
__device__ void blackbody(float T, float* rgb){
    float t = clampf(T, 1000.0f, 40000.0f) / 100.0f;
    float r, g, b;
    if (t <= 66){ r = 255; g = 99.4708025861f * logf(t) - 161.1195681661f; }
    else { r = 329.698727446f  * powf(t - 60, -0.1332047592f);
           g = 288.1221695283f * powf(t - 60, -0.0755148492f); }
    if (t >= 66) b = 255;
    else if (t <= 19) b = 0;
    else b = 138.5177312231f * logf(t - 10) - 305.0447927307f;
    rgb[0] = powf(clampf(r, 0, 255) / 255.0f, 2.2f);
    rgb[1] = powf(clampf(g, 0, 255) / 255.0f, 2.2f);
    rgb[2] = powf(clampf(b, 0, 255) / 255.0f, 2.2f);
}
__device__ __forceinline__ float aces(float x){
    return clampf(x*(2.51f*x + 0.03f)/(x*(2.43f*x + 0.59f) + 0.14f), 0.0f, 1.0f);
}
// AgX minimal (Sobotka/Wrensch): inset matrix -> log2 encode -> 6th-order
// contrast approx -> outset matrix. Output is display-referred (sRGB-encoded).
__device__ __forceinline__ float agx_c(float x){
    float x2 = x*x, x4 = x2*x2;
    return 15.5f*x4*x2 - 40.14f*x4*x + 31.96f*x4 - 6.868f*x2*x
         + 0.4298f*x2 + 0.1191f*x - 0.00232f;
}
__device__ void agx(float& r, float& g, float& b){
    float ir = 0.842479062253094f*r + 0.0784335999999992f*g + 0.0792237451477643f*b;
    float ig = 0.0423282422610123f*r + 0.878468636469772f *g + 0.0791661274605434f*b;
    float ib = 0.0423756549057051f*r + 0.0784336f          *g + 0.879142973793104f *b;
    const float lo = -12.47393f, hi = 4.026069f;
    ir = clampf((log2f(fmaxf(ir, 1e-10f)) - lo) / (hi - lo), 0.0f, 1.0f);
    ig = clampf((log2f(fmaxf(ig, 1e-10f)) - lo) / (hi - lo), 0.0f, 1.0f);
    ib = clampf((log2f(fmaxf(ib, 1e-10f)) - lo) / (hi - lo), 0.0f, 1.0f);
    ir = agx_c(ir); ig = agx_c(ig); ib = agx_c(ib);
    r = clampf( 1.19687900512017f*ir - 0.0980208811401368f*ig - 0.0990297440797205f*ib, 0.0f, 1.0f);
    g = clampf(-0.0528968517574562f*ir + 1.15190312990417f*ig - 0.0989611768448433f*ib, 0.0f, 1.0f);
    b = clampf(-0.0529716355144438f*ir - 0.0980434501171241f*ig + 1.15107367264116f*ib, 0.0f, 1.0f);
}

// ----------------------------------------------------------------------------
// Kernels: init, tick, splat, histogram, bloom, composite, HUD
// ----------------------------------------------------------------------------
// analytic enclosed-mass fraction of the galaxy profile (disk 68% / bulge 12% / halo 20%)
__device__ float encMassFrac(float r){
    float fd = 0.0f;
    if (r > 14.0f){ float x = (r - 14.0f)/45.0f; fd = 1.0f - (1.0f + x)*expf(-x); }
    float t = r*0.1f;
    float fb = erff(t*0.70710678f) - 0.7978845f*t*expf(-0.5f*t*t);
    float fh = (r <= 40.0f) ? 0.0f : fminf(logf(r/40.0f)/1.7917595f, 1.0f);
    return 0.68f*fd + 0.12f*fb + 0.20f*fh;
}
__device__ __forceinline__ float vcirc(float r, float Mtot){
    return sqrtf(G_DIAL*Mtot*encMassFrac(r)/fmaxf(r, 2.0f));
}

__global__ void kInit(float4* pos, float4* mom, float* tau, unsigned* regime,
                      int N, unsigned long long seed, float Mtot){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= N) return;
    using namespace orrery;
    double u0 = counter_uniform(seed, i, 0, 0);
    float x, y, z, vx, vy, vz, T;
    if (u0 < 0.68){                                   // spiral disk
        double u1 = counter_uniform(seed, i, 1, 0);
        float r = -45.0f * logf((float)fmax(1.0 - u1, 1e-9));
        r = clampf(r + 14.0f, 14.0f, 230.0f);         // inner cutoff: bulge owns the core
        int arm = (counter_uniform(seed, i, 2, 0) < 0.5) ? 0 : 1;
        float th = logf(r / 6.0f) / tanf(0.24f) + arm*PI_F
                 + 0.22f*(float)counter_gauss(seed, i, 3, 0);
        float zz = (float)counter_gauss(seed, i, 4, 0) * (3.5f*expf(-r/70.0f) + 0.5f);
        x = r*cosf(th); y = r*sinf(th); z = zz;
        float vc = vcirc(r, Mtot);
        vx = -vc*sinf(th) + 0.08f*vc*(float)counter_gauss(seed, i, 5, 0);
        vy =  vc*cosf(th) + 0.08f*vc*(float)counter_gauss(seed, i, 6, 0);
        vz =  0.05f*vc*(float)counter_gauss(seed, i, 7, 0);
        double ut = counter_uniform(seed, i, 8, 0);
        if (ut < 0.07) T = 9500.0f + 5500.0f*(float)counter_uniform(seed, i, 9, 0);  // blue giants
        else           T = 3300.0f + 3200.0f*(float)(ut*ut);
    } else if (u0 < 0.80){                            // bulge
        x = 10.0f*(float)counter_gauss(seed, i, 10, 0);
        y = 10.0f*(float)counter_gauss(seed, i, 11, 0);
        z = 7.0f*(float)counter_gauss(seed, i, 12, 0);
        float rr = sqrtf(x*x + y*y + z*z) + 2.0f;
        float vc = vcirc(rr, Mtot);
        float ir = 1.0f / (sqrtf(x*x + y*y) + 2.0f);
        vx = -0.6f*vc*y*ir + 0.35f*vc*(float)counter_gauss(seed, i, 22, 0);
        vy =  0.6f*vc*x*ir + 0.35f*vc*(float)counter_gauss(seed, i, 13, 0);
        vz =  0.35f*vc*(float)counter_gauss(seed, i, 23, 0);
        T = 3100.0f + 2100.0f*(float)counter_uniform(seed, i, 14, 0);
    } else {                                          // halo (log-uniform shell)
        float r = 40.0f * powf(240.0f/40.0f, (float)counter_uniform(seed, i, 15, 0));
        float ct = 2.0f*(float)counter_uniform(seed, i, 16, 0) - 1.0f;
        float st = sqrtf(fmaxf(1.0f - ct*ct, 0.0f));
        float ph = 2.0f*PI_F*(float)counter_uniform(seed, i, 17, 0);
        x = r*st*cosf(ph); y = r*st*sinf(ph); z = r*ct;
        float sig = 0.35f*vcirc(r, Mtot);
        vx = sig*(float)counter_gauss(seed, i, 18, 0);
        vy = sig*(float)counter_gauss(seed, i, 19, 0);
        vz = sig*(float)counter_gauss(seed, i, 20, 0);
        T = 2800.0f + 1800.0f*(float)counter_uniform(seed, i, 21, 0);
    }
    pos[i] = make_float4(x, y, z, T);
    mom[i] = make_float4(vx, vy, vz, 1.0f);          // mass 1 (p = mv, v << c)
    tau[i] = 0.0f;
    regime[i] = 0x02u;                                // CLASSICAL (derived; inert scene)
}

__global__ void kInitCloud(float4* pos, float4* mom, float* tau, unsigned* regime,
                           int N, unsigned long long seed){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= N) return;
    using namespace orrery;
    float gx = (float)counter_gauss(seed, i, 0, 1);
    float gy = (float)counter_gauss(seed, i, 1, 1);
    float gz = (float)counter_gauss(seed, i, 2, 1);
    float il = rsqrtf(gx*gx + gy*gy + gz*gz + 1e-12f);
    float r = 100.0f * cbrtf((float)counter_uniform(seed, i, 3, 1));
    pos[i] = make_float4(r*gx*il, r*gy*il, r*gz*il, 3800.0f);
    // sub-virial support (0.35 v_vir): a truly cold pressureless collapse is
    // singular — unresolvable at any fixed dt (contract changelog, D-014)
    float sigv = 0.35f*sqrtf(G_DIAL*(float)N/100.0f);
    mom[i] = make_float4(sigv*(float)counter_gauss(seed, i, 4, 1),
                         sigv*(float)counter_gauss(seed, i, 5, 1),
                         sigv*(float)counter_gauss(seed, i, 6, 1), 1.0f);
    tau[i] = 0.0f; regime[i] = 0x02u;
}

// --- KDK pieces (frame contract: ping-pong publish; Kahan drift compensation) ---
__global__ void kKick(float4* mom, const float4* acc, int N, float hdt){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= N) return;
    float4 m = mom[i], a = acc[i];
    m.x += m.w*a.x*hdt; m.y += m.w*a.y*hdt; m.z += m.w*a.z*hdt;
    mom[i] = m;
}
__global__ void kDriftK(const float4* pin, float4* pout, float4* perr,
                        const float4* mom, float* tau, unsigned* regime,
                        int N, float dt, float inv2c2){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= N) return;
    float4 p = pin[i], m = mom[i], e = perr[i];
    float im = 1.0f/m.w;
    float3 d = make_float3(m.x*im*dt, m.y*im*dt, m.z*im*dt);
    float y, t;
    y = d.x - e.x; t = p.x + y; e.x = (t - p.x) - y; p.x = t;
    y = d.y - e.y; t = p.y + y; e.y = (t - p.y) - y; p.y = t;
    y = d.z - e.z; t = p.z + y; e.z = (t - p.z) - y; p.z = t;
    pout[i] = p; perr[i] = e;
    float v2 = (m.x*m.x + m.y*m.y + m.z*m.z)*im*im;
    tau[i] += dt * (1.0f - v2*inv2c2);
    regime[i] = (v2 > REL_V2) ? 0x04u : 0x02u;
}

// --- PM pipeline: fixed-point CIC deposit -> cuFFT -> Green -> force grid -> gather ---
__global__ void kZeroFix(unsigned long long* g, int n){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < n) g[i] = 0ull;
}
__device__ __forceinline__ int wrapc(int c){ return c & (PMN - 1); }
__global__ void kDeposit(const float4* pos, const float4* mom, int N, unsigned long long* g){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= N) return;
    float4 p = pos[i];
    float m = mom[i].w;
    float gx = (p.x + BOXL*0.5f)/CELL, gy = (p.y + BOXL*0.5f)/CELL, gz = (p.z + BOXL*0.5f)/CELL;
    int x0 = (int)floorf(gx), y0 = (int)floorf(gy), z0 = (int)floorf(gz);
    float fx = gx - x0, fy = gy - y0, fz = gz - z0;
    for (int dz = 0; dz < 2; dz++)
        for (int dy = 0; dy < 2; dy++)
            for (int dx = 0; dx < 2; dx++){
                float w = (dx ? fx : 1.0f - fx)*(dy ? fy : 1.0f - fy)*(dz ? fz : 1.0f - fz);
                int c = (wrapc(x0 + dx)*PMN + wrapc(y0 + dy))*PMN + wrapc(z0 + dz);
                orrery::fixed_atomic_add(&g[c], (double)(w*m));   // order-invariant (Invariant 4)
            }
}
__global__ void kFixToReal(const unsigned long long* g, float* r, int n){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < n) r[i] = (float)orrery::fixed_decode((long long)g[i]);
}
__global__ void kGreen(cufftComplex* s){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    int total = PMN*PMN*PMNZC;
    if (i >= total) return;
    int z = i % PMNZC, y = (i/PMNZC) % PMN, x = i/(PMNZC*PMN);
    int fx = (x <= PMN/2) ? x : x - PMN;
    int fy = (y <= PMN/2) ? y : y - PMN;
    int fz = z;
    float kf = 2.0f*PI_F/BOXL;
    float k2 = kf*kf*(float)(fx*fx + fy*fy + fz*fz);
    float f = (k2 > 0.0f)
        ? -4.0f*PI_F*G_DIAL/(k2 * (CELL*CELL*CELL) * (float)(PMN*PMN*PMN))
        : 0.0f;                                       // k=0: mean density does not gravitate
    s[i].x *= f; s[i].y *= f;
}
__global__ void kForceGrid(const float* phi, float4* F){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    int n = PMN*PMN*PMN;
    if (i >= n) return;
    int z = i % PMN, y = (i/PMN) % PMN, x = i/(PMN*PMN);
    float inv2h = 1.0f/(2.0f*CELL);
    #define PHI(X,Y,Z) phi[(wrapc(X)*PMN + wrapc(Y))*PMN + wrapc(Z)]
    float fx = -(PHI(x+1,y,z) - PHI(x-1,y,z))*inv2h;
    float fy = -(PHI(x,y+1,z) - PHI(x,y-1,z))*inv2h;
    float fz = -(PHI(x,y,z+1) - PHI(x,y,z-1))*inv2h;
    #undef PHI
    F[i] = make_float4(fx, fy, fz, 0.0f);
}
__global__ void kGather(const float4* pos, const float4* F, float4* acc, int N){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= N) return;
    float4 p = pos[i];
    float gx = (p.x + BOXL*0.5f)/CELL, gy = (p.y + BOXL*0.5f)/CELL, gz = (p.z + BOXL*0.5f)/CELL;
    int x0 = (int)floorf(gx), y0 = (int)floorf(gy), z0 = (int)floorf(gz);
    float fx = gx - x0, fy = gy - y0, fz = gz - z0;
    float3 a = make_float3(0, 0, 0);
    for (int dz = 0; dz < 2; dz++)
        for (int dy = 0; dy < 2; dy++)
            for (int dx = 0; dx < 2; dx++){
                float w = (dx ? fx : 1.0f - fx)*(dy ? fy : 1.0f - fy)*(dz ? fz : 1.0f - fz);
                float4 f = F[(wrapc(x0 + dx)*PMN + wrapc(y0 + dy))*PMN + wrapc(z0 + dz)];
                a.x += w*f.x; a.y += w*f.y; a.z += w*f.z;
            }
    acc[i] = make_float4(a.x, a.y, a.z, 0.0f);
}
__device__ float phiAt(const float* phi, float4 p){
    float gx = (p.x + BOXL*0.5f)/CELL, gy = (p.y + BOXL*0.5f)/CELL, gz = (p.z + BOXL*0.5f)/CELL;
    int x0 = (int)floorf(gx), y0 = (int)floorf(gy), z0 = (int)floorf(gz);
    float fx = gx - x0, fy = gy - y0, fz = gz - z0;
    float s = 0;
    for (int dz = 0; dz < 2; dz++)
        for (int dy = 0; dy < 2; dy++)
            for (int dx = 0; dx < 2; dx++){
                float w = (dx ? fx : 1.0f - fx)*(dy ? fy : 1.0f - fy)*(dz ? fz : 1.0f - fz);
                s += w*phi[(wrapc(x0 + dx)*PMN + wrapc(y0 + dy))*PMN + wrapc(z0 + dz)];
            }
    return s;
}

// --- direct solver (N <= 4096): fixed order, register accumulation, no atomics ---
__global__ void kDirect(const float4* pos, const float4* mom, float4* acc, int N){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= N) return;
    float4 pi = pos[i];
    float3 a = make_float3(0, 0, 0);
    for (int j = 0; j < N; j++){
        if (j == i) continue;
        float4 pj = pos[j];
        float dx = pj.x - pi.x, dy = pj.y - pi.y, dz = pj.z - pi.z;
        float r2 = dx*dx + dy*dy + dz*dz + 1e-6f;
        float ir = rsqrtf(r2);
        float f = G_DIAL*mom[j].w*ir*ir*ir;
        a.x += f*dx; a.y += f*dy; a.z += f*dz;
    }
    acc[i] = make_float4(a.x, a.y, a.z, 0.0f);
}

// --- tiny solver (N <= 32): one block, many ticks per launch, fp64 internal ---
// (D-014: tiny scenarios are the sim's precision references; declared state
//  round-trips through the fp32 frame-contract buffers only at batch edges)
__device__ double3 tinyAcc(const double3* sp, const double* sm, int N, int i){
    double3 a = make_double3(0, 0, 0);
    for (int j = 0; j < N; j++){
        if (j == i) continue;
        double dx = sp[j].x - sp[i].x, dy = sp[j].y - sp[i].y, dz = sp[j].z - sp[i].z;
        double r2 = dx*dx + dy*dy + dz*dz + 1e-6;
        double ir = rsqrt(r2);
        double f = (double)G_DIAL*sm[j]*ir*ir*ir;
        a.x += f*dx; a.y += f*dy; a.z += f*dz;
    }
    return a;
}
__global__ void kTiny(const float4* pin, float4* pout, float4* mom, float* tau,
                      unsigned* regime, int N, int steps, float dtf, float inv2c2f){
    __shared__ double3 sp[32], sv[32], sa[32];
    __shared__ double  sm[32];
    int i = threadIdx.x;
    if (i >= N) return;
    double dt = (double)dtf, inv2c2 = (double)inv2c2f;
    float4 P = pin[i], M = mom[i];
    sp[i] = make_double3(P.x, P.y, P.z);
    sm[i] = M.w;
    sv[i] = make_double3(M.x/(double)M.w, M.y/(double)M.w, M.z/(double)M.w);
    double tu = tau[i];
    __syncthreads();
    sa[i] = tinyAcc(sp, sm, N, i);
    __syncthreads();
    double hdt = 0.5*dt;
    for (int s = 0; s < steps; s++){
        sv[i].x += sa[i].x*hdt; sv[i].y += sa[i].y*hdt; sv[i].z += sa[i].z*hdt;
        sp[i].x += sv[i].x*dt;  sp[i].y += sv[i].y*dt;  sp[i].z += sv[i].z*dt;
        __syncthreads();
        sa[i] = tinyAcc(sp, sm, N, i);
        __syncthreads();
        sv[i].x += sa[i].x*hdt; sv[i].y += sa[i].y*hdt; sv[i].z += sa[i].z*hdt;
        double v2 = sv[i].x*sv[i].x + sv[i].y*sv[i].y + sv[i].z*sv[i].z;
        tu += dt*(1.0 - v2*inv2c2);
    }
    pout[i] = make_float4((float)sp[i].x, (float)sp[i].y, (float)sp[i].z, P.w);
    mom[i]  = make_float4((float)(sv[i].x*sm[i]), (float)(sv[i].y*sm[i]),
                          (float)(sv[i].z*sm[i]), (float)sm[i]);
    tau[i]  = (float)tu;
    double v2 = sv[i].x*sv[i].x + sv[i].y*sv[i].y + sv[i].z*sv[i].z;
    regime[i] = (v2 > (double)REL_V2) ? 0x04u : 0x02u;
}

// --- conservation meters (fixed-point accumulators; PE via PM potential when available)
// acc slots: 0 KE · 1 PE · 2-4 Px,Py,Pz · 5 Lz · 6 nRel · 7 nClassical
__global__ void kMeters(const float4* pos, const float4* mom, const float* phi,
                        int N, unsigned long long* a8){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= N) return;
    float4 p = pos[i], m = mom[i];
    float im = 1.0f/m.w;
    double ke = 0.5*(double)(m.x*m.x + m.y*m.y + m.z*m.z)*im;
    orrery::fixed_atomic_add(&a8[0], ke);
    if (phi) orrery::fixed_atomic_add(&a8[1], 0.5*(double)m.w*(double)phiAt(phi, p));
    orrery::fixed_atomic_add(&a8[2], (double)m.x);
    orrery::fixed_atomic_add(&a8[3], (double)m.y);
    orrery::fixed_atomic_add(&a8[4], (double)m.z);
    orrery::fixed_atomic_add(&a8[5], (double)(p.x*m.y - p.y*m.x));
    atomicAdd(&a8[(m.x*m.x + m.y*m.y + m.z*m.z)*im*im > REL_V2 ? 6 : 7], 1ull);
}

__global__ void kClear(float4* hdr, int n){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < n) hdr[i] = make_float4(0, 0, 0, 0);
}

__global__ void kSplat(const float4* pos, int N, float4* hdr, RP rp){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= N) return;
    float4 p = pos[i];
    float rx = p.x - rp.cpos[0], ry = p.y - rp.cpos[1], rz = p.z - rp.cpos[2];
    float cz = rx*rp.fwd[0] + ry*rp.fwd[1] + rz*rp.fwd[2];
    if (cz < 1.0f) return;
    float cx = rx*rp.rgt[0] + ry*rp.rgt[1] + rz*rp.rgt[2];
    float cy = rx*rp.upv[0] + ry*rp.upv[1] + rz*rp.upv[2];
    float ndx = cx / (cz * rp.tf * rp.aspect);
    float ndy = cy / (cz * rp.tf);
    float px = (ndx*0.5f + 0.5f) * rp.RW - 0.5f;
    float py = (0.5f - ndy*0.5f) * rp.RH - 0.5f;
    if (px < 0 || py < 0 || px >= rp.RW - 1 || py >= rp.RH - 1) return;
    float T = p.w;
    float L = rp.mode ? powf(T/5800.0f, 4.0f)         // physical: L ~ T^4
                      : powf(T/5800.0f, 2.2f);        // cinematic: softened
    float flux = 4200.0f * L / (cz*cz + 25.0f);
    float bb[3]; blackbody(T, bb);
    int x0 = (int)px, y0 = (int)py;
    float fx = px - x0, fy = py - y0;
    float w00 = (1-fx)*(1-fy), w10 = fx*(1-fy), w01 = (1-fx)*fy, w11 = fx*fy;
    float4* h;
    h = &hdr[y0*rp.RW + x0];
    atomicAdd(&h->x, bb[0]*flux*w00); atomicAdd(&h->y, bb[1]*flux*w00); atomicAdd(&h->z, bb[2]*flux*w00);
    h = &hdr[y0*rp.RW + x0 + 1];
    atomicAdd(&h->x, bb[0]*flux*w10); atomicAdd(&h->y, bb[1]*flux*w10); atomicAdd(&h->z, bb[2]*flux*w10);
    h = &hdr[(y0+1)*rp.RW + x0];
    atomicAdd(&h->x, bb[0]*flux*w01); atomicAdd(&h->y, bb[1]*flux*w01); atomicAdd(&h->z, bb[2]*flux*w01);
    h = &hdr[(y0+1)*rp.RW + x0 + 1];
    atomicAdd(&h->x, bb[0]*flux*w11); atomicAdd(&h->y, bb[1]*flux*w11); atomicAdd(&h->z, bb[2]*flux*w11);
}

__global__ void kHist(const float4* hdr, int n, unsigned* hist){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= n) return;
    float4 c = hdr[i];
    float lum = 0.2126f*c.x + 0.7152f*c.y + 0.0722f*c.z;
    int bin = (int)((log2f(lum + 1e-8f) + 16.0f) * 8.0f);
    bin = bin < 0 ? 0 : (bin > 255 ? 255 : bin);
    atomicAdd(&hist[bin], 1u);
}

// 13-tap CoD-style downsample (threshold-free)
__global__ void kDown13(const float4* src, int sw, int sh, float4* dst, int dw, int dh){
    int x = blockIdx.x*blockDim.x + threadIdx.x;
    int y = blockIdx.y*blockDim.y + threadIdx.y;
    if (x >= dw || y >= dh) return;
    int cx = x*2, cy = y*2;
    #define S(ox, oy) src[min(max(cy + (oy), 0), sh - 1)*sw + min(max(cx + (ox), 0), sw - 1)]
    float4 acc = make_float4(0, 0, 0, 0);
    float4 s;
    // inner quad (w = 0.125 each)
    s = S(-1,-1); acc.x += 0.125f*s.x; acc.y += 0.125f*s.y; acc.z += 0.125f*s.z;
    s = S( 1,-1); acc.x += 0.125f*s.x; acc.y += 0.125f*s.y; acc.z += 0.125f*s.z;
    s = S(-1, 1); acc.x += 0.125f*s.x; acc.y += 0.125f*s.y; acc.z += 0.125f*s.z;
    s = S( 1, 1); acc.x += 0.125f*s.x; acc.y += 0.125f*s.y; acc.z += 0.125f*s.z;
    // center (0.125) + edges (0.0625) + corners (0.03125)
    s = S( 0, 0); acc.x += 0.125f*s.x; acc.y += 0.125f*s.y; acc.z += 0.125f*s.z;
    s = S(-2, 0); acc.x += 0.0625f*s.x; acc.y += 0.0625f*s.y; acc.z += 0.0625f*s.z;
    s = S( 2, 0); acc.x += 0.0625f*s.x; acc.y += 0.0625f*s.y; acc.z += 0.0625f*s.z;
    s = S( 0,-2); acc.x += 0.0625f*s.x; acc.y += 0.0625f*s.y; acc.z += 0.0625f*s.z;
    s = S( 0, 2); acc.x += 0.0625f*s.x; acc.y += 0.0625f*s.y; acc.z += 0.0625f*s.z;
    s = S(-2,-2); acc.x += 0.03125f*s.x; acc.y += 0.03125f*s.y; acc.z += 0.03125f*s.z;
    s = S( 2,-2); acc.x += 0.03125f*s.x; acc.y += 0.03125f*s.y; acc.z += 0.03125f*s.z;
    s = S(-2, 2); acc.x += 0.03125f*s.x; acc.y += 0.03125f*s.y; acc.z += 0.03125f*s.z;
    s = S( 2, 2); acc.x += 0.03125f*s.x; acc.y += 0.03125f*s.y; acc.z += 0.03125f*s.z;
    #undef S
    dst[y*dw + x] = acc;
}

// 3x3 tent upsample of lo, added into hi
__global__ void kUpTent(const float4* lo, int lw, int lh, float4* hi, int hw, int hh){
    int x = blockIdx.x*blockDim.x + threadIdx.x;
    int y = blockIdx.y*blockDim.y + threadIdx.y;
    if (x >= hw || y >= hh) return;
    int lx = x >> 1, ly = y >> 1;
    float4 acc = make_float4(0, 0, 0, 0);
    const float w[3] = {0.25f, 0.5f, 0.25f};
    for (int dy = -1; dy <= 1; dy++)
        for (int dx = -1; dx <= 1; dx++){
            int sx = min(max(lx + dx, 0), lw - 1);
            int sy = min(max(ly + dy, 0), lh - 1);
            float ww = w[dx+1]*w[dy+1];
            float4 s = lo[sy*lw + sx];
            acc.x += ww*s.x; acc.y += ww*s.y; acc.z += ww*s.z;
        }
    float4* d = &hi[y*hw + x];
    d->x += acc.x; d->y += acc.y; d->z += acc.z;
}

__global__ void kComposite(const float4* hdr, const float4* bloom, int bw, int bh,
                           uchar4* out, RP rp){
    int ox = blockIdx.x*blockDim.x + threadIdx.x;
    int oy = blockIdx.y*blockDim.y + threadIdx.y;
    if (ox >= rp.OW || oy >= rp.OH) return;
    float r = 0, g = 0, b = 0;
    for (int sy = 0; sy < rp.ss; sy++)
        for (int sx = 0; sx < rp.ss; sx++){
            int x = ox*rp.ss + sx, y = oy*rp.ss + sy;
            float4 c = hdr[y*rp.RW + x];
            float cr = c.x, cg = c.y, cb = c.z;
            if (rp.bloomOn){
                int lx = min(x >> 1, bw - 1), ly = min(y >> 1, bh - 1);
                float4 bl = bloom[ly*bw + lx];
                cr += rp.bloomAmt*bl.x; cg += rp.bloomAmt*bl.y; cb += rp.bloomAmt*bl.z;
            }
            r += cr; g += cg; b += cb;
        }
    float inv = 1.0f / (rp.ss*rp.ss);
    r *= inv*rp.exposure; g *= inv*rp.exposure; b *= inv*rp.exposure;
    if (!rp.mode){                                    // cinematic (declared deviations, HUD):
        float lum = 0.2126f*r + 0.7152f*g + 0.0722f*b;
        if (lum > 1.0f){                              // astro stretch: highlights only (L>1),
            const float Wp = 40.0f;                   // extended Reinhard on (L-1); mids untouched
            float e = lum - 1.0f;
            float ld = 1.0f + e*(1.0f + e/(Wp*Wp))/(1.0f + e);
            float s = ld/lum;
            r *= s; g *= s; b *= s;
            lum = ld;
        }
        r = lum + (r - lum)*rp.satCine;               // +10% saturation
        g = lum + (g - lum)*rp.satCine;
        b = lum + (b - lum)*rp.satCine;
        r = fmaxf(r, 0.0f); g = fmaxf(g, 0.0f); b = fmaxf(b, 0.0f);
    }
    if (rp.tonemap == 0){ agx(r, g, b); }
    else { r = powf(aces(r), 1.0f/2.2f); g = powf(aces(g), 1.0f/2.2f); b = powf(aces(b), 1.0f/2.2f); }
    // triangular dither before quantization (deterministic per pixel+frame)
    using namespace orrery;
    float d1 = (float)u01(hash4(rp.frame, ox, oy, 1));
    float d2 = (float)u01(hash4(rp.frame, ox, oy, 2));
    float dith = (d1 + d2 - 1.0f) * (1.0f/255.0f);
    r = clampf(r + dith, 0.0f, 1.0f);
    g = clampf(g + dith, 0.0f, 1.0f);
    b = clampf(b + dith, 0.0f, 1.0f);
    out[oy*rp.OW + ox] = make_uchar4((unsigned char)(r*255.0f + 0.5f),
                                     (unsigned char)(g*255.0f + 0.5f),
                                     (unsigned char)(b*255.0f + 0.5f), 255);
}

// ----------------------------------------------------------------------------
// HUD: 5x7 bitmap font, drawn by CUDA (typed panels: physics cyan, render amber)
// ----------------------------------------------------------------------------
struct Glyph { short x, y; unsigned char ch, col; };

__device__ __constant__ unsigned char dFont[44][7] = {
    {0x0E,0x11,0x13,0x15,0x19,0x11,0x0E},{0x04,0x0C,0x04,0x04,0x04,0x04,0x0E}, // 0 1
    {0x0E,0x11,0x01,0x06,0x08,0x10,0x1F},{0x1F,0x02,0x04,0x02,0x01,0x11,0x0E}, // 2 3
    {0x02,0x06,0x0A,0x12,0x1F,0x02,0x02},{0x1F,0x10,0x1E,0x01,0x01,0x11,0x0E}, // 4 5
    {0x06,0x08,0x10,0x1E,0x11,0x11,0x0E},{0x1F,0x01,0x02,0x04,0x08,0x08,0x08}, // 6 7
    {0x0E,0x11,0x11,0x0E,0x11,0x11,0x0E},{0x0E,0x11,0x11,0x0F,0x01,0x02,0x0C}, // 8 9
    {0x0E,0x11,0x11,0x1F,0x11,0x11,0x11},{0x1E,0x11,0x11,0x1E,0x11,0x11,0x1E}, // A B
    {0x0E,0x11,0x10,0x10,0x10,0x11,0x0E},{0x1C,0x12,0x11,0x11,0x11,0x12,0x1C}, // C D
    {0x1F,0x10,0x10,0x1E,0x10,0x10,0x1F},{0x1F,0x10,0x10,0x1E,0x10,0x10,0x10}, // E F
    {0x0E,0x11,0x10,0x17,0x11,0x11,0x0F},{0x11,0x11,0x11,0x1F,0x11,0x11,0x11}, // G H
    {0x0E,0x04,0x04,0x04,0x04,0x04,0x0E},{0x07,0x02,0x02,0x02,0x02,0x12,0x0C}, // I J
    {0x11,0x12,0x14,0x18,0x14,0x12,0x11},{0x10,0x10,0x10,0x10,0x10,0x10,0x1F}, // K L
    {0x11,0x1B,0x15,0x15,0x11,0x11,0x11},{0x11,0x19,0x15,0x13,0x11,0x11,0x11}, // M N
    {0x0E,0x11,0x11,0x11,0x11,0x11,0x0E},{0x1E,0x11,0x11,0x1E,0x10,0x10,0x10}, // O P
    {0x0E,0x11,0x11,0x11,0x15,0x12,0x0D},{0x1E,0x11,0x11,0x1E,0x14,0x12,0x11}, // Q R
    {0x0F,0x10,0x10,0x0E,0x01,0x01,0x1E},{0x1F,0x04,0x04,0x04,0x04,0x04,0x04}, // S T
    {0x11,0x11,0x11,0x11,0x11,0x11,0x0E},{0x11,0x11,0x11,0x11,0x11,0x0A,0x04}, // U V
    {0x11,0x11,0x11,0x15,0x15,0x1B,0x11},{0x11,0x11,0x0A,0x04,0x0A,0x11,0x11}, // W X
    {0x11,0x11,0x0A,0x04,0x04,0x04,0x04},{0x1F,0x01,0x02,0x04,0x08,0x10,0x1F}, // Y Z
    {0x00,0x00,0x00,0x00,0x00,0x0C,0x0C},{0x00,0x0C,0x0C,0x00,0x0C,0x0C,0x00}, // . :
    {0x00,0x04,0x04,0x1F,0x04,0x04,0x00},{0x00,0x00,0x00,0x1F,0x00,0x00,0x00}, // + -
    {0x01,0x02,0x04,0x04,0x08,0x10,0x00},{0x00,0x00,0x1F,0x00,0x1F,0x00,0x00}, // / =
    {0x02,0x04,0x08,0x08,0x08,0x04,0x02},{0x08,0x04,0x02,0x02,0x02,0x04,0x08}, // ( )
};

__host__ __device__ inline int glyphIndex(char c){
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'A' && c <= 'Z') return 10 + (c - 'A');
    switch (c){ case '.': return 36; case ':': return 37; case '+': return 38;
                case '-': return 39; case '/': return 40; case '=': return 41;
                case '(': return 42; case ')': return 43; }
    return -1;
}

__global__ void kHud(uchar4* out, int W, int H, const Glyph* gl, int n){
    int g = blockIdx.x;
    int t = threadIdx.x;                              // 0..69 (5x7 at 2x scale)
    if (g >= n || t >= 70) return;
    Glyph q = gl[g];
    int idx = glyphIndex((char)q.ch);
    if (idx < 0) return;
    int gx = t % 10, gy = t / 10;                     // 10x14 block, 2x scale
    int fx = gx >> 1, fy = gy >> 1;
    if (!((dFont[idx][fy] >> (4 - fx)) & 1)) return;
    unsigned char cr, cg, cb;
    if (q.col == 0){ cr = 120; cg = 225; cb = 255; }   // physics: cyan
    else            { cr = 255; cg = 200; cb = 110; }  // render: amber
    int px = q.x + gx, py = q.y + gy;
    if (px < 0 || py < 0 || px >= W - 1 || py >= H - 1) return;
    out[(py+1)*W + px + 1] = make_uchar4(0, 0, 0, 255);         // drop shadow
    out[py*W + px] = make_uchar4(cr, cg, cb, 255);
}

// ----------------------------------------------------------------------------
// Host state
// ----------------------------------------------------------------------------
static int   gN = 1000000;
static unsigned long long gSeed = 20260711ull;
static int   gOW = 1920, gOH = 1080, gSS = 1;
static const float DT = 1.0f/240.0f;                  // dial
static const float C_LIGHT = 20.0f;                   // dial

static float4 *dPos[2] = {0,0}, *dMom = 0;
static float  *dTau = 0;
static unsigned *dRegime = 0, *dHist = 0;
static float4 *dHdr = 0, *dMip[7] = {0};
static int    mipW[7], mipH[7];
static uchar4 *dOut = 0;
static Glyph  *dGlyphs = 0;
static cudaStream_t simS, prsS;
static cudaEvent_t  tickDone;
static int    gPub = 0;                               // published pos buffer

// --- M2 gravity state ---
static Scenario gScenario = SC_GALAXY;
static Solver   gSolver = SOLV_PM;
static float4 *dAcc = 0, *dPosErr = 0, *dFor = 0;
static unsigned long long *dFix = 0, *dMet = 0;
static float  *dReal = 0;
static cufftComplex *dSpec = 0;
static cufftHandle planF, planI;
static double gMTot = 1e6;
struct Met { double KE, PE, E, Px, Py, Pz, Lz; unsigned long long nRel, nCls; };
static Met  gM0 = {};                                  // meters at t=0
static double gDE = 0, gDP = 0;                        // HUD drift readouts
#define CUFFT_CHECK(call) do { cufftResult _r = (call); if (_r != CUFFT_SUCCESS){ \
    fprintf(stderr, "cuFFT error %d at %s:%d\n", (int)_r, __FILE__, __LINE__); exit(2);} } while (0)

// camera (damped inertial orbit — CINEMATIC §5)
struct Spring { float cur, tgt, vel; };
static Spring cAz = {35, 35, 0}, cEl = {24, 24, 0}, cDist = {5.94f, 5.94f, 0}; // dist is log
static int   gDrag = 0, gLastX = 0, gLastY = 0, gAutoOrbit = 1, gPaused = 0;
static int   gMode = 0, gTonemap = 0, gBloomOn = 1, gAutoExp = 1, gHud = 1, gVsync = 1;
static float gEvOffset = 0.0f, gExposure = 1.0f, gLavgSmooth = -1.0f;
static double gSimTime = 0.0; static unsigned long long gTick = 0;
static unsigned gFrame = 0;
static bool gRunning = true;
static float gFps = 0, gMs = 0;

static void springStep(Spring& s, float dt){
    const float k = 90.0f, c = 19.0f;                 // ~critically damped, ~9.5/s
    s.vel += (k*(s.tgt - s.cur) - c*s.vel)*dt;
    s.cur += s.vel*dt;
}
static void applyPreset(int i){
    const float D0[4] = {380.0f, 650.0f, 9.0f, 420.0f};   // per-scenario base distance
    float d = D0[gScenario];
    if (i == 0){ cAz.tgt = 35;  cEl.tgt = 24; cDist.tgt = logf(d); }
    if (i == 1){ cAz.tgt = 120; cEl.tgt = 10; cDist.tgt = logf(0.45f*d); }
    if (i == 2){ cAz.tgt = 80;  cEl.tgt = 86; cDist.tgt = logf(0.87f*d); }
}

static void buildCamera(RP& rp){
    float az = cAz.cur*PI_F/180.0f, el = cEl.cur*PI_F/180.0f, d = expf(cDist.cur);
    float ce = cosf(el), se = sinf(el);
    rp.cpos[0] = d*ce*cosf(az); rp.cpos[1] = d*ce*sinf(az); rp.cpos[2] = d*se;
    float f[3] = {-rp.cpos[0]/d, -rp.cpos[1]/d, -rp.cpos[2]/d};
    float up[3] = {0, 0, 1};
    if (fabsf(f[2]) > 0.985f){ up[0] = 0; up[1] = 1; up[2] = 0; }
    float r[3] = { f[1]*up[2]-f[2]*up[1], f[2]*up[0]-f[0]*up[2], f[0]*up[1]-f[1]*up[0] };
    float rl = sqrtf(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    r[0]/=rl; r[1]/=rl; r[2]/=rl;
    float u[3] = { r[1]*f[2]-r[2]*f[1], r[2]*f[0]-r[0]*f[2], r[0]*f[1]-r[1]*f[0] };
    memcpy(rp.fwd, f, sizeof f); memcpy(rp.rgt, r, sizeof r); memcpy(rp.upv, u, sizeof u);
    rp.tf = tanf(55.0f*0.5f*PI_F/180.0f);
    rp.aspect = (float)gOW/(float)gOH;
}

// ----------------------------------------------------------------------------
// Render chain (presentStream)
// ----------------------------------------------------------------------------
static unsigned hHist[256];

static void renderFrame(RP& rp){
    int RW = rp.RW, RH = rp.RH, npx = RW*RH;
    dim3 b1(256), g1((npx + 255)/256);
    kClear<<<g1, b1, 0, prsS>>>(dHdr, npx);
    kSplat<<<(gN + 255)/256, b1, 0, prsS>>>(dPos[gPub], gN, dHdr, rp);
    // auto-exposure histogram
    CUDA_CHECK(cudaMemsetAsync(dHist, 0, 256*sizeof(unsigned), prsS));
    kHist<<<g1, b1, 0, prsS>>>(dHdr, npx, dHist);
    CUDA_CHECK(cudaMemcpyAsync(hHist, dHist, 256*sizeof(unsigned),
                               cudaMemcpyDeviceToHost, prsS));
    // bloom mip chain (threshold-free)
    if (gBloomOn){
        dim3 b2(16, 16);
        const float4* src = dHdr; int sw = RW, sh = RH;
        for (int m = 1; m <= 6; m++){
            dim3 g2((mipW[m] + 15)/16, (mipH[m] + 15)/16);
            kDown13<<<g2, b2, 0, prsS>>>(src, sw, sh, dMip[m], mipW[m], mipH[m]);
            src = dMip[m]; sw = mipW[m]; sh = mipH[m];
        }
        for (int m = 5; m >= 1; m--){
            dim3 g2((mipW[m] + 15)/16, (mipH[m] + 15)/16);
            kUpTent<<<g2, b2, 0, prsS>>>(dMip[m+1], mipW[m+1], mipH[m+1],
                                         dMip[m], mipW[m], mipH[m]);
        }
    }
    CUDA_CHECK(cudaStreamSynchronize(prsS));          // hHist ready; EV on host
    double tot = 0; for (int i = 0; i < 256; i++) tot += hHist[i];
    double lo = tot*0.02, hi = tot*0.98, cum = 0, sum = 0, cnt = 0;
    for (int i = 0; i < 256; i++){
        double c0 = cum; cum += hHist[i];
        double take = fmax(0.0, fmin(cum, hi) - fmax(c0, lo));
        double lum = (i + 0.5)/8.0 - 16.0;            // log2 lum bin center
        sum += take*lum; cnt += take;
    }
    if (cnt > 0){
        float Lavg = exp2f((float)(sum/cnt));
        if (gLavgSmooth < 0) gLavgSmooth = Lavg;
        float a = 1.0f - expf(-1.0f/(60.0f*0.8f));    // ~0.8s adaptation
        gLavgSmooth += (Lavg - gLavgSmooth)*a;
    }
    float expo = gAutoExp && gLavgSmooth > 0 ? 0.06f/gLavgSmooth : 1.0f;
    gExposure = expo * exp2f(gEvOffset);
    rp.exposure = gExposure;
    dim3 b3(16, 16), g3((gOW + 15)/16, (gOH + 15)/16);
    kComposite<<<g3, b3, 0, prsS>>>(dHdr, dMip[1], mipW[1], mipH[1], dOut, rp);
}

static void drawHud(const RP& rp){
    if (!gHud) return;
    std::vector<Glyph> gs;
    char line[128];
    auto emit = [&](int x, int y, const char* s, unsigned char col){
        for (int i = 0; s[i]; i++){
            if (s[i] != ' ') gs.push_back({(short)(x + i*12), (short)y, (unsigned char)s[i], col});
        }
    };
    // physics panel (cyan, left)
    emit(16, 14, "TINY UNIVERSE M1 FIRST LIGHT", 0);
    snprintf(line, sizeof line, "N=%d DT=1/240 TICK=%llu", gN, (unsigned long long)gTick);
    emit(16, 32, line, 0);
    snprintf(line, sizeof line, "SIM T=%.1FS  C=20 HBAR=0.5 G=0.002", gSimTime);
    emit(16, 50, line, 0);
    if (gMode) emit(16, 68, "PHYSICAL: L=T4 SAT=1.0 NO CLIP", 0);
    else       emit(16, 68, "CINEMATIC: L=T2.2 SAT=1.2 STRETCH", 0);
    static const char* SC_HUD[4] = {"GALAXY", "KEPLER", "THREEBODY", "CLOUD"};
    snprintf(line, sizeof line, "%s %s  DE %+.1E  DP %.1E", SC_HUD[gScenario],
             gSolver == SOLV_PM ? "PM" : (gSolver == SOLV_TINY ? "TINY" : "DIRECT"), gDE, gDP);
    emit(16, 86, line, 0);
    // render panel (amber, right)
    int rx = gOW - 12*30 - 16;
    snprintf(line, sizeof line, "FPS %.0F  %.1F MS", gFps, gMs);
    emit(rx, 14, line, 1);
    snprintf(line, sizeof line, "EV %+.1F %s  EXPO %.3G", gEvOffset, gAutoExp ? "AUTO" : "MAN", gExposure);
    emit(rx, 32, line, 1);
    snprintf(line, sizeof line, "%s  BLOOM %s  SSAA %dX", gTonemap ? "ACES" : "AGX",
             gBloomOn ? "0.06" : "OFF", gSS);
    emit(rx, 50, line, 1);
    snprintf(line, sizeof line, "%s %s", gPaused ? "PAUSED" : "RUNNING", gVsync ? "VSYNC" : "");
    emit(rx, 68, line, 1);
    int n = (int)gs.size(); if (!n) return;
    if (n > 512) n = 512;
    CUDA_CHECK(cudaMemcpyAsync(dGlyphs, gs.data(), n*sizeof(Glyph),
                               cudaMemcpyHostToDevice, prsS));
    kHud<<<n, 70, 0, prsS>>>(dOut, gOW, gOH, dGlyphs, n);
}

// ----------------------------------------------------------------------------
// Sim ticks (simStream, ping-pong publish)
// ----------------------------------------------------------------------------
// invariant: dAcc = a(dPos[gPub]) after every full tick and after init
static void forcePass(){
    dim3 b(256);
    if (gSolver == SOLV_PM){
        int nc = PMN*PMN*PMN, ns = PMN*PMN*PMNZC;
        kZeroFix<<<(nc + 255)/256, b, 0, simS>>>(dFix, nc);
        kDeposit<<<(gN + 255)/256, b, 0, simS>>>(dPos[gPub], dMom, gN, dFix);
        kFixToReal<<<(nc + 255)/256, b, 0, simS>>>(dFix, dReal, nc);
        CUFFT_CHECK(cufftExecR2C(planF, dReal, dSpec));
        kGreen<<<(ns + 255)/256, b, 0, simS>>>(dSpec);
        CUFFT_CHECK(cufftExecC2R(planI, dSpec, dReal));   // dReal now holds Phi
        kForceGrid<<<(nc + 255)/256, b, 0, simS>>>(dReal, dFor);
        kGather<<<(gN + 255)/256, b, 0, simS>>>(dPos[gPub], dFor, dAcc, gN);
    } else if (gSolver == SOLV_DIRECT){
        kDirect<<<(gN + 255)/256, b, 0, simS>>>(dPos[gPub], dMom, dAcc, gN);
    }
}

static void runTicks(int nt){
    const float inv2c2 = 1.0f/(2.0f*C_LIGHT*C_LIGHT), hdt = 0.5f*DT;
    if (gSolver == SOLV_TINY){
        int src = gPub, dst = 1 - gPub;
        kTiny<<<1, 32, 0, simS>>>(dPos[src], dPos[dst], dMom, dTau, dRegime,
                                  gN, nt, DT, inv2c2);
        gPub = dst; gTick += nt; gSimTime += (double)nt*DT;
    } else {
        dim3 b(256), g((gN + 255)/256);
        for (int t = 0; t < nt; t++){
            kKick<<<g, b, 0, simS>>>(dMom, dAcc, gN, hdt);
            int src = gPub, dst = 1 - gPub;
            kDriftK<<<g, b, 0, simS>>>(dPos[src], dPos[dst], dPosErr, dMom, dTau,
                                       dRegime, gN, DT, inv2c2);
            gPub = dst;
            forcePass();
            kKick<<<g, b, 0, simS>>>(dMom, dAcc, gN, hdt);
            gTick++; gSimTime += DT;
        }
    }
    CUDA_CHECK(cudaEventRecord(tickDone, simS));
    CUDA_CHECK(cudaStreamWaitEvent(prsS, tickDone, 0));
}

// conservation meters. Big-N: device fixed-point accumulators (order-invariant,
// 2^-32 quanta declared). Tiny-N: host fp64 with exact pairwise PE.
static Met metersNow(){
    Met M = {};
    if (gSolver == SOLV_TINY || gN <= 4096){
        std::vector<float4> hp(gN), hm(gN);
        CUDA_CHECK(cudaStreamSynchronize(simS));
        CUDA_CHECK(cudaMemcpy(hp.data(), dPos[gPub], gN*sizeof(float4), cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(hm.data(), dMom, gN*sizeof(float4), cudaMemcpyDeviceToHost));
        for (int i = 0; i < gN; i++){
            double m = hm[i].w;
            M.KE += 0.5*((double)hm[i].x*hm[i].x + (double)hm[i].y*hm[i].y + (double)hm[i].z*hm[i].z)/m;
            M.Px += hm[i].x; M.Py += hm[i].y; M.Pz += hm[i].z;
            M.Lz += (double)hp[i].x*hm[i].y - (double)hp[i].y*hm[i].x;
            double v2 = ((double)hm[i].x*hm[i].x + (double)hm[i].y*hm[i].y + (double)hm[i].z*hm[i].z)/(m*m);
            if (v2 > REL_V2) M.nRel++; else M.nCls++;
            for (int j = i + 1; j < gN; j++){
                double dx = hp[j].x - hp[i].x, dy = hp[j].y - hp[i].y, dz = hp[j].z - hp[i].z;
                M.PE -= (double)G_DIAL*m*hm[j].w/sqrt(dx*dx + dy*dy + dz*dz + 1e-6);
            }
        }
    } else {
        CUDA_CHECK(cudaMemsetAsync(dMet, 0, 8*sizeof(unsigned long long), simS));
        kMeters<<<(gN + 255)/256, 256, 0, simS>>>(dPos[gPub], dMom,
                gSolver == SOLV_PM ? dReal : nullptr, gN, dMet);
        unsigned long long h[8];
        CUDA_CHECK(cudaStreamSynchronize(simS));
        CUDA_CHECK(cudaMemcpy(h, dMet, sizeof h, cudaMemcpyDeviceToHost));
        M.KE = orrery::fixed_decode((long long)h[0]);
        M.PE = orrery::fixed_decode((long long)h[1]);
        M.Px = orrery::fixed_decode((long long)h[2]);
        M.Py = orrery::fixed_decode((long long)h[3]);
        M.Pz = orrery::fixed_decode((long long)h[4]);
        M.Lz = orrery::fixed_decode((long long)h[5]);
        M.nRel = h[6]; M.nCls = h[7];
    }
    M.E = M.KE + M.PE;
    return M;
}

// ----------------------------------------------------------------------------
// BMP writer (zero-dependency screenshot)
// ----------------------------------------------------------------------------
static bool writeBMP(const char* path, const uchar4* px, int W, int H){
    FILE* f = fopen(path, "wb");
    if (!f) return false;
    int rowBytes = W*3, pad = (4 - rowBytes % 4) % 4, imgSize = (rowBytes + pad)*H;
    unsigned char hdr[54] = {0};
    hdr[0]='B'; hdr[1]='M';
    unsigned fileSize = 54 + imgSize;
    memcpy(hdr+2, &fileSize, 4); unsigned off = 54; memcpy(hdr+10, &off, 4);
    unsigned ihs = 40; memcpy(hdr+14, &ihs, 4);
    memcpy(hdr+18, &W, 4); memcpy(hdr+22, &H, 4);
    unsigned short planes = 1, bpp = 24; memcpy(hdr+26, &planes, 2); memcpy(hdr+28, &bpp, 2);
    memcpy(hdr+34, &imgSize, 4);
    fwrite(hdr, 1, 54, f);
    std::vector<unsigned char> row(rowBytes + pad, 0);
    for (int y = H - 1; y >= 0; y--){
        for (int x = 0; x < W; x++){
            uchar4 c = px[y*W + x];
            row[x*3+0] = c.z; row[x*3+1] = c.y; row[x*3+2] = c.x;
        }
        fwrite(row.data(), 1, rowBytes + pad, f);
    }
    fclose(f);
    return true;
}

// ----------------------------------------------------------------------------
// Allocation + init
// ----------------------------------------------------------------------------
static void allocAll(){
    int RW = gOW*gSS, RH = gOH*gSS;
    CUDA_CHECK(cudaMalloc(&dPos[0], (size_t)gN*sizeof(float4)));
    CUDA_CHECK(cudaMalloc(&dPos[1], (size_t)gN*sizeof(float4)));
    CUDA_CHECK(cudaMalloc(&dMom,    (size_t)gN*sizeof(float4)));
    CUDA_CHECK(cudaMalloc(&dTau,    (size_t)gN*sizeof(float)));
    CUDA_CHECK(cudaMalloc(&dRegime, (size_t)gN*sizeof(unsigned)));
    CUDA_CHECK(cudaMalloc(&dHdr,    (size_t)RW*RH*sizeof(float4)));
    mipW[0] = RW; mipH[0] = RH;
    for (int m = 1; m <= 6; m++){
        mipW[m] = (mipW[m-1] + 1)/2; mipH[m] = (mipH[m-1] + 1)/2;
        CUDA_CHECK(cudaMalloc(&dMip[m], (size_t)mipW[m]*mipH[m]*sizeof(float4)));
    }
    CUDA_CHECK(cudaMalloc(&dOut,  (size_t)gOW*gOH*sizeof(uchar4)));
    CUDA_CHECK(cudaMalloc(&dHist, 256*sizeof(unsigned)));
    CUDA_CHECK(cudaMalloc(&dGlyphs, 512*sizeof(Glyph)));
    int lo, hi;
    CUDA_CHECK(cudaDeviceGetStreamPriorityRange(&lo, &hi));
    CUDA_CHECK(cudaStreamCreateWithPriority(&simS, cudaStreamNonBlocking, lo));
    CUDA_CHECK(cudaStreamCreateWithPriority(&prsS, cudaStreamNonBlocking, hi));
    CUDA_CHECK(cudaEventCreateWithFlags(&tickDone, cudaEventDisableTiming));
    // gravity buffers + plans
    CUDA_CHECK(cudaMalloc(&dAcc,    (size_t)gN*sizeof(float4)));
    CUDA_CHECK(cudaMalloc(&dPosErr, (size_t)gN*sizeof(float4)));
    CUDA_CHECK(cudaMemset(dPosErr, 0, (size_t)gN*sizeof(float4)));
    CUDA_CHECK(cudaMemset(dAcc, 0, (size_t)gN*sizeof(float4)));
    CUDA_CHECK(cudaMalloc(&dFix,  (size_t)PMN*PMN*PMN*sizeof(unsigned long long)));
    CUDA_CHECK(cudaMalloc(&dReal, (size_t)PMN*PMN*PMN*sizeof(float)));
    CUDA_CHECK(cudaMalloc(&dSpec, (size_t)PMN*PMN*PMNZC*sizeof(cufftComplex)));
    CUDA_CHECK(cudaMalloc(&dFor,  (size_t)PMN*PMN*PMN*sizeof(float4)));
    CUDA_CHECK(cudaMalloc(&dMet,  8*sizeof(unsigned long long)));
    CUFFT_CHECK(cufftPlan3d(&planF, PMN, PMN, PMN, CUFFT_R2C));
    CUFFT_CHECK(cufftPlan3d(&planI, PMN, PMN, PMN, CUFFT_C2R));
    CUFFT_CHECK(cufftSetStream(planF, simS));
    CUFFT_CHECK(cufftSetStream(planI, simS));

    // scenario init (contracts/newton.contract.md)
    switch (gScenario){
    case SC_GALAXY:
        gMTot = (double)gN * 1.0;
        kInit<<<(gN + 255)/256, 256>>>(dPos[0], dMom, dTau, dRegime, gN, gSeed, (float)gMTot);
        gSolver = SOLV_PM;
        break;
    case SC_CLOUD:
        gMTot = (double)gN * 1.0;
        kInitCloud<<<(gN + 255)/256, 256>>>(dPos[0], dMom, dTau, dRegime, gN, gSeed);
        gSolver = SOLV_PM;
        break;
    case SC_KEPLER: {
        gN = 2; gMTot = 10001.0;
        double mu = (double)G_DIAL*10001.0, rap = 240.0, a = 150.0;
        double vap = sqrt(mu*(2.0/rap - 1.0/a));
        float4 hp[2] = { make_float4((float)(-rap/10001.0), 0, 0, 6500.0f),
                         make_float4((float)( rap*10000.0/10001.0), 0, 0, 9000.0f) };
        float4 hm[2] = { make_float4(0, (float)(-vap/10001.0*1e4), 0, 1e4f),
                         make_float4(0, (float)( vap*10000.0/10001.0), 0, 1.0f) };
        float ht[2] = {0, 0}; unsigned hr[2] = {0x02u, 0x02u};
        CUDA_CHECK(cudaMemcpy(dPos[0], hp, sizeof hp, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(dMom, hm, sizeof hm, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(dTau, ht, sizeof ht, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(dRegime, hr, sizeof hr, cudaMemcpyHostToDevice));
        gSolver = SOLV_TINY;
    } break;
    case SC_THREEBODY: {
        gN = 3; gMTot = 3.0;
        double sf = sqrt((double)G_DIAL);              // velocity rescale for G != 1
        float4 hp[3] = { make_float4( 0.97000436f, -0.24308753f, 0, 6000.0f),
                         make_float4(-0.97000436f,  0.24308753f, 0, 8000.0f),
                         make_float4(0, 0, 0, 11000.0f) };
        float4 hm[3] = { make_float4((float)(0.4662036850*sf), (float)(0.4323657300*sf), 0, 1.0f),
                         make_float4((float)(0.4662036850*sf), (float)(0.4323657300*sf), 0, 1.0f),
                         make_float4((float)(-0.93240737*sf),  (float)(-0.86473146*sf),  0, 1.0f) };
        float ht[3] = {0, 0, 0}; unsigned hr[3] = {0x02u, 0x02u, 0x02u};
        CUDA_CHECK(cudaMemcpy(dPos[0], hp, sizeof hp, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(dMom, hm, sizeof hm, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(dTau, ht, sizeof ht, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(dRegime, hr, sizeof hr, cudaMemcpyHostToDevice));
        gSolver = SOLV_TINY;
    } break;
    }
    CUDA_CHECK(cudaMemcpy(dPos[1], dPos[0], (size_t)gN*sizeof(float4),
                          cudaMemcpyDeviceToDevice));
    CUDA_CHECK(cudaDeviceSynchronize());
    forcePass();                                      // establish dAcc = a(x0); Phi for meters
    CUDA_CHECK(cudaStreamSynchronize(simS));
    gM0 = metersNow();
}

// ----------------------------------------------------------------------------
// Win32 + GL presentation (D-012 P0: thin blit of the finished image)
// ----------------------------------------------------------------------------
static HWND gWnd; static HDC gDC; static HGLRC gRC; static GLuint gTex;
static cudaGraphicsResource* gRes = 0;
typedef BOOL (WINAPI* PFNSWAPINT)(int);
static PFNSWAPINT pSwapInterval = 0;

static LRESULT CALLBACK WndProc(HWND h, UINT msg, WPARAM wp, LPARAM lp){
    switch (msg){
    case WM_CLOSE: case WM_DESTROY: gRunning = false; return 0;
    case WM_LBUTTONDOWN:
        gDrag = 1; gLastX = GET_X_LPARAM(lp); gLastY = GET_Y_LPARAM(lp);
        SetCapture(h); break;
    case WM_LBUTTONUP: gDrag = 0; ReleaseCapture(); break;
    case WM_MOUSEMOVE:
        if (gDrag){
            int mx = GET_X_LPARAM(lp), my = GET_Y_LPARAM(lp);
            cAz.tgt += (mx - gLastX)*0.28f;
            cEl.tgt = fminf(fmaxf(cEl.tgt + (my - gLastY)*0.20f, -88.0f), 88.0f);
            gLastX = mx; gLastY = my;
        }
        break;
    case WM_MOUSEWHEEL: {
        int z = GET_WHEEL_DELTA_WPARAM(wp);
        cDist.tgt += -0.105f*(z/120.0f);              // exp zoom, x0.90 per notch
        cDist.tgt = fminf(fmaxf(cDist.tgt, logf(40.0f)), logf(1200.0f));
    } break;
    case WM_KEYDOWN: {
        WPARAM k = wp;
        if (k == VK_ESCAPE) gRunning = false;
        else if (k >= '1' && k <= '3') applyPreset((int)(k - '1'));
        else if (k == 'C') gMode ^= 1;
        else if (k == 'T') gTonemap ^= 1;
        else if (k == 'B') gBloomOn ^= 1;
        else if (k == 'A') gAutoExp ^= 1;
        else if (k == 'O') gAutoOrbit ^= 1;
        else if (k == 'P') gPaused ^= 1;
        else if (k == VK_F1) gHud ^= 1;
        else if (k == 'V'){ gVsync ^= 1; if (pSwapInterval) pSwapInterval(gVsync); }
        else if (k == VK_OEM_PLUS  || k == VK_ADD)      gEvOffset = fminf(gEvOffset + 0.25f,  8.0f);
        else if (k == VK_OEM_MINUS || k == VK_SUBTRACT) gEvOffset = fmaxf(gEvOffset - 0.25f, -8.0f);
    } break;
    }
    return DefWindowProcA(h, msg, wp, lp);
}

static void initWindowGL(){
    HINSTANCE hi = GetModuleHandleA(NULL);
    WNDCLASSA wc = {};
    wc.style = CS_OWNDC; wc.lpfnWndProc = WndProc; wc.hInstance = hi;
    wc.hCursor = LoadCursorA(NULL, MAKEINTRESOURCEA(32512));
    wc.lpszClassName = "TinyUniverseWnd";
    RegisterClassA(&wc);
    DWORD style = WS_OVERLAPPED | WS_CAPTION | WS_SYSMENU | WS_MINIMIZEBOX;
    RECT r = {0, 0, gOW, gOH};
    AdjustWindowRect(&r, style, FALSE);
    gWnd = CreateWindowA("TinyUniverseWnd", "TINY UNIVERSE - first light",
                         style | WS_VISIBLE, CW_USEDEFAULT, CW_USEDEFAULT,
                         r.right - r.left, r.bottom - r.top, NULL, NULL, hi, NULL);
    gDC = GetDC(gWnd);
    PIXELFORMATDESCRIPTOR pfd = {};
    pfd.nSize = sizeof(pfd); pfd.nVersion = 1;
    pfd.dwFlags = PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL | PFD_DOUBLEBUFFER;
    pfd.iPixelType = PFD_TYPE_RGBA; pfd.cColorBits = 32;
    SetPixelFormat(gDC, ChoosePixelFormat(gDC, &pfd), &pfd);
    gRC = wglCreateContext(gDC);
    wglMakeCurrent(gDC, gRC);
    pSwapInterval = (PFNSWAPINT)wglGetProcAddress("wglSwapIntervalEXT");
    if (pSwapInterval) pSwapInterval(gVsync);
    glGenTextures(1, &gTex);
    glBindTexture(GL_TEXTURE_2D, gTex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, gOW, gOH, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
    CUDA_CHECK(cudaGraphicsGLRegisterImage(&gRes, gTex, GL_TEXTURE_2D,
               cudaGraphicsRegisterFlagsWriteDiscard));
}

static void present(){
    CUDA_CHECK(cudaGraphicsMapResources(1, &gRes, prsS));
    cudaArray_t arr;
    CUDA_CHECK(cudaGraphicsSubResourceGetMappedArray(&arr, gRes, 0, 0));
    CUDA_CHECK(cudaMemcpy2DToArrayAsync(arr, 0, 0, dOut, gOW*sizeof(uchar4),
               gOW*sizeof(uchar4), gOH, cudaMemcpyDeviceToDevice, prsS));
    CUDA_CHECK(cudaGraphicsUnmapResources(1, &gRes, prsS));
    CUDA_CHECK(cudaStreamSynchronize(prsS));
    glClear(GL_COLOR_BUFFER_BIT);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, gTex);
    glBegin(GL_QUADS);
    glTexCoord2f(0, 1); glVertex2f(-1, -1);
    glTexCoord2f(1, 1); glVertex2f( 1, -1);
    glTexCoord2f(1, 0); glVertex2f( 1,  1);
    glTexCoord2f(0, 0); glVertex2f(-1,  1);
    glEnd();
    SwapBuffers(gDC);
}

// ----------------------------------------------------------------------------
// main
// ----------------------------------------------------------------------------
int main(int argc, char** argv){
    SetProcessDPIAware();
    const char* shotPath = 0;
    int benchFrames = 0, shotFrames = 180;
    long long warmTicks = 0, runSteps = 0;
    bool jsonMode = false, goldenMode = false;
    for (int i = 1; i < argc; i++){
        std::string a = argv[i];
        if      (a == "--n" && i+1 < argc)      gN = atoi(argv[++i]);
        else if (a == "--seed" && i+1 < argc)   gSeed = strtoull(argv[++i], 0, 10);
        else if (a == "--w" && i+1 < argc)      gOW = atoi(argv[++i]);
        else if (a == "--h" && i+1 < argc)      gOH = atoi(argv[++i]);
        else if (a == "--ssaa")                 gSS = 2;
        else if (a == "--phys")                 gMode = 1;
        else if (a == "--aces")                 gTonemap = 1;
        else if (a == "--bench" && i+1 < argc)  benchFrames = atoi(argv[++i]);
        else if (a == "--shot" && i+1 < argc)   shotPath = argv[++i];
        else if (a == "--frames" && i+1 < argc) shotFrames = atoi(argv[++i]);
        else if (a == "--warm" && i+1 < argc)   warmTicks = atoll(argv[++i]);
        else if (a == "--steps" && i+1 < argc)  runSteps = atoll(argv[++i]);
        else if (a == "--json")                 jsonMode = true;
        else if (a == "--golden")               goldenMode = true;
        else if (a == "--scenario" && i+1 < argc){
            std::string s = argv[++i];
            if      (s == "galaxy")    gScenario = SC_GALAXY;
            else if (s == "kepler")    gScenario = SC_KEPLER;
            else if (s == "threebody") gScenario = SC_THREEBODY;
            else if (s == "cloud")     gScenario = SC_CLOUD;
            else { fprintf(stderr, "error: unknown scenario\n"); return 2; }
        }
        else { fprintf(stderr, "usage: tinyuniverse [--scenario galaxy|kepler|threebody|cloud] "
                       "[--w W --h H] [--n N] [--seed S] [--ssaa] [--phys] [--aces] [--bench F] "
                       "[--shot PATH.bmp [--frames N] [--warm T]] [--steps S --json] [--golden]\n");
               return 2; }
    }
    gOW &= ~1; gOH &= ~1;
    if (gN < 1000) gN = 1000;
    if (goldenMode) gSeed = 20260711ull;              // frozen seed at init time

    CUDA_CHECK(cudaSetDevice(0));
    cudaDeviceProp prop;
    CUDA_CHECK(cudaGetDeviceProperties(&prop, 0));
    printf("TINY UNIVERSE first light | CUDA: %s (%d SMs) | N=%d | %dx%d ss=%d\n",
           prop.name, prop.multiProcessorCount, gN, gOW, gOH, gSS);

    allocAll();
    RP rp = {};
    rp.RW = gOW*gSS; rp.RH = gOH*gSS; rp.OW = gOW; rp.OH = gOH; rp.ss = gSS;
    rp.bloomAmt = 0.035f; rp.satCine = 1.20f;

    if (jsonMode || goldenMode){                      // headless instrument face (newton contract)
        if (goldenMode){                              // frozen golden params (seed forced pre-init)
            const long long GS[4] = {10000, 1000000, 1000000, 12000};
            runSteps = GS[gScenario];
        }
        if (runSteps <= 0){
            const long long GS[4] = {10000, 1000000, 1000000, 12000};
            runSteps = GS[gScenario];
        }
        long long done = 0, nextPct = 10;
        const int batch = (gSolver == SOLV_TINY) ? 20000 : 250;
        while (done < runSteps){
            int nt = (int)((runSteps - done) < batch ? (runSteps - done) : batch);
            runTicks(nt);
            done += nt;
            if (done*100 >= runSteps*nextPct){
                CUDA_CHECK(cudaStreamSynchronize(simS));
                fprintf(stderr, "  %lld%%\r", nextPct); nextPct += 10;
            }
        }
        CUDA_CHECK(cudaStreamSynchronize(simS));
        Met M1 = metersNow();

        // declared state hash: raw pos/mom/tau/regime bytes at the final tick
        std::string bytes;
        bytes.resize((size_t)gN*(16 + 16 + 4 + 4));
        char* w = &bytes[0];
        CUDA_CHECK(cudaMemcpy(w, dPos[gPub], (size_t)gN*16, cudaMemcpyDeviceToHost)); w += (size_t)gN*16;
        CUDA_CHECK(cudaMemcpy(w, dMom,       (size_t)gN*16, cudaMemcpyDeviceToHost)); w += (size_t)gN*16;
        CUDA_CHECK(cudaMemcpy(w, dTau,       (size_t)gN*4,  cudaMemcpyDeviceToHost)); w += (size_t)gN*4;
        CUDA_CHECK(cudaMemcpy(w, dRegime,    (size_t)gN*4,  cudaMemcpyDeviceToHost));
        std::string shash = orrery::blake2b_hex(bytes, 32);

        double pnorm = sqrt(2.0*gMTot*(gM0.KE > 1e-12 ? gM0.KE : fabs(gM0.E) + 1e-12));
        double de = fabs((M1.E - gM0.E)/(fabs(gM0.E) > 1e-12 ? fabs(gM0.E) : 1.0));
        double dp = sqrt((M1.Px-gM0.Px)*(M1.Px-gM0.Px) + (M1.Py-gM0.Py)*(M1.Py-gM0.Py)
                       + (M1.Pz-gM0.Pz)*(M1.Pz-gM0.Pz)) / (pnorm > 1e-12 ? pnorm : 1.0);
        double lnorm = fabs(gM0.Lz) > pnorm*BOXL/8.0 ? fabs(gM0.Lz) : pnorm*BOXL/8.0;
        double dl = fabs(M1.Lz - gM0.Lz)/lnorm;
        double maxr = 0;
        if (gSolver == SOLV_TINY){
            std::vector<float4> hp(gN);
            CUDA_CHECK(cudaMemcpy(hp.data(), dPos[gPub], gN*sizeof(float4), cudaMemcpyDeviceToHost));
            for (int q = 0; q < gN; q++){
                double r = sqrt((double)hp[q].x*hp[q].x + (double)hp[q].y*hp[q].y + (double)hp[q].z*hp[q].z);
                if (r > maxr) maxr = r;
            }
        }
        bool pass = false;
        std::string gates;
        using namespace orrery;
        auto B = [](bool b){ return std::string(b ? "true" : "false"); };
        switch (gScenario){
        case SC_KEPLER:
            pass = (de < 2e-3) && (dl < 2e-3);
            gates = "\"de_lt_2e-3\":" + B(de < 2e-3) + ",\"dl_lt_2e-3\":" + B(dl < 2e-3);
            break;
        case SC_THREEBODY:
            pass = (de < 2e-3) && (maxr < 5.0);
            gates = "\"de_lt_2e-3\":" + B(de < 2e-3) + ",\"max_r\":" + fmt6(maxr)
                  + ",\"bounded_r_lt_5\":" + B(maxr < 5.0);
            break;
        case SC_CLOUD:
            pass = (de < 0.08) && (dp < 1e-3);
            gates = "\"de_lt_8e-2\":" + B(de < 0.08) + ",\"dp_lt_1e-3\":" + B(dp < 1e-3);
            break;
        case SC_GALAXY:
            pass = (de < 0.02) && (dp < 1e-3);
            gates = "\"de_lt_2e-2\":" + B(de < 0.02) + ",\"dp_lt_1e-3\":" + B(dp < 1e-3);
            break;
        }
        const char* solvName = gSolver == SOLV_PM ? "pm" : (gSolver == SOLV_TINY ? "tiny" : "direct");
        std::string body =
              std::string("\"seed\":") + fmti((long long)gSeed)
            + ",\"params\":{\"scenario\":\"" + SC_NAME[gScenario] + "\",\"n\":" + fmti(gN)
            + ",\"steps\":" + fmti(runSteps) + ",\"solver\":\"" + solvName + "\""
            + ",\"c\":" + fmt6(20.0) + ",\"hbar\":" + fmt6(0.5) + ",\"G\":" + fmt6(0.002)
            + ",\"dt\":" + fmt6(1.0/240.0) + ",\"L_box\":" + fmt6(512.0) + "}"
            + ",\"result\":{\"state_b2b\":\"" + shash + "\""
            + ",\"e0\":" + fmt6(gM0.E) + ",\"e1\":" + fmt6(M1.E)
            + ",\"de_rel\":" + fmt6(de) + ",\"p_drift\":" + fmt6(dp) + ",\"l_drift\":" + fmt6(dl)
            + ",\"n_rel\":" + fmti((long long)M1.nRel) + ",\"n_classical\":" + fmti((long long)M1.nCls) + "}"
            + ",\"gates\":{" + gates + "}"
            + ",\"verdict\":\"" + (pass ? "pass" : "fail") + "\"";
        std::string env = full_envelope("tinyuniverse", "0.2.0", body,
                                        "newton scenario run; visuals non-declared");
        if (goldenMode)
            return golden_check(SC_NAME[gScenario], declared_object(body), env);
        printf("%s\n", env.c_str());
        fprintf(stderr, "declared[%s] %s\n", SC_NAME[gScenario],
                blake2b_hex(declared_object(body)).c_str());
        return pass ? 0 : 1;
    }

    if (shotPath){                                    // headless render-to-BMP
        applyPreset(0);
        cAz.cur = cAz.tgt; cEl.cur = cEl.tgt; cDist.cur = cDist.tgt;
        while (warmTicks > 0){                        // evolve before the beauty pass
            int nt = (int)(warmTicks < 250 ? warmTicks : 250);
            runTicks(nt);
            warmTicks -= nt;
        }
        CUDA_CHECK(cudaStreamSynchronize(simS));
        for (int f = 0; f < shotFrames; f++){
            if (!gPaused) runTicks(4);
            rp.mode = gMode; rp.tonemap = gTonemap; rp.bloomOn = gBloomOn; rp.frame = gFrame++;
            buildCamera(rp);
            renderFrame(rp);
            CUDA_CHECK(cudaStreamSynchronize(prsS));
        }
        drawHud(rp);
        CUDA_CHECK(cudaStreamSynchronize(prsS));
        std::vector<uchar4> hOut((size_t)gOW*gOH);
        CUDA_CHECK(cudaMemcpy(hOut.data(), dOut, (size_t)gOW*gOH*sizeof(uchar4),
                              cudaMemcpyDeviceToHost));
        if (!writeBMP(shotPath, hOut.data(), gOW, gOH)){
            fprintf(stderr, "error: cannot write %s\n", shotPath); return 2;
        }
        printf("shot written: %s (%d frames settled, expo %.4g)\n", shotPath, shotFrames, gExposure);
        return 0;
    }

    initWindowGL();
    if (benchFrames){ gVsync = 0; if (pSwapInterval) pSwapInterval(0); }
    applyPreset(0);
    cAz.cur = cAz.tgt; cEl.cur = cEl.tgt; cDist.cur = cDist.tgt;

    LARGE_INTEGER freq, tPrev, tFps;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&tPrev); tFps = tPrev;
    int fpsN = 0, benchN = 0, warm = 0;
    double benchSum = 0, benchMin = 1e9;

    printf("keys: drag orbit | wheel zoom | 1-3 presets | O orbit | P pause | C cine/phys\n"
           "      T AgX/ACES | B bloom | A auto-exp | +/- EV | F1 HUD | V vsync | Esc\n");

    while (gRunning){
        MSG msg;
        while (PeekMessageA(&msg, NULL, 0, 0, PM_REMOVE)){
            if (msg.message == WM_QUIT) gRunning = false;
            TranslateMessage(&msg); DispatchMessageA(&msg);
        }
        if (!gRunning) break;

        LARGE_INTEGER now; QueryPerformanceCounter(&now);
        float dtF = (float)((double)(now.QuadPart - tPrev.QuadPart)/freq.QuadPart);
        tPrev = now;
        if (dtF > 0.1f) dtF = 0.1f;

        if (gAutoOrbit && !gDrag) cAz.tgt += dtF*2.2f;
        springStep(cAz, dtF); springStep(cEl, dtF); springStep(cDist, dtF);

        if (!gPaused){
            int nt = (int)(dtF/DT + 0.5f); nt = nt < 1 ? 1 : (nt > 8 ? 8 : nt);
            runTicks(nt);
        }
        if (gFrame % 120 == 60){                      // HUD conservation meters (~2 s cadence)
            Met m = metersNow();
            gDE = (m.E - gM0.E)/(fabs(gM0.E) > 1e-12 ? fabs(gM0.E) : 1.0);
            double pn = sqrt(2.0*gMTot*(gM0.KE > 1e-12 ? gM0.KE : 1.0));
            gDP = sqrt((m.Px-gM0.Px)*(m.Px-gM0.Px) + (m.Py-gM0.Py)*(m.Py-gM0.Py)
                     + (m.Pz-gM0.Pz)*(m.Pz-gM0.Pz))/pn;
        }
        rp.mode = gMode; rp.tonemap = gTonemap; rp.bloomOn = gBloomOn; rp.frame = gFrame++;
        buildCamera(rp);
        renderFrame(rp);
        drawHud(rp);
        present();

        fpsN++;
        double dtf = (double)(now.QuadPart - tFps.QuadPart)/freq.QuadPart;
        if (dtf >= 0.5){
            gFps = (float)(fpsN/dtf); gMs = (float)(1000.0*dtf/fpsN);
            fpsN = 0; tFps = now;
            char title[256];
            snprintf(title, sizeof title,
                     "TINY UNIVERSE - first light | N=%d | %s %s | %.0f fps",
                     gN, gMode ? "PHYSICAL" : "CINEMATIC", gTonemap ? "ACES" : "AGX", gFps);
            SetWindowTextA(gWnd, title);
        }
        if (benchFrames){
            if (warm < 60) warm++;
            else {
                benchSum += dtF; benchN++;
                if (dtF > 1e-9 && 1.0/dtF < benchMin) benchMin = 1.0/dtF;
                if (benchN >= benchFrames) gRunning = false;
            }
        }
    }

    if (benchFrames && benchN){
        double avg = benchN/benchSum;
        printf("BENCH: %d frames | avg %.1f fps (%.2f ms) | min %.1f fps | N=%d %dx%d ss=%d\n",
               benchN, avg, 1000.0*benchSum/benchN, benchMin, gN, gOW, gOH, gSS);
    }
    cudaDeviceSynchronize();
    return 0;
}
