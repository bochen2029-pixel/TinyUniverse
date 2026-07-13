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
#include <cuda_fp16.h>
#include <complex>
#include "../core/lib/rng.cuh"
#include "../core/lib/reduce.cuh"
#include "../core/lib/envelope.h"

#define PI_F 3.14159265358979f

// --- M2 newton: dials + PM constants (contracts/newton.contract.md) ---
#define G_DIAL   2.0e-3f
#define C_LIGHT  20.0f
#define PMN      128
#define PMNZC    (PMN/2 + 1)
#define CELL     4.0f
#define BOXL     512.0f
#define REL_V2   1.0f            // (0.05c)^2 = 1 su^2/s^2 — regime threshold
enum Scenario { SC_GALAXY = 0, SC_KEPLER, SC_THREEBODY, SC_CLOUD,
                SC_MERGER, SC_ECHO, SC_RATCHET, SC_DETECTOR,
                SC_KEPREL, SC_CLOCKS, SC_PHOTONS,
                SC_COLLAPSE, SC_ISCO, SC_HAWKING,
                SC_DOUBLESLIT, SC_TUNNELING, SC_SHOQ,
                SC_CIRCUMNAV, SC_EXPAND, SC_BUBBLES };
enum Solver   { SOLV_PM = 0, SOLV_DIRECT, SOLV_TINY, SOLV_NONE };
static const char* SC_NAME[20] = {"galaxy", "kepler", "threebody", "cloud",
                                  "merger", "echo", "ratchet", "detector",
                                  "keprel", "clocks", "photons",
                                  "collapse", "isco", "hawking",
                                  "doubleslit", "tunneling", "shoq",
                                  "circumnav", "expand", "bubbles"};
// planck contract: psi engine (2D 256^2 over a 64x64 su window)
#define QN   256
#define QL   64.0f
#define QDX  (QL/QN)
// cosmos contract: roaming 3D bubbles (64^3 over a 16^3 su window, dx = 0.25)
#define QB    64
#define QB3   (QB*QB*QB)
#define QBL   16.0f
#define QBDX  (QBL/QB)
#define NBUB  16
// gargantua contract: BH constants (dial-derived)
#define NBH_MAX  8
#define RS_PER_M 1.0e-5f                   // 2G/c^2
#define M_FORM   1.0e5f                    // r_s = CELL/4 (the N1 N_BH crossover)
#define M_POP    50.0f
#define REG_DEAD 0x100u
// arrow contract: ratchet constants (N6/ORRERY), record thresholds
#define RATCHET_PDOWN 0.4166666666666667   // p/(p+(1-p)rho), p=0.3 rho=0.6
#define REC_THRESH 16u
#define NRECS 3000000

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
    // M5: screen-space BH lenses (render-res pixels; visuals non-declared)
    int   nLens;
    float lensX[4], lensY[4], lensTE2[4], lensShadow[4];
    // M7 cosmos: projection + light-history (visuals non-declared)
    int   proj;            // 0 perspective, 1 stereographic little planet
    float lpR;             // little-planet horizon radius in render px
    int   histN, histHead, histDepth;   // valid snapshots / ring head / ring size
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

// M7 cosmos: Zel'dovich growing-mode lattice (expand scenario; cosmos contract).
// 100^3 lattice, single mode k1 along x, displacement psi = Ad*sin(k1*q),
// growing-mode peculiar velocity v = H0*psi at a = 1 (EdS f = 1). No RNG.
__global__ void kInitZeldovich(float4* pos, float4* mom, float* tau, unsigned* regime,
                               float Ad, float k1, float H0){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= 1000000) return;
    int ix = i % 100, iy = (i/100) % 100, iz = i/10000;
    float qx = -256.0f + (ix + 0.5f)*5.12f;
    float qy = -256.0f + (iy + 0.5f)*5.12f;
    float qz = -256.0f + (iz + 0.5f)*5.12f;
    float psi = Ad*sinf(k1*qx);
    pos[i] = make_float4(qx + psi, qy, qz, 5600.0f);   // uniform fog; structure will heat it visually
    mom[i] = make_float4(H0*psi, 0, 0, 1.0f);          // p = a^2 x_dot = v at a=1, m=1
    tau[i] = 0.0f; regime[i] = 0x02u;
}

__global__ void kInitCloud(float4* pos, float4* mom, float* tau, unsigned* regime,
                           int N, unsigned long long seed, float R, float sigFrac){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= N) return;
    using namespace orrery;
    float gx = (float)counter_gauss(seed, i, 0, 1);
    float gy = (float)counter_gauss(seed, i, 1, 1);
    float gz = (float)counter_gauss(seed, i, 2, 1);
    float il = rsqrtf(gx*gx + gy*gy + gz*gz + 1e-12f);
    float r = R * cbrtf((float)counter_uniform(seed, i, 3, 1));
    pos[i] = make_float4(r*gx*il, r*gy*il, r*gz*il, 3800.0f);
    // sub-virial support: a truly cold pressureless collapse is singular —
    // unresolvable at any fixed dt (D-014); collapse scenario runs colder (D-017)
    float sigv = sigFrac*sqrtf(G_DIAL*(float)N/R);
    mom[i] = make_float4(sigv*(float)counter_gauss(seed, i, 4, 1),
                         sigv*(float)counter_gauss(seed, i, 5, 1),
                         sigv*(float)counter_gauss(seed, i, 6, 1), 1.0f);
    tau[i] = 0.0f; regime[i] = 0x02u;
}

__device__ float phiAt(const float* phi, float4 p);   // defined with the PM pipeline below

// --- KDK pieces, relativistic (einstein contract): p = gamma*m*v state,
//     v = p/sqrt(m^2 + |p|^2/c^2) (cancellation-free), photons v = c*p_hat ---
__device__ __forceinline__ float3 velOf(float4 m){
    float p2 = m.x*m.x + m.y*m.y + m.z*m.z;
    if (m.w > 0.0f){
        float ir = rsqrtf(m.w*m.w + p2/(C_LIGHT*C_LIGHT));
        return make_float3(m.x*ir, m.y*ir, m.z*ir);
    }
    float ir = C_LIGHT*rsqrtf(fmaxf(p2, 1e-30f));
    return make_float3(m.x*ir, m.y*ir, m.z*ir);
}
__global__ void kKick(float4* mom, const float4* acc, const unsigned* regime,
                      int N, float hdt){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= N) return;
    if (regime[i] & REG_DEAD) return;                 // absorbed: frozen
    float4 m = mom[i], a = acc[i];                    // a.xyz = g, a.w = Phi
    const float ic2 = 1.0f/(C_LIGHT*C_LIGHT);
    if (m.w > 0.0f){
        // dp = m*g*dt (D-016: the 1PN field term was withdrawn after measurement
        // contradicted the 7-pi superposition claim; strong-field precession
        // returns with M5's pseudo-potential oracles. Q-006 holds the derivation.)
        (void)ic2;
        m.x += m.w*a.x*hdt; m.y += m.w*a.y*hdt; m.z += m.w*a.z*hdt;
    } else {
        float pm = sqrtf(m.x*m.x + m.y*m.y + m.z*m.z);
        float f = 2.0f*pm/C_LIGHT*hdt;                // weak-field photon bending (factor 2)
        m.x += f*a.x; m.y += f*a.y; m.z += f*a.z;
    }
    mom[i] = m;
}
__global__ void kDriftK(const float4* pin, float4* pout, float4* perr,
                        const float4* mom, float* tau, unsigned* regime,
                        const float* phi, int N, float dt, float invc2, float ds,
                        const unsigned* partBub){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= N) return;
    if (regime[i] & REG_DEAD){                        // absorbed: park out of the box
        pout[i] = make_float4(1.0e6f + (float)i*1.0e-3f, 0, 0, pin[i].w);
        return;
    }
    if (partBub && partBub[i]){                       // bubbled: classical drift + tau suspend
        pout[i] = pin[i];                             // (declared tau-gap, cosmos contract)
        return;
    }
    float4 p = pin[i], m = mom[i], e = perr[i];
    float3 v = velOf(m);
    float ddt = dt*ds;                                // ds = 1/a^2 in expand mode, else exactly 1
    float3 d = make_float3(v.x*ddt, v.y*ddt, v.z*ddt);
    float y, t;
    y = d.x - e.x; t = p.x + y; e.x = (t - p.x) - y; p.x = t;
    y = d.y - e.y; t = p.y + y; e.y = (t - p.y) - y; p.y = t;
    y = d.z - e.z; t = p.z + y; e.z = (t - p.z) - y; p.z = t;
    // M7: canonical torus coordinates [-256, 256). +/-512 is exact fp32 near the
    // seam (power of two), so wrapping is bit-reversible; max speed c => at most
    // one crossing per tick. Kahan compensation survives an exact-op shift.
    if      (p.x >=  256.0f) p.x -= 512.0f; else if (p.x < -256.0f) p.x += 512.0f;
    if      (p.y >=  256.0f) p.y -= 512.0f; else if (p.y < -256.0f) p.y += 512.0f;
    if      (p.z >=  256.0f) p.z -= 512.0f; else if (p.z < -256.0f) p.z += 512.0f;
    pout[i] = p; perr[i] = e;
    float v2 = v.x*v.x + v.y*v.y + v.z*v.z;
    float ph = phi ? phiAt(phi, p) : 0.0f;            // tick-lagged Phi (declared)
    if (m.w > 0.0f)                                   // photons structurally ageless (einstein
        tau[i] += dt * sqrtf(fmaxf(1.0f + (2.0f*ph - v2)*invc2, 0.0f));   // contract; the fp32
                                                      // clamp leaks ~2e-4/tick at v = c — D-019
    unsigned sp = (m.w == 0.0f) ? (0x20u | 0x04u) : ((v2 > REL_V2) ? 0x04u : 0x02u);
    regime[i] = (regime[i] & 0xC0u) | sp;
}

// --- PM pipeline: fixed-point CIC deposit -> cuFFT -> Green -> force grid -> gather ---
__global__ void kZeroFix(unsigned long long* g, int n){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < n) g[i] = 0ull;
}
__device__ __forceinline__ int wrapc(int c){ return c & (PMN - 1); }
// M7 cosmos: minimum-image separation on the 3-torus. 512 = 2^9, so the scale,
// round, and shift are all exact fp32 ops; for |d| < 256 the result is bit-identical d.
__host__ __device__ __forceinline__ float mimg(float d){ return d - 512.0f*roundf(d*(1.0f/512.0f)); }
// M7 cosmos: shared camera projection. proj 0 = perspective (fd = view depth),
// proj 1 = stereographic little planet from the nadir antipode (fd = true distance).
__device__ __forceinline__ bool projCam(const RP& rp, float rx, float ry, float rz,
                                        float& px, float& py, float& fd){
    float cx = rx*rp.rgt[0] + ry*rp.rgt[1] + rz*rp.rgt[2];
    float cy = rx*rp.upv[0] + ry*rp.upv[1] + rz*rp.upv[2];
    float cz = rx*rp.fwd[0] + ry*rp.fwd[1] + rz*rp.fwd[2];
    if (rp.proj == 0){
        if (cz < 1.0f) return false;
        px = ((cx/(cz*rp.tf*rp.aspect))*0.5f + 0.5f)*rp.RW - 0.5f;
        py = (0.5f - (cy/(cz*rp.tf))*0.5f)*rp.RH - 0.5f;
        fd = cz;
        return true;
    }
    float D = sqrtf(rx*rx + ry*ry + rz*rz + 1e-12f);
    float id = 1.0f/D;
    float cd = cz*id;                                 // cos(angle from nadir)
    if (cd < -0.985f) return false;                   // zenith blowup guard
    float s = 1.0f/(1.0f + cd);
    px = rp.RW*0.5f + (cx*id)*s*rp.lpR - 0.5f;
    py = rp.RH*0.5f - (cy*id)*s*rp.lpR - 0.5f;
    fd = D;
    return true;
}
__global__ void kDeposit(const float4* pos, const float4* mom, const unsigned* regime,
                         int N, unsigned long long* g){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= N) return;
    if (regime[i] & REG_DEAD) return;                 // absorbed mass lives in the BH now
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
__global__ void kGreen(cufftComplex* s, float pref){
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
    f *= pref;                                        // comoving Poisson x1/a (expand); 1 otherwise
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
__global__ void kGather(const float4* pos, const float4* F, const float* phi,
                        float4* acc, int N){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= N) return;
    float4 p = pos[i];
    float gx = (p.x + BOXL*0.5f)/CELL, gy = (p.y + BOXL*0.5f)/CELL, gz = (p.z + BOXL*0.5f)/CELL;
    int x0 = (int)floorf(gx), y0 = (int)floorf(gy), z0 = (int)floorf(gz);
    float fx = gx - x0, fy = gy - y0, fz = gz - z0;
    float3 a = make_float3(0, 0, 0);
    float pw = 0;
    for (int dz = 0; dz < 2; dz++)
        for (int dy = 0; dy < 2; dy++)
            for (int dx = 0; dx < 2; dx++){
                float w = (dx ? fx : 1.0f - fx)*(dy ? fy : 1.0f - fy)*(dz ? fz : 1.0f - fz);
                int c = (wrapc(x0 + dx)*PMN + wrapc(y0 + dy))*PMN + wrapc(z0 + dz);
                float4 f = F[c];
                a.x += w*f.x; a.y += w*f.y; a.z += w*f.z;
                pw += w*phi[c];
            }
    acc[i] = make_float4(a.x, a.y, a.z, pw);
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
    float pw = 0;
    for (int j = 0; j < N; j++){
        if (j == i) continue;
        float4 pj = pos[j];
        float dx = mimg(pj.x - pi.x), dy = mimg(pj.y - pi.y), dz = mimg(pj.z - pi.z);
        float r2 = dx*dx + dy*dy + dz*dz + 1e-6f;
        float ir = rsqrtf(r2);
        float f = G_DIAL*mom[j].w*ir*ir*ir;
        a.x += f*dx; a.y += f*dy; a.z += f*dz;
        pw -= G_DIAL*mom[j].w*ir;
    }
    acc[i] = make_float4(a.x, a.y, a.z, pw);
}

// --- tiny solver (N <= 32): one block, many ticks per launch, fp64 internal ---
// (D-014: tiny scenarios are the sim's precision references; declared state
//  round-trips through the fp32 frame-contract buffers only at batch edges)
__device__ double3 tinyAcc(const double3* sp, const double* sm, int N, int i, double* phiOut){
    double3 a = make_double3(0, 0, 0);
    double ph = 0;
    for (int j = 0; j < N; j++){
        if (j == i) continue;
        double dx = sp[j].x - sp[i].x, dy = sp[j].y - sp[i].y, dz = sp[j].z - sp[i].z;
        double r2 = dx*dx + dy*dy + dz*dz + 1e-6;
        double ir = rsqrt(r2);
        double f = (double)G_DIAL*sm[j]*ir*ir*ir;
        a.x += f*dx; a.y += f*dy; a.z += f*dz;
        ph -= (double)G_DIAL*sm[j]*ir;
    }
    *phiOut = ph;
    return a;
}
// 1PN field term + SR inertia in u = p/m form (einstein contract, fp64)
__device__ double3 tinyVel(double3 u, double c2){
    double ig = rsqrt(1.0 + (u.x*u.x + u.y*u.y + u.z*u.z)/c2);
    return make_double3(u.x*ig, u.y*ig, u.z*ig);
}
__device__ double3 tinyForce(double3 g, double ph, double3 v, double c2){
    double v2 = v.x*v.x + v.y*v.y + v.z*v.z;
    double gv = g.x*v.x + g.y*v.y + g.z*v.z;
    double c1 = 1.0 + (v2 + 4.0*ph)/c2;
    return make_double3(c1*g.x - 4.0/c2*gv*v.x,
                        c1*g.y - 4.0/c2*gv*v.y,
                        c1*g.z - 4.0/c2*gv*v.z);
}
__global__ void kTiny(const float4* pin, float4* pout, float4* mom, float* tau,
                      unsigned* regime, int N, int steps, float dtf, float invc2f){
    __shared__ double3 sp[32], su[32], sa[32];
    __shared__ double  sm[32], sph[32];
    int i = threadIdx.x;
    if (i >= N) return;
    double dt = (double)dtf, c2 = 1.0/(double)invc2f;
    float4 P = pin[i], M = mom[i];
    sp[i] = make_double3(P.x, P.y, P.z);
    sm[i] = M.w;
    su[i] = make_double3(M.x/(double)M.w, M.y/(double)M.w, M.z/(double)M.w);  // u = p/m
    double tu = tau[i];
    __syncthreads();
    sa[i] = tinyAcc(sp, sm, N, i, &sph[i]);
    __syncthreads();
    double hdt = 0.5*dt;
    for (int s = 0; s < steps; s++){
        su[i].x += sa[i].x*hdt; su[i].y += sa[i].y*hdt; su[i].z += sa[i].z*hdt;
        double3 v = tinyVel(su[i], c2);
        sp[i].x += v.x*dt;  sp[i].y += v.y*dt;  sp[i].z += v.z*dt;
        __syncthreads();
        sa[i] = tinyAcc(sp, sm, N, i, &sph[i]);
        __syncthreads();
        su[i].x += sa[i].x*hdt; su[i].y += sa[i].y*hdt; su[i].z += sa[i].z*hdt;
        v = tinyVel(su[i], c2);
        double v2 = v.x*v.x + v.y*v.y + v.z*v.z;
        tu += dt*sqrt(fmax(1.0 + (2.0*sph[i] - v2)/c2, 0.0));
    }
    pout[i] = make_float4((float)sp[i].x, (float)sp[i].y, (float)sp[i].z, P.w);
    mom[i]  = make_float4((float)(su[i].x*sm[i]), (float)(su[i].y*sm[i]),
                          (float)(su[i].z*sm[i]), (float)sm[i]);
    tau[i]  = (float)tu;
    double3 vf = tinyVel(su[i], c2);
    double v2f = vf.x*vf.x + vf.y*vf.y + vf.z*vf.z;
    regime[i] = (v2f > (double)REL_V2) ? 0x04u : 0x02u;
}

// --- conservation meters (fixed-point accumulators; PE via PM potential when available)
// acc slots: 0 KE · 1 PE · 2-4 Px,Py,Pz · 5 Lz · 6 nRel · 7 nClassical
__global__ void kMeters(const float4* pos, const float4* mom, const unsigned* regime,
                        const float* phi, int N, unsigned long long* a8){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= N) return;
    if (regime[i] & REG_DEAD) return;                 // ledgered with the BH (declared)
    float4 p = pos[i], m = mom[i];
    double p2 = (double)m.x*m.x + (double)m.y*m.y + (double)m.z*m.z;
    double c2 = (double)C_LIGHT*C_LIGHT;
    double ke = (m.w > 0.0f)
        ? p2*c2/(sqrt((double)m.w*m.w*c2*c2 + p2*c2) + (double)m.w*c2)   // (gamma-1)mc^2, stable
        : sqrt(p2)*C_LIGHT;                                              // photon |p|c
    orrery::fixed_atomic_add(&a8[0], ke);
    if (phi && m.w > 0.0f)
        orrery::fixed_atomic_add(&a8[1], 0.5*(double)m.w*(double)phiAt(phi, p));
    float im = (m.w > 0.0f) ? 1.0f/m.w : 0.0f;
    orrery::fixed_atomic_add(&a8[2], (double)m.x);
    orrery::fixed_atomic_add(&a8[3], (double)m.y);
    orrery::fixed_atomic_add(&a8[4], (double)m.z);
    orrery::fixed_atomic_add(&a8[5], (double)(p.x*m.y - p.y*m.x));
    atomicAdd(&a8[(m.x*m.x + m.y*m.y + m.z*m.z)*im*im > REL_V2 ? 6 : 7], 1ull);
}

// --- M3 arrow: entropy histograms, momentum flip, scenarios, ratchet engine ---
__global__ void kCount(const float4* pos, const float4* mom, int N,
                       unsigned* cx, unsigned* cv){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= N) return;
    float4 p = pos[i], m = mom[i];
    float im = 1.0f/m.w;
    int xi = min(max((int)((p.x + 256.0f)/16.0f), 0), 31);
    int yi = min(max((int)((p.y + 256.0f)/16.0f), 0), 31);
    int zi = min(max((int)((p.z + 256.0f)/16.0f), 0), 31);
    atomicAdd(&cx[(xi*32 + yi)*32 + zi], 1u);
    int vx = min(max((int)((m.x*im + 8.0f)/0.5f), 0), 31);
    int vy = min(max((int)((m.y*im + 8.0f)/0.5f), 0), 31);
    int vz = min(max((int)((m.z*im + 8.0f)/0.5f), 0), 31);
    atomicAdd(&cv[(vx*32 + vy)*32 + vz], 1u);
}
__global__ void kFlipMom(float4* mom, int N){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= N) return;
    float4 m = mom[i];
    mom[i] = make_float4(-m.x, -m.y, -m.z, m.w);
}
__global__ void kInitMerger(float4* pos, float4* mom, float* tau, unsigned* regime,
                            int N, unsigned long long seed){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= N) return;
    using namespace orrery;
    int half = N/2;
    int a = (i < half);
    float gx = (float)counter_gauss(seed, i, 0, 3);
    float gy = (float)counter_gauss(seed, i, 1, 3);
    float gz = (float)counter_gauss(seed, i, 2, 3);
    float il = rsqrtf(gx*gx + gy*gy + gz*gz + 1e-12f);
    float r = 60.0f * cbrtf((float)counter_uniform(seed, i, 3, 3));
    float cxr = a ? -80.0f : 80.0f, cz = a ? -10.0f : 10.0f, vx0 = a ? 1.5f : -1.5f;
    pos[i] = make_float4(cxr + r*gx*il, r*gy*il, cz + r*gz*il, a ? 4500.0f : 9500.0f);
    float sigv = 0.35f*sqrtf(G_DIAL*(float)(N/2)/60.0f);
    mom[i] = make_float4(vx0 + sigv*(float)counter_gauss(seed, i, 4, 3),
                         sigv*(float)counter_gauss(seed, i, 5, 3),
                         sigv*(float)counter_gauss(seed, i, 6, 3), 1.0f);
    tau[i] = 0.0f; regime[i] = 0x02u;
}
__global__ void kInitStream(float4* pos, float4* mom, float* tau, unsigned* regime,
                            int N, unsigned long long seed){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= N) return;
    using namespace orrery;
    float x = -12.0f - 20.0f*(float)counter_uniform(seed, i, 0, 4);
    float y = 20.0f*(float)counter_gauss(seed, i, 1, 4);
    float z = 20.0f*(float)counter_gauss(seed, i, 2, 4);
    pos[i] = make_float4(x, y, z, 5200.0f);
    mom[i] = make_float4(2.0f, 0, 0, 1.0f);
    tau[i] = 0.0f; regime[i] = 0x02u;
}
__global__ void kRatchetInitRecs(unsigned* rec, int nr){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= nr) return;
    rec[i] = (i < nr/3) ? 1u : (i < 2*(nr/3) ? 4u : 16u);
}
// one walk event per unresolved record per tick; absorb at 0 or R0+30
__global__ void kRatchetStep(unsigned* rec, int nr, unsigned long long tick,
                             unsigned long long seed){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= nr) return;
    unsigned n = rec[i];
    unsigned R0 = (i < nr/3) ? 1u : (i < 2*(nr/3) ? 4u : 16u);
    unsigned cap = R0 + 30u;
    if (n == 0u || n >= cap) return;
    double u = orrery::counter_uniform(seed ^ 0xA77CE7ull, i, tick, 5);
    rec[i] = n + ((u < RATCHET_PDOWN) ? -1 : +1);
}
// detector slab x in [0,8]: writes +1/tick in-slab (cap 64); ratchet walk outside;
// bit 0x80 = ever-inscribed, bit 0x40 = RECORDED (n >= 16)
__global__ void kDetectorStep(const float4* pos, unsigned* regime, unsigned* rec,
                              int N, unsigned long long tick, unsigned long long seed){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= N) return;
    float x = pos[i].x;
    unsigned n = rec[i];
    bool in = (x >= 0.0f && x <= 8.0f);
    if (in){
        n = (n < 64u) ? n + 1u : 64u;
        regime[i] |= 0x80u;
    } else if (n > 0u && n < 64u){
        double u = orrery::counter_uniform(seed ^ 0xDE7EC7ull, i, tick, 6);
        n += (u < RATCHET_PDOWN) ? -1 : +1;
    }
    rec[i] = n;
    if (n >= REC_THRESH) regime[i] |= 0x40u;
}

// --- M5 gargantua: BH force, absorption, formation, step (gargantua contract) ---
__global__ void kBHForce(const float4* pos, unsigned* regime, float4* acc, int N,
                         const float4* bhpos, const float4* bhmom, const unsigned* bhn){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= N) return;
    if (regime[i] & REG_DEAD) return;
    int nbh = (int)bhn[0];
    float4 p = pos[i], a = acc[i];
    for (int b = 0; b < nbh; b++){
        float4 bp = bhpos[b];
        float rs = bhmom[b].w;
        float dx = mimg(bp.x - p.x), dy = mimg(bp.y - p.y), dz = mimg(bp.z - p.z);
        float r = sqrtf(dx*dx + dy*dy + dz*dz + 1e-12f);
        float d = fmaxf(r - rs, 0.05f);
        float g = G_DIAL*bp.w/(d*d);
        a.x += g*dx/r; a.y += g*dy/r; a.z += g*dz/r;
        a.w -= G_DIAL*bp.w/d;                          // Phi_PW into tau/1PN channel
        if (r < 10.0f*rs) regime[i] |= 0x08u;          // COMPACT
    }
    acc[i] = a;
}
__global__ void kAbsorb(const float4* pos, const float4* mom, unsigned* regime, int N,
                        const float4* bhpos, const float4* bhmom, const unsigned* bhn,
                        unsigned long long* bhacc, unsigned long long* ledger){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= N) return;
    if (regime[i] & REG_DEAD) return;
    int nbh = (int)bhn[0];
    float4 p = pos[i], m = mom[i];
    for (int b = 0; b < nbh; b++){
        float4 bp = bhpos[b];
        float rs = bhmom[b].w;
        float cap = fmaxf(1.2f*rs, (bp.w < M_FORM) ? 0.75f*CELL : 0.3f);
        float dx = mimg(p.x - bp.x), dy = mimg(p.y - bp.y), dz = mimg(p.z - bp.z);
        if (dx*dx + dy*dy + dz*dz < cap*cap){
            regime[i] |= REG_DEAD;
            orrery::fixed_atomic_add(&bhacc[b*4 + 0], (double)m.w);
            orrery::fixed_atomic_add(&bhacc[b*4 + 1], (double)m.x);
            orrery::fixed_atomic_add(&bhacc[b*4 + 2], (double)m.y);
            orrery::fixed_atomic_add(&bhacc[b*4 + 3], (double)m.z);
            orrery::fixed_atomic_add(&ledger[0], (double)m.w);
            return;
        }
    }
}
__global__ void kBHDetect(const unsigned long long* fix, int n, unsigned long long* argmax){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= n) return;
    double m = orrery::fixed_decode((long long)fix[i]);
    if (m < 1000.0) return;                            // probe floor (argmax candidates only)
    unsigned long long pack = ((unsigned long long)llrint(m*8.0) << 24) | (unsigned long long)i;
    atomicMax(argmax, pack);
}
__global__ void kBHSpawn(unsigned long long* argmax, float4* bhpos, float4* bhmom,
                         unsigned* bhn){
    if (threadIdx.x != 0 || blockIdx.x != 0) return;
    unsigned long long pk = *argmax;
    *argmax = 0ull;
    double m = (double)(pk >> 24)/8.0;
    if ((unsigned)m > bhn[2]) bhn[2] = (unsigned)m;    // peak-cell-mass probe (declared)
    if (m < (double)M_FORM) return;
    unsigned n = bhn[0];
    if (n >= NBH_MAX) return;
    int cell = (int)(pk & 0xFFFFFFull);
    int cx = cell/(PMN*PMN), cy = (cell/PMN) % PMN, cz = cell % PMN;
    float x = -BOXL*0.5f + (cx + 0.5f)*CELL;
    float y = -BOXL*0.5f + (cy + 0.5f)*CELL;
    float z = -BOXL*0.5f + (cz + 0.5f)*CELL;
    for (unsigned b = 0; b < n; b++){
        float dx = mimg(bhpos[b].x - x), dy = mimg(bhpos[b].y - y), dz = mimg(bhpos[b].z - z);
        if (dx*dx + dy*dy + dz*dz < 64.0f) return;     // existing BH owns this core
    }
    bhpos[n] = make_float4(x, y, z, 0.0f);             // seed mass 0: formation-tick capture binds the core
    bhmom[n] = make_float4(0, 0, 0, 0.0f);
    bhn[0] = n + 1;
}
// one-block housekeeping: apply absorption deltas, Hawking cube-step, BH motion
__global__ void kBHStep(float4* bhpos, float4* bhmom, unsigned long long* bhacc,
                        unsigned* bhn, unsigned long long* ledger,
                        const float4* For, const float* phi, float dt, double kHawk,
                        unsigned long long tick){
    if (blockIdx.x != 0 || threadIdx.x != 0) return;
    unsigned n = bhn[0];
    double c2 = (double)C_LIGHT*C_LIGHT;
    for (unsigned b = 0; b < n; b++){
        float4 bp = bhpos[b], bm = bhmom[b];
        // absorption deltas (accumulated during the previous tick)
        double dm = orrery::fixed_decode((long long)bhacc[b*4 + 0]);
        bp.w += (float)dm;
        bm.x += (float)orrery::fixed_decode((long long)bhacc[b*4 + 1]);
        bm.y += (float)orrery::fixed_decode((long long)bhacc[b*4 + 2]);
        bm.z += (float)orrery::fixed_decode((long long)bhacc[b*4 + 3]);
        bhacc[b*4+0] = bhacc[b*4+1] = bhacc[b*4+2] = bhacc[b*4+3] = 0ull;
        // Hawking: exact cube step (gargantua contract)
        if (bp.w > 0.0f){
            double M = bp.w;
            double m3 = M*M*M - 3.0*kHawk*(double)dt;
            double Mn = (m3 > 0.0) ? cbrt(m3) : 0.0;
            if (Mn <= (double)M_POP){
                orrery::fixed_atomic_add(&ledger[1], M*c2);        // pop: remainder radiates
                if (ledger[2] == 0ull) ledger[2] = tick;           // first-pop tick (declared)
                bp.w = -1.0f;                                      // popped sentinel (seeds are 0)
            } else {
                orrery::fixed_atomic_add(&ledger[1], (M - Mn)*c2);
                bp.w = (float)Mn;
            }
        }
        // motion: PM field at BH position + Newtonian BH-BH (declared)
        if (bp.w > 0.0f){
            float3 g = make_float3(0, 0, 0);
            if (For) {
                float gx = (bp.x + BOXL*0.5f)/CELL, gy = (bp.y + BOXL*0.5f)/CELL, gz = (bp.z + BOXL*0.5f)/CELL;
                int x0 = (int)floorf(gx), y0 = (int)floorf(gy), z0 = (int)floorf(gz);
                float fx = gx - x0, fy = gy - y0, fz = gz - z0;
                for (int dz2 = 0; dz2 < 2; dz2++)
                    for (int dy2 = 0; dy2 < 2; dy2++)
                        for (int dx2 = 0; dx2 < 2; dx2++){
                            float w = (dx2 ? fx : 1.0f-fx)*(dy2 ? fy : 1.0f-fy)*(dz2 ? fz : 1.0f-fz);
                            float4 f = For[(wrapc(x0+dx2)*PMN + wrapc(y0+dy2))*PMN + wrapc(z0+dz2)];
                            g.x += w*f.x; g.y += w*f.y; g.z += w*f.z;
                        }
            }
            for (unsigned o = 0; o < n; o++){
                if (o == b || bhpos[o].w <= 0.0f) continue;
                float dx = mimg(bhpos[o].x - bp.x), dy = mimg(bhpos[o].y - bp.y), dz = mimg(bhpos[o].z - bp.z);
                float r2 = dx*dx + dy*dy + dz*dz + 1.0f;
                float ir = rsqrtf(r2);
                float gg = G_DIAL*bhpos[o].w*ir*ir*ir;
                g.x += gg*dx; g.y += gg*dy; g.z += gg*dz;
            }
            bm.x += bp.w*g.x*dt; bm.y += bp.w*g.y*dt; bm.z += bp.w*g.z*dt;
            float im = 1.0f/bp.w;
            bp.x += bm.x*im*dt; bp.y += bm.y*im*dt; bp.z += bm.z*im*dt;
            if      (bp.x >=  256.0f) bp.x -= 512.0f; else if (bp.x < -256.0f) bp.x += 512.0f;
            if      (bp.y >=  256.0f) bp.y -= 512.0f; else if (bp.y < -256.0f) bp.y += 512.0f;
            if      (bp.z >=  256.0f) bp.z -= 512.0f; else if (bp.z < -256.0f) bp.z += 512.0f;
        }
        bm.w = fmaxf(RS_PER_M*bp.w, 0.0f);             // r_s cache
        bhpos[b] = bp; bhmom[b] = bm;
    }
    // compact popped BHs only — mass-0 SEEDS survive to their first meal (D-017 bug fix)
    unsigned w = 0;
    for (unsigned b = 0; b < n; b++){
        if (bhpos[b].w >= 0.0f){
            if (w != b){ bhpos[w] = bhpos[b]; bhmom[w] = bhmom[b]; }
            w++;
        }
    }
    for (unsigned b = w; b < n; b++){                  // zero vacated slots (hash hygiene)
        bhpos[b] = make_float4(0, 0, 0, 0);
        bhmom[b] = make_float4(0, 0, 0, 0);
    }
    bhn[0] = w;
}
// screen-space point-lens warp (visuals, non-declared; gargantua contract)
__global__ void kLens(const float4* src, float4* dst, int W, int H, RP rp){
    int x = blockIdx.x*blockDim.x + threadIdx.x;
    int y = blockIdx.y*blockDim.y + threadIdx.y;
    if (x >= W || y >= H) return;
    float bestI = 0; int bb = -1;
    for (int b = 0; b < rp.nLens; b++){
        float dx = x - rp.lensX[b], dy = y - rp.lensY[b];
        float r2 = dx*dx + dy*dy + 1e-6f;
        float I = rp.lensTE2[b]/r2;
        if (I > bestI){ bestI = I; bb = b; }
    }
    if (bb < 0){ dst[y*W + x] = src[y*W + x]; return; }
    float dx = x - rp.lensX[bb], dy = y - rp.lensY[bb];
    float r = sqrtf(dx*dx + dy*dy + 1e-6f);
    if (r < rp.lensShadow[bb]){ dst[y*W + x] = make_float4(0, 0, 0, 1); return; }
    float s = 1.0f - rp.lensTE2[bb]/(r*r);             // beta = theta(1 - thetaE^2/theta^2)
    float sx = rp.lensX[bb] + dx*s, sy = rp.lensY[bb] + dy*s;
    int x0 = min(max((int)sx, 0), W - 2), y0 = min(max((int)sy, 0), H - 2);
    float fx = fminf(fmaxf(sx - x0, 0.0f), 1.0f), fy = fminf(fmaxf(sy - y0, 0.0f), 1.0f);
    float4 a = src[y0*W + x0], b2 = src[y0*W + x0 + 1];
    float4 c = src[(y0+1)*W + x0], d = src[(y0+1)*W + x0 + 1];
    float4 o;
    o.x = (a.x*(1-fx) + b2.x*fx)*(1-fy) + (c.x*(1-fx) + d.x*fx)*fy;
    o.y = (a.y*(1-fx) + b2.y*fx)*(1-fy) + (c.y*(1-fx) + d.y*fx)*fy;
    o.z = (a.z*(1-fx) + b2.z*fx)*(1-fy) + (c.z*(1-fx) + d.z*fx)*fy;
    o.w = 1;
    dst[y*W + x] = o;
}
// Hawking glow: BHs splat as blackbody points at their horizon temperature
// (M7: all 27 periodic images; no history — the hole itself is bright "now")
__global__ void kSplatBH(const float4* bhpos, int nbh, float4* hdr, RP rp){
    int b = blockIdx.x*blockDim.x + threadIdx.x;
    if (b >= nbh) return;
    float4 p = bhpos[b];
    if (p.w <= 0.0f) return;
    float bx = mimg(p.x - rp.cpos[0]);
    float by = mimg(p.y - rp.cpos[1]);
    float bz = mimg(p.z - rp.cpos[2]);
    float TH = 79577.0f/p.w;                           // horizon temperature (energy units)
    float Tdisp = clampf(TH*40.0f, 1200.0f, 39000.0f);
    float bb[3]; blackbody(Tdisp, bb);
    for (int oz = -1; oz <= 1; oz++)
    for (int oy2 = -1; oy2 <= 1; oy2++)
    for (int ox2 = -1; ox2 <= 1; ox2++){
        float px, py, fd;
        if (!projCam(rp, bx + 512.0f*ox2, by + 512.0f*oy2, bz + 512.0f*oz, px, py, fd)) continue;
        if (px < 2 || py < 2 || px >= rp.RW - 3 || py >= rp.RH - 3) continue;
        float flux = 40.0f*(250.0f/p.w)*(250.0f/p.w)*4200.0f/(fd*fd + 25.0f);
        int x0 = (int)px, y0 = (int)py;
        for (int oy = -1; oy <= 1; oy++)
            for (int ox = -1; ox <= 1; ox++){
                float w = expf(-0.5f*(ox*ox + oy*oy));
                float4* h = &hdr[(y0+oy)*rp.RW + x0 + ox];
                atomicAdd(&h->x, bb[0]*flux*w);
                atomicAdd(&h->y, bb[1]*flux*w);
                atomicAdd(&h->z, bb[2]*flux*w);
            }
    }
}

// --- M6 planck: psi engine (2D 256^2 split-step; planck contract) ---
__global__ void kQGauss(cufftComplex* q, float x0, float y0, float sx, float sy, float kx0){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if (idx >= QN*QN) return;
    int ix = idx % QN, iy = idx/QN;
    float x = -QL*0.5f + (ix + 0.5f)*QDX;
    float y = -QL*0.5f + (iy + 0.5f)*QDX;
    float dx = (x - x0)/sx, dy = (y - y0)/sy;
    float env = expf(-0.25f*(dx*dx + dy*dy));
    float ph = kx0*x;
    q[idx].x = env*cosf(ph);
    q[idx].y = env*sinf(ph);
}
__global__ void kQPhaseK(cufftComplex* q, float coef){       // e^{-i coef k^2}
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if (idx >= QN*QN) return;
    int ix = idx % QN, iy = idx/QN;
    int fx = (ix <= QN/2) ? ix : ix - QN;
    int fy = (iy <= QN/2) ? iy : iy - QN;
    float kf = 2.0f*PI_F/QL;
    float k2 = kf*kf*(float)(fx*fx + fy*fy);
    float ph = -coef*k2;
    float c = cosf(ph), s = sinf(ph);
    cufftComplex v = q[idx];
    q[idx].x = v.x*c - v.y*s;
    q[idx].y = v.x*s + v.y*c;
}
__global__ void kQDecayK(cufftComplex* q, float coef){       // e^{-coef k^2}
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if (idx >= QN*QN) return;
    int ix = idx % QN, iy = idx/QN;
    int fx = (ix <= QN/2) ? ix : ix - QN;
    int fy = (iy <= QN/2) ? iy : iy - QN;
    float kf = 2.0f*PI_F/QL;
    float d = expf(-coef*kf*kf*(float)(fx*fx + fy*fy));
    q[idx].x *= d; q[idx].y *= d;
}
__global__ void kQPhaseV(cufftComplex* q, const float* V, float coef){   // e^{-i V coef}
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if (idx >= QN*QN) return;
    float ph = -V[idx]*coef;
    float c = cosf(ph), s = sinf(ph);
    cufftComplex v = q[idx];
    q[idx].x = v.x*c - v.y*s;
    q[idx].y = v.x*s + v.y*c;
}
__global__ void kQDecayV(cufftComplex* q, const float* V, float coef){   // e^{-V coef}
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if (idx >= QN*QN) return;
    float d = expf(-V[idx]*coef);
    q[idx].x *= d; q[idx].y *= d;
}
__global__ void kQMul(cufftComplex* q, const float* m){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if (idx >= QN*QN) return;
    q[idx].x *= m[idx]; q[idx].y *= m[idx];
}
__global__ void kQScale(cufftComplex* q, float s){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if (idx >= QN*QN) return;
    q[idx].x *= s; q[idx].y *= s;
}
__global__ void kQNormAcc(const cufftComplex* q, unsigned long long* acc){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if (idx >= QN*QN) return;
    orrery::fixed_atomic_add(acc, (double)q[idx].x*q[idx].x + (double)q[idx].y*q[idx].y);
}

// --- M7 cosmos: 3D bubble kernels (cosmos contract; M6 split-step lineage) ---
// grid coords are bubble-local: r = -QBL/2 + (i + 0.5)*QBDX, center = owner at spawn
__global__ void kQ3Gauss(cufftComplex* q, float kx, float ky, float kz, float s0){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if (idx >= QB3) return;
    int ix = idx % QB, iy = (idx/QB) % QB, iz = idx/(QB*QB);
    float x = -QBL*0.5f + (ix + 0.5f)*QBDX;
    float y = -QBL*0.5f + (iy + 0.5f)*QBDX;
    float z = -QBL*0.5f + (iz + 0.5f)*QBDX;
    float env = expf(-0.25f*(x*x + y*y + z*z)/(s0*s0));
    float ph = kx*x + ky*y + kz*z;
    q[idx].x = env*cosf(ph);
    q[idx].y = env*sinf(ph);
}
__global__ void kQ3PhaseK(cufftComplex* q, int total, float coef){   // e^{-i coef k^2}, all slots
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if (idx >= total) return;
    int c = idx % QB3;
    int ix = c % QB, iy = (c/QB) % QB, iz = c/(QB*QB);
    int fx = (ix <= QB/2) ? ix : ix - QB;
    int fy = (iy <= QB/2) ? iy : iy - QB;
    int fz = (iz <= QB/2) ? iz : iz - QB;
    float kf = 2.0f*PI_F/QBL;
    float k2 = kf*kf*(float)(fx*fx + fy*fy + fz*fz);
    float ph = -coef*k2;
    float cs = cosf(ph), sn = sinf(ph);
    cufftComplex v = q[idx];
    q[idx].x = v.x*cs - v.y*sn;
    q[idx].y = v.x*sn + v.y*cs;
}
__global__ void kQ3Edge(cufftComplex* q, int total, const float* edge){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if (idx >= total) return;
    float m = edge[idx % QB3];
    q[idx].x *= m; q[idx].y *= m;
}
__global__ void kQ3Scale(cufftComplex* q, int total, float s){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if (idx >= total) return;
    q[idx].x *= s; q[idx].y *= s;
}
__global__ void kQ3NormAcc(const cufftComplex* q, unsigned long long* acc){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if (idx >= QB3) return;
    orrery::fixed_atomic_add(acc, (double)q[idx].x*q[idx].x + (double)q[idx].y*q[idx].y);
}
// isolation: a massive, live, non-bubbled particle whose 3^3 PM-cell neighborhood
// holds no deposited mass but its own (fixed-point exact compare; cosmos contract)
__global__ void kBubIso(const float4* pos, const float4* mom, const unsigned* regime,
                        const unsigned* partBub, int N, const unsigned long long* fix,
                        unsigned char* iso){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= N) return;
    iso[i] = 0;
    if (mom[i].w <= 0.0f) return;
    if (regime[i] & REG_DEAD) return;
    if (regime[i] & 0x40u) return;    // RECORDED stays classical: its which-path is
                                      // inscribed in the environment (M3 semantics)
    if (partBub[i]) return;
    float4 p = pos[i];
    int cx = (int)floorf((p.x + BOXL*0.5f)/CELL);
    int cy = (int)floorf((p.y + BOXL*0.5f)/CELL);
    int cz = (int)floorf((p.z + BOXL*0.5f)/CELL);
    double msum = 0.0;
    for (int dz = -1; dz <= 1; dz++)
        for (int dy = -1; dy <= 1; dy++)
            for (int dx = -1; dx <= 1; dx++)
                msum += orrery::fixed_decode((long long)
                    fix[(wrapc(cx + dx)*PMN + wrapc(cy + dy))*PMN + wrapc(cz + dz)]);
    if (msum - (double)mom[i].w < 0.5) iso[i] = 1;
}
// single-thread slot assignment in ascending particle id (deterministic)
__global__ void kBubSpawnScan(const unsigned char* iso, const float4* pos,
                              unsigned* regime, unsigned* partBub, int N,
                              int* bubOwner, float4* bubCtr, unsigned* bubRec,
                              int* newList, int* newCount){
    if (threadIdx.x || blockIdx.x) return;
    int nn = 0;
    for (int i = 0; i < N && nn < NBUB; i++){
        if (!iso[i] || partBub[i]) continue;
        int slot = -1;
        for (int s = 0; s < NBUB; s++) if (bubOwner[s] < 0){ slot = s; break; }
        if (slot < 0) break;
        bubOwner[slot] = i;
        bubCtr[slot] = pos[i];                        // xyz = grid center; w carries emitT
        bubRec[slot] = 0u;
        partBub[i] = (unsigned)(slot + 1);
        regime[i] = (regime[i] & 0xC0u) | 0x01u;      // QUANTUM (latches preserved)
        newList[nn++] = slot;
    }
    *newCount = nn;
}
// intrusion: any foreign massive live particle inside a bubble's window this tick
__global__ void kBubIntrude(const float4* pos, const float4* mom, const unsigned* regime,
                            const unsigned* partBub, int N, const int* bubOwner,
                            const float4* bubCtr, unsigned* bubFlag){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= N) return;
    if (mom[i].w <= 0.0f) return;
    if (regime[i] & REG_DEAD) return;
    if (partBub[i]) return;
    float4 p = pos[i];
    for (int s = 0; s < NBUB; s++){
        if (bubOwner[s] < 0) continue;
        float4 c = bubCtr[s];
        if (fabsf(mimg(p.x - c.x)) < QBL*0.5f &&
            fabsf(mimg(p.y - c.y)) < QBL*0.5f &&
            fabsf(mimg(p.z - c.z)) < QBL*0.5f) atomicOr(&bubFlag[s], 1u);
    }
}
// per-tick record walk: strong write (+1 per intruded tick); collapse past threshold
__global__ void kBubTick(const int* bubOwner, unsigned* bubRec, unsigned* bubFlag,
                         unsigned* bubPend){
    if (threadIdx.x || blockIdx.x) return;
    unsigned pend = 0;
    for (int s = 0; s < NBUB; s++){
        if (bubOwner[s] >= 0 && bubFlag[s]){
            bubRec[s]++;
            if (bubRec[s] >= REC_THRESH) pend |= (1u << s);
        }
        bubFlag[s] = 0u;
    }
    *bubPend = pend;
}

__global__ void kClear(float4* hdr, int n){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < n) hdr[i] = make_float4(0, 0, 0, 0);
}

// M7 cosmos: torus-aware splat — 27 periodic images around the minimum image,
// each sampled from the light-history ring at retarded time t - D/c (one
// refinement pass; visuals non-declared, cosmos contract).
__global__ void kSplat(const float4* pos, const unsigned* reg, int N, float4* hdr, RP rp,
                       const ushort4* hist){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= N) return;
    float4 p = pos[i];
    if (reg && (reg[i] & REG_DEAD)) return;                    // absorbed: not rendered
    if (reg && (reg[i] & 0x01u)) return;                       // bubbled: the cloud renders, not the dot
    bool rec = reg && (reg[i] & 0x40u);
    float bx = mimg(p.x - rp.cpos[0]);
    float by = mimg(p.y - rp.cpos[1]);
    float bz = mimg(p.z - rp.cpos[2]);
    for (int oz = -1; oz <= 1; oz++)
    for (int oy = -1; oy <= 1; oy++)
    for (int ox = -1; ox <= 1; ox++){
        float rx = bx + 512.0f*ox, ry = by + 512.0f*oy, rz = bz + 512.0f*oz;
        float px, py, fd;
        if (!projCam(rp, rx, ry, rz, px, py, fd)) continue;    // pre-cull on current pos
        float T = p.w;
        if (hist && rp.histN > 0){
            float D = sqrtf(rx*rx + ry*ry + rz*rz);
            int back = (int)(D*0.5f + 0.5f);          // 24-tick snapshots = 2 su light travel
            if (back > 0){
                if (back > rp.histN) back = rp.histN;
                int slot = rp.histHead - (back - 1); if (slot < 0) slot += rp.histDepth;
                ushort4 h = hist[(size_t)slot*N + i];
                float hx = __half2float(__ushort_as_half(h.x));
                float hy = __half2float(__ushort_as_half(h.y));
                float hz = __half2float(__ushort_as_half(h.z));
                if (fabsf(hx) > 300.0f) continue;     // parked at that epoch (absorbed)
                rx = mimg(hx - rp.cpos[0]) + 512.0f*ox;
                ry = mimg(hy - rp.cpos[1]) + 512.0f*oy;
                rz = mimg(hz - rp.cpos[2]) + 512.0f*oz;
                float D2 = sqrtf(rx*rx + ry*ry + rz*rz);
                int back2 = (int)(D2*0.5f + 0.5f);    // one retarded-time refinement
                if (back2 > rp.histN) back2 = rp.histN;
                if (back2 > 0 && back2 != back){
                    slot = rp.histHead - (back2 - 1); if (slot < 0) slot += rp.histDepth;
                    h = hist[(size_t)slot*N + i];
                    hx = __half2float(__ushort_as_half(h.x));
                    hy = __half2float(__ushort_as_half(h.y));
                    hz = __half2float(__ushort_as_half(h.z));
                    if (fabsf(hx) > 300.0f) continue;
                    rx = mimg(hx - rp.cpos[0]) + 512.0f*ox;
                    ry = mimg(hy - rp.cpos[1]) + 512.0f*oy;
                    rz = mimg(hz - rp.cpos[2]) + 512.0f*oz;
                }
                T = __half2float(__ushort_as_half(h.w));
                if (!projCam(rp, rx, ry, rz, px, py, fd)) continue;
            }
        }
        if (rec) T = fmaxf(T, 13000.0f);              // RECORDED: blue-white tint
        if (px < 0 || py < 0 || px >= rp.RW - 1 || py >= rp.RH - 1) continue;
        float L = rp.mode ? powf(T/5800.0f, 4.0f)     // physical: L ~ T^4
                          : powf(T/5800.0f, 2.2f);    // cinematic: softened
        float flux = 4200.0f * L / (fd*fd + 25.0f);
        float bb[3]; blackbody(T, bb);
        int x0 = (int)px, y0 = (int)py;
        float fx = px - x0, fy = py - y0;
        float w00 = (1-fx)*(1-fy), w10 = fx*(1-fy), w01 = (1-fx)*fy, w11 = fx*fy;
        float4* h4;
        h4 = &hdr[y0*rp.RW + x0];
        atomicAdd(&h4->x, bb[0]*flux*w00); atomicAdd(&h4->y, bb[1]*flux*w00); atomicAdd(&h4->z, bb[2]*flux*w00);
        h4 = &hdr[y0*rp.RW + x0 + 1];
        atomicAdd(&h4->x, bb[0]*flux*w10); atomicAdd(&h4->y, bb[1]*flux*w10); atomicAdd(&h4->z, bb[2]*flux*w10);
        h4 = &hdr[(y0+1)*rp.RW + x0];
        atomicAdd(&h4->x, bb[0]*flux*w01); atomicAdd(&h4->y, bb[1]*flux*w01); atomicAdd(&h4->z, bb[2]*flux*w01);
        h4 = &hdr[(y0+1)*rp.RW + x0 + 1];
        atomicAdd(&h4->x, bb[0]*flux*w11); atomicAdd(&h4->y, bb[1]*flux*w11); atomicAdd(&h4->z, bb[2]*flux*w11);
    }
}

// M7 cosmos: bubbleVis — |psi|^2 volume splat as a hot quantum glow (nearest image;
// psi is read live off simS: shimmer is non-declared)
__global__ void kSplatPsi(const cufftComplex* bub, const int* bubOwner, const float4* bubCtr,
                          float4* hdr, RP rp){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if (idx >= NBUB*QB3) return;
    int s = idx/QB3, c = idx % QB3;
    if (bubOwner[s] < 0) return;
    cufftComplex a = bub[idx];
    float p2 = a.x*a.x + a.y*a.y;
    if (p2 < 1e-7f) return;
    float4 ctr = bubCtr[s];
    int ix = c % QB, iy = (c/QB) % QB, iz = c/(QB*QB);
    float wx = ctr.x + (-QBL*0.5f + (ix + 0.5f)*QBDX);
    float wy = ctr.y + (-QBL*0.5f + (iy + 0.5f)*QBDX);
    float wz = ctr.z + (-QBL*0.5f + (iz + 0.5f)*QBDX);
    float px, py, fd;
    if (!projCam(rp, mimg(wx - rp.cpos[0]), mimg(wy - rp.cpos[1]), mimg(wz - rp.cpos[2]),
                 px, py, fd)) return;
    if (px < 1 || py < 1 || px >= rp.RW - 2 || py >= rp.RH - 2) return;
    float flux = 26000.0f*p2*(QBDX*QBDX*QBDX)/(fd*fd + 25.0f);
    float bb[3]; blackbody(15000.0f, bb);
    float4* h = &hdr[(int)py*rp.RW + (int)px];
    atomicAdd(&h->x, bb[0]*flux);
    atomicAdd(&h->y, bb[1]*flux);
    atomicAdd(&h->z, bb[2]*flux);
}

// M7 cosmos: light-history snapshot (half4; parked DEAD positions overflow fp16
// to +inf and are skipped at sample time)
__global__ void kSnap(const float4* pos, ushort4* dst, int N){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= N) return;
    float4 p = pos[i];
    dst[i] = make_ushort4(__half_as_ushort(__float2half(p.x)),
                          __half_as_ushort(__float2half(p.y)),
                          __half_as_ushort(__float2half(p.z)),
                          __half_as_ushort(__float2half(p.w)));
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
static const float DT = 1.0f/240.0f;                  // dial (C_LIGHT defined above)

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
static unsigned *dCntX = 0, *dCntV = 0, *dRec = 0;    // M3: entropy histograms + records
static double gS = 0, gdS = 0;                        // HUD entropy readouts
// M5 gargantua: BH entities (device-resident, deterministic)
static float4 *dBHpos = 0, *dBHmom = 0;               // pos.w = mass, mom.w = r_s cache
static unsigned long long *dBHacc = 0;                // per-BH deltas [b*4+{dm,dpx,dpy,dpz}]
static unsigned long long *dLedger = 0;               // [0] absorbed mass, [1] radiated energy
static unsigned long long *dArgmax = 0;               // packed formation argmax
static unsigned *dBHn = 0;                            // [0] count
static bool  gBHFormEnabled = false, gBHActive = false;
static const double gKHawk = 80000.0/(15360.0*3.14159265358979*4.0e-6);  // hbar c^4/(15360 pi G^2)
static unsigned hBHn = 0;
static float4 hBHpos[NBH_MAX];
static float4 *dHdrW = 0;                             // lensed HDR (render res)
// M6 planck: psi engine state
static cufftComplex *dQ = 0;
static float *dQWall = 0, *dQEdge = 0, *dQV = 0;
static unsigned long long *dQNorm = 0;
static cufftHandle planQ;
static bool gQAlloc = false;
static std::string gPsiBytes;                         // final psi (joins declared hash)
struct QMet { double c_psi, c_dots, c_det, kpk_rel, T, Tref, E, sigx; long long nrec; };
static QMet gQM = {};
static double gMTot = 1e6;
// M7 cosmos: EdS expansion state (cosmos contract). a(t) closed-form, a(0) = 1;
// H0 from the actual box density: H0 = sqrt(8 pi G rho_bar / 3).
static bool   gExpand = false;
static double gH0 = 0.0;
static double aOfT(double t){ return gExpand ? pow(1.0 + 1.5*gH0*t, 2.0/3.0) : 1.0; }
// M7 cosmos: roaming 3D bubble state (cosmos contract)
static cufftComplex *dBub = 0;                        // NBUB x QB^3 grids
static float    *dBubEdge = 0;                        // shared cos^2 absorber
static int      *dBubOwner = 0, *dBubNew = 0, *dBubNewCount = 0;
static float4   *dBubCtr = 0;
static unsigned *dBubRec = 0, *dBubFlag = 0, *dBubPend = 0;
static unsigned *dPartBub = 0;                        // per-particle handle (slot+1, 0 = none)
static unsigned char *dIso = 0;
static unsigned long long *dBubNorm = 0;
static cufftHandle planB;
static bool gBubActive = false, gBubAlloc = false, gBubWant = false;
static int  hBubOwner[NBUB];
static unsigned gBubSpawned = 0, gBubCollapsed = 0;
static int  gBubLive = 0;
static double gBubSig = 0, gBubSigExp = 0;
// M7 cosmos: projection + light-history host state (visuals non-declared)
static int  gProj = 0;                                // 0 perspective, 1 little planet
static float gShotAz = -999, gShotEl = -999, gShotDist = -999;   // shot framing overrides
static ushort4 *dHistBuf = 0;
static int  gHistDepth = 0, gHistN = 0, gHistHead = -1;
static bool gHistOn = false, gWantHist = true;        // false in the headless json/golden face
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
    const float D0[20] = {380.0f, 650.0f, 9.0f, 420.0f, 420.0f, 420.0f,
                          300.0f, 120.0f, 45.0f, 260.0f, 260.0f,
                          90.0f, 25.0f, 30.0f, 50.0f, 50.0f, 50.0f,
                          300.0f, 700.0f, 260.0f};
    float d = D0[gScenario];
    if (gScenario == SC_DOUBLESLIT && i == 0){        // face the screen plane
        cAz.tgt = 0; cEl.tgt = 6; cDist.tgt = logf(d);
        return;
    }
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
    // BH mirror + screen-space lens parameters (visuals, non-declared)
    rp.proj = gProj;
    rp.lpR = 0.42f*(float)rp.RH;
    rp.histN = gHistOn ? gHistN : 0;
    rp.histHead = gHistHead; rp.histDepth = gHistDepth;
    rp.nLens = 0;
    if (gBHActive){
        CUDA_CHECK(cudaMemcpy(&hBHn, dBHn, sizeof hBHn, cudaMemcpyDeviceToHost));
        if (hBHn > 0)
            CUDA_CHECK(cudaMemcpy(hBHpos, dBHpos, hBHn*sizeof(float4), cudaMemcpyDeviceToHost));
        for (unsigned b = 0; b < hBHn && rp.nLens < 4; b++){
            float4 bp = hBHpos[b];
            if (bp.w <= 0.0f) continue;
            float rx = mimg(bp.x - rp.cpos[0]), ry = mimg(bp.y - rp.cpos[1]), rz = mimg(bp.z - rp.cpos[2]);
            float cz = rx*rp.fwd[0] + ry*rp.fwd[1] + rz*rp.fwd[2];
            float cx = rx*rp.rgt[0] + ry*rp.rgt[1] + rz*rp.rgt[2];
            float cy = rx*rp.upv[0] + ry*rp.upv[1] + rz*rp.upv[2];
            float D = sqrtf(rx*rx + ry*ry + rz*rz + 1e-12f);
            float thE = sqrtf(4.0f*G_DIAL*bp.w/((C_LIGHT*C_LIGHT)*D));
            float lx, ly, pxpr;
            if (rp.proj == 0){
                if (cz < 2.0f) continue;
                lx = ((cx/(cz*rp.tf*rp.aspect))*0.5f + 0.5f)*rp.RW;
                ly = (0.5f - (cy/(cz*rp.tf))*0.5f)*rp.RH;
                pxpr = rp.RH/(2.0f*rp.tf);            // px per radian (vertical, center)
            } else {
                float cd = cz/D;
                if (cd < -0.9f) continue;
                float s = 1.0f/(1.0f + cd);
                lx = rp.RW*0.5f + (cx/D)*s*rp.lpR;
                ly = rp.RH*0.5f - (cy/D)*s*rp.lpR;
                float sig2 = ((cx*cx + cy*cy)/(D*D))*s*s;
                pxpr = rp.lpR*0.5f*(1.0f + sig2);     // stereographic d(sigma)/d(theta)
            }
            float tE_px = thE*pxpr;
            int L = rp.nLens++;
            rp.lensX[L] = lx;
            rp.lensY[L] = ly;
            rp.lensTE2[L] = tE_px*tE_px;
            rp.lensShadow[L] = 2.6f*(RS_PER_M*bp.w)/D*pxpr;
        }
    }
    dim3 b1(256), g1((npx + 255)/256);
    kClear<<<g1, b1, 0, prsS>>>(dHdr, npx);
    kSplat<<<(gN + 255)/256, b1, 0, prsS>>>(dPos[gPub], dRegime, gN, dHdr, rp,
                                            gHistOn ? dHistBuf : nullptr);
    if (gBHActive && hBHn > 0)
        kSplatBH<<<1, NBH_MAX, 0, prsS>>>(dBHpos, (int)hBHn, dHdr, rp);
    if (gBubActive && gBubLive > 0)
        kSplatPsi<<<(NBUB*QB3 + 255)/256, b1, 0, prsS>>>(dBub, dBubOwner, dBubCtr, dHdr, rp);
    const float4* rsrc = dHdr;
    if (rp.nLens > 0){
        dim3 b2l(16, 16), g2l((RW + 15)/16, (RH + 15)/16);
        kLens<<<g2l, b2l, 0, prsS>>>(dHdr, dHdrW, RW, RH, rp);
        rsrc = dHdrW;
    }
    // auto-exposure histogram
    CUDA_CHECK(cudaMemsetAsync(dHist, 0, 256*sizeof(unsigned), prsS));
    kHist<<<g1, b1, 0, prsS>>>(rsrc, npx, dHist);
    CUDA_CHECK(cudaMemcpyAsync(hHist, dHist, 256*sizeof(unsigned),
                               cudaMemcpyDeviceToHost, prsS));
    // bloom mip chain (threshold-free)
    if (gBloomOn){
        dim3 b2(16, 16);
        const float4* src = rsrc; int sw = RW, sh = RH;
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
    kComposite<<<g3, b3, 0, prsS>>>(rsrc, dMip[1], mipW[1], mipH[1], dOut, rp);
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
    emit(16, 14, "TINY UNIVERSE M7 COSMOS", 0);
    snprintf(line, sizeof line, "N=%d DT=1/240 TICK=%llu", gN, (unsigned long long)gTick);
    emit(16, 32, line, 0);
    snprintf(line, sizeof line, "SIM T=%.1FS  C=20 HBAR=0.5 G=0.002", gSimTime);
    emit(16, 50, line, 0);
    if (gMode) emit(16, 68, "PHYSICAL: L=T4 SAT=1.0 NO CLIP", 0);
    else       emit(16, 68, "CINEMATIC: L=T2.2 SAT=1.2 STRETCH", 0);
    static const char* SC_HUD[20] = {"GALAXY", "KEPLER", "THREEBODY", "CLOUD",
                                     "MERGER", "ECHO", "RATCHET", "DETECTOR",
                                     "KEPREL", "CLOCKS", "PHOTONS",
                                     "COLLAPSE", "ISCO", "HAWKING",
                                     "DOUBLESLIT", "TUNNELING", "SHOQ",
                                     "CIRCUMNAV", "EXPAND", "BUBBLES"};
    snprintf(line, sizeof line, "%s %s  DE %+.1E  DP %.1E", SC_HUD[gScenario],
             gSolver == SOLV_PM ? "PM" : (gSolver == SOLV_TINY ? "TINY" :
             (gSolver == SOLV_NONE ? "NONE" : "DIRECT")), gDE, gDP);
    emit(16, 86, line, 0);
    if (gScenario == SC_MERGER || gScenario == SC_ECHO){
        snprintf(line, sizeof line, "S %.3F  DS/DT %+.1E", gS, gdS);
        emit(16, 104, line, 0);
    }
    if (gBHActive){
        if (hBHn > 0)
            snprintf(line, sizeof line, "BH N=%u  M=%.3E  TH=%.2F", hBHn,
                     (double)hBHpos[0].w, hBHpos[0].w > 0 ? 79577.0/hBHpos[0].w : 0.0);
        else
            snprintf(line, sizeof line, "BH N=0 (EVAPORATED)");
        emit(16, 104, line, 0);
    }
    if (gExpand){
        snprintf(line, sizeof line, "A(T) %.3F  H0 %.4F  COMOVING", aOfT((double)gTick/240.0), gH0);
        emit(16, 122, line, 0);
    }
    if (gBubActive){
        snprintf(line, sizeof line, "BUBBLES %d LIVE  %u SPAWNED  %u COLLAPSED",
                 gBubLive, gBubSpawned, gBubCollapsed);
        emit(16, 122, line, 0);
    }
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
    snprintf(line, sizeof line, "%s  HIST %.1FS", gProj ? "PLANET" : "PERSP",
             gHistOn ? gHistN*0.1 : 0.0);
    emit(rx, 86, line, 1);
    int n = (int)gs.size(); if (!n) return;
    if (n > 512) n = 512;
    CUDA_CHECK(cudaMemcpyAsync(dGlyphs, gs.data(), n*sizeof(Glyph),
                               cudaMemcpyHostToDevice, prsS));
    kHud<<<n, 70, 0, prsS>>>(dOut, gOW, gOH, dGlyphs, n);
}

// ----------------------------------------------------------------------------
// Sim ticks (simStream, ping-pong publish)
// ----------------------------------------------------------------------------
// coarse-grained entropy S = S_x + S_v (arrow contract: declared machinery)
static double entropyNow(){
    CUDA_CHECK(cudaMemsetAsync(dCntX, 0, 32768*sizeof(unsigned), simS));
    CUDA_CHECK(cudaMemsetAsync(dCntV, 0, 32768*sizeof(unsigned), simS));
    kCount<<<(gN + 255)/256, 256, 0, simS>>>(dPos[gPub], dMom, gN, dCntX, dCntV);
    static std::vector<unsigned> hx(32768), hv(32768);
    CUDA_CHECK(cudaStreamSynchronize(simS));
    CUDA_CHECK(cudaMemcpy(hx.data(), dCntX, 32768*4, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(hv.data(), dCntV, 32768*4, cudaMemcpyDeviceToHost));
    double S = 0;
    for (int pass = 0; pass < 2; pass++){
        const std::vector<unsigned>& h = pass ? hv : hx;
        for (int i = 0; i < 32768; i++){
            if (!h[i]) continue;
            double q = (double)h[i]/(double)gN;
            S -= q*log(q);
        }
    }
    return S;
}

// --- M7 cosmos: bubble host lifecycle ---
static void bubAllocOnce(){
    if (gBubAlloc) return;
    CUDA_CHECK(cudaMalloc(&dBub, (size_t)NBUB*QB3*sizeof(cufftComplex)));
    CUDA_CHECK(cudaMemset(dBub, 0, (size_t)NBUB*QB3*sizeof(cufftComplex)));
    CUDA_CHECK(cudaMalloc(&dBubEdge, (size_t)QB3*sizeof(float)));
    CUDA_CHECK(cudaMalloc(&dBubOwner, NBUB*sizeof(int)));
    CUDA_CHECK(cudaMalloc(&dBubCtr, NBUB*sizeof(float4)));
    CUDA_CHECK(cudaMemset(dBubCtr, 0, NBUB*sizeof(float4)));
    CUDA_CHECK(cudaMalloc(&dBubRec, NBUB*sizeof(unsigned)));
    CUDA_CHECK(cudaMemset(dBubRec, 0, NBUB*sizeof(unsigned)));
    CUDA_CHECK(cudaMalloc(&dBubFlag, NBUB*sizeof(unsigned)));
    CUDA_CHECK(cudaMemset(dBubFlag, 0, NBUB*sizeof(unsigned)));
    CUDA_CHECK(cudaMalloc(&dBubPend, sizeof(unsigned)));
    CUDA_CHECK(cudaMalloc(&dBubNew, NBUB*sizeof(int)));
    CUDA_CHECK(cudaMalloc(&dBubNewCount, sizeof(int)));
    CUDA_CHECK(cudaMalloc(&dBubNorm, sizeof(unsigned long long)));
    CUDA_CHECK(cudaMalloc(&dIso, (size_t)gN));
    int minus1[NBUB];
    for (int s = 0; s < NBUB; s++){ minus1[s] = -1; hBubOwner[s] = -1; }
    CUDA_CHECK(cudaMemcpy(dBubOwner, minus1, sizeof minus1, cudaMemcpyHostToDevice));
    // cos^2 absorbing border, 8 cells per face (M6 idiom, 3D product)
    std::vector<float> m((size_t)QB3);
    auto s1 = [](int i){
        const int B = 8;
        float t = 1.0f;
        if (i < B) t = sinf(PI_F*0.5f*(i + 0.5f)/B);
        else if (i >= QB - B) t = sinf(PI_F*0.5f*(QB - 0.5f - i)/B);
        return t*t;
    };
    for (int iz = 0; iz < QB; iz++)
        for (int iy = 0; iy < QB; iy++)
            for (int ix = 0; ix < QB; ix++)
                m[((size_t)iz*QB + iy)*QB + ix] = s1(ix)*s1(iy)*s1(iz);
    CUDA_CHECK(cudaMemcpy(dBubEdge, m.data(), m.size()*4, cudaMemcpyHostToDevice));
    int dims[3] = {QB, QB, QB};
    CUFFT_CHECK(cufftPlanMany(&planB, 3, dims, NULL, 1, QB3, NULL, 1, QB3,
                              CUFFT_C2C, NBUB));
    CUFFT_CHECK(cufftSetStream(planB, simS));
    gBubAlloc = true;
}
// one split-step tick for all slots (V = 0 free evolution: single kinetic step is
// exact; dead slots are zero grids). Runs on simS, per tick, when any bubble lives.
static void bubEvolveTick(){
    if (gBubLive == 0) return;
    const int total = NBUB*QB3;
    dim3 b(256), g((total + 255)/256);
    const float coefK = 0.25f*DT;                     // hbar/(2m)*dt
    CUFFT_CHECK(cufftExecC2C(planB, dBub, dBub, CUFFT_FORWARD));
    kQ3PhaseK<<<g, b, 0, simS>>>(dBub, total, coefK);
    CUFFT_CHECK(cufftExecC2C(planB, dBub, dBub, CUFFT_INVERSE));
    kQ3Scale<<<g, b, 0, simS>>>(dBub, total, 1.0f/(float)QB3);
    kQ3Edge<<<g, b, 0, simS>>>(dBub, total, dBubEdge);
}
// collapse by inscription: sample the 3D |psi|^2 CDF (host fp64, counter-keyed),
// re-pointify the owner there, free the slot (cosmos contract)
static void bubCollapse(int slot, int owner){
    std::vector<cufftComplex> h((size_t)QB3);
    CUDA_CHECK(cudaStreamSynchronize(simS));
    CUDA_CHECK(cudaMemcpy(h.data(), dBub + (size_t)slot*QB3, (size_t)QB3*sizeof(cufftComplex),
                          cudaMemcpyDeviceToHost));
    std::vector<double> cdf((size_t)QB3);
    double s = 0;
    for (int c = 0; c < QB3; c++){
        s += (double)h[c].x*h[c].x + (double)h[c].y*h[c].y;
        cdf[c] = s;
    }
    double u = orrery::counter_uniform(gSeed ^ 0xB0BB1Eull, owner, gTick, 10) * s;
    int lo = 0, hi = QB3 - 1;
    while (lo < hi){ int mid = (lo + hi)/2; if (cdf[mid] < u) lo = mid + 1; else hi = mid; }
    int ix = lo % QB, iy = (lo/QB) % QB, iz = lo/(QB*QB);
    double u1 = orrery::counter_uniform(gSeed ^ 0xB0BB1Eull, owner, gTick, 11);
    double u2 = orrery::counter_uniform(gSeed ^ 0xB0BB1Eull, owner, gTick, 12);
    double u3 = orrery::counter_uniform(gSeed ^ 0xB0BB1Eull, owner, gTick, 13);
    float4 ctr;
    CUDA_CHECK(cudaMemcpy(&ctr, dBubCtr + slot, sizeof ctr, cudaMemcpyDeviceToHost));
    float4 np = make_float4(
        ctr.x + (float)(-QBL*0.5 + (ix + u1)*QBDX),
        ctr.y + (float)(-QBL*0.5 + (iy + u2)*QBDX),
        ctr.z + (float)(-QBL*0.5 + (iz + u3)*QBDX), ctr.w);
    if      (np.x >=  256.0f) np.x -= 512.0f; else if (np.x < -256.0f) np.x += 512.0f;
    if      (np.y >=  256.0f) np.y -= 512.0f; else if (np.y < -256.0f) np.y += 512.0f;
    if      (np.z >=  256.0f) np.z -= 512.0f; else if (np.z < -256.0f) np.z += 512.0f;
    CUDA_CHECK(cudaMemcpy(dPos[gPub] + owner, &np, sizeof np, cudaMemcpyHostToDevice));
    float4 zero4 = make_float4(0, 0, 0, 0);
    CUDA_CHECK(cudaMemcpy(dPosErr + owner, &zero4, sizeof zero4, cudaMemcpyHostToDevice));
    unsigned reg = 0xC2u;                             // RECORDED|INSCRIBED latched + classical
    CUDA_CHECK(cudaMemcpy(dRegime + owner, &reg, sizeof reg, cudaMemcpyHostToDevice));
    unsigned zu = 0;
    CUDA_CHECK(cudaMemcpy(dPartBub + owner, &zu, sizeof zu, cudaMemcpyHostToDevice));
    int m1 = -1;
    CUDA_CHECK(cudaMemcpy(dBubOwner + slot, &m1, sizeof m1, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemset(dBub + (size_t)slot*QB3, 0, (size_t)QB3*sizeof(cufftComplex)));
    hBubOwner[slot] = -1;
    gBubCollapsed++; gBubLive--;
}
// per-axis 2nd central moment of |psi|^2, averaged over axes (host fp64 observable)
static double bubSigma(int slot){
    std::vector<cufftComplex> h((size_t)QB3);
    CUDA_CHECK(cudaStreamSynchronize(simS));
    CUDA_CHECK(cudaMemcpy(h.data(), dBub + (size_t)slot*QB3, (size_t)QB3*sizeof(cufftComplex),
                          cudaMemcpyDeviceToHost));
    double n = 0, mx[3] = {0,0,0}, mxx[3] = {0,0,0};
    for (int c = 0; c < QB3; c++){
        double p = (double)h[c].x*h[c].x + (double)h[c].y*h[c].y;
        double r[3] = { -QBL*0.5 + ((c % QB) + 0.5)*QBDX,
                        -QBL*0.5 + (((c/QB) % QB) + 0.5)*QBDX,
                        -QBL*0.5 + ((c/(QB*QB)) + 0.5)*QBDX };
        n += p;
        for (int a = 0; a < 3; a++){ mx[a] += p*r[a]; mxx[a] += p*r[a]*r[a]; }
    }
    double sig = 0;
    for (int a = 0; a < 3; a++){
        double m1 = mx[a]/n;
        sig += sqrt(mxx[a]/n - m1*m1);
    }
    return sig/3.0;
}
// per-tick bubble hooks: evolve -> intrude -> records/collapse; spawn checks every 24
static void bubTickHost(){
    dim3 b(256), gP((gN + 255)/256);
    bubEvolveTick();
    if (gBubLive > 0){
        kBubIntrude<<<gP, b, 0, simS>>>(dPos[gPub], dMom, dRegime, dPartBub, gN,
                                        dBubOwner, dBubCtr, dBubFlag);
        kBubTick<<<1, 1, 0, simS>>>(dBubOwner, dBubRec, dBubFlag, dBubPend);
        unsigned pend = 0;
        CUDA_CHECK(cudaStreamSynchronize(simS));
        CUDA_CHECK(cudaMemcpy(&pend, dBubPend, sizeof pend, cudaMemcpyDeviceToHost));
        for (int s = 0; s < NBUB; s++)
            if (pend & (1u << s)) bubCollapse(s, hBubOwner[s]);
    }
    if ((gTick % 24) == 0 && gBubLive < NBUB){        // isolation check cadence
        if (gSolver != SOLV_PM){                      // no PM force pass: deposit explicitly
            int nc = PMN*PMN*PMN;
            kZeroFix<<<(nc + 255)/256, b, 0, simS>>>(dFix, nc);
            kDeposit<<<gP, b, 0, simS>>>(dPos[gPub], dMom, dRegime, gN, dFix);
        }
        kBubIso<<<gP, b, 0, simS>>>(dPos[gPub], dMom, dRegime, dPartBub, gN, dFix, dIso);
        kBubSpawnScan<<<1, 1, 0, simS>>>(dIso, dPos[gPub], dRegime, dPartBub, gN,
                                         dBubOwner, dBubCtr, dBubRec, dBubNew, dBubNewCount);
        int nn = 0, slots[NBUB];
        CUDA_CHECK(cudaStreamSynchronize(simS));
        CUDA_CHECK(cudaMemcpy(&nn, dBubNewCount, sizeof nn, cudaMemcpyDeviceToHost));
        if (nn > 0){
            CUDA_CHECK(cudaMemcpy(slots, dBubNew, nn*sizeof(int), cudaMemcpyDeviceToHost));
            int owners[NBUB];
            CUDA_CHECK(cudaMemcpy(owners, dBubOwner, sizeof owners, cudaMemcpyDeviceToHost));
            dim3 gB((QB3 + 255)/256);
            for (int q = 0; q < nn; q++){
                int slot = slots[q], owner = owners[slot];
                hBubOwner[slot] = owner;
                float4 pm;
                CUDA_CHECK(cudaMemcpy(&pm, dMom + owner, sizeof pm, cudaMemcpyDeviceToHost));
                // psi0: Gaussian sigma0 = 1 su with phase k = p/hbar = 2p (hbar = 0.5, m = 1)
                kQ3Gauss<<<gB, b, 0, simS>>>(dBub + (size_t)slot*QB3,
                                             2.0f*pm.x, 2.0f*pm.y, 2.0f*pm.z, 1.0f);
                CUDA_CHECK(cudaMemsetAsync(dBubNorm, 0, 8, simS));
                kQ3NormAcc<<<gB, b, 0, simS>>>(dBub + (size_t)slot*QB3, dBubNorm);
                unsigned long long hn;
                CUDA_CHECK(cudaStreamSynchronize(simS));
                CUDA_CHECK(cudaMemcpy(&hn, dBubNorm, 8, cudaMemcpyDeviceToHost));
                double nrm = orrery::fixed_decode((long long)hn)*(double)QBDX*QBDX*QBDX;
                kQ3Scale<<<gB, b, 0, simS>>>(dBub + (size_t)slot*QB3, QB3,
                                             (float)(1.0/sqrt(nrm)));
                gBubSpawned++; gBubLive++;
            }
        }
    }
}

// M7 cosmos: |rho_hat(k1)| growth observable (expand scenario). Measurement-only
// deposit + forward FFT; it clobbers dReal/dSpec, so the caller must re-run
// forcePass() afterwards to restore the dAcc/Phi invariant.
static double expandAmp(){
    int nc = PMN*PMN*PMN;
    dim3 b(256);
    kZeroFix<<<(nc + 255)/256, b, 0, simS>>>(dFix, nc);
    kDeposit<<<(gN + 255)/256, b, 0, simS>>>(dPos[gPub], dMom, dRegime, gN, dFix);
    kFixToReal<<<(nc + 255)/256, b, 0, simS>>>(dFix, dReal, nc);
    CUFFT_CHECK(cufftExecR2C(planF, dReal, dSpec));
    cufftComplex h;
    CUDA_CHECK(cudaStreamSynchronize(simS));
    CUDA_CHECK(cudaMemcpy(&h, dSpec + (size_t)(4*PMN + 0)*PMNZC + 0, sizeof h,
                          cudaMemcpyDeviceToHost));
    return sqrt((double)h.x*h.x + (double)h.y*h.y);
}

// --- M6 planck: host side of the psi engine ---
static void hfft(std::vector<std::complex<double>>& a, bool inv){   // radix-2, in-place
    const size_t n = a.size();
    for (size_t i = 1, j = 0; i < n; i++){
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }
    for (size_t len = 2; len <= n; len <<= 1){
        double ang = 2.0*3.14159265358979323846/(double)len*(inv ? 1.0 : -1.0);
        std::complex<double> wl(cos(ang), sin(ang));
        for (size_t i = 0; i < n; i += len){
            std::complex<double> w(1, 0);
            for (size_t k = 0; k < len/2; k++){
                std::complex<double> u = a[i+k], v = a[i+k+len/2]*w;
                a[i+k] = u + v; a[i+k+len/2] = u - v;
                w *= wl;
            }
        }
    }
    if (inv) for (auto& x : a) x /= (double)n;
}
static void qAllocOnce(){
    if (gQAlloc) return;
    CUDA_CHECK(cudaMalloc(&dQ, (size_t)QN*QN*sizeof(cufftComplex)));
    CUDA_CHECK(cudaMalloc(&dQWall, (size_t)QN*QN*sizeof(float)));
    CUDA_CHECK(cudaMalloc(&dQEdge, (size_t)QN*QN*sizeof(float)));
    CUDA_CHECK(cudaMalloc(&dQV,    (size_t)QN*QN*sizeof(float)));
    CUDA_CHECK(cudaMalloc(&dQNorm, sizeof(unsigned long long)));
    CUFFT_CHECK(cufftPlan2d(&planQ, QN, QN, CUFFT_C2C));
    CUFFT_CHECK(cufftSetStream(planQ, simS));
    gQAlloc = true;
}
static double qNormHost(){
    CUDA_CHECK(cudaMemsetAsync(dQNorm, 0, 8, simS));
    kQNormAcc<<<(QN*QN + 255)/256, 256, 0, simS>>>(dQ, dQNorm);
    unsigned long long h;
    CUDA_CHECK(cudaStreamSynchronize(simS));
    CUDA_CHECK(cudaMemcpy(&h, dQNorm, 8, cudaMemcpyDeviceToHost));
    return orrery::fixed_decode((long long)h)*(double)QDX*QDX;
}
static void qNormalize(){
    double n = qNormHost();
    kQScale<<<(QN*QN + 255)/256, 256, 0, simS>>>(dQ, (float)(1.0/sqrt(n)));
}
static void qDownload(std::vector<cufftComplex>& h){
    h.resize((size_t)QN*QN);
    CUDA_CHECK(cudaStreamSynchronize(simS));
    CUDA_CHECK(cudaMemcpy(h.data(), dQ, (size_t)QN*QN*sizeof(cufftComplex),
                          cudaMemcpyDeviceToHost));
}
// Strang split-step, real time (dt = dial): V half, wall, K, V half, wall, edge
static void qStepReal(int steps, bool useV, bool useWall){
    dim3 b(256), g((QN*QN + 255)/256);
    const float coefK = 0.25f*DT;                     // hbar/(2m)*dt, hbar=0.5 m=1
    const float halfV = DT;                            // V*dt/(2*hbar) = V*dt (1/(2hbar)=1)
    for (int s = 0; s < steps; s++){
        if (useV)    kQPhaseV<<<g, b, 0, simS>>>(dQ, dQV, halfV);
        if (useWall) kQMul<<<g, b, 0, simS>>>(dQ, dQWall);
        CUFFT_CHECK(cufftExecC2C(planQ, dQ, dQ, CUFFT_FORWARD));
        kQPhaseK<<<g, b, 0, simS>>>(dQ, coefK);
        CUFFT_CHECK(cufftExecC2C(planQ, dQ, dQ, CUFFT_INVERSE));
        kQScale<<<g, b, 0, simS>>>(dQ, 1.0f/(QN*QN));
        if (useV)    kQPhaseV<<<g, b, 0, simS>>>(dQ, dQV, halfV);
        if (useWall) kQMul<<<g, b, 0, simS>>>(dQ, dQWall);
        kQMul<<<g, b, 0, simS>>>(dQ, dQEdge);
    }
}
static void qStepImag(int iters, float tau){
    dim3 b(256), g((QN*QN + 255)/256);
    const float coefK = 0.25f*tau;                    // hbar/(2m)*tau
    const float halfV = tau;                          // V*tau/(2*hbar)
    for (int s = 0; s < iters; s++){
        kQDecayV<<<g, b, 0, simS>>>(dQ, dQV, halfV);
        CUFFT_CHECK(cufftExecC2C(planQ, dQ, dQ, CUFFT_FORWARD));
        kQDecayK<<<g, b, 0, simS>>>(dQ, coefK);
        CUFFT_CHECK(cufftExecC2C(planQ, dQ, dQ, CUFFT_INVERSE));
        kQScale<<<g, b, 0, simS>>>(dQ, 1.0f/(QN*QN));
        kQDecayV<<<g, b, 0, simS>>>(dQ, dQV, halfV);
        qNormalize();
    }
}
// screen-window y-marginal, normalized to sum 1
static std::vector<double> qMarginalY(const std::vector<cufftComplex>& h,
                                      double xlo, double xhi){
    std::vector<double> P(QN, 0.0);
    double tot = 0;
    for (int iy = 0; iy < QN; iy++)
        for (int ix = 0; ix < QN; ix++){
            double x = -QL*0.5 + (ix + 0.5)*QDX;
            if (x < xlo || x > xhi) continue;
            double p = (double)h[iy*QN + ix].x*h[iy*QN + ix].x
                     + (double)h[iy*QN + ix].y*h[iy*QN + ix].y;
            P[iy] += p; tot += p;
        }
    if (tot > 0) for (auto& p : P) p /= tot;
    return P;
}
static double qYof(int iy){ return -QL*0.5 + (iy + 0.5)*QDX; }
static double qContrast(const std::vector<double>& yv, double kstar){   // 2|sum e^{-iky}|/N
    double re = 0, im = 0;
    for (double y : yv){ re += cos(kstar*y); im -= sin(kstar*y); }
    return 2.0*sqrt(re*re + im*im)/(double)yv.size();
}
static double qContrastP(const std::vector<double>& P, double kstar){
    double re = 0, im = 0;
    for (int iy = 0; iy < QN; iy++){ re += P[iy]*cos(kstar*qYof(iy)); im -= P[iy]*sin(kstar*qYof(iy)); }
    return 2.0*sqrt(re*re + im*im);
}
static double qSampleY(const std::vector<double>& cdf, double u){       // inverse CDF
    int lo = 0, hi = QN - 1;
    while (lo < hi){ int mid = (lo + hi)/2; if (cdf[mid] < u) lo = mid + 1; else hi = mid; }
    double c1 = cdf[lo], c0 = lo ? cdf[lo-1] : 0.0;
    double f = (c1 > c0) ? (u - c0)/(c1 - c0) : 0.5;
    return -QL*0.5 + (lo + f)*QDX;
}
static void qBuildEdgeMask(){
    std::vector<float> m((size_t)QN*QN);
    auto s1 = [](int i){
        const int B = 8;
        float t = 1.0f;
        if (i < B) t = sinf(PI_F*0.5f*(i + 0.5f)/B);
        else if (i >= QN - B) t = sinf(PI_F*0.5f*(QN - 0.5f - i)/B);
        return t*t;
    };
    for (int iy = 0; iy < QN; iy++)
        for (int ix = 0; ix < QN; ix++) m[(size_t)iy*QN + ix] = s1(ix)*s1(iy);
    CUDA_CHECK(cudaMemcpy(dQEdge, m.data(), m.size()*4, cudaMemcpyHostToDevice));
}
static void qBuildWall(int mode){                     // 0 both slits · 1 A only · 2 B only
    std::vector<float> w((size_t)QN*QN, 1.0f);
    for (int iy = 0; iy < QN; iy++)
        for (int ix = 0; ix < QN; ix++){
            double x = -QL*0.5 + (ix + 0.5)*QDX;
            if (x < 0.0 || x > 1.0) continue;
            double y = qYof(iy);
            bool inA = (y >  1.25 && y <  2.75);
            bool inB = (y > -2.75 && y < -1.25);
            bool open = (mode == 0) ? (inA || inB) : (mode == 1 ? inA : inB);
            if (!open) w[(size_t)iy*QN + ix] = 0.0f;
        }
    CUDA_CHECK(cudaMemcpy(dQWall, w.data(), w.size()*4, cudaMemcpyHostToDevice));
}
static void qRunDoubleslit(std::vector<double>& dotsY){
    qAllocOnce(); qBuildEdgeMask();
    const double p0 = 2.0, hb = 0.5, d = 4.0, Lsc = 14.0;
    const double kstar = d*p0/(hb*Lsc);
    const int NDOTS = 4096, STEPS = 3600;
    dim3 b(256), g((QN*QN + 255)/256);
    auto evolve = [&](int wallMode, std::vector<cufftComplex>& out){
        qBuildWall(wallMode);
        kQGauss<<<g, b, 0, simS>>>(dQ, -15.0f, 0.0f, 1.5f, 4.0f, (float)(p0/hb));
        qNormalize();
        qStepReal(STEPS, false, true);
        qDownload(out);
    };
    std::vector<cufftComplex> hBoth, hA, hB;
    evolve(0, hBoth); evolve(1, hA); evolve(2, hB);
    auto Pb = qMarginalY(hBoth, 12.0, 16.0);
    auto Pa = qMarginalY(hA, 12.0, 16.0);
    auto Pc = qMarginalY(hB, 12.0, 16.0);
    gQM.c_psi = qContrastP(Pb, kstar);
    double bestk = 0, bestv = -1;
    for (int q = 0; q <= 100; q++){
        double k = kstar*(0.5 + q*0.01);
        double v = qContrastP(Pb, k);
        if (v > bestv){ bestv = v; bestk = k; }
    }
    gQM.kpk_rel = fabs(bestk/kstar - 1.0);
    auto mkcdf = [](const std::vector<double>& P){
        std::vector<double> c(QN); double s = 0;
        for (int i = 0; i < QN; i++){ s += P[i]; c[i] = s; }
        for (int i = 0; i < QN; i++) c[i] /= (s > 0 ? s : 1.0);
        return c;
    };
    auto cb = mkcdf(Pb), ca = mkcdf(Pa), cc = mkcdf(Pc);
    dotsY.resize(NDOTS);
    std::vector<double> detY(NDOTS);
    for (int j = 0; j < NDOTS; j++){
        dotsY[j] = qSampleY(cb, orrery::counter_uniform(gSeed, j, 0, 8));
        bool useA = orrery::counter_uniform(gSeed, j, 1, 8) < 0.5;   // symmetric slits (declared)
        detY[j] = qSampleY(useA ? ca : cc, orrery::counter_uniform(gSeed, j, 2, 8));
    }
    gQM.c_dots = qContrast(dotsY, kstar);
    gQM.c_det  = qContrast(detY, kstar);
    gQM.nrec = NDOTS;
    gPsiBytes.assign((const char*)hBoth.data(), hBoth.size()*sizeof(cufftComplex));
}
static void qRunTunneling(){
    qAllocOnce(); qBuildEdgeMask();
    const double V0 = 1.8, p0 = 1.5;
    const int STEPS = 2400;
    std::vector<float> V((size_t)QN*QN, 0.0f);
    for (int iy = 0; iy < QN; iy++)
        for (int ix = 0; ix < QN; ix++){
            double x = -QL*0.5 + (ix + 0.5)*QDX;
            if (x >= 0.0 && x <= 1.0) V[(size_t)iy*QN + ix] = (float)V0;
        }
    CUDA_CHECK(cudaMemcpy(dQV, V.data(), V.size()*4, cudaMemcpyHostToDevice));
    dim3 b(256), g((QN*QN + 255)/256);
    kQGauss<<<g, b, 0, simS>>>(dQ, -8.0f, 0.0f, 2.0f, 6.0f, (float)(p0/0.5));
    qNormalize();
    qStepReal(STEPS, true, false);
    std::vector<cufftComplex> h;
    qDownload(h);
    double num = 0, tot = 0;
    for (int iy = 0; iy < QN; iy++)
        for (int ix = 0; ix < QN; ix++){
            double x = -QL*0.5 + (ix + 0.5)*QDX;
            double p = (double)h[(size_t)iy*QN+ix].x*h[(size_t)iy*QN+ix].x
                     + (double)h[(size_t)iy*QN+ix].y*h[(size_t)iy*QN+ix].y;
            tot += p; if (x > 1.5) num += p;
        }
    gQM.T = num/tot;
    // in-scenario oracle: host fp64 1D split-step, SAME dx as the GPU grid —
    // tunneling is exponentially sensitive to barrier-edge discretization, so the
    // oracle isolates implementation (fp32 vs fp64, 2D vs 1D), not grid resolution
    const int n1 = QN; const double dx1 = (double)QL/n1;
    std::vector<std::complex<double>> psi(n1);
    for (int i = 0; i < n1; i++){
        double x = -QL*0.5 + (i + 0.5)*dx1;
        double dxp = (x + 8.0)/2.0;
        psi[i] = std::polar(exp(-0.25*dxp*dxp), (p0/0.5)*x);
    }
    { double s = 0; for (auto& c : psi) s += std::norm(c); s = sqrt(s*dx1); for (auto& c : psi) c /= s; }
    std::vector<double> V1(n1, 0.0), E1(n1);
    const int B1 = (int)(2.0/dx1);
    for (int i = 0; i < n1; i++){
        double x = -QL*0.5 + (i + 0.5)*dx1;
        if (x >= 0.0 && x <= 1.0) V1[i] = V0;
        double t = 1.0;
        if (i < B1) t = sin(3.14159265358979324*0.5*(i + 0.5)/B1);
        else if (i >= n1 - B1) t = sin(3.14159265358979324*0.5*(n1 - 0.5 - i)/B1);
        E1[i] = t*t;
    }
    const double dtq = 1.0/240.0, coefK = 0.25*dtq, halfV = dtq;
    for (int s = 0; s < STEPS; s++){
        for (int i = 0; i < n1; i++) psi[i] *= std::polar(1.0, -V1[i]*halfV);
        hfft(psi, false);
        for (int i = 0; i < n1; i++){
            int f = (i <= n1/2) ? i : i - n1;
            double k = 2.0*3.14159265358979324*f/QL;
            psi[i] *= std::polar(1.0, -coefK*k*k);
        }
        hfft(psi, true);
        for (int i = 0; i < n1; i++) psi[i] *= std::polar(1.0, -V1[i]*halfV);
        for (int i = 0; i < n1; i++) psi[i] *= E1[i];
    }
    double nrm = 0, trn = 0;
    for (int i = 0; i < n1; i++){
        double x = -QL*0.5 + (i + 0.5)*dx1;
        double p = std::norm(psi[i]);
        nrm += p; if (x > 1.5) trn += p;
    }
    gQM.Tref = trn/nrm;
    gPsiBytes.assign((const char*)h.data(), h.size()*sizeof(cufftComplex));
}
static void qRunShoq(){
    qAllocOnce();
    std::vector<float> V((size_t)QN*QN);
    for (int iy = 0; iy < QN; iy++)
        for (int ix = 0; ix < QN; ix++){
            double x = -QL*0.5 + (ix + 0.5)*QDX, y = qYof(iy);
            V[(size_t)iy*QN + ix] = (float)(0.5*(x*x + y*y));
        }
    CUDA_CHECK(cudaMemcpy(dQV, V.data(), V.size()*4, cudaMemcpyHostToDevice));
    dim3 b(256), g((QN*QN + 255)/256);
    kQGauss<<<g, b, 0, simS>>>(dQ, 0.4f, -0.3f, 1.5f, 1.2f, 0.0f);
    qNormalize();
    qStepImag(15000, 0.002f);
    std::vector<cufftComplex> h;
    qDownload(h);
    std::vector<std::complex<double>> a((size_t)QN*QN);
    double nn = 0, vv = 0, mx = 0, mxx = 0;
    for (int iy = 0; iy < QN; iy++)
        for (int ix = 0; ix < QN; ix++){
            std::complex<double> c(h[(size_t)iy*QN+ix].x, h[(size_t)iy*QN+ix].y);
            a[(size_t)iy*QN + ix] = c;
            double p = std::norm(c);
            double x = -QL*0.5 + (ix + 0.5)*QDX;
            nn += p; vv += p*V[(size_t)iy*QN+ix]; mx += p*x; mxx += p*x*x;
        }
    vv /= nn; mx /= nn; mxx /= nn;
    gQM.sigx = sqrt(mxx - mx*mx);
    std::vector<std::complex<double>> row(QN), col(QN);
    for (int iy = 0; iy < QN; iy++){
        for (int ix = 0; ix < QN; ix++) row[ix] = a[(size_t)iy*QN + ix];
        hfft(row, false);
        for (int ix = 0; ix < QN; ix++) a[(size_t)iy*QN + ix] = row[ix];
    }
    for (int ix = 0; ix < QN; ix++){
        for (int iy = 0; iy < QN; iy++) col[iy] = a[(size_t)iy*QN + ix];
        hfft(col, false);
        for (int iy = 0; iy < QN; iy++) a[(size_t)iy*QN + ix] = col[iy];
    }
    double kn = 0, tt = 0;
    for (int iy = 0; iy < QN; iy++)
        for (int ix = 0; ix < QN; ix++){
            int fx = (ix <= QN/2) ? ix : ix - QN;
            int fy = (iy <= QN/2) ? iy : iy - QN;
            double kf = 2.0*3.14159265358979324/QL;
            double k2 = kf*kf*((double)fx*fx + (double)fy*fy);
            double p = std::norm(a[(size_t)iy*QN + ix]);
            kn += p; tt += p*0.125*k2;                 // hbar^2 k^2 / 2m = 0.25*k^2/2
        }
    gQM.E = tt/kn + vv;
    gPsiBytes.assign((const char*)h.data(), h.size()*sizeof(cufftComplex));
}

// invariant: dAcc = a(dPos[gPub]) after every full tick and after init.
// tState = sim time of the position state being differentiated (drives a(t) in expand).
static void forcePass(double tState){
    dim3 b(256);
    if (gSolver == SOLV_NONE){
        if (!gBHActive) return;                       // pure drift (dAcc stays zero)
        CUDA_CHECK(cudaMemsetAsync(dAcc, 0, (size_t)gN*sizeof(float4), simS));
        kBHForce<<<(gN + 255)/256, b, 0, simS>>>(dPos[gPub], dRegime, dAcc, gN,
                                                 dBHpos, dBHmom, dBHn);
        return;
    }
    if (gSolver == SOLV_PM){
        int nc = PMN*PMN*PMN, ns = PMN*PMN*PMNZC;
        kZeroFix<<<(nc + 255)/256, b, 0, simS>>>(dFix, nc);
        kDeposit<<<(gN + 255)/256, b, 0, simS>>>(dPos[gPub], dMom, dRegime, gN, dFix);
        kFixToReal<<<(nc + 255)/256, b, 0, simS>>>(dFix, dReal, nc);
        CUFFT_CHECK(cufftExecR2C(planF, dReal, dSpec));
        kGreen<<<(ns + 255)/256, b, 0, simS>>>(dSpec, (float)(1.0/aOfT(tState)));
        CUFFT_CHECK(cufftExecC2R(planI, dSpec, dReal));   // dReal now holds Phi
        kForceGrid<<<(nc + 255)/256, b, 0, simS>>>(dReal, dFor);
        kGather<<<(gN + 255)/256, b, 0, simS>>>(dPos[gPub], dFor, dReal, dAcc, gN);
    } else if (gSolver == SOLV_DIRECT){
        kDirect<<<(gN + 255)/256, b, 0, simS>>>(dPos[gPub], dMom, dAcc, gN);
    }
    if (gBHActive)
        kBHForce<<<(gN + 255)/256, b, 0, simS>>>(dPos[gPub], dRegime, dAcc, gN,
                                                 dBHpos, dBHmom, dBHn);
}

static void runTicks(int nt){
    const float invc2 = 1.0f/(C_LIGHT*C_LIGHT), hdt = 0.5f*DT;
    if (gSolver == SOLV_TINY){
        int src = gPub, dst = 1 - gPub;
        kTiny<<<1, 32, 0, simS>>>(dPos[src], dPos[dst], dMom, dTau, dRegime,
                                  gN, nt, DT, invc2);
        gPub = dst; gTick += nt; gSimTime += (double)nt*DT;
    } else {
        dim3 b(256), g((gN + 255)/256);
        const float* phiP = (gSolver == SOLV_PM) ? dReal : nullptr;
        for (int t = 0; t < nt; t++){
            kKick<<<g, b, 0, simS>>>(dMom, dAcc, dRegime, gN, hdt);
            int src = gPub, dst = 1 - gPub;
            float ds = 1.0f;                          // expand: drift x_dot = p/a^2 at mid-step
            if (gExpand){
                double am = aOfT(((double)gTick + 0.5)*(1.0/240.0));
                ds = (float)(1.0/(am*am));
            }
            kDriftK<<<g, b, 0, simS>>>(dPos[src], dPos[dst], dPosErr, dMom, dTau,
                                       dRegime, phiP, gN, DT, invc2, ds, dPartBub);
            gPub = dst;
            forcePass(((double)gTick + 1.0)*(1.0/240.0));
            kKick<<<g, b, 0, simS>>>(dMom, dAcc, dRegime, gN, hdt);
            if (gBHActive){
                kAbsorb<<<g, b, 0, simS>>>(dPos[gPub], dMom, dRegime, gN,
                                           dBHpos, dBHmom, dBHn, dBHacc, dLedger);
                if (gBHFormEnabled && (gTick % 24) == 0){
                    int nc = PMN*PMN*PMN;
                    kBHDetect<<<(nc + 255)/256, b, 0, simS>>>(dFix, nc, dArgmax);
                    kBHSpawn<<<1, 1, 0, simS>>>(dArgmax, dBHpos, dBHmom, dBHn);
                }
                kBHStep<<<1, 1, 0, simS>>>(dBHpos, dBHmom, dBHacc, dBHn, dLedger,
                        gSolver == SOLV_PM ? dFor : nullptr,
                        gSolver == SOLV_PM ? dReal : nullptr, DT, gKHawk, gTick);
            }
            if (gScenario == SC_DETECTOR)
                kDetectorStep<<<g, b, 0, simS>>>(dPos[gPub], dRegime, dRec, gN, gTick, gSeed);
            if (gScenario == SC_RATCHET)
                kRatchetStep<<<(NRECS + 255)/256, b, 0, simS>>>(dRec, NRECS, gTick, gSeed);
            if (gBubActive)
                bubTickHost();                        // evolve -> inscribe/collapse -> spawn
            gTick++; gSimTime += DT;
            if (gHistOn && (gTick % 24) == 0){        // light-history snapshot (0.1 s cadence)
                gHistHead = (gHistHead + 1) % gHistDepth;
                if (gHistN < gHistDepth) gHistN++;
                kSnap<<<g, b, 0, simS>>>(dPos[gPub], dHistBuf + (size_t)gHistHead*gN, gN);
            }
            if (gScenario == SC_ECHO && gTick == 6000){   // Loschmidt flip (declared)
                kFlipMom<<<g, b, 0, simS>>>(dMom, gN);
                CUDA_CHECK(cudaMemsetAsync(dPosErr, 0, (size_t)gN*sizeof(float4), simS));
                forcePass((double)gTick*(1.0/240.0));     // re-establish dAcc invariant
            }
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
        std::vector<unsigned> hgr(gN);
        CUDA_CHECK(cudaStreamSynchronize(simS));
        CUDA_CHECK(cudaMemcpy(hp.data(), dPos[gPub], gN*sizeof(float4), cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(hm.data(), dMom, gN*sizeof(float4), cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(hgr.data(), dRegime, gN*sizeof(unsigned), cudaMemcpyDeviceToHost));
        for (int i = 0; i < gN; i++){
            if (hgr[i] & REG_DEAD) continue;
            double m = hm[i].w;
            double p2 = (double)hm[i].x*hm[i].x + (double)hm[i].y*hm[i].y + (double)hm[i].z*hm[i].z;
            double c2 = (double)C_LIGHT*C_LIGHT;
            M.KE += (m > 0) ? p2*c2/(sqrt(m*m*c2*c2 + p2*c2) + m*c2) : sqrt(p2)*C_LIGHT;
            M.Px += hm[i].x; M.Py += hm[i].y; M.Pz += hm[i].z;
            M.Lz += (double)hp[i].x*hm[i].y - (double)hp[i].y*hm[i].x;
            if (m <= 0){ M.nRel++; }
            else {
                double u2 = p2/(m*m), v2 = u2/(1.0 + u2/c2);
                if (v2 > REL_V2) M.nRel++; else M.nCls++;
            }
            for (int j = i + 1; j < gN; j++){
                if (hgr[j] & REG_DEAD) continue;
                double dx = hp[j].x - hp[i].x, dy = hp[j].y - hp[i].y, dz = hp[j].z - hp[i].z;
                M.PE -= (double)G_DIAL*(double)hm[i].w*hm[j].w/sqrt(dx*dx + dy*dy + dz*dz + 1e-6);
            }
        }
    } else {
        CUDA_CHECK(cudaMemsetAsync(dMet, 0, 8*sizeof(unsigned long long), simS));
        kMeters<<<(gN + 255)/256, 256, 0, simS>>>(dPos[gPub], dMom, dRegime,
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
    CUDA_CHECK(cudaMalloc(&dCntX, 32768*sizeof(unsigned)));
    CUDA_CHECK(cudaMalloc(&dCntV, 32768*sizeof(unsigned)));
    CUDA_CHECK(cudaMalloc(&dRec,  (size_t)NRECS*sizeof(unsigned)));
    CUDA_CHECK(cudaMemset(dRec, 0, (size_t)NRECS*sizeof(unsigned)));
    CUDA_CHECK(cudaMalloc(&dPartBub, (size_t)gN*sizeof(unsigned)));
    CUDA_CHECK(cudaMemset(dPartBub, 0, (size_t)gN*sizeof(unsigned)));
    // M5 gargantua: BH state + lensed HDR
    CUDA_CHECK(cudaMalloc(&dBHpos, NBH_MAX*sizeof(float4)));
    CUDA_CHECK(cudaMalloc(&dBHmom, NBH_MAX*sizeof(float4)));
    CUDA_CHECK(cudaMemset(dBHpos, 0, NBH_MAX*sizeof(float4)));
    CUDA_CHECK(cudaMemset(dBHmom, 0, NBH_MAX*sizeof(float4)));
    CUDA_CHECK(cudaMalloc(&dBHacc, NBH_MAX*4*sizeof(unsigned long long)));
    CUDA_CHECK(cudaMemset(dBHacc, 0, NBH_MAX*4*sizeof(unsigned long long)));
    CUDA_CHECK(cudaMalloc(&dLedger, 4*sizeof(unsigned long long)));
    CUDA_CHECK(cudaMemset(dLedger, 0, 4*sizeof(unsigned long long)));
    CUDA_CHECK(cudaMalloc(&dArgmax, sizeof(unsigned long long)));
    CUDA_CHECK(cudaMemset(dArgmax, 0, sizeof(unsigned long long)));
    CUDA_CHECK(cudaMalloc(&dBHn, 4*sizeof(unsigned)));
    CUDA_CHECK(cudaMemset(dBHn, 0, 4*sizeof(unsigned)));
    CUDA_CHECK(cudaMalloc(&dHdrW, (size_t)RW*RH*sizeof(float4)));
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
        kInitCloud<<<(gN + 255)/256, 256>>>(dPos[0], dMom, dTau, dRegime, gN, gSeed, 100.0f, 0.35f);
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
    case SC_MERGER:
    case SC_ECHO:
        gMTot = (double)gN * 1.0;
        kInitMerger<<<(gN + 255)/256, 256>>>(dPos[0], dMom, dTau, dRegime, gN, gSeed);
        gSolver = SOLV_PM;
        break;
    case SC_RATCHET:
        gN = 1000; gMTot = 1000.0;
        kInitCloud<<<(gN + 255)/256, 256>>>(dPos[0], dMom, dTau, dRegime, gN, gSeed, 100.0f, 0.35f);
        kRatchetInitRecs<<<(NRECS + 255)/256, 256>>>(dRec, NRECS);
        gSolver = SOLV_NONE;
        break;
    case SC_DETECTOR:
        gN = 200000; gMTot = 200000.0;
        kInitStream<<<(gN + 255)/256, 256>>>(dPos[0], dMom, dTau, dRegime, gN, gSeed);
        gSolver = SOLV_NONE;
        break;
    case SC_KEPREL: {                                 // einstein contract: eps = 5e-3 orbit
        gN = 2; gMTot = 10001.0;
        double mu = (double)G_DIAL*10001.0, rap = 16.0, a = 10.0;
        double vap = sqrt(mu*(2.0/rap - 1.0/a));
        float4 hp[2] = { make_float4((float)(-rap/10001.0), 0, 0, 6500.0f),
                         make_float4((float)( rap*10000.0/10001.0), 0, 0, 9000.0f) };
        float4 hm[2] = { make_float4(0, (float)(-vap/10001.0*1e4), 0, 1e4f),
                         make_float4(0, (float)( vap*10000.0/10001.0), 0, 1.0f) };
        float ht[2] = {0, 0}; unsigned hr[2] = {0x04u, 0x04u};
        CUDA_CHECK(cudaMemcpy(dPos[0], hp, sizeof hp, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(dMom, hm, sizeof hm, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(dTau, ht, sizeof ht, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(dRegime, hr, sizeof hr, cudaMemcpyHostToDevice));
        gSolver = SOLV_TINY;
    } break;
    case SC_CLOCKS: {
        gN = 3; gMTot = 10002.0;
        double GM = (double)G_DIAL*1e4;
        double vA = sqrt(GM/2.0), vB = sqrt(GM/100.0);
        float4 hp[3] = { make_float4(0, 0, 0, 6500.0f),
                         make_float4(2.0f, 0, 0, 12000.0f),
                         make_float4(100.0f, 0, 0, 4000.0f) };
        float4 hm[3] = { make_float4(0, 0, 0, 1e4f),
                         make_float4(0, (float)vA, 0, 1.0f),
                         make_float4(0, (float)vB, 0, 1.0f) };
        float ht[3] = {0, 0, 0}; unsigned hr[3] = {0x02u, 0x04u, 0x02u};
        CUDA_CHECK(cudaMemcpy(dPos[0], hp, sizeof hp, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(dMom, hm, sizeof hm, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(dTau, ht, sizeof ht, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(dRegime, hr, sizeof hr, cudaMemcpyHostToDevice));
        gSolver = SOLV_TINY;
    } break;
    case SC_PHOTONS: {
        gN = 65; gMTot = 1e4;
        std::vector<float4> hp(65), hm(65);
        std::vector<float> ht(65, 0.0f);
        std::vector<unsigned> hr(65, 0x24u);
        hp[0] = make_float4(0, 0, 0, 6500.0f);
        hm[0] = make_float4(0, 0, 0, 1e4f);
        hr[0] = 0x02u;
        for (int q = 1; q <= 64; q++){
            hp[q] = make_float4(-100.0f, 2.0f + 0.25f*(q - 1), 0, 20000.0f);
            hm[q] = make_float4(1.0f, 0, 0, 0.0f);    // photon: |p| = 1, m = 0
        }
        CUDA_CHECK(cudaMemcpy(dPos[0], hp.data(), 65*sizeof(float4), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(dMom, hm.data(), 65*sizeof(float4), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(dTau, ht.data(), 65*sizeof(float), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(dRegime, hr.data(), 65*sizeof(unsigned), cudaMemcpyHostToDevice));
        gSolver = SOLV_DIRECT;
    } break;
    case SC_COLLAPSE:                                 // gargantua contract: unscripted formation
        gMTot = (double)gN * 1.0;
        kInitCloud<<<(gN + 255)/256, 256>>>(dPos[0], dMom, dTau, dRegime, gN, gSeed, 40.0f, 0.15f);
        gSolver = SOLV_PM; gBHActive = true; gBHFormEnabled = true;
        break;
    case SC_ISCO: {
        gN = 8; gMTot = 1e5 + 8;
        const float rr[8] = {2.4f, 2.55f, 2.7f, 4.5f, 5.0f, 5.5f, 6.5f, 8.0f};
        float4 hp[8], hm[8]; float ht[8] = {0};
        unsigned hr8[8];
        double GM = (double)G_DIAL*1e5, rs = 1.0, c2 = (double)C_LIGHT*C_LIGHT;
        for (int q = 0; q < 8; q++){
            double r = rr[q], d = r - rs;
            double gacc = GM/(d*d);
            double v2 = gacc*r;                       // SR-circular: v^2 = g r / gamma, iterate
            for (int it = 0; it < 40; it++) v2 = gacc*r*sqrt(1.0 - v2/c2);
            double u = sqrt(v2)/sqrt(1.0 - v2/c2);
            double ph = 2.0*3.14159265358979*q/8.0;
            hp[q] = make_float4((float)(r*cos(ph)), (float)(r*sin(ph)), 0, 8000.0f + 500.0f*q);
            // -1e-3 u_phi inward kick: unstable circulars must roll off (nexus N4 idiom)
            hm[q] = make_float4((float)(-u*sin(ph) - 1e-3*u*cos(ph)),
                                (float)( u*cos(ph) - 1e-3*u*sin(ph)), 0, 1.0f);
            ht[q] = 0; hr8[q] = 0x04u;
        }
        CUDA_CHECK(cudaMemcpy(dPos[0], hp, sizeof hp, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(dMom, hm, sizeof hm, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(dTau, ht, sizeof ht, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(dRegime, hr8, sizeof hr8, cudaMemcpyHostToDevice));
        float4 bp0 = make_float4(0, 0, 0, 1e5f), bm0 = make_float4(0, 0, 0, RS_PER_M*1e5f);
        CUDA_CHECK(cudaMemcpy(dBHpos, &bp0, sizeof bp0, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(dBHmom, &bm0, sizeof bm0, cudaMemcpyHostToDevice));
        unsigned one = 1;
        CUDA_CHECK(cudaMemcpy(dBHn, &one, sizeof one, cudaMemcpyHostToDevice));
        gSolver = SOLV_NONE; gBHActive = true;
    } break;
    case SC_HAWKING: {
        gN = 32; gMTot = 250.0 + 32.0;
        float4 hp[32], hm[32]; float ht[32] = {0};
        unsigned hr32[32];
        double GM = (double)G_DIAL*250.0;
        double v = sqrt(GM/5.0);
        for (int q = 0; q < 32; q++){
            double ph = 2.0*3.14159265358979*q/32.0;
            hp[q] = make_float4((float)(5.0*cos(ph)), (float)(5.0*sin(ph)), 0, 4200.0f + 90.0f*q);
            hm[q] = make_float4((float)(-v*sin(ph)), (float)(v*cos(ph)), 0, 1.0f);
            hr32[q] = 0x02u;
        }
        CUDA_CHECK(cudaMemcpy(dPos[0], hp, sizeof hp, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(dMom, hm, sizeof hm, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(dTau, ht, sizeof ht, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(dRegime, hr32, sizeof hr32, cudaMemcpyHostToDevice));
        float4 bp0 = make_float4(0, 0, 0, 250.0f), bm0 = make_float4(0, 0, 0, RS_PER_M*250.0f);
        CUDA_CHECK(cudaMemcpy(dBHpos, &bp0, sizeof bp0, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(dBHmom, &bm0, sizeof bm0, cudaMemcpyHostToDevice));
        unsigned one = 1;
        CUDA_CHECK(cudaMemcpy(dBHn, &one, sizeof one, cudaMemcpyHostToDevice));
        gSolver = SOLV_NONE; gBHActive = true;
    } break;
    case SC_DOUBLESLIT: {                             // planck contract: the experiment
        gN = 4096; gMTot = 4096.0;
        std::vector<double> dotsY;
        qRunDoubleslit(dotsY);
        std::vector<float4> hp(gN), hm(gN, make_float4(0, 0, 0, 1.0f));
        std::vector<float> ht(gN, 0.0f);
        std::vector<unsigned> hr(gN, 0x02u);
        for (int j = 0; j < gN; j++){
            float zj = (float)(4.0*orrery::counter_gauss(gSeed, j, 3, 8));  // curtain depth, deterministic
            hp[j] = make_float4(14.5f, (float)dotsY[j], zj, 12000.0f);
        }
        CUDA_CHECK(cudaMemcpy(dPos[0], hp.data(), (size_t)gN*16, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(dMom, hm.data(), (size_t)gN*16, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(dTau, ht.data(), (size_t)gN*4, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(dRegime, hr.data(), (size_t)gN*4, cudaMemcpyHostToDevice));
        gSolver = SOLV_NONE;
    } break;
    case SC_TUNNELING:
        gN = 1000; gMTot = 1000.0;
        kInitCloud<<<(gN + 255)/256, 256>>>(dPos[0], dMom, dTau, dRegime, gN, gSeed, 100.0f, 0.35f);
        qRunTunneling();
        gSolver = SOLV_NONE;
        break;
    case SC_SHOQ:
        gN = 1000; gMTot = 1000.0;
        kInitCloud<<<(gN + 255)/256, 256>>>(dPos[0], dMom, dTau, dRegime, gN, gSeed, 100.0f, 0.35f);
        qRunShoq();
        gSolver = SOLV_NONE;
        break;
    case SC_CIRCUMNAV: {                              // cosmos contract: light laps the universe
        gN = 2; gMTot = 1.0;
        double u05 = 10.0/sqrt(0.75);                 // u = v/sqrt(1-v^2/c^2) at v = c/2
        float4 hp[2] = { make_float4(0,  2.0f, 0, 20000.0f),
                         make_float4(0, -2.0f, 0,  8000.0f) };
        float4 hm[2] = { make_float4(1.0f, 0, 0, 0.0f),          // photon: |p| = 1, m = 0
                         make_float4((float)u05, 0, 0, 1.0f) };  // massive at 0.5c
        float ht[2] = {0, 0}; unsigned hr[2] = {0x24u, 0x04u};
        CUDA_CHECK(cudaMemcpy(dPos[0], hp, sizeof hp, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(dMom, hm, sizeof hm, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(dTau, ht, sizeof ht, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(dRegime, hr, sizeof hr, cudaMemcpyHostToDevice));
        gSolver = SOLV_NONE;
    } break;
    case SC_EXPAND: {                                 // cosmos contract: EdS comoving cosmology
        gN = 1000000; gMTot = 1e6;                    // 100^3 Zel'dovich lattice (forced in main)
        gExpand = true;
        gH0 = sqrt(8.0*3.14159265358979323846*(double)G_DIAL
                   * ((double)gN/(512.0*512.0*512.0)) / 3.0);
        const float k1 = 2.0f*PI_F*4.0f/512.0f;       // mode (4,0,0): 32 PM cells per wavelength
        const float Ad = 0.02f/k1;                    // delta_0 = Ad*k1 = 0.02
        kInitZeldovich<<<(gN + 255)/256, 256>>>(dPos[0], dMom, dTau, dRegime,
                                                Ad, k1, (float)gH0);
        gSolver = SOLV_PM;
    } break;
    case SC_BUBBLES: {                                // cosmos contract: roaming 3D bubbles
        gN = 8 + 4*512; gMTot = (double)gN;
        std::vector<float4> hp(gN), hm(gN);
        std::vector<float> ht(gN, 0.0f);
        std::vector<unsigned> hr(gN, 0x02u);
        for (int q = 0; q < 8; q++){                  // loners: isolated ring, p = 0.2 tangential
            double ph = 2.0*3.14159265358979*q/8.0;
            hp[q] = make_float4((float)(150.0*cos(ph)), (float)(150.0*sin(ph)), 0, 12000.0f);
            hm[q] = make_float4((float)(-0.2*sin(ph)), (float)(0.2*cos(ph)), 0, 1.0f);
        }
        for (int c = 0; c < 4; c++){                  // rigid intruder clusters aimed at loners 0-3
            double ph = 2.0*3.14159265358979*c/8.0;
            float ccx = (float)(210.0*cos(ph)), ccy = (float)(210.0*sin(ph));
            float cvx = (float)(-6.0*cos(ph)), cvy = (float)(-6.0*sin(ph));
            for (int j = 0; j < 512; j++){
                int i = 8 + c*512 + j;
                using namespace orrery;
                float gx = (float)counter_gauss(gSeed, i, 0, 9);
                float gy = (float)counter_gauss(gSeed, i, 1, 9);
                float gz = (float)counter_gauss(gSeed, i, 2, 9);
                float il = 1.0f/sqrtf(gx*gx + gy*gy + gz*gz + 1e-12f);
                float r = 4.0f*cbrtf((float)counter_uniform(gSeed, i, 3, 9));
                hp[i] = make_float4(ccx + r*gx*il, ccy + r*gy*il, r*gz*il, 4500.0f);
                hm[i] = make_float4(cvx, cvy, 0, 1.0f);
            }
        }
        CUDA_CHECK(cudaMemcpy(dPos[0], hp.data(), (size_t)gN*16, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(dMom, hm.data(), (size_t)gN*16, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(dTau, ht.data(), (size_t)gN*4, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(dRegime, hr.data(), (size_t)gN*4, cudaMemcpyHostToDevice));
        gSolver = SOLV_NONE;
    } break;
    }
    if (gScenario == SC_BUBBLES || gBubWant) gBubActive = true;
    if (gBubActive) bubAllocOnce();
    if (gWantHist){                                   // light-history ring (visuals only):
        size_t freeB = 0, totB = 0;                   // adaptive depth, <= 40% of free VRAM
        cudaMemGetInfo(&freeB, &totB);
        size_t per = (size_t)gN*sizeof(ushort4);
        int depth = (int)((double)freeB*0.4/(double)per);
        if (depth > 384) depth = 384;
        if (depth >= 48){
            CUDA_CHECK(cudaMalloc(&dHistBuf, (size_t)depth*per));
            gHistDepth = depth; gHistOn = true;
        }
    }
    CUDA_CHECK(cudaMemcpy(dPos[1], dPos[0], (size_t)gN*sizeof(float4),
                          cudaMemcpyDeviceToDevice));
    CUDA_CHECK(cudaDeviceSynchronize());
    forcePass((double)gTick*(1.0/240.0));             // establish dAcc = a(x0); Phi for meters
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
        else if (k == 'L') gProj ^= 1;                // little planet <-> perspective
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
        else if (a == "--bubbles")              gBubWant = true;   // live-universe bubbles (declared)
        else if (a == "--planet")               gProj = 1;         // little-planet projection
        else if (a == "--az" && i+1 < argc)     gShotAz = (float)atof(argv[++i]);   // shot framing
        else if (a == "--el" && i+1 < argc)     gShotEl = (float)atof(argv[++i]);
        else if (a == "--dist" && i+1 < argc)   gShotDist = (float)atof(argv[++i]);
        else if (a == "--ev" && i+1 < argc)     gEvOffset = (float)atof(argv[++i]);
        else if (a == "--nohud")                gHud = 0;          // clean render (non-declared)
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
            else if (s == "merger")    gScenario = SC_MERGER;
            else if (s == "echo")      gScenario = SC_ECHO;
            else if (s == "ratchet")   gScenario = SC_RATCHET;
            else if (s == "detector")  gScenario = SC_DETECTOR;
            else if (s == "keprel")    gScenario = SC_KEPREL;
            else if (s == "clocks")    gScenario = SC_CLOCKS;
            else if (s == "photons")   gScenario = SC_PHOTONS;
            else if (s == "collapse")  gScenario = SC_COLLAPSE;
            else if (s == "isco")      gScenario = SC_ISCO;
            else if (s == "hawking")   gScenario = SC_HAWKING;
            else if (s == "doubleslit") gScenario = SC_DOUBLESLIT;
            else if (s == "tunneling")  gScenario = SC_TUNNELING;
            else if (s == "shoq")       gScenario = SC_SHOQ;
            else if (s == "circumnav")  gScenario = SC_CIRCUMNAV;
            else if (s == "expand")     gScenario = SC_EXPAND;
            else if (s == "bubbles")    gScenario = SC_BUBBLES;
            else { fprintf(stderr, "error: unknown scenario\n"); return 2; }
        }
        else { fprintf(stderr, "usage: tinyuniverse [--scenario galaxy|kepler|threebody|cloud] "
                       "[--w W --h H] [--n N] [--seed S] [--ssaa] [--phys] [--aces] [--bench F] "
                       "[--shot PATH.bmp [--frames N] [--warm T]] [--steps S --json] [--golden]\n");
               return 2; }
    }
    gOW &= ~1; gOH &= ~1;
    if (gN < 1000) gN = 1000;
    if (gScenario == SC_EXPAND) gN = 1000000;         // lattice is exactly 100^3 (contract)
    if (goldenMode) gSeed = 20260711ull;              // frozen seed at init time
    if (jsonMode || goldenMode) gWantHist = false;    // history is render-only, never declared

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
            const long long GS[20] = {10000, 1000000, 1000000, 12000, 12000, 12000, 4000, 4000,
                                      320000, 24000, 3000, 12000, 6000, 4000, 0, 0, 0,
                                      12288, 24000, 4800};
            runSteps = GS[gScenario];
        }
        if (runSteps <= 0){
            const long long GS[20] = {10000, 1000000, 1000000, 12000, 12000, 12000, 4000, 4000,
                                      320000, 24000, 3000, 12000, 6000, 4000, 0, 0, 0,
                                      12288, 24000, 4800};
            runSteps = GS[gScenario];
        }
        long long done = 0, nextPct = 10;
        const int batch = (gSolver == SOLV_TINY)
                        ? ((gScenario == SC_KEPREL) ? 1000 : 20000)
                        : ((gScenario == SC_BUBBLES) ? 240 : 250);
        const bool wantS = (gScenario == SC_MERGER || gScenario == SC_ECHO);
        double S0 = 0, Smid = 0, Sprev = 0;
        int monoUp = 0, monoTot = 0;
        double prcPrev = 0, prcTot = 0, prcCum = 0; bool prcInit = false;
        std::vector<double> angs;
        std::vector<double> amps, avals;              // expand: growth checkpoints
        if (wantS){ S0 = entropyNow(); Sprev = S0; }
        if (gScenario == SC_EXPAND){
            amps.push_back(expandAmp()); avals.push_back(1.0);
            forcePass(0.0);                           // restore Phi/dAcc after measurement
        }
        while (done < runSteps){
            int nt = (int)((runSteps - done) < batch ? (runSteps - done) : batch);
            runTicks(nt);
            done += nt;
            if (gScenario == SC_EXPAND && done % 2000 == 0){
                amps.push_back(expandAmp()); avals.push_back(aOfT((double)done/240.0));
                forcePass((double)done/240.0);
            }
            if (gScenario == SC_BUBBLES && done == 1440){
                // sigma probe: slot 0 (lowest owner id), untouched until first intrusion
                // at ~tick 1920; bubbles spawned at check tick 0 have evolved done-1 steps
                double tEv = (double)(done - 1)/240.0;
                gBubSig = bubSigma(0);
                gBubSigExp = sqrt(1.0 + 0.0625*tEv*tEv);   // sigma0 sqrt(1+(hbar t/2m sigma0^2)^2)
            }
            if (gScenario == SC_KEPREL){              // LRL angle, batch-sampled + unwrapped
                float4 sp2[2], sm2[2];
                CUDA_CHECK(cudaStreamSynchronize(simS));
                CUDA_CHECK(cudaMemcpy(sp2, dPos[gPub], sizeof sp2, cudaMemcpyDeviceToHost));
                CUDA_CHECK(cudaMemcpy(sm2, dMom, sizeof sm2, cudaMemcpyDeviceToHost));
                double c2 = (double)C_LIGHT*C_LIGHT;
                double rx = sp2[1].x - sp2[0].x, ry = sp2[1].y - sp2[0].y;
                double u1x = sm2[1].x/sm2[1].w, u1y = sm2[1].y/sm2[1].w;
                double u0x = sm2[0].x/sm2[0].w, u0y = sm2[0].y/sm2[0].w;
                double g1 = sqrt(1.0 + (u1x*u1x + u1y*u1y)/c2), g0 = sqrt(1.0 + (u0x*u0x + u0y*u0y)/c2);
                double vx = u1x/g1 - u0x/g0, vy = u1y/g1 - u0y/g0;
                double mu = (double)G_DIAL*10001.0;
                double L = rx*vy - ry*vx, rn = sqrt(rx*rx + ry*ry);
                double ex = vy*L/mu - rx/rn, ey = -vx*L/mu - ry/rn;
                double ang = atan2(ey, ex);
                if (!prcInit){ prcPrev = ang; prcInit = true; angs.push_back(0.0); }
                else {
                    double dgl = ang - prcPrev;
                    while (dgl >  PI_F) dgl -= 2.0*PI_F;
                    while (dgl < -PI_F) dgl += 2.0*PI_F;
                    prcCum += dgl; prcPrev = ang;
                    angs.push_back(prcCum);
                }
            }
            if (wantS){
                double S = entropyNow();
                if (S >= Sprev - 0.01) monoUp++;      // fluctuation-scale tolerance (D-015)
                monoTot++;
                if (gScenario == SC_ECHO && done == 6000) Smid = S;
                Sprev = S;
            }
            if (done*100 >= runSteps*nextPct){
                CUDA_CHECK(cudaStreamSynchronize(simS));
                fprintf(stderr, "  %lld%%\r", nextPct); nextPct += 10;
            }
        }
        double Send = Sprev;
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
        unsigned nbFinal = 0;
        float4 bhP[NBH_MAX] = {}, bhM[NBH_MAX] = {};
        unsigned long long lg[4] = {};
        if (gBHActive){                               // BH block joins the declared state
            CUDA_CHECK(cudaMemcpy(&nbFinal, dBHn, 4, cudaMemcpyDeviceToHost));
            CUDA_CHECK(cudaMemcpy(bhP, dBHpos, sizeof bhP, cudaMemcpyDeviceToHost));
            CUDA_CHECK(cudaMemcpy(bhM, dBHmom, sizeof bhM, cudaMemcpyDeviceToHost));
            CUDA_CHECK(cudaMemcpy(lg, dLedger, sizeof lg, cudaMemcpyDeviceToHost));
            bytes.append((const char*)&nbFinal, 4);
            bytes.append((const char*)bhP, sizeof bhP);
            bytes.append((const char*)bhM, sizeof bhM);
            bytes.append((const char*)lg, sizeof lg);
        }
        if (!gPsiBytes.empty()) bytes.append(gPsiBytes);   // psi joins the declared state (planck)
        int hOwn[NBUB] = {}; unsigned hRecs[NBUB] = {};
        if (gBubActive){                              // bubble block joins the declared state (cosmos)
            float4 hCtrs[NBUB];
            CUDA_CHECK(cudaMemcpy(hOwn, dBubOwner, sizeof hOwn, cudaMemcpyDeviceToHost));
            CUDA_CHECK(cudaMemcpy(hRecs, dBubRec, sizeof hRecs, cudaMemcpyDeviceToHost));
            CUDA_CHECK(cudaMemcpy(hCtrs, dBubCtr, sizeof hCtrs, cudaMemcpyDeviceToHost));
            unsigned bs[3] = {gBubSpawned, gBubCollapsed, (unsigned)gBubLive};
            bytes.append((const char*)bs, sizeof bs);
            bytes.append((const char*)hOwn, sizeof hOwn);
            bytes.append((const char*)hRecs, sizeof hRecs);
            bytes.append((const char*)hCtrs, sizeof hCtrs);
            std::vector<cufftComplex> hg((size_t)QB3);
            for (int s = 0; s < NBUB; s++){           // surviving psi grids, slot order
                if (hOwn[s] < 0) continue;
                CUDA_CHECK(cudaMemcpy(hg.data(), dBub + (size_t)s*QB3,
                                      (size_t)QB3*sizeof(cufftComplex), cudaMemcpyDeviceToHost));
                bytes.append((const char*)hg.data(), (size_t)QB3*sizeof(cufftComplex));
            }
        }
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
        std::string gates, resExtra;
        using namespace orrery;
        auto B = [](bool b){ return std::string(b ? "true" : "false"); };
        switch (gScenario){
        case SC_DOUBLESLIT: {
            bool g1 = gQM.c_psi > 0.35, g2 = gQM.c_dots > 0.30, g3 = gQM.c_det < 0.12;
            bool g4 = gQM.kpk_rel < 0.08, g5 = fabs(gQM.c_dots - gQM.c_psi) < 0.1;  // near-field (D-018)
            bool g6 = (gQM.nrec == 4096);
            pass = g1 && g2 && g3 && g4 && g5 && g6;
            gates = "\"c_psi_gt_0.35\":" + B(g1) + ",\"c_dots_gt_0.30\":" + B(g2)
                  + ",\"c_det_lt_0.12\":" + B(g3) + ",\"k_peak_8pc\":" + B(g4)
                  + ",\"dots_match_psi\":" + B(g5) + ",\"recorded_4096\":" + B(g6);
            resExtra = ",\"c_psi\":" + fmt6(gQM.c_psi) + ",\"c_dots\":" + fmt6(gQM.c_dots)
                     + ",\"c_det\":" + fmt6(gQM.c_det) + ",\"kpk_rel\":" + fmt6(gQM.kpk_rel)
                     + ",\"n_recorded\":" + fmti(gQM.nrec);
        } break;
        case SC_TUNNELING: {
            double dT = fabs(gQM.T - gQM.Tref);
            bool g1 = dT < 1e-3, g2 = (gQM.T > 0.005 && gQM.T < 0.2);
            pass = g1 && g2;
            gates = "\"t_vs_fp64_lt_1e-3\":" + B(g1) + ",\"t_band\":" + B(g2);
            resExtra = ",\"t_2d\":" + fmt6(gQM.T) + ",\"t_ref_fp64\":" + fmt6(gQM.Tref);
        } break;
        case SC_SHOQ: {
            double erel = fabs(gQM.E - 0.5)/0.5, srel = fabs(gQM.sigx - 0.5)/0.5;
            bool g1 = erel < 1e-3, g2 = srel < 1e-2;
            pass = g1 && g2;
            gates = "\"e0_lt_1e-3\":" + B(g1) + ",\"sigma_lt_1e-2\":" + B(g2);
            resExtra = ",\"e0\":" + fmt6(gQM.E) + ",\"e0_exp\":" + fmt6(0.5)
                     + ",\"sigma_x\":" + fmt6(gQM.sigx);
        } break;
        case SC_COLLAPSE: {
            std::vector<unsigned> hg(gN);
            CUDA_CHECK(cudaMemcpy(hg.data(), dRegime, (size_t)gN*4, cudaMemcpyDeviceToHost));
            long long nDead = 0;
            for (int q = 0; q < gN; q++) if (hg[q] & REG_DEAD) nDead++;
            double mBH = (nbFinal > 0) ? (double)bhP[0].w : 0.0;
            double absorbed = orrery::fixed_decode((long long)lg[0]);
            double ledgErr = fabs(absorbed - (double)nDead);
            pass = (nbFinal >= 1) && (mBH >= 1.5e5) && (ledgErr < 1e-6);
            gates = "\"formed_ge_1\":" + B(nbFinal >= 1) + ",\"m_bh_ge_1.5e5\":" + B(mBH >= 1.5e5)
                  + ",\"ledger_exact\":" + B(ledgErr < 1e-6);
            unsigned pk4[4];
            CUDA_CHECK(cudaMemcpy(pk4, dBHn, sizeof pk4, cudaMemcpyDeviceToHost));
            resExtra = ",\"n_bh\":" + fmti((long long)nbFinal) + ",\"m_bh\":" + fmt6(mBH)
                     + ",\"n_absorbed\":" + fmti(nDead) + ",\"absorbed_mass\":" + fmt6(absorbed)
                     + ",\"peak_cell_mass\":" + fmti((long long)pk4[2]);
        } break;
        case SC_ISCO: {
            unsigned hg8[8];
            CUDA_CHECK(cudaMemcpy(hg8, dRegime, sizeof hg8, cudaMemcpyDeviceToHost));
            int plunged = 0, survived = 0;
            for (int q = 0; q < 3; q++) if (hg8[q] & REG_DEAD) plunged++;
            for (int q = 3; q < 8; q++) if (!(hg8[q] & REG_DEAD)) survived++;
            pass = (plunged == 3) && (survived == 5);
            gates = "\"inner3_plunged\":" + B(plunged == 3) + ",\"outer5_alive\":" + B(survived == 5);
            resExtra = ",\"n_plunged\":" + fmti(plunged) + ",\"n_survived\":" + fmti(survived)
                     + ",\"m_bh\":" + fmt6((double)bhP[0].w);
        } break;
        case SC_HAWKING: {
            double radiated = orrery::fixed_decode((long long)lg[1]);
            double E0 = 250.0*(double)C_LIGHT*C_LIGHT;
            double ledgRel = fabs(radiated - E0)/E0;
            long long popTick = (long long)lg[2];
            double nExact = (250.0*250.0*250.0 - (double)M_POP*M_POP*M_POP)/(3.0*gKHawk*(1.0/240.0));
            long long popExp = (long long)ceil(nExact - 1e-9) - 1;   // tick index of the popping step
            bool popOk = fabs((double)(popTick - popExp)) <= 0.005*nExact;
            pass = (nbFinal == 0) && popOk && (ledgRel < 1e-6);   // fp32 BH mass round-trip (D-017)
            gates = "\"evaporated\":" + B(nbFinal == 0) + ",\"pop_tick_0.5pc\":" + B(popOk)
                  + ",\"ledger_lt_1e-6\":" + B(ledgRel < 1e-6);
            resExtra = ",\"pop_tick\":" + fmti(popTick) + ",\"pop_exp\":" + fmti(popExp)
                     + ",\"radiated\":" + fmt6(radiated) + ",\"e0c2\":" + fmt6(E0);
        } break;
        case SC_KEPREL: {
            // exact Sommerfeld precession for relativistic Kepler (dp/dt = -mu m r_hat/r^2):
            // per orbit 2*pi*(1/Lambda - 1), Lambda = sqrt(1 - (mu/(c*l))^2), l = |r x u| at init
            double mu = (double)G_DIAL*10001.0, aa = 10.0, rap = 16.0;
            double uap = sqrt(mu*(2.0/rap - 1.0/aa));
            double l = rap*uap;
            double lam = sqrt(1.0 - (mu/((double)C_LIGHT*l))*(mu/((double)C_LIGHT*l)));
            double rate = 2.0*3.14159265358979*(1.0/lam - 1.0);
            double Tn = 2.0*3.14159265358979*sqrt(aa*aa*aa/mu);
            double prcExp = rate*((double)runSteps/240.0)/Tn;
            // least-squares slope over dense samples averages out the intra-orbit
            // LRL oscillation (endpoint sampling aliases it — measured, D-016)
            {
                size_t n = angs.size();
                double tb = (n - 1)/2.0, num = 0, den = 0, ab = 0;
                for (size_t q = 0; q < n; q++) ab += angs[q];
                ab /= n;
                for (size_t q = 0; q < n; q++){ num += (q - tb)*(angs[q] - ab); den += (q - tb)*(q - tb); }
                prcTot = (num/den)*(double)(n - 1);
            }
            double rel = fabs(prcTot/prcExp - 1.0);
            pass = (rel < 0.05) && (maxr < 30.0);
            gates = "\"prec_lt_5pc\":" + B(rel < 0.05) + ",\"bounded_r_lt_30\":" + B(maxr < 30.0);
            resExtra = ",\"prec_meas\":" + fmt6(prcTot) + ",\"prec_exp\":" + fmt6(prcExp)
                     + ",\"prec_rel_err\":" + fmt6(rel);
        } break;
        case SC_CLOCKS: {
            float ht3[3];
            CUDA_CHECK(cudaMemcpy(ht3, dTau, sizeof ht3, cudaMemcpyDeviceToHost));
            double GM = (double)G_DIAL*1e4, c2 = (double)C_LIGHT*C_LIGHT;
            auto aOf = [&](double r0){
                double u0 = sqrt(GM/r0), v0 = u0/sqrt(1.0 + u0*u0/c2);
                return -GM/(2.0*(0.5*v0*v0 - GM/r0));
            };
            double aA = aOf(2.0), aB = aOf(100.0);
            double expd = sqrt(1.0 - 3.0*GM/(aA*c2))/sqrt(1.0 - 3.0*GM/(aB*c2));
            double ratio = (double)ht3[1]/(double)ht3[2];
            pass = fabs(ratio - expd) < 2e-3;
            gates = "\"tau_ratio_lt_2e-3\":" + B(pass);
            resExtra = ",\"tau_a\":" + fmt6(ht3[1]) + ",\"tau_b\":" + fmt6(ht3[2])
                     + ",\"ratio\":" + fmt6(ratio) + ",\"ratio_exp\":" + fmt6(expd);
        } break;
        case SC_PHOTONS: {
            std::vector<float4> hm65(65);
            CUDA_CHECK(cudaMemcpy(hm65.data(), dMom, 65*sizeof(float4), cudaMemcpyDeviceToHost));
            double GM = (double)G_DIAL*1e4, c2 = (double)C_LIGHT*C_LIGHT;
            double sum = 0; int cnt = 0;
            for (int q = 19; q <= 64; q++){
                double bq = 2.0 + 0.25*(q - 1);
                double angm = atan2((double)hm65[q].y, (double)hm65[q].x);
                double ange = -4.0*GM/(c2*bq);
                sum += fabs(angm/ange - 1.0); cnt++;
            }
            double meanrel = sum/cnt;
            pass = (meanrel < 0.05);
            gates = "\"defl_mean_lt_5pc\":" + B(pass);
            resExtra = ",\"defl_mean_rel\":" + fmt6(meanrel) + ",\"photons_gated\":" + fmti(cnt);
        } break;
        case SC_MERGER: {
            double rise = Send - S0, mf = monoTot ? (double)monoUp/monoTot : 0;
            pass = (rise > 0.3) && (mf >= 0.75) && (de < 0.08);
            gates = "\"ds_gt_0.3\":" + B(rise > 0.3) + ",\"mono_ge_0.75\":" + B(mf >= 0.75)
                  + ",\"de_lt_8e-2\":" + B(de < 0.08);
            resExtra = ",\"s0\":" + fmt6(S0) + ",\"s_end\":" + fmt6(Send)
                     + ",\"mono_frac\":" + fmt6(mf);
        } break;
        case SC_ECHO: {
            double rise = Smid - S0, ret = Send - S0;
            pass = (rise > 0.25) && (ret < 0.35*rise);
            gates = "\"rise_gt_0.25\":" + B(rise > 0.25)
                  + ",\"echo_return_lt_0.35rise\":" + B(ret < 0.35*rise);
            resExtra = ",\"s0\":" + fmt6(S0) + ",\"s_mid\":" + fmt6(Smid)
                     + ",\"s_end\":" + fmt6(Send);
        } break;
        case SC_RATCHET: {
            std::vector<unsigned> hr(NRECS);
            CUDA_CHECK(cudaMemcpy(hr.data(), dRec, (size_t)NRECS*4, cudaMemcpyDeviceToHost));
            long long unw[3] = {0, 0, 0};
            const int cls = NRECS/3;
            for (int q = 0; q < NRECS; q++)
                if (hr[q] == 0u) unw[q/cls]++;         // unresolved counts as survived (declared)
            const double lam = 0.3/(0.7*0.6);
            double fr[3], ex[3], rel[3];
            bool ok = true;
            const double tol[3] = {0.01, 0.01, 0.05};
            const int Rv[3] = {1, 4, 16};
            for (int q = 0; q < 3; q++){
                fr[q] = (double)unw[q]/cls;
                ex[q] = pow(lam < 1.0 ? lam : 1.0, Rv[q]);
                rel[q] = fabs(fr[q]/ex[q] - 1.0);
                ok &= (rel[q] < tol[q]);
            }
            pass = ok;
            gates = "\"r1_lt_1pc\":" + B(rel[0] < 0.01) + ",\"r4_lt_1pc\":" + B(rel[1] < 0.01)
                  + ",\"r16_lt_5pc\":" + B(rel[2] < 0.05);
            resExtra = ",\"unwrite_r1\":" + fmt6(fr[0]) + ",\"unwrite_r4\":" + fmt6(fr[1])
                     + ",\"unwrite_r16\":" + fmt6(fr[2]);
        } break;
        case SC_DETECTOR: {
            std::vector<unsigned> hg(gN);
            CUDA_CHECK(cudaMemcpy(hg.data(), dRegime, (size_t)gN*4, cudaMemcpyDeviceToHost));
            long long crossed = 0, recorded = 0;
            for (int q = 0; q < gN; q++){
                if (hg[q] & 0x80u) crossed++;
                if (hg[q] & 0x40u) recorded++;
            }
            double frac = crossed ? (double)recorded/crossed : 0;
            pass = (frac > 0.95) && (crossed > 50000);
            gates = "\"rec_frac_gt_0.95\":" + B(frac > 0.95)
                  + ",\"crossed_gt_5e4\":" + B(crossed > 50000);
            resExtra = ",\"n_crossed\":" + fmti(crossed) + ",\"n_recorded\":" + fmti(recorded);
        } break;
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
        case SC_CIRCUMNAV: {
            float4 hp2[2]; float ht2[2];
            CUDA_CHECK(cudaMemcpy(hp2, dPos[gPub], sizeof hp2, cudaMemcpyDeviceToHost));
            CUDA_CHECK(cudaMemcpy(ht2, dTau, sizeof ht2, cudaMemcpyDeviceToHost));
            double dph = sqrt((double)hp2[0].x*hp2[0].x + ((double)hp2[0].y - 2.0)*((double)hp2[0].y - 2.0)
                            + (double)hp2[0].z*hp2[0].z);
            double dma = sqrt((double)hp2[1].x*hp2[1].x + ((double)hp2[1].y + 2.0)*((double)hp2[1].y + 2.0)
                            + (double)hp2[1].z*hp2[1].z);
            double tExp = ((double)runSteps/240.0)*sqrt(0.75);
            double tRel = fabs((double)ht2[1]/tExp - 1.0);
            bool g1 = dph < 0.02, g2 = dma < 0.02;
            bool g3 = (ht2[0] == 0.0f), g4 = tRel < 2e-3;
            pass = g1 && g2 && g3 && g4;
            gates = "\"photon_return_lt_0.02\":" + B(g1) + ",\"massive_return_lt_0.02\":" + B(g2)
                  + ",\"photon_tau_zero\":" + B(g3) + ",\"tau_sr_lt_2e-3\":" + B(g4);
            resExtra = ",\"photon_return\":" + fmt6(dph) + ",\"massive_return\":" + fmt6(dma)
                     + ",\"tau_massive\":" + fmt6((double)ht2[1]) + ",\"tau_exp\":" + fmt6(tExp)
                     + ",\"tau_rel_err\":" + fmt6(tRel);
        } break;
        case SC_EXPAND: {
            // linear growth: ln A vs ln a LS slope must be ~1 (D proportional to a, EdS)
            size_t na = amps.size();
            double sx = 0, sy = 0, sxx = 0, sxy = 0;
            for (size_t q = 0; q < na; q++){
                double lx = log(avals[q]), ly = log(amps[q]);
                sx += lx; sy += ly; sxx += lx*lx; sxy += lx*ly;
            }
            double slope = (na*sxy - sx*sy)/(na*sxx - sx*sx);
            double ratRel = fabs((amps.back()/amps.front())/(avals.back()/avals.front()) - 1.0);
            bool g1 = (slope > 0.92 && slope < 1.08), g2 = ratRel < 0.08, g3 = dp < 1e-3;
            pass = g1 && g2 && g3;
            gates = "\"growth_slope_1pm0.08\":" + B(g1) + ",\"growth_ratio_lt_8pc\":" + B(g2)
                  + ",\"dp_lt_1e-3\":" + B(g3);
            resExtra = ",\"a_end\":" + fmt6(avals.back()) + ",\"growth_slope\":" + fmt6(slope)
                     + ",\"growth_ratio_rel\":" + fmt6(ratRel) + ",\"h0\":" + fmt6(gH0)
                     + ",\"amp0\":" + fmt6(amps.front()) + ",\"amp_end\":" + fmt6(amps.back());
        } break;
        case SC_BUBBLES: {
            double srel = (gBubSigExp > 0) ? fabs(gBubSig/gBubSigExp - 1.0) : 1.0;
            bool g1 = (gBubSpawned == 8), g2 = srel < 0.02;
            bool g3 = (gBubCollapsed == 4), g4 = (gBubLive == 4);
            pass = g1 && g2 && g3 && g4;
            gates = "\"spawned_8\":" + B(g1) + ",\"sigma_lt_2pc\":" + B(g2)
                  + ",\"collapsed_4\":" + B(g3) + ",\"alive_4\":" + B(g4);
            resExtra = ",\"n_spawned\":" + fmti((long long)gBubSpawned)
                     + ",\"n_collapsed\":" + fmti((long long)gBubCollapsed)
                     + ",\"n_alive\":" + fmti((long long)gBubLive)
                     + ",\"sigma_meas\":" + fmt6(gBubSig) + ",\"sigma_exp\":" + fmt6(gBubSigExp)
                     + ",\"sigma_rel\":" + fmt6(srel);
        } break;
        default: break;
        }
        const char* solvName = gSolver == SOLV_PM ? "pm" : (gSolver == SOLV_TINY ? "tiny"
                             : (gSolver == SOLV_NONE ? "none" : "direct"));
        std::string body =
              std::string("\"seed\":") + fmti((long long)gSeed)
            + ",\"params\":{\"scenario\":\"" + SC_NAME[gScenario] + "\",\"n\":" + fmti(gN)
            + ",\"steps\":" + fmti(runSteps) + ",\"solver\":\"" + solvName + "\""
            + ",\"c\":" + fmt6(20.0) + ",\"hbar\":" + fmt6(0.5) + ",\"G\":" + fmt6(0.002)
            + ",\"dt\":" + fmt6(1.0/240.0) + ",\"L_box\":" + fmt6(512.0) + "}"
            + ",\"result\":{\"state_b2b\":\"" + shash + "\""
            + ",\"e0\":" + fmt6(gM0.E) + ",\"e1\":" + fmt6(M1.E)
            + ",\"de_rel\":" + fmt6(de) + ",\"p_drift\":" + fmt6(dp) + ",\"l_drift\":" + fmt6(dl)
            + ",\"n_rel\":" + fmti((long long)M1.nRel) + ",\"n_classical\":" + fmti((long long)M1.nCls) + resExtra + "}"
            + ",\"gates\":{" + gates + "}"
            + ",\"verdict\":\"" + (pass ? "pass" : "fail") + "\"";
        std::string env = full_envelope("tinyuniverse", "0.5.0", body,
                                        "scenario run; visuals non-declared");
        if (goldenMode)
            return golden_check(SC_NAME[gScenario], declared_object(body), env);
        printf("%s\n", env.c_str());
        fprintf(stderr, "declared[%s] %s\n", SC_NAME[gScenario],
                blake2b_hex(declared_object(body)).c_str());
        return pass ? 0 : 1;
    }

    if (shotPath){                                    // headless render-to-BMP
        applyPreset(0);
        if (gShotAz   > -998) cAz.tgt = gShotAz;
        if (gShotEl   > -998) cEl.tgt = gShotEl;
        if (gShotDist > -998) cDist.tgt = logf(gShotDist);
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
            if (gScenario == SC_MERGER || gScenario == SC_ECHO){
                double Sn = entropyNow();
                if (gS > 0) gdS = (Sn - gS)/2.0;      // ~2 s cadence
                gS = Sn;
            }
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
