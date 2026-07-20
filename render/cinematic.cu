// render/cinematic.cu — R1: the CINEMATIC post-chain on the R0 interop path.
// Contract: contracts/cinematic.contract.md v0.1.0 APPROVED (operator: "okay try
// supernova first ... proceed"). Scene = SUPERNOVA: ~2500 counter-hashed blackbody
// stars + progenitor star 0 whose light curve flashes x1e5 (hot blue) then decays —
// the in-frame dynamic range that makes bloom and auto-exposure DO something (§2).
// Pipeline (CINEMATIC §1, exact order): fp32 linear accumulation (≥RGBA16F floor) →
// mip-chain threshold-free bloom → log2-luminance histogram auto-exposure (EV,
// tick-smoothed) → AgX (ACES-Narkowicz parity mode) → [cinematic: astro-stretch W=40,
// grain 0.3%] → triangular dither → the ONLY sRGB encode → 8-bit BGRA pack (fixes the
// R0 R/B swizzle) → the shared buffer → swapchain. Splats are analytic gaussians
// (σ≥1.35 px, band-limited — no shimmer; the R1 AA resolution). NO float atomics in
// the declared path (per-pixel fixed-order gather; histogram = integer atomics).
// Faces: windowed (damped orbit, C/T/B/D/E keys, HUD) · --headless --frames N
// [--golden|--json] (fixed camera, NO HUD in declared bytes) · --shot PATH.bmp ·
// --selftest (KATs).
// Build (BUILD.md): nvcc -O3 -arch=sm_89 -o build\cinematic.exe render\cinematic.cu
//                   -I"C:\VulkanSDK\1.4.350.0\Include" "C:\VulkanSDK\1.4.350.0\Lib\vulkan-1.lib" user32.lib
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#define VK_USE_PLATFORM_WIN32_KHR
#include <vulkan/vulkan.h>
#include <cuda_runtime.h>
#include <cstdio>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <string>
#include <vector>

// ----------------------------------------------------------------- blake2b-256 (house idiom)
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
static std::string hash256_hex(const uint8_t* data, size_t n){
    uint64_t h[8]; for (int i=0;i<8;i++) h[i]=IV[i]; h[0]^=0x01010000ull^32ull;
    size_t off=0; uint8_t block[128];
    while (n-off>128){ memcpy(block,data+off,128); off+=128; compress(h,block,(uint64_t)off,false); }
    size_t rem=n-off; memset(block,0,128); if (rem) memcpy(block,data+off,rem);
    compress(h,block,(uint64_t)n,true);
    char hex[65]; for (int i=0;i<32;i++){ unsigned b=(unsigned)((h[i/8]>>(8*(i%8)))&0xFF); snprintf(hex+2*i,3,"%02x",b); }
    return std::string(hex,64);
}
static std::string hash256_hex(const std::string& s){ return hash256_hex((const uint8_t*)s.data(), s.size()); }
}

#define CU_CHECK(x) do{ cudaError_t e=(x); if(e!=cudaSuccess){ \
    fprintf(stderr,"CUDA ERROR %s:%d %s\n",__FILE__,__LINE__,cudaGetErrorString(e)); exit(2);} }while(0)
#define VK_CHECK(x) do{ VkResult r_=(x); if(r_!=VK_SUCCESS){ \
    fprintf(stderr,"VK ERROR %s:%d code=%d\n",__FILE__,__LINE__,(int)r_); exit(2);} }while(0)

static const uint32_t W = 1280, H = 720;
static const size_t   BYTES = (size_t)W*H*4;
static const int NSTAR = 2500;
static const int NBIN = 256;

// ----------------------------------------------------------------- shared device/host math
__host__ __device__ static inline float clampf(float x, float a, float b){ return x<a?a:(x>b?b:x); }

// Tanner–Helland blackbody fit (GARGANTUA pattern) + the 2.2 linearization
__host__ __device__ static void blackbody(float T, float* r, float* g, float* b){
    float t = T*0.01f;
    float R,G,B;
    R = (t<=66.0f)? 255.0f : 329.698727446f*powf(t-60.0f,-0.1332047592f);
    G = (t<=66.0f)? (99.4708025861f*logf(t)-161.1195681661f)
                  : 288.1221695283f*powf(t-60.0f,-0.0755148492f);
    B = (t>=66.0f)? 255.0f : (t<=19.0f? 0.0f : 138.5177312231f*logf(t-10.0f)-305.0447927307f);
    R=clampf(R,0.0f,255.0f)/255.0f; G=clampf(G,0.0f,255.0f)/255.0f; B=clampf(B,0.0f,255.0f)/255.0f;
    *r=powf(R,2.2f); *g=powf(G,2.2f); *b=powf(B,2.2f);
}

// minimal AgX (inset matrix -> log2 window -> 6th-order sigmoid -> outset -> 2.2 EOTF)
__host__ __device__ static void agx3(float r, float g, float b, float* orr, float* og, float* ob){
    // inset
    float ir = 0.842479062253094f*r + 0.0423282422610123f*g + 0.0423756549057051f*b;
    float ig = 0.0784335999999992f*r + 0.878468636469772f *g + 0.0784336f        *b;
    float ib = 0.0792237451477643f*r + 0.0791661274605434f*g + 0.879142973793104f*b;
    // log2 window
    const float LO=-12.47393f, HI=4.026069f;
    float v[3]={ir,ig,ib};
    for (int i=0;i<3;i++){
        float x = clampf((log2f(fmaxf(v[i],1e-10f)) - LO)/(HI-LO), 0.0f, 1.0f);
        float x2=x*x, x4=x2*x2;
        v[i] = 15.5f*x4*x2 - 40.14f*x4*x + 31.96f*x4 - 6.868f*x2*x + 0.4298f*x2 + 0.1191f*x - 0.00232f;
        v[i] = clampf(v[i], 0.0f, 1.0f);
    }
    // outset
    float rr = 1.19687900512017f  *v[0] - 0.0528968517574562f*v[1] - 0.0529716355144438f*v[2];
    float gg = -0.0980208811401368f*v[0] + 1.15190312990417f  *v[1] - 0.0980434501171241f*v[2];
    float bb = -0.0990297440797205f*v[0] - 0.0989611768448433f*v[1] + 1.15107367264116f  *v[2];
    // AgX EOTF (2.2) -> LINEAR out; the pipeline's single sRGB encode follows later
    *orr = powf(clampf(rr,0.0f,1.0f), 2.2f);
    *og  = powf(clampf(gg,0.0f,1.0f), 2.2f);
    *ob  = powf(clampf(bb,0.0f,1.0f), 2.2f);
}
__host__ __device__ static float aces1(float x){       // Narkowicz (GARGANTUA parity), per channel
    return clampf(x*(2.51f*x+0.03f)/(x*(2.43f*x+0.59f)+0.14f), 0.0f, 1.0f);
}
__host__ __device__ static float srgb1(float c){
    c = clampf(c, 0.0f, 1.0f);
    return (c<=0.0031308f)? 12.92f*c : 1.055f*powf(c,1.0f/2.4f)-0.055f;
}
__host__ __device__ static uint32_t hash32(uint32_t x){
    x ^= x>>16; x *= 0x7feb352dU; x ^= x>>15; x *= 0x846ca68bU; x ^= x>>16; return x;
}
// supernova light curve — pure function of the tick (loop period 1200)
__host__ __device__ static void snCurve(uint32_t n, float* L, float* T){
    float t = (float)(n % 1200u);
    const float L0=2.0f, T0=5800.0f;
    if (t < 60.0f){ *L=L0; *T=T0; return; }
    if (t < 150.0f){
        float u=(t-60.0f)/90.0f;
        *L = L0*powf(10.0f, 5.0f*u*u*(3.0f-2.0f*u));       // smoothstep in log-L
        *T = T0 + (15000.0f-T0)*u;
        return;
    }
    float d=(t-150.0f);
    *L = L0*1e5f*expf(-d/180.0f) + L0;
    *T = 3500.0f + (15000.0f-3500.0f)*expf(-d/120.0f);
}

// ----------------------------------------------------------------- kernels
__global__ void kStarPrep(const float4* dirT, const float* lum, float4* sxy, float4* srgb,
                          float4 Rx, float4 Ry, float4 Rz, float fpx, uint32_t n){
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= NSTAR) return;
    float4 d = dirT[i];
    float L = lum[i], T = d.w;
    if (i == 0) snCurve(n, &L, &T);
    float cx = Rx.x*d.x + Rx.y*d.y + Rx.z*d.z;
    float cy = Ry.x*d.x + Ry.y*d.y + Ry.z*d.z;
    float cz = Rz.x*d.x + Rz.y*d.y + Rz.z*d.z;
    float4 o; o.z = 0.0f; o.w = 0.0f;
    if (cz < 0.05f){ o.x = -1e6f; o.y = -1e6f; }
    else { o.x = 0.5f*W + fpx*cx/cz; o.y = 0.5f*H + fpx*cy/cz; }
    sxy[i] = o;
    float r,g,b; blackbody(T,&r,&g,&b);
    srgb[i] = make_float4(L*r, L*g, L*b, 0.0f);
}

__global__ void kGather(float4* acc, const float4* sxy, const float4* srgb){
    int x = blockIdx.x*blockDim.x + threadIdx.x;
    int y = blockIdx.y*blockDim.y + threadIdx.y;
    if (x >= (int)W || y >= (int)H) return;
    const float SIG = 1.35f, R = 7.0f, inv2s2 = 1.0f/(2.0f*SIG*SIG);
    float r=0.0f, g=0.0f, b=0.0f;
    for (int i=0;i<NSTAR;i++){
        float dx = (float)x + 0.5f - sxy[i].x;
        float dy = (float)y + 0.5f - sxy[i].y;
        if (dx>R || dx<-R || dy>R || dy<-R) continue;
        float w = expf(-(dx*dx+dy*dy)*inv2s2);
        float4 c = srgb[i];
        r += w*c.x; g += w*c.y; b += w*c.z;
    }
    acc[y*W+x] = make_float4(0.008f*r, 0.008f*g, 0.008f*b, 1.0f);
}

__global__ void kHist(const float4* acc, unsigned int* hist){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if (idx >= (int)(W*H)) return;
    float4 c = acc[idx];
    float lum = 0.2126f*c.x + 0.7152f*c.y + 0.0722f*c.z;
    float l2 = log2f(fmaxf(lum, 1e-8f));
    int bin = (int)clampf((l2+16.0f)/32.0f*(float)NBIN, 0.0f, (float)(NBIN-1));
    atomicAdd(&hist[bin], 1u);                       // integer atomics: order-independent, house-legal
}

__global__ void kDown13(const float4* src, int sw, int sh, float4* dst, int dw, int dh){
    int x = blockIdx.x*blockDim.x + threadIdx.x;
    int y = blockIdx.y*blockDim.y + threadIdx.y;
    if (x >= dw || y >= dh) return;
    // CoD-style: center 2x2 box (w 0.5) + four corner 2x2 boxes (w 0.125 each)
    auto box = [&](int cx, int cy){
        int x0=min(max(cx,0),sw-1), x1=min(cx+1,sw-1);
        int y0=min(max(cy,0),sh-1), y1=min(cy+1,sh-1);
        float4 a=src[y0*sw+x0], b=src[y0*sw+x1], c=src[y1*sw+x0], d=src[y1*sw+x1];
        return make_float4(0.25f*(a.x+b.x+c.x+d.x),0.25f*(a.y+b.y+c.y+d.y),0.25f*(a.z+b.z+c.z+d.z),1.0f);
    };
    int sx=2*x, sy=2*y;
    float4 C=box(sx,sy), TL=box(sx-2,sy-2), TR=box(sx+2,sy-2), BL=box(sx-2,sy+2), BR=box(sx+2,sy+2);
    dst[y*dw+x] = make_float4(
        0.5f*C.x + 0.125f*(TL.x+TR.x+BL.x+BR.x),
        0.5f*C.y + 0.125f*(TL.y+TR.y+BL.y+BR.y),
        0.5f*C.z + 0.125f*(TL.z+TR.z+BL.z+BR.z), 1.0f);
}

__global__ void kUpTent(const float4* coarse, int cw, int ch, float4* fine, int fw, int fh){
    int x = blockIdx.x*blockDim.x + threadIdx.x;
    int y = blockIdx.y*blockDim.y + threadIdx.y;
    if (x >= fw || y >= fh) return;
    int cx = x/2, cy = y/2;
    auto S = [&](int a, int b){
        a=min(max(a,0),cw-1); b=min(max(b,0),ch-1);
        return coarse[b*cw+a];
    };
    const float w0=4.0f/16.0f, w1=2.0f/16.0f, w2=1.0f/16.0f;
    float4 s0=S(cx,cy);
    float4 sl=S(cx-1,cy), sr=S(cx+1,cy), su=S(cx,cy-1), sd=S(cx,cy+1);
    float4 a=S(cx-1,cy-1), b=S(cx+1,cy-1), c=S(cx-1,cy+1), d=S(cx+1,cy+1);
    float4 o = fine[y*fw+x];
    o.x += w0*s0.x + w1*(sl.x+sr.x+su.x+sd.x) + w2*(a.x+b.x+c.x+d.x);
    o.y += w0*s0.y + w1*(sl.y+sr.y+su.y+sd.y) + w2*(a.y+b.y+c.y+d.y);
    o.z += w0*s0.z + w1*(sl.z+sr.z+su.z+sd.z) + w2*(a.z+b.z+c.z+d.z);
    fine[y*fw+x] = o;
}

// tone modes: 0 AgX, 1 ACES, 2 raw clamp (the didactic shocker)
__global__ void kCompose(const float4* acc, const float4* bloom, uchar4* out,
                         float exposure, float bloomMix, int tone, int cinematic,
                         int dither, uint32_t frame){
    int x = blockIdx.x*blockDim.x + threadIdx.x;
    int y = blockIdx.y*blockDim.y + threadIdx.y;
    if (x >= (int)W || y >= (int)H) return;
    int idx = y*W+x;
    float4 a = acc[idx], bl = bloom[idx];
    float r = (a.x + bloomMix*bl.x)*exposure;
    float g = (a.y + bloomMix*bl.y)*exposure;
    float b = (a.z + bloomMix*bl.z)*exposure;
    if (cinematic){
        // astro stretch (CINEMATIC §4a): identity for L<=1, extended Reinhard on L-1, W=40
        float L = 0.2126f*r + 0.7152f*g + 0.0722f*b;
        if (L > 1.0f){
            float Le = L-1.0f, Wp = 40.0f;
            float Ls = 1.0f + Le*(1.0f+Le/(Wp*Wp))/(1.0f+Le);
            float s = Ls/L;
            r*=s; g*=s; b*=s;
        }
    }
    float tr,tg,tb;
    if (tone==0) agx3(r,g,b,&tr,&tg,&tb);
    else if (tone==1){ tr=aces1(r); tg=aces1(g); tb=aces1(b); }
    else { tr=clampf(r,0.0f,1.0f); tg=clampf(g,0.0f,1.0f); tb=clampf(b,0.0f,1.0f); }
    if (cinematic){
        // grain 0.3% (§1 bound 0.5%)
        uint32_t hg = hash32((uint32_t)idx ^ (frame*0x9e3779b9u) ^ 0x51ed270bu);
        float gn = ((float)(hg & 0xFFFF)/65535.0f - 0.5f)*0.006f;
        tr = clampf(tr+gn,0.0f,1.0f); tg = clampf(tg+gn,0.0f,1.0f); tb = clampf(tb+gn,0.0f,1.0f);
    }
    // sRGB encode — the ONLY gamma
    tr = srgb1(tr); tg = srgb1(tg); tb = srgb1(tb);
    // triangular dither BEFORE quantization (mandatory)
    float dr=0.0f;
    if (dither){
        uint32_t h1 = hash32((uint32_t)idx ^ frame*0x85ebca6bu);
        uint32_t h2 = hash32(h1 ^ 0xc2b2ae35u);
        dr = (((float)(h1&0xFFFF)+(float)(h2&0xFFFF))/65535.0f - 1.0f)*(1.0f/255.0f);
    }
    unsigned char R8 = (unsigned char)clampf(255.0f*(tr+dr)+0.5f, 0.0f, 255.0f);
    unsigned char G8 = (unsigned char)clampf(255.0f*(tg+dr)+0.5f, 0.0f, 255.0f);
    unsigned char B8 = (unsigned char)clampf(255.0f*(tb+dr)+0.5f, 0.0f, 255.0f);
    out[idx] = make_uchar4(B8, G8, R8, 255);          // BGRA (the R0 swizzle, fixed here)
}

// minimal 5x7 HUD font (windowed only; NOT in the declared bytes) — render-amber
__constant__ unsigned char FONT[16*7] = {
    // 0-9 A C E G S V X M N P  (16 glyph slots x 7 rows, 5-bit rows)
    0x0E,0x11,0x13,0x15,0x19,0x11,0x0E,  0x04,0x0C,0x04,0x04,0x04,0x04,0x0E,
    0x0E,0x11,0x01,0x06,0x08,0x10,0x1F,  0x0E,0x11,0x01,0x06,0x01,0x11,0x0E,
    0x02,0x06,0x0A,0x12,0x1F,0x02,0x02,  0x1F,0x10,0x1E,0x01,0x01,0x11,0x0E,
    0x06,0x08,0x10,0x1E,0x11,0x11,0x0E,  0x1F,0x01,0x02,0x04,0x08,0x08,0x08,
    0x0E,0x11,0x11,0x0E,0x11,0x11,0x0E,  0x0E,0x11,0x11,0x0F,0x01,0x02,0x0C,
    0x0E,0x11,0x11,0x1F,0x11,0x11,0x11,  0x0E,0x11,0x10,0x10,0x10,0x11,0x0E,
    0x1F,0x10,0x10,0x1E,0x10,0x10,0x1F,  0x0F,0x10,0x10,0x13,0x11,0x11,0x0F,
    0x0F,0x10,0x10,0x0E,0x01,0x01,0x1E,  0x11,0x11,0x11,0x15,0x15,0x1B,0x11,
};
// glyph map: '0'-'9' -> 0..9, 'A'=10 'C'=11 'E'=12 'G'=13 'S'=14 'M'=15
__global__ void kHud(uchar4* out, const unsigned char* text, int len, int ox, int oy){
    int i = blockIdx.x;                     // one block per char
    int px = threadIdx.x, py = threadIdx.y; // 5x7
    if (i>=len || px>=5 || py>=7) return;
    unsigned char ch = text[i];
    if (ch == 255) return;                  // space
    unsigned char row = FONT[ch*7+py];
    if (!(row & (0x10>>px))) return;
    int X = ox + i*12 + px*2, Y = oy + py*2;
    for (int dy=0;dy<2;dy++) for (int dx=0;dx<2;dx++){
        int xx=X+dx, yy=Y+dy;
        if (xx>=0 && xx<(int)W && yy>=0 && yy<(int)H)
            out[yy*W+xx] = make_uchar4(30, 190, 255, 255);   // render-amber (BGRA)
    }
}

// ----------------------------------------------------------------- Vulkan plumbing (R0 pattern)
struct Ctx {
    bool validation=false, windowed=false;
    VkInstance inst=VK_NULL_HANDLE;
    VkDebugUtilsMessengerEXT dbg=VK_NULL_HANDLE;
    VkPhysicalDevice phys=VK_NULL_HANDLE;
    VkDevice dev=VK_NULL_HANDLE;
    uint32_t qfam=0;
    VkQueue q=VK_NULL_HANDLE;
    VkCommandPool pool=VK_NULL_HANDLE;
    char devName[256]={0};
    uint8_t cuLuid[8]={0};
    bool luidMatch=false;
    VkBuffer sharedBuf=VK_NULL_HANDLE;
    VkDeviceMemory sharedMem=VK_NULL_HANDLE;
    cudaExternalMemory_t extMem=nullptr;
    uchar4* dOut=nullptr;
    VkSemaphore tlCuda=VK_NULL_HANDLE, tlVk=VK_NULL_HANDLE;
    cudaExternalSemaphore_t extCuda=nullptr, extVk=nullptr;
    uint64_t vErrors=0, vWarnings=0;
};
static Ctx G;

static VKAPI_ATTR VkBool32 VKAPI_CALL dbgCb(VkDebugUtilsMessageSeverityFlagBitsEXT sev,
        VkDebugUtilsMessageTypeFlagsEXT, const VkDebugUtilsMessengerCallbackDataEXT* data, void*){
    if (sev & VK_DEBUG_UTILS_MESSAGE_SEVERITY_ERROR_BIT_EXT){ G.vErrors++; fprintf(stderr,"[VK-ERROR] %s\n",data->pMessage); }
    else if (sev & VK_DEBUG_UTILS_MESSAGE_SEVERITY_WARNING_BIT_EXT){ G.vWarnings++; fprintf(stderr,"[VK-WARN ] %s\n",data->pMessage); }
    return VK_FALSE;
}
static void createInstance(){
    VkApplicationInfo ai{VK_STRUCTURE_TYPE_APPLICATION_INFO};
    ai.pApplicationName="tinyuniverse-cinematic-r1"; ai.apiVersion=VK_API_VERSION_1_3;
    std::vector<const char*> exts, layers;
    if (G.windowed){ exts.push_back("VK_KHR_surface"); exts.push_back("VK_KHR_win32_surface"); }
    if (G.validation){ exts.push_back(VK_EXT_DEBUG_UTILS_EXTENSION_NAME); layers.push_back("VK_LAYER_KHRONOS_validation"); }
    VkInstanceCreateInfo ci{VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO};
    ci.pApplicationInfo=&ai;
    ci.enabledExtensionCount=(uint32_t)exts.size(); ci.ppEnabledExtensionNames=exts.data();
    ci.enabledLayerCount=(uint32_t)layers.size();   ci.ppEnabledLayerNames=layers.data();
    VK_CHECK(vkCreateInstance(&ci,nullptr,&G.inst));
    if (G.validation){
        auto f=(PFN_vkCreateDebugUtilsMessengerEXT)vkGetInstanceProcAddr(G.inst,"vkCreateDebugUtilsMessengerEXT");
        VkDebugUtilsMessengerCreateInfoEXT di{VK_STRUCTURE_TYPE_DEBUG_UTILS_MESSENGER_CREATE_INFO_EXT};
        di.messageSeverity=VK_DEBUG_UTILS_MESSAGE_SEVERITY_WARNING_BIT_EXT|VK_DEBUG_UTILS_MESSAGE_SEVERITY_ERROR_BIT_EXT;
        di.messageType=VK_DEBUG_UTILS_MESSAGE_TYPE_GENERAL_BIT_EXT|VK_DEBUG_UTILS_MESSAGE_TYPE_VALIDATION_BIT_EXT|
                       VK_DEBUG_UTILS_MESSAGE_TYPE_PERFORMANCE_BIT_EXT;
        di.pfnUserCallback=dbgCb;
        if (f) VK_CHECK(f(G.inst,&di,nullptr,&G.dbg));
    }
}
static void pickDeviceByLuid(){
    CU_CHECK(cudaSetDevice(0));
    cudaDeviceProp p{}; CU_CHECK(cudaGetDeviceProperties(&p,0));
    memcpy(G.cuLuid,p.luid,8);
    uint32_t n=0; VK_CHECK(vkEnumeratePhysicalDevices(G.inst,&n,nullptr));
    std::vector<VkPhysicalDevice> devs(n);
    VK_CHECK(vkEnumeratePhysicalDevices(G.inst,&n,devs.data()));
    for (auto d:devs){
        VkPhysicalDeviceIDProperties idp{VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_ID_PROPERTIES};
        VkPhysicalDeviceProperties2 p2{VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_PROPERTIES_2};
        p2.pNext=&idp;
        vkGetPhysicalDeviceProperties2(d,&p2);
        if (idp.deviceLUIDValid && memcmp(idp.deviceLUID,G.cuLuid,8)==0){
            G.phys=d; snprintf(G.devName,sizeof(G.devName),"%s",p2.properties.deviceName);
            G.luidMatch=true; return;
        }
    }
    fprintf(stderr,"G-DEVICE FAIL: no LUID match\n"); exit(2);
}
static void createDevice(){
    uint32_t n=0; vkGetPhysicalDeviceQueueFamilyProperties(G.phys,&n,nullptr);
    std::vector<VkQueueFamilyProperties> qf(n);
    vkGetPhysicalDeviceQueueFamilyProperties(G.phys,&n,qf.data());
    bool found=false;
    for (uint32_t i=0;i<n;i++) if (qf[i].queueFlags & VK_QUEUE_GRAPHICS_BIT){ G.qfam=i; found=true; break; }
    if (!found){ fprintf(stderr,"no graphics queue\n"); exit(2); }
    float prio=1.0f;
    VkDeviceQueueCreateInfo qi{VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO};
    qi.queueFamilyIndex=G.qfam; qi.queueCount=1; qi.pQueuePriorities=&prio;
    std::vector<const char*> exts={
        VK_KHR_EXTERNAL_MEMORY_EXTENSION_NAME, VK_KHR_EXTERNAL_MEMORY_WIN32_EXTENSION_NAME,
        VK_KHR_EXTERNAL_SEMAPHORE_EXTENSION_NAME, VK_KHR_EXTERNAL_SEMAPHORE_WIN32_EXTENSION_NAME };
    if (G.windowed) exts.push_back(VK_KHR_SWAPCHAIN_EXTENSION_NAME);
    VkPhysicalDeviceVulkan12Features f12{VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_VULKAN_1_2_FEATURES};
    f12.timelineSemaphore=VK_TRUE;
    VkDeviceCreateInfo ci{VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO};
    ci.pNext=&f12; ci.queueCreateInfoCount=1; ci.pQueueCreateInfos=&qi;
    ci.enabledExtensionCount=(uint32_t)exts.size(); ci.ppEnabledExtensionNames=exts.data();
    VK_CHECK(vkCreateDevice(G.phys,&ci,nullptr,&G.dev));
    vkGetDeviceQueue(G.dev,G.qfam,0,&G.q);
    VkCommandPoolCreateInfo pi{VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO};
    pi.queueFamilyIndex=G.qfam; pi.flags=VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT;
    VK_CHECK(vkCreateCommandPool(G.dev,&pi,nullptr,&G.pool));
}
static uint32_t memType(uint32_t bits, VkMemoryPropertyFlags props){
    VkPhysicalDeviceMemoryProperties mp;
    vkGetPhysicalDeviceMemoryProperties(G.phys,&mp);
    for (uint32_t i=0;i<mp.memoryTypeCount;i++)
        if ((bits&(1u<<i)) && (mp.memoryTypes[i].propertyFlags&props)==props) return i;
    fprintf(stderr,"no memory type\n"); exit(2);
}
static void createSharedBuffer(){
    VkExternalMemoryBufferCreateInfo ext{VK_STRUCTURE_TYPE_EXTERNAL_MEMORY_BUFFER_CREATE_INFO};
    ext.handleTypes=VK_EXTERNAL_MEMORY_HANDLE_TYPE_OPAQUE_WIN32_BIT;
    VkBufferCreateInfo bi{VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO};
    bi.pNext=&ext; bi.size=BYTES;
    bi.usage=VK_BUFFER_USAGE_TRANSFER_SRC_BIT|VK_BUFFER_USAGE_TRANSFER_DST_BIT;
    bi.sharingMode=VK_SHARING_MODE_EXCLUSIVE;
    VK_CHECK(vkCreateBuffer(G.dev,&bi,nullptr,&G.sharedBuf));
    VkMemoryRequirements mr; vkGetBufferMemoryRequirements(G.dev,G.sharedBuf,&mr);
    VkExportMemoryAllocateInfo exp{VK_STRUCTURE_TYPE_EXPORT_MEMORY_ALLOCATE_INFO};
    exp.handleTypes=VK_EXTERNAL_MEMORY_HANDLE_TYPE_OPAQUE_WIN32_BIT;
    VkMemoryAllocateInfo mai{VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO};
    mai.pNext=&exp; mai.allocationSize=mr.size;
    mai.memoryTypeIndex=memType(mr.memoryTypeBits,VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
    VK_CHECK(vkAllocateMemory(G.dev,&mai,nullptr,&G.sharedMem));
    VK_CHECK(vkBindBufferMemory(G.dev,G.sharedBuf,G.sharedMem,0));
    auto getH=(PFN_vkGetMemoryWin32HandleKHR)vkGetDeviceProcAddr(G.dev,"vkGetMemoryWin32HandleKHR");
    if (!getH){ fprintf(stderr,"no vkGetMemoryWin32HandleKHR\n"); exit(2); }
    VkMemoryGetWin32HandleInfoKHR gi{VK_STRUCTURE_TYPE_MEMORY_GET_WIN32_HANDLE_INFO_KHR};
    gi.memory=G.sharedMem; gi.handleType=VK_EXTERNAL_MEMORY_HANDLE_TYPE_OPAQUE_WIN32_BIT;
    HANDLE h=nullptr; VK_CHECK(getH(G.dev,&gi,&h));
    cudaExternalMemoryHandleDesc hd{};
    hd.type=cudaExternalMemoryHandleTypeOpaqueWin32; hd.handle.win32.handle=h; hd.size=mr.size;
    CU_CHECK(cudaImportExternalMemory(&G.extMem,&hd));
    CloseHandle(h);
    cudaExternalMemoryBufferDesc bd{}; bd.offset=0; bd.size=BYTES;
    void* dp=nullptr; CU_CHECK(cudaExternalMemoryGetMappedBuffer(&dp,G.extMem,&bd));
    G.dOut=(uchar4*)dp;
}
static VkSemaphore createTimeline(cudaExternalSemaphore_t* cuOut){
    VkExportSemaphoreCreateInfo exp{VK_STRUCTURE_TYPE_EXPORT_SEMAPHORE_CREATE_INFO};
    exp.handleTypes=VK_EXTERNAL_SEMAPHORE_HANDLE_TYPE_OPAQUE_WIN32_BIT;
    VkSemaphoreTypeCreateInfo ty{VK_STRUCTURE_TYPE_SEMAPHORE_TYPE_CREATE_INFO};
    ty.pNext=&exp; ty.semaphoreType=VK_SEMAPHORE_TYPE_TIMELINE; ty.initialValue=0;
    VkSemaphoreCreateInfo ci{VK_STRUCTURE_TYPE_SEMAPHORE_CREATE_INFO};
    ci.pNext=&ty;
    VkSemaphore s; VK_CHECK(vkCreateSemaphore(G.dev,&ci,nullptr,&s));
    auto getH=(PFN_vkGetSemaphoreWin32HandleKHR)vkGetDeviceProcAddr(G.dev,"vkGetSemaphoreWin32HandleKHR");
    if (!getH){ fprintf(stderr,"no vkGetSemaphoreWin32HandleKHR\n"); exit(2); }
    VkSemaphoreGetWin32HandleInfoKHR gi{VK_STRUCTURE_TYPE_SEMAPHORE_GET_WIN32_HANDLE_INFO_KHR};
    gi.semaphore=s; gi.handleType=VK_EXTERNAL_SEMAPHORE_HANDLE_TYPE_OPAQUE_WIN32_BIT;
    HANDLE h=nullptr; VK_CHECK(getH(G.dev,&gi,&h));
    cudaExternalSemaphoreHandleDesc sd{};
    sd.type=cudaExternalSemaphoreHandleTypeTimelineSemaphoreWin32; sd.handle.win32.handle=h;
    CU_CHECK(cudaImportExternalSemaphore(cuOut,&sd));
    CloseHandle(h);
    return s;
}
static void cuWait(cudaExternalSemaphore_t s, uint64_t v){
    cudaExternalSemaphoreWaitParams p{}; p.params.fence.value=v;
    CU_CHECK(cudaWaitExternalSemaphoresAsync(&s,&p,1,0));
}
static void cuSignal(cudaExternalSemaphore_t s, uint64_t v){
    cudaExternalSemaphoreSignalParams p{}; p.params.fence.value=v;
    CU_CHECK(cudaSignalExternalSemaphoresAsync(&s,&p,1,0));
}
static void submitTimeline(VkCommandBuffer cb, uint64_t v,
                           VkSemaphore extraWaitBin=VK_NULL_HANDLE, VkSemaphore extraSigBin=VK_NULL_HANDLE){
    uint64_t waitVals[2]={v,0}, sigVals[2]={v,0};
    VkSemaphore waits[2]={G.tlCuda,extraWaitBin}, sigs[2]={G.tlVk,extraSigBin};
    VkPipelineStageFlags st[2]={VK_PIPELINE_STAGE_TRANSFER_BIT,VK_PIPELINE_STAGE_TRANSFER_BIT};
    uint32_t nw=extraWaitBin?2u:1u, ns=extraSigBin?2u:1u;
    VkTimelineSemaphoreSubmitInfo tsi{VK_STRUCTURE_TYPE_TIMELINE_SEMAPHORE_SUBMIT_INFO};
    tsi.waitSemaphoreValueCount=nw; tsi.pWaitSemaphoreValues=waitVals;
    tsi.signalSemaphoreValueCount=ns; tsi.pSignalSemaphoreValues=sigVals;
    VkSubmitInfo si{VK_STRUCTURE_TYPE_SUBMIT_INFO};
    si.pNext=&tsi;
    si.waitSemaphoreCount=nw; si.pWaitSemaphores=waits; si.pWaitDstStageMask=st;
    si.commandBufferCount=1; si.pCommandBuffers=&cb;
    si.signalSemaphoreCount=ns; si.pSignalSemaphores=sigs;
    VK_CHECK(vkQueueSubmit(G.q,1,&si,VK_NULL_HANDLE));
}
static void hostWaitVkDone(uint64_t v){
    VkSemaphoreWaitInfo wi{VK_STRUCTURE_TYPE_SEMAPHORE_WAIT_INFO};
    wi.semaphoreCount=1; wi.pSemaphores=&G.tlVk; wi.pValues=&v;
    VK_CHECK(vkWaitSemaphores(G.dev,&wi,~0ull));
}
static void setupCommon(){
    createInstance(); pickDeviceByLuid(); createDevice(); createSharedBuffer();
    G.tlCuda=createTimeline(&G.extCuda); G.tlVk=createTimeline(&G.extVk);
}

// ----------------------------------------------------------------- the renderer core
struct Pipe {
    float4 *acc=nullptr, *mip[6]={nullptr};
    int mw[6], mh[6];
    unsigned int* dHist=nullptr;
    float4 *dSxy=nullptr, *dSrgb=nullptr, *dDirT=nullptr;
    float *dLum=nullptr;
    unsigned char* dText=nullptr;
    float ev=0.0f;                 // smoothed log2 scene key
    bool evInit=false;
};
static Pipe P;

static uint64_t splitmix(uint64_t& s){
    s += 0x9E3779B97f4A7C15ull;
    uint64_t z=s;
    z=(z^(z>>30))*0xBF58476D1CE4E5B9ull;
    z=(z^(z>>27))*0x94D049BB133111EBull;
    return z^(z>>31);
}
static void initScene(){
    std::vector<float4> dirT(NSTAR);
    std::vector<float> lum(NSTAR);
    uint64_t s=42ull;
    for (int i=0;i<NSTAR;i++){
        double u1=(double)(splitmix(s)>>11)/9007199254740992.0;
        double u2=(double)(splitmix(s)>>11)/9007199254740992.0;
        double u3=(double)(splitmix(s)>>11)/9007199254740992.0;
        double u4=(double)(splitmix(s)>>11)/9007199254740992.0;
        float z=(float)(2.0*u1-1.0), ph=(float)(6.28318530717958648*u2);
        float rxy=sqrtf(fmaxf(0.0f,1.0f-z*z));
        dirT[i]=make_float4(rxy*cosf(ph), rxy*sinf(ph), z, 3200.0f+7800.0f*(float)(u3*u3*u3));
        float L=0.02f+0.6f*(float)(u4*u4*u4);
        if ((splitmix(s)&0xFF) < 3) L*=40.0f;          // a few bright foreground stars
        lum[i]=L;
    }
    // star 0 = the progenitor: fixed direction near the view axis
    float3 d0=make_float3(0.10f,0.04f,1.0f);
    float n0=sqrtf(d0.x*d0.x+d0.y*d0.y+d0.z*d0.z);
    dirT[0]=make_float4(d0.x/n0,d0.y/n0,d0.z/n0,5800.0f);
    lum[0]=2.0f;
    CU_CHECK(cudaMalloc(&P.dDirT,NSTAR*sizeof(float4)));
    CU_CHECK(cudaMalloc(&P.dLum,NSTAR*sizeof(float)));
    CU_CHECK(cudaMalloc(&P.dSxy,NSTAR*sizeof(float4)));
    CU_CHECK(cudaMalloc(&P.dSrgb,NSTAR*sizeof(float4)));
    CU_CHECK(cudaMemcpy(P.dDirT,dirT.data(),NSTAR*sizeof(float4),cudaMemcpyHostToDevice));
    CU_CHECK(cudaMemcpy(P.dLum,lum.data(),NSTAR*sizeof(float),cudaMemcpyHostToDevice));
    CU_CHECK(cudaMalloc(&P.acc,W*H*sizeof(float4)));
    int w=W,h=H;
    for (int m=0;m<6;m++){ w=(w+1)/2; h=(h+1)/2; P.mw[m]=w; P.mh[m]=h; CU_CHECK(cudaMalloc(&P.mip[m],(size_t)w*h*sizeof(float4))); }
    CU_CHECK(cudaMalloc(&P.dHist,NBIN*sizeof(unsigned int)));
    CU_CHECK(cudaMalloc(&P.dText,64));
}

// camera basis from az/el (orbit looking at origin-forward)
static void camBasis(float az, float el, float4* Rx, float4* Ry, float4* Rz){
    float ca=cosf(az), sa=sinf(az), ce=cosf(el), se=sinf(el);
    // forward
    float fx=ce*sa, fy=se, fz=ce*ca;
    // right = normalize(cross(up,f)), up=(0,1,0)
    float rx=fz, ry=0.0f, rz=-fx;
    float rn=sqrtf(rx*rx+rz*rz)+1e-9f; rx/=rn; rz/=rn;
    // true up = cross(f,right)
    float ux=fy*rz-fz*ry, uy=fz*rx-fx*rz, uz=fx*ry-fy*rx;
    *Rx=make_float4(rx,ry,rz,0.0f);
    *Ry=make_float4(ux,uy,uz,0.0f);
    *Rz=make_float4(fx,fy,fz,0.0f);
}

struct FrameStats { float evUsed; float rangeLog2; };

// render one frame into G.dOut; tone 0/1/2, cinematic flag, dither flag; hud text optional
static FrameStats renderFrame(uint32_t n, float az, float el, int tone, int cinematic,
                              int dither, const char* hud){
    const float FOVdeg=58.0f;
    float fpx = 0.5f*(float)W/tanf(0.5f*FOVdeg*3.14159265358979f/180.0f);
    float4 Rx,Ry,Rz; camBasis(az,el,&Rx,&Ry,&Rz);
    kStarPrep<<<(NSTAR+127)/128,128>>>(P.dDirT,P.dLum,P.dSxy,P.dSrgb,Rx,Ry,Rz,fpx,n);
    dim3 bs(16,16), gs((W+15)/16,(H+15)/16);
    kGather<<<gs,bs>>>(P.acc,P.dSxy,P.dSrgb);
    // exposure meter (2%-98% log2-luminance histogram)
    CU_CHECK(cudaMemset(P.dHist,0,NBIN*sizeof(unsigned int)));
    kHist<<<(W*H+255)/256,256>>>(P.acc,P.dHist);
    unsigned int hist[NBIN];
    CU_CHECK(cudaMemcpy(hist,P.dHist,sizeof(hist),cudaMemcpyDeviceToHost));
    // astro metering, ENERGY-WEIGHTED: key = sum(E*log2L)/sum(E) over lit bins
    // (measured, two steps: with >96% true-black pixels the 2-98% percentile band sat
    // in bin 0, EV pinned at -15.94; excluding bin 0 the percentile MEAN still moved
    // only 0.03 EV at the flash — percentile statistics are robust against small
    // bright regions BY DESIGN, and a supernova IS one. A camera is blinded by a
    // point flash because the flash carries the frame's ENERGY: weight by luminance.)
    double sum=0.0, wsum=0.0;
    int minB=NBIN, maxB=-1;
    for (int i=1;i<NBIN;i++){
        if (!hist[i]) continue;
        if (i<minB) minB=i;
        if (i>maxB) maxB=i;
        double l2=-16.0+32.0*((double)i+0.5)/NBIN;
        double E=(double)hist[i]*exp2(l2);
        sum+=E*l2; wsum+=E;
    }
    float key = (wsum>0.0)? (float)(sum/wsum) : -8.0f;
    const float TAU=90.0f;
    float alpha = 1.0f-expf(-1.0f/TAU);
    if (!P.evInit){ P.ev=key; P.evInit=true; } else P.ev += (key-P.ev)*alpha;
    float exposure = 0.18f/exp2f(P.ev);
    float rangeLog2 = (maxB>=minB)? 32.0f*(float)(maxB-minB)/NBIN : 0.0f;
    // bloom pyramid
    kDown13<<<dim3((P.mw[0]+15)/16,(P.mh[0]+15)/16),bs>>>(P.acc,W,H,P.mip[0],P.mw[0],P.mh[0]);
    for (int m=1;m<6;m++)
        kDown13<<<dim3((P.mw[m]+15)/16,(P.mh[m]+15)/16),bs>>>(P.mip[m-1],P.mw[m-1],P.mh[m-1],P.mip[m],P.mw[m],P.mh[m]);
    for (int m=4;m>=0;m--)
        kUpTent<<<dim3((P.mw[m]+15)/16,(P.mh[m]+15)/16),bs>>>(P.mip[m+1],P.mw[m+1],P.mh[m+1],P.mip[m],P.mw[m],P.mh[m]);
    // compose into a base-res bloom buffer: upsample mip0 -> full res happens inside kCompose?
    // simplest: tent-upsample mip[0] into acc-sized bloom = reuse mip storage trick — allocate once:
    static float4* bloomFull=nullptr;
    if (!bloomFull) CU_CHECK(cudaMalloc(&bloomFull,W*H*sizeof(float4)));
    CU_CHECK(cudaMemset(bloomFull,0,W*H*sizeof(float4)));
    kUpTent<<<gs,bs>>>(P.mip[0],P.mw[0],P.mh[0],bloomFull,W,H);
    kCompose<<<gs,bs>>>(P.acc,bloomFull,G.dOut,exposure,0.06f,tone,cinematic,dither,n);
    if (hud){
        unsigned char txt[64]; int len=0;
        for (const char* c=hud; *c && len<64; c++){
            unsigned char g=255;
            if (*c>='0'&&*c<='9') g=(unsigned char)(*c-'0');
            else if (*c=='A') g=10; else if (*c=='C') g=11; else if (*c=='E') g=12;
            else if (*c=='G') g=13; else if (*c=='S') g=14; else if (*c=='M') g=15;
            txt[len++]=g;
        }
        CU_CHECK(cudaMemcpy(P.dText,txt,len,cudaMemcpyHostToDevice));
        kHud<<<len,dim3(5,7)>>>(G.dOut,P.dText,len,12,12);
    }
    FrameStats fs; fs.evUsed=P.ev; fs.rangeLog2=rangeLog2;
    return fs;
}

// ----------------------------------------------------------------- BMP writer (--shot)
static void writeBMP(const char* path, const uint8_t* bgra){
    FILE* f=fopen(path,"wb");
    if (!f){ fprintf(stderr,"cannot open %s\n",path); exit(2); }
    uint32_t rowB=W*3, pad=(4-(rowB&3))&3, img=(rowB+pad)*H, off=54, size=off+img;
    uint8_t hdr[54]={0};
    hdr[0]='B'; hdr[1]='M';
    memcpy(hdr+2,&size,4); memcpy(hdr+10,&off,4);
    uint32_t ih=40; memcpy(hdr+14,&ih,4);
    int32_t w=(int32_t)W, h=(int32_t)H;
    memcpy(hdr+18,&w,4); memcpy(hdr+22,&h,4);
    uint16_t planes=1, bpp=24; memcpy(hdr+26,&planes,2); memcpy(hdr+28,&bpp,2);
    memcpy(hdr+34,&img,4);
    fwrite(hdr,1,54,f);
    std::vector<uint8_t> row(rowB+pad,0);
    for (int y=H-1;y>=0;y--){
        for (uint32_t x=0;x<W;x++){
            const uint8_t* p=bgra+((size_t)y*W+x)*4;
            row[x*3+0]=p[0]; row[x*3+1]=p[1]; row[x*3+2]=p[2];
        }
        fwrite(row.data(),1,rowB+pad,f);
    }
    fclose(f);
}

// ----------------------------------------------------------------- faces
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

static int headless(int frames, bool golden, bool json, const char* shot){
    G.validation=true;
    setupCommon();
    initScene();
    // readback buffer
    VkBuffer rb; VkDeviceMemory rbMem;
    {
        VkBufferCreateInfo bi{VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO};
        bi.size=BYTES; bi.usage=VK_BUFFER_USAGE_TRANSFER_DST_BIT; bi.sharingMode=VK_SHARING_MODE_EXCLUSIVE;
        VK_CHECK(vkCreateBuffer(G.dev,&bi,nullptr,&rb));
        VkMemoryRequirements mr; vkGetBufferMemoryRequirements(G.dev,rb,&mr);
        VkMemoryAllocateInfo mai{VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO};
        mai.allocationSize=mr.size;
        mai.memoryTypeIndex=memType(mr.memoryTypeBits,VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT|VK_MEMORY_PROPERTY_HOST_COHERENT_BIT);
        VK_CHECK(vkAllocateMemory(G.dev,&mai,nullptr,&rbMem));
        VK_CHECK(vkBindBufferMemory(G.dev,rb,rbMem,0));
    }
    void* rbMap=nullptr; VK_CHECK(vkMapMemory(G.dev,rbMem,0,BYTES,0,&rbMap));
    VkCommandBuffer cb;
    {
        VkCommandBufferAllocateInfo ai{VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO};
        ai.commandPool=G.pool; ai.level=VK_COMMAND_BUFFER_LEVEL_PRIMARY; ai.commandBufferCount=1;
        VK_CHECK(vkAllocateCommandBuffers(G.dev,&ai,&cb));
        VkCommandBufferBeginInfo bi{VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO};
        VK_CHECK(vkBeginCommandBuffer(cb,&bi));
        VkBufferCopy cp{0,0,BYTES};
        vkCmdCopyBuffer(cb,G.sharedBuf,rb,1,&cp);
        VK_CHECK(vkEndCommandBuffer(cb));
    }
    std::string h1,hMid,hLast;
    float ev50=0, ev130=0, evLast=0, rangeMax=0;
    const float AZ=0.0f, EL=0.05f;                      // fixed headless camera
    for (int n=1;n<=frames;n++){
        cuWait(G.extVk,(uint64_t)(n-1));
        FrameStats fs = renderFrame((uint32_t)n, AZ, EL, 0/*AgX*/, 1/*cinematic*/, 1/*dither*/, nullptr);
        cuSignal(G.extCuda,(uint64_t)n);
        submitTimeline(cb,(uint64_t)n);
        hostWaitVkDone((uint64_t)n);
        if (fs.rangeLog2>rangeMax) rangeMax=fs.rangeLog2;
        if (n==50)  ev50=fs.evUsed;
        if (n==130) ev130=fs.evUsed;
        if (n==frames) evLast=fs.evUsed;
        const uint8_t* r=(const uint8_t*)rbMap;
        if (n==1)        h1  =blake2b::hash256_hex(r,BYTES);
        if (n==frames/2) hMid=blake2b::hash256_hex(r,BYTES);
        if (n==frames){  hLast=blake2b::hash256_hex(r,BYTES);
            if (shot) writeBMP(shot,r);
        }
        if (shot && n==130){ char p2[512]; snprintf(p2,sizeof(p2),"%s.flash.bmp",shot); writeBMP(p2,r); }
    }
    VK_CHECK(vkQueueWaitIdle(G.q));
    bool gVal=(G.vErrors==0 && G.vWarnings==0);
    bool gRange=(rangeMax>=10.0f);                       // >= 2^10 = 1024x in-frame
    bool gExpo=(fabsf(ev130-ev50)>=1.5f);                // adaptation demonstrably moved
    bool verdict=gVal&&gRange&&gExpo&&G.luidMatch;
    char buf[1024];
    snprintf(buf,sizeof(buf),
        "{\"module\":\"cinematic\",\"ver\":\"0.1.0\",\"face\":\"headless\",\"W\":%u,\"H\":%u,"
        "\"frames\":%d,\"scene\":\"supernova\",\"tone\":\"agx\",\"h1\":\"%.16s\",\"hmid\":\"%.16s\","
        "\"hlast\":\"%.16s\",\"ev50\":%.6f,\"ev130\":%.6f,\"evlast\":%.6f,\"range_log2\":%.3f,"
        "\"gates\":{\"range\":%d,\"exposure\":%d,\"device_luid\":%d,\"validation_clean\":%d},"
        "\"verdict\":%d}",
        W,H,frames,h1.c_str(),hMid.c_str(),hLast.c_str(),ev50,ev130,evLast,rangeMax,
        gRange?1:0,gExpo?1:0,G.luidMatch?1:0,gVal?1:0,verdict?1:0);
    std::string declared(buf);
    if (golden) return goldenFace(declared,"goldens/cinematic/golden.hash");
    if (json){ printf("%s\n",declared.c_str()); return verdict?0:1; }
    printf("cinematic v0.1.0 - R1 HEADLESS (supernova scene, full CINEMATIC chain)\n");
    printf("-------------------------------------------------------\n");
    printf("  device            %s (LUID %s)\n",G.devName,G.luidMatch?"PASS":"FAIL");
    printf("  frames            %d @ %ux%u  AgX, cinematic mode, dithered\n",frames,W,H);
    printf("  G-RANGE           %s   (max in-frame range 2^%.2f; gate >= 2^10)\n",gRange?"PASS":"FAIL",rangeMax);
    printf("  G-EXPOSURE        %s   (EV 50->130 moved %.2f; gate >= 1.5)\n",gExpo?"PASS":"FAIL",fabsf(ev130-ev50));
    printf("  G-VALIDATION      %s   (%llu err, %llu warn)\n",gVal?"PASS":"FAIL",
           (unsigned long long)G.vErrors,(unsigned long long)G.vWarnings);
    printf("  EV(50/130/last)   %.3f / %.3f / %.3f\n",ev50,ev130,evLast);
    printf("  hashes            %.16s %.16s %.16s\n",h1.c_str(),hMid.c_str(),hLast.c_str());
    printf("-------------------------------------------------------\n");
    printf("VERDICT: %s\n",verdict?"PASS":"FAIL");
    printf("declared hash: %s\n",blake2b::hash256_hex(declared).c_str());
    return verdict?0:1;
}

static int selftestFace(){
    // KATs, CPU-side (host __host__ paths of the same functions)
    printf("cinematic v0.1.0 --selftest (KATs)\n");
    float r,g,b; bool ok=true;
    blackbody(3200.0f,&r,&g,&b);
    bool k1=(r>g && g>b && b>=0.0f);
    printf("[KAT] blackbody(3200K)  = (%.4f, %.4f, %.4f)  r>g>b %s\n",r,g,b,k1?"PASS":"FAIL");
    blackbody(6500.0f,&r,&g,&b);
    bool k2=(fabsf(r/g-1.0f)<0.25f && fabsf(b/g-1.0f)<0.35f);
    printf("[KAT] blackbody(6500K)  = (%.4f, %.4f, %.4f)  near-white %s\n",r,g,b,k2?"PASS":"FAIL");
    blackbody(11000.0f,&r,&g,&b);
    bool k3=(b>r);
    printf("[KAT] blackbody(11000K) = (%.4f, %.4f, %.4f)  blue>red %s\n",r,g,b,k3?"PASS":"FAIL");
    // AgX monotone on a gray ramp + range
    bool k4=true; float prev=-1.0f;
    for (int i=0;i<=40;i++){
        float x=powf(10.0f,-4.0f+6.0f*(float)i/40.0f);
        float orr,og,ob; agx3(x,x,x,&orr,&og,&ob);
        if (orr<prev-1e-6f || orr<0.0f || orr>1.0f) k4=false;
        prev=orr;
    }
    printf("[KAT] AgX gray ramp monotone in [0,1]: %s\n",k4?"PASS":"FAIL");
    // sRGB round trip
    bool k5=true;
    for (int i=0;i<=20;i++){
        float c=(float)i/20.0f, e=srgb1(c);
        float d=(e<=0.04045f)? e/12.92f : powf((e+0.055f)/1.055f,2.4f);
        if (fabsf(d-c)>2e-6f) k5=false;
    }
    printf("[KAT] sRGB encode/decode roundtrip <= 2e-6: %s\n",k5?"PASS":"FAIL");
    // triangular dither statistics
    double mean=0.0, am=0.0; const int ND=1000000;
    for (int i=0;i<ND;i++){
        uint32_t h1=hash32((uint32_t)i^0x85ebca6bu), h2=hash32(h1^0xc2b2ae35u);
        double d=(((double)(h1&0xFFFF)+(double)(h2&0xFFFF))/65535.0-1.0)/255.0;
        mean+=d; am+=fabs(d);
    }
    mean/=ND; am/=ND;
    // E|tri| = 1/3 LSB; mean gate at ~4 sigma for 1e6 draws (SE ~ 1.6e-6)
    bool k6=(fabs(mean)<6.5e-6 && fabs(am-1.0/(3.0*255.0))<3e-5);
    printf("[KAT] triangular dither: mean=%.2e  E|d|=%.6f (expect %.6f)  %s\n",
           mean,am,1.0/(3.0*255.0),k6?"PASS":"FAIL");
    ok=k1&&k2&&k3&&k4&&k5&&k6;
    printf("VERDICT: %s\n",ok?"PASS":"FAIL");
    return ok?0:1;
}

static LRESULT CALLBACK wndProc(HWND h, UINT m, WPARAM w, LPARAM l){
    if (m==WM_DESTROY || (m==WM_KEYDOWN && w==VK_ESCAPE)){ PostQuitMessage(0); return 0; }
    return DefWindowProcA(h,m,w,l);
}

static int windowed(){
    G.windowed=true;
    setupCommon();
    initScene();
    HINSTANCE hi=GetModuleHandleA(nullptr);
    WNDCLASSA wc{}; wc.lpfnWndProc=wndProc; wc.hInstance=hi; wc.lpszClassName="TU_CINEMATIC";
    wc.hCursor=LoadCursorA(nullptr,IDC_ARROW);
    RegisterClassA(&wc);
    RECT rc{0,0,(LONG)W,(LONG)H};
    AdjustWindowRect(&rc,WS_OVERLAPPEDWINDOW,FALSE);
    HWND hwnd=CreateWindowA("TU_CINEMATIC","TINY UNIVERSE - R1 cinematic (SUPERNOVA)",
        WS_OVERLAPPEDWINDOW|WS_VISIBLE,CW_USEDEFAULT,CW_USEDEFAULT,
        rc.right-rc.left,rc.bottom-rc.top,nullptr,nullptr,hi,nullptr);
    VkSurfaceKHR surf;
    VkWin32SurfaceCreateInfoKHR sci{VK_STRUCTURE_TYPE_WIN32_SURFACE_CREATE_INFO_KHR};
    sci.hinstance=hi; sci.hwnd=hwnd;
    VK_CHECK(vkCreateWin32SurfaceKHR(G.inst,&sci,nullptr,&surf));
    VkBool32 sup=VK_FALSE;
    VK_CHECK(vkGetPhysicalDeviceSurfaceSupportKHR(G.phys,G.qfam,surf,&sup));
    if (!sup){ fprintf(stderr,"present unsupported\n"); return 2; }
    VkSurfaceCapabilitiesKHR caps;
    VK_CHECK(vkGetPhysicalDeviceSurfaceCapabilitiesKHR(G.phys,surf,&caps));
    VkSwapchainCreateInfoKHR sc{VK_STRUCTURE_TYPE_SWAPCHAIN_CREATE_INFO_KHR};
    sc.surface=surf; sc.minImageCount=caps.minImageCount+1;
    if (caps.maxImageCount && sc.minImageCount>caps.maxImageCount) sc.minImageCount=caps.maxImageCount;
    sc.imageFormat=VK_FORMAT_B8G8R8A8_UNORM;
    sc.imageColorSpace=VK_COLOR_SPACE_SRGB_NONLINEAR_KHR;
    sc.imageExtent={W,H}; sc.imageArrayLayers=1;
    sc.imageUsage=VK_IMAGE_USAGE_TRANSFER_DST_BIT|VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT;
    sc.imageSharingMode=VK_SHARING_MODE_EXCLUSIVE;
    sc.preTransform=caps.currentTransform;
    sc.compositeAlpha=VK_COMPOSITE_ALPHA_OPAQUE_BIT_KHR;
    sc.presentMode=VK_PRESENT_MODE_FIFO_KHR;
    sc.clipped=VK_TRUE;
    VkSwapchainKHR swap; VK_CHECK(vkCreateSwapchainKHR(G.dev,&sc,nullptr,&swap));
    uint32_t nImg=0; VK_CHECK(vkGetSwapchainImagesKHR(G.dev,swap,&nImg,nullptr));
    std::vector<VkImage> imgs(nImg);
    VK_CHECK(vkGetSwapchainImagesKHR(G.dev,swap,&nImg,imgs.data()));
    std::vector<VkCommandBuffer> cbs(nImg);
    {
        VkCommandBufferAllocateInfo ai{VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO};
        ai.commandPool=G.pool; ai.level=VK_COMMAND_BUFFER_LEVEL_PRIMARY; ai.commandBufferCount=nImg;
        VK_CHECK(vkAllocateCommandBuffers(G.dev,&ai,cbs.data()));
        for (uint32_t i=0;i<nImg;i++){
            VkCommandBufferBeginInfo bi{VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO};
            VK_CHECK(vkBeginCommandBuffer(cbs[i],&bi));
            VkImageMemoryBarrier b1{VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER};
            b1.oldLayout=VK_IMAGE_LAYOUT_UNDEFINED; b1.newLayout=VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL;
            b1.srcQueueFamilyIndex=VK_QUEUE_FAMILY_IGNORED; b1.dstQueueFamilyIndex=VK_QUEUE_FAMILY_IGNORED;
            b1.image=imgs[i]; b1.subresourceRange={VK_IMAGE_ASPECT_COLOR_BIT,0,1,0,1};
            b1.srcAccessMask=0; b1.dstAccessMask=VK_ACCESS_TRANSFER_WRITE_BIT;
            vkCmdPipelineBarrier(cbs[i],VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT,VK_PIPELINE_STAGE_TRANSFER_BIT,
                                 0,0,nullptr,0,nullptr,1,&b1);
            VkBufferImageCopy cp{};
            cp.imageSubresource={VK_IMAGE_ASPECT_COLOR_BIT,0,0,1};
            cp.imageExtent={W,H,1};
            vkCmdCopyBufferToImage(cbs[i],G.sharedBuf,imgs[i],VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,1,&cp);
            VkImageMemoryBarrier b2=b1;
            b2.oldLayout=VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL; b2.newLayout=VK_IMAGE_LAYOUT_PRESENT_SRC_KHR;
            b2.srcAccessMask=VK_ACCESS_TRANSFER_WRITE_BIT; b2.dstAccessMask=0;
            vkCmdPipelineBarrier(cbs[i],VK_PIPELINE_STAGE_TRANSFER_BIT,VK_PIPELINE_STAGE_BOTTOM_OF_PIPE_BIT,
                                 0,0,nullptr,0,nullptr,1,&b2);
            VK_CHECK(vkEndCommandBuffer(cbs[i]));
        }
    }
    VkSemaphore acq,ren;
    VkSemaphoreCreateInfo bsi{VK_STRUCTURE_TYPE_SEMAPHORE_CREATE_INFO};
    VK_CHECK(vkCreateSemaphore(G.dev,&bsi,nullptr,&acq));
    VK_CHECK(vkCreateSemaphore(G.dev,&bsi,nullptr,&ren));
    // damped inertial orbit (§5): spring on az/el toward mouse-set targets
    float az=0.0f, el=0.05f, vaz=0.0f, vel=0.0f, tAz=0.0f, tEl=0.05f;
    int tone=0, cinematic=1, dither=1;
    bool dragging=false; POINT lastMouse{};
    uint64_t n=0; MSG msg{};
    ULONGLONG t0=GetTickCount64(); uint64_t f0=0;
    for (;;){
        while (PeekMessageA(&msg,nullptr,0,0,PM_REMOVE)){
            if (msg.message==WM_QUIT) goto done;
            if (msg.message==WM_KEYDOWN){
                if (msg.wParam=='T') tone=(tone+1)%3;
                else if (msg.wParam=='C') cinematic^=1;
                else if (msg.wParam=='D') dither^=1;
            }
            if (msg.message==WM_LBUTTONDOWN){ dragging=true; GetCursorPos(&lastMouse); }
            if (msg.message==WM_LBUTTONUP) dragging=false;
            TranslateMessage(&msg); DispatchMessageA(&msg);
        }
        if (dragging){
            POINT mp; GetCursorPos(&mp);
            tAz += 0.005f*(float)(mp.x-lastMouse.x);
            tEl += 0.005f*(float)(mp.y-lastMouse.y);
            tEl = clampf(tEl,-1.35f,1.35f);
            lastMouse=mp;
        }
        // critically-damped spring, k=10 s^-1 at ~60fps -> per-frame
        const float dt=1.0f/60.0f, k=10.0f;
        vaz += (k*k*(tAz-az) - 2.0f*k*vaz)*dt; az += vaz*dt;
        vel += (k*k*(tEl-el) - 2.0f*k*vel)*dt; el += vel*dt;
        n++;
        cuWait(G.extVk,n-1);
        char hud[64];
        FrameStats fs; fs.evUsed=P.ev;
        {
            const char* tn=(tone==0)?"AGX":((tone==1)?"ACES":"CS");
            snprintf(hud,sizeof(hud),"%s E%d %s",tn,(int)lroundf(P.ev),(cinematic?"C":"P"));
        }
        fs=renderFrame((uint32_t)n,az,el,tone,cinematic,dither,hud);
        cuSignal(G.extCuda,n);
        uint32_t idx=0;
        VkResult ar=vkAcquireNextImageKHR(G.dev,swap,~0ull,acq,VK_NULL_HANDLE,&idx);
        if (ar!=VK_SUCCESS && ar!=VK_SUBOPTIMAL_KHR){ fprintf(stderr,"acquire %d\n",(int)ar); break; }
        submitTimeline(cbs[idx],n,acq,ren);
        VkPresentInfoKHR pi{VK_STRUCTURE_TYPE_PRESENT_INFO_KHR};
        pi.waitSemaphoreCount=1; pi.pWaitSemaphores=&ren;
        pi.swapchainCount=1; pi.pSwapchains=&swap; pi.pImageIndices=&idx;
        VkResult pr=vkQueuePresentKHR(G.q,&pi);
        if (pr!=VK_SUCCESS && pr!=VK_SUBOPTIMAL_KHR){ fprintf(stderr,"present %d\n",(int)pr); break; }
        ULONGLONG t=GetTickCount64();
        if (t-t0>=1000){
            char title[160];
            uint32_t cyc=(uint32_t)(n%1200u);
            const char* phase=(cyc<60)?"quiescent":((cyc<150)?"RISING":"decay");
            snprintf(title,sizeof(title),
                "TINY UNIVERSE - R1 cinematic | SUPERNOVA %s | %s | EV %.1f | %.0f fps  [T tone, C cin/phy, D dither, drag orbit]",
                phase,(tone==0?"AgX":(tone==1?"ACES":"RAW-CLAMP")),P.ev,
                1000.0*(double)(n-f0)/(double)(t-t0));
            SetWindowTextA(hwnd,title);
            t0=t; f0=n;
        }
    }
done:
    VK_CHECK(vkQueueWaitIdle(G.q));
    printf("[windowed] clean exit after %llu frames\n",(unsigned long long)n);
    return 0;
}

int main(int argc, char** argv){
    bool json=false, golden=false, self=false, head=false;
    int frames=240;
    const char* shot=nullptr;
    for (int i=1;i<argc;i++){
        std::string a=argv[i];
        if      (a=="--json") json=true;
        else if (a=="--golden") golden=true;
        else if (a=="--selftest") self=true;
        else if (a=="--headless") head=true;
        else if (a=="--frames" && i+1<argc) frames=atoi(argv[++i]);
        else if (a=="--shot" && i+1<argc){ shot=argv[++i]; head=true; }
        else { fprintf(stderr,"usage: cinematic [--headless [--frames N] [--golden|--json]] [--shot out.bmp] [--selftest]\n"); return 2; }
    }
    if (self) return selftestFace();
    if (head) return headless(frames,golden,json,shot);
    return windowed();
}
