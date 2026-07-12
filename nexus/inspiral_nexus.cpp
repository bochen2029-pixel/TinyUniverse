// ============================================================================
//  inspiral_nexus — TINY UNIVERSE v1 polish: 2.5PN gravitational-wave inspiral
//  Contract: contracts/inspiral.contract.md v1.0.0
//
//  Single-file CPU fp64 (no GPU). The Peters (1964) orbit-averaged 2.5PN radiation
//  reaction: a binary radiates gravitational waves, so its orbit SHRINKS and
//  CIRCULARIZES, spiralling to merger (the physics behind every LIGO chirp).
//    da/dt = -(64/5)(G^3/c^5)(m1 m2 M/a^3)(1-e^2)^-7/2 (1 + 73/24 e^2 + 37/96 e^4)
//    de/dt = -(304/15)(G^3/c^5)(m1 m2 M/a^4) e (1-e^2)^-5/2 (1 + 121/304 e^2)
//  Oracles (exact GR): circular merger time T_c = (5/256) c^5 a0^4/(G^3 m1 m2 M),
//  and the circularization curve a(e) = c0 e^(12/19)/(1-e^2) [1+121/304 e^2]^(870/2299).
//
//  Envelope face (--scenario X --json|--golden|--selftest). Exit 0/1/2.
//  Build: cl /std:c++17 /EHsc /O2 /W4 nexus\inspiral_nexus.cpp /Fe:build\inspiral_nexus.exe
//  Dials: c=20, G=2e-3, m1=m2=1e4 (r_s=2GM/c^2=0.2, a0=10 -> a/r_s=50, weak-field 2.5PN).
// ============================================================================

#define _CRT_SECURE_NO_WARNINGS
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <string>

struct Dials { double c = 20.0; double G = 2.0e-3; double m1 = 1.0e4; double m2 = 1.0e4; };

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
    for (int i = 0; i < 16; i++){ uint64_t w=0; for (int j=7;j>=0;j--) w=(w<<8)|block[i*8+j]; m[i]=w; }
    for (int r = 0; r < 12; r++){
        const uint8_t* s = SIGMA[r % 10];
        Gm(v,0,4, 8,12,m[s[ 0]],m[s[ 1]]); Gm(v,1,5, 9,13,m[s[ 2]],m[s[ 3]]);
        Gm(v,2,6,10,14,m[s[ 4]],m[s[ 5]]); Gm(v,3,7,11,15,m[s[ 6]],m[s[ 7]]);
        Gm(v,0,5,10,15,m[s[ 8]],m[s[ 9]]); Gm(v,1,6,11,12,m[s[10]],m[s[11]]);
        Gm(v,2,7, 8,13,m[s[12]],m[s[13]]); Gm(v,3,4, 9,14,m[s[14]],m[s[15]]);
    }
    for (int i = 0; i < 8; i++) h[i] ^= v[i] ^ v[i + 8];
}
static std::string hash256_hex(const std::string& in){
    uint64_t h[8]; for (int i=0;i<8;i++) h[i]=IV[i];
    h[0] ^= 0x01010000ull ^ 32ull;
    size_t n=in.size(), off=0; uint8_t block[128];
    while (n-off > 128){ memcpy(block,in.data()+off,128); off+=128; compress(h,block,(uint64_t)off,false); }
    size_t rem=n-off; memset(block,0,128); if (rem) memcpy(block,in.data()+off,rem);
    compress(h,block,(uint64_t)n,true);
    char hex[65]; for (int i=0;i<32;i++){ unsigned b=(unsigned)((h[i/8]>>(8*(i%8)))&0xFF); snprintf(hex+2*i,3,"%02x",b); }
    return std::string(hex,64);
}
} // namespace blake2b

static std::string fmt6(double x){ if(!std::isfinite(x)) return "9999999.999999"; if(x==0.0)x=0.0; char b[48]; snprintf(b,sizeof(b),"%.6f",x); return std::string(b); }
static std::string fmt9(double x){ if(!std::isfinite(x)) return "9.999999999e+99"; if(x==0.0)x=0.0; char b[48]; snprintf(b,sizeof(b),"%.9e",x); return std::string(b); }

// ============================================================================
//  PETERS (1964) INTEGRATORS (fp64 RK4)
// ============================================================================

// Peters circularization function g(e): a(e) = c0 * g(e).
static double petersG(double e){
    return std::pow(e, 12.0/19.0) / (1.0 - e*e)
         * std::pow(1.0 + (121.0/304.0)*e*e, 870.0/2299.0);
}

// Circular (e=0) inspiral: integrate a(t) via da/dt = -(64/5) K / a^3 from a0 until
// a <= a_end (K = G^3 m1 m2 M / c^5). Returns the merger time (interp at the crossing).
static double circularMergerTime(double K, double a0, double a_end, long steps){
    double a=a0, t=0.0;
    double T_c = (5.0/256.0)*std::pow(a0,4.0)/K;      // full merger time (for step sizing)
    double dt = T_c/(double)steps;
    auto dadt=[&](double a){ return -(64.0/5.0)*K/(a*a*a); };
    for (long i=0;i<steps*2;i++){
        double aprev=a, tprev=t;
        double k1=dadt(a), k2=dadt(a+0.5*dt*k1), k3=dadt(a+0.5*dt*k2), k4=dadt(a+dt*k3);
        a += dt/6.0*(k1+2*k2+2*k3+k4); t += dt;
        if (a <= a_end){ double frac=(aprev-a_end)/(aprev-a); return tprev + frac*dt; }
    }
    return t;
}

// Eccentric: integrate a(e) via da/de = (12/19)(a/e)(1+73/24 e^2+37/96 e^4)/
// [(1-e^2)(1+121/304 e^2)] from e0 down to the checkpoints; record a at each.
static void eccentricCurve(double a0, double e0, const double* echk, int nchk, double* aout, long steps){
    double a=a0, e=e0, emin=echk[nchk-1];
    double de=(emin-e0)/(double)steps;               // negative: e decreases
    auto dade=[&](double a,double e){
        double num=1.0+(73.0/24.0)*e*e+(37.0/96.0)*e*e*e*e;
        double den=(1.0-e*e)*(1.0+(121.0/304.0)*e*e);
        return (12.0/19.0)*(a/e)*num/den;
    };
    int ci=0;
    for (long i=0;i<steps;i++){
        double eprev=e, aprev=a;
        double k1=dade(a,e);
        double k2=dade(a+0.5*de*k1, e+0.5*de);
        double k3=dade(a+0.5*de*k2, e+0.5*de);
        double k4=dade(a+de*k3, e+de);
        a += de/6.0*(k1+2*k2+2*k3+k4); e += de;
        while (ci<nchk && e <= echk[ci]){
            double frac=(eprev-echk[ci])/(eprev-e);
            aout[ci]=aprev+frac*(a-aprev); ci++;
        }
    }
    while (ci<nchk) aout[ci++]=a;
}

// ============================================================================
//  Result + JSON
// ============================================================================
struct Result {
    std::string scenario;
    double c=0,G=0,m1=0,m2=0,a0=0,e0=0,a_end=0; long steps=0;
    // circular
    double T_merge=0, T_target=0, T_c=0, T_rel=0;
    // eccentric
    double echk[5]={0}, a_num[5]={0}, a_pet[5]={0}; double a_rel_max=0;
    bool nan_free=true, gate_primary=false, gate_nan=false, verdict=false;
};

static std::string declaredJson(const Dials& D, uint64_t seed, const Result& R){
    std::string s; s.reserve(1536);
    double M=D.m1+D.m2, r_s=2.0*D.G*M/(D.c*D.c);
    s += "{\"tool\":\"inspiral_nexus\",\"version\":\"1.0.0\",\"seed\":" + std::to_string(seed);
    s += ",\"params\":{\"scenario\":\"" + R.scenario + "\"";
    s += ",\"c\":" + fmt6(D.c) + ",\"G\":" + fmt9(D.G) + ",\"m1\":" + fmt6(D.m1) + ",\"m2\":" + fmt6(D.m2);
    s += ",\"r_s\":" + fmt6(r_s) + ",\"a0\":" + fmt6(R.a0) + ",\"e0\":" + fmt6(R.e0);
    s += ",\"steps\":" + std::to_string(R.steps);
    if (R.scenario == "circular"){
        s += ",\"a_end\":" + fmt6(R.a_end) + "},\"result\":{";
        s += "\"T_merge\":" + fmt9(R.T_merge) + ",\"T_target\":" + fmt9(R.T_target);
        s += ",\"T_c_full\":" + fmt9(R.T_c) + ",\"T_rel\":" + fmt9(R.T_rel);
    } else {
        s += "},\"result\":{";
        s += "\"a_num\":[" + fmt6(R.a_num[0]);
        for (int i=1;i<5;i++) s += "," + fmt6(R.a_num[i]);
        s += "],\"a_peters\":[" + fmt6(R.a_pet[0]);
        for (int i=1;i<5;i++) s += "," + fmt6(R.a_pet[i]);
        s += "],\"e_checkpoints\":[" + fmt6(R.echk[0]);
        for (int i=1;i<5;i++) s += "," + fmt6(R.echk[i]);
        s += "],\"a_rel_max\":" + fmt9(R.a_rel_max);
    }
    s += ",\"nan_free\":" + std::string(R.nan_free?"1":"0");
    s += "},\"gates\":{\"primary\":" + std::string(R.gate_primary?"true":"false");
    s += ",\"nan\":" + std::string(R.gate_nan?"true":"false");
    s += "},\"verdict\":\"" + std::string(R.verdict?"pass":"fail") + "\"";
    return s;
}

// ============================================================================
//  SCENARIOS
// ============================================================================
static Result runCircular(const Dials& D){
    Result R; R.scenario="circular";
    R.c=D.c;R.G=D.G;R.m1=D.m1;R.m2=D.m2;R.a0=10.0;R.e0=0.0;R.a_end=1.0;R.steps=1000000;
    double M=D.m1+D.m2, K=std::pow(D.G,3.0)*D.m1*D.m2*M/std::pow(D.c,5.0);
    R.T_merge = circularMergerTime(K, R.a0, R.a_end, R.steps);
    R.T_c = (5.0/256.0)*std::pow(R.a0,4.0)/K;                           // full T_c (to a=0)
    R.T_target = (5.0/256.0)*(std::pow(R.a0,4.0)-std::pow(R.a_end,4.0))/K;  // a0 -> a_end
    R.T_rel = std::fabs(R.T_merge/R.T_target - 1.0);
    R.nan_free = std::isfinite(R.T_merge);
    R.gate_primary = (R.T_rel < 1e-3);
    R.gate_nan = R.nan_free;
    R.verdict = R.gate_primary && R.gate_nan;
    return R;
}
static Result runEccentric(const Dials& D){
    Result R; R.scenario="eccentric";
    R.c=D.c;R.G=D.G;R.m1=D.m1;R.m2=D.m2;R.a0=10.0;R.e0=0.7;R.steps=1000000;
    const double echk[5]={0.6,0.5,0.4,0.3,0.2};
    for (int i=0;i<5;i++) R.echk[i]=echk[i];
    eccentricCurve(R.a0, R.e0, echk, 5, R.a_num, R.steps);
    double g0=petersG(R.e0), maxrel=0;
    for (int i=0;i<5;i++){
        R.a_pet[i] = R.a0 * petersG(echk[i]) / g0;
        double rel=std::fabs(R.a_num[i]/R.a_pet[i]-1.0);
        if (rel>maxrel) maxrel=rel;
    }
    R.a_rel_max=maxrel;
    R.nan_free = std::isfinite(maxrel);
    R.gate_primary = (maxrel < 1e-3);
    R.gate_nan = R.nan_free;
    R.verdict = R.gate_primary && R.gate_nan;
    return R;
}

static void printHuman(const Result& R){
    double M=R.m1+R.m2, r_s=2.0*R.G*M/(R.c*R.c);
    printf("inspiral_nexus v1.0.0 - TINY UNIVERSE v1 polish: 2.5PN GW inspiral (Peters 1964)\n");
    printf("scenario: %s   m1=%.0f m2=%.0f  r_s=%.4f  a0=%.1f (a0/r_s=%.0f)\n",
           R.scenario.c_str(), R.m1, R.m2, r_s, R.a0, R.a0/r_s);
    printf("-------------------------------------------------------\n");
    if (R.scenario == "circular"){
        printf("  circular inspiral -> merger (gravitational radiation shrinks the orbit):\n");
        printf("    T_merge (a0->a_end=%.1f)  %.6e s\n", R.a_end, R.T_merge);
        printf("    exact GR T_target        %.6e s  (Peters T_c full = %.6e)\n", R.T_target, R.T_c);
        printf("    rel error %.3e   [gate < 1e-3]  %s\n", R.T_rel, R.gate_primary?"PASS":"FAIL");
    } else {
        printf("  eccentric inspiral circularizes along the Peters a(e) curve:\n");
        printf("      e      a_numerical    a_Peters      rel\n");
        for (int i=0;i<5;i++)
            printf("    %.2f   %.6f    %.6f   %.2e\n", R.echk[i], R.a_num[i], R.a_pet[i],
                   std::fabs(R.a_num[i]/R.a_pet[i]-1.0));
        printf("    max rel error %.3e   [gate < 1e-3]  %s\n", R.a_rel_max, R.gate_primary?"PASS":"FAIL");
        printf("    (a shrinks 10.00 -> %.2f as e: 0.70 -> 0.20 -- the orbit circularizes)\n", R.a_num[4]);
    }
    printf("  nan_free      %s\n", R.nan_free?"yes":"NO");
    printf("-------------------------------------------------------\n");
    printf("VERDICT: %s\n", R.verdict?"PASS":"FAIL");
}

int main(int argc, char** argv){
    Dials D; uint64_t seed=20260711ull;
    std::string scenario="circular";
    bool json=false, selftest=false, golden=false;
    for (int i=1;i<argc;i++){
        std::string a=argv[i];
        if      (a=="--json") json=true;
        else if (a=="--selftest") selftest=true;
        else if (a=="--golden") golden=true;
        else if (a=="--scenario" && i+1<argc) scenario=argv[++i];
        else if (a=="--seed" && i+1<argc) seed=strtoull(argv[++i],nullptr,10);
        else { fprintf(stderr,"usage: inspiral_nexus --scenario circular|eccentric [--seed N] [--json] [--golden] [--selftest]\n"); return 2; }
    }

    if (selftest){
        // K=0 (radiation off): the circular orbit does NOT inspiral -- a stays a0.
        double K=0.0, a=10.0, dt=1.0;
        auto dadt=[&](double aa){ return -(64.0/5.0)*K/(aa*aa*aa); };
        for (int i=0;i<100000;i++){
            double k1=dadt(a),k2=dadt(a+0.5*dt*k1),k3=dadt(a+0.5*dt*k2),k4=dadt(a+dt*k3);
            a += dt/6.0*(k1+2*k2+2*k3+k4);
        }
        bool ok1 = (std::fabs(a-10.0) < 1e-12);
        bool ok2 = (petersG(0.6) > petersG(0.3)) && (petersG(0.3) > petersG(0.1));  // g(e) monotone -> circularization
        printf("[selftest] K=0 (no radiation): circular a unchanged  a=%.9f (expect 10)  [%s]\n", a, ok1?"PASS":"FAIL");
        printf("[selftest] Peters g(e) monotone (orbit circularizes as e drops)  [%s]\n", ok2?"PASS":"FAIL");
        printf("VERDICT: %s\n", (ok1&&ok2)?"PASS":"FAIL");
        return (ok1&&ok2)?0:1;
    }

    if (golden){ D = Dials(); seed=20260711ull; }

    Result R = (scenario=="circular") ? runCircular(D) : (scenario=="eccentric") ? runEccentric(D)
             : (fprintf(stderr,"error: unknown scenario '%s' (circular|eccentric)\n", scenario.c_str()), std::exit(2), Result());
    std::string declared = declaredJson(D, seed, R);
    std::string hash = blake2b::hash256_hex(declared);

    if (golden){
        std::string path = "goldens/inspiral_" + R.scenario + "/golden.hash";
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
