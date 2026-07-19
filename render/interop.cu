// render/interop.cu — R0: the Vulkan⇄CUDA presentation rung (contract:
// contracts/interop.contract.md v0.1.0 APPROVED; D-002 Path A, D-012 P1, D-034 SDK pin).
// CUDA writes a shared device-local VkBuffer (exported OPAQUE_WIN32, imported via
// cudaImportExternalMemory); Vulkan presents it by vkCmdCopyBufferToImage; ordering by
// two exported TIMELINE semaphores (cudaDone / vkDone). Zero CPU pixel round-trips in
// the steady loop; no vkQueueWaitIdle. Device selection: LUID match vs the CUDA device
// (the Intel iGPU also enumerates on this machine — D-034 measured; index-0 is a defect).
// Faces: windowed (default, "the sim breathes"), --headless --frames N (the gated
// declared face; validation layers ON), --selftest (capability probe).
// Build (BUILD.md): nvcc -O3 -arch=sm_89 -o build\interop.exe render\interop.cu
//                   -I"%VULKAN_SDK%\Include" "%VULKAN_SDK%\Lib\vulkan-1.lib" user32.lib
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#define VK_USE_PLATFORM_WIN32_KHR
#include <vulkan/vulkan.h>
#include <cuda_runtime.h>
#include <cstdio>
#include <cstring>
#include <cstdint>
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
} // namespace blake2b

// ----------------------------------------------------------------- checks + config
#define CU_CHECK(x) do{ cudaError_t e=(x); if(e!=cudaSuccess){ \
    fprintf(stderr,"CUDA ERROR %s:%d %s\n",__FILE__,__LINE__,cudaGetErrorString(e)); exit(2);} }while(0)
#define VK_CHECK(x) do{ VkResult r_=(x); if(r_!=VK_SUCCESS){ \
    fprintf(stderr,"VK ERROR %s:%d code=%d\n",__FILE__,__LINE__,(int)r_); exit(2);} }while(0)

static const uint32_t W = 1280, H = 720;
static const size_t   BYTES = (size_t)W*H*4;

// ----------------------------------------------------------------- the CUDA pattern
// Deterministic, tick-parameterized (NO wall-clock): three drifting interference waves
// -> RGBA8. Pixel (0,0)..(0,3) carries the frame counter bytes (the G-SYNC probe).
__global__ void kPattern(uchar4* buf, uint32_t n){
    uint32_t x = blockIdx.x*blockDim.x + threadIdx.x;
    uint32_t y = blockIdx.y*blockDim.y + threadIdx.y;
    if (x >= W || y >= H) return;
    float fx = (float)x/(float)W, fy = (float)y/(float)H, t = 0.02f*(float)n;
    float a = 0.5f + 0.5f*sinf(12.0f*fx + 9.0f*fy + t);
    float b = 0.5f + 0.5f*sinf(17.0f*fx - 6.0f*fy - 1.7f*t + 2.1f);
    float c = 0.5f + 0.5f*sinf(7.0f*(fx-0.5f)*(fx-0.5f)*24.0f + 21.0f*fy*fy + 0.9f*t);
    float rr = a*0.8f + c*0.2f, gg = 0.25f*a + 0.75f*b, bb = 0.6f*b + 0.4f*c;
    buf[y*W + x] = make_uchar4((unsigned char)(255.0f*rr*rr),
                               (unsigned char)(255.0f*gg),
                               (unsigned char)(255.0f*sqrtf(bb)), 255);
    if (y == 0 && x < 4){
        unsigned char byte = (unsigned char)((n >> (8*x)) & 0xFF);
        buf[x] = make_uchar4(byte, byte, byte, 255);
    }
}

// ----------------------------------------------------------------- Vulkan context
struct Ctx {
    bool validation = false, windowed = false;
    VkInstance inst = VK_NULL_HANDLE;
    VkDebugUtilsMessengerEXT dbg = VK_NULL_HANDLE;
    VkPhysicalDevice phys = VK_NULL_HANDLE;
    VkDevice dev = VK_NULL_HANDLE;
    uint32_t qfam = 0;
    VkQueue q = VK_NULL_HANDLE;
    VkCommandPool pool = VK_NULL_HANDLE;
    char devName[256] = {0};
    uint8_t vkLuid[8] = {0};
    int cudaDev = 0;
    uint8_t cuLuid[8] = {0};
    bool luidMatch = false;
    // shared buffer
    VkBuffer sharedBuf = VK_NULL_HANDLE;
    VkDeviceMemory sharedMem = VK_NULL_HANDLE;
    cudaExternalMemory_t extMem = nullptr;
    uchar4* dPix = nullptr;                 // CUDA view of the shared buffer
    // timelines
    VkSemaphore tlCuda = VK_NULL_HANDLE, tlVk = VK_NULL_HANDLE;   // cudaDone / vkDone
    cudaExternalSemaphore_t extCuda = nullptr, extVk = nullptr;
    // validation counters
    uint64_t vErrors = 0, vWarnings = 0;
};
static Ctx G;

static VKAPI_ATTR VkBool32 VKAPI_CALL dbgCb(VkDebugUtilsMessageSeverityFlagBitsEXT sev,
        VkDebugUtilsMessageTypeFlagsEXT, const VkDebugUtilsMessengerCallbackDataEXT* data, void*){
    if (sev & VK_DEBUG_UTILS_MESSAGE_SEVERITY_ERROR_BIT_EXT){ G.vErrors++;   fprintf(stderr,"[VK-ERROR] %s\n", data->pMessage); }
    else if (sev & VK_DEBUG_UTILS_MESSAGE_SEVERITY_WARNING_BIT_EXT){ G.vWarnings++; fprintf(stderr,"[VK-WARN ] %s\n", data->pMessage); }
    return VK_FALSE;
}

static void createInstance(){
    VkApplicationInfo ai{VK_STRUCTURE_TYPE_APPLICATION_INFO};
    ai.pApplicationName = "tinyuniverse-interop-r0";
    ai.apiVersion = VK_API_VERSION_1_3;
    std::vector<const char*> exts, layers;
    if (G.windowed){ exts.push_back("VK_KHR_surface"); exts.push_back("VK_KHR_win32_surface"); }
    if (G.validation){ exts.push_back(VK_EXT_DEBUG_UTILS_EXTENSION_NAME); layers.push_back("VK_LAYER_KHRONOS_validation"); }
    VkInstanceCreateInfo ci{VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO};
    ci.pApplicationInfo = &ai;
    ci.enabledExtensionCount = (uint32_t)exts.size(); ci.ppEnabledExtensionNames = exts.data();
    ci.enabledLayerCount = (uint32_t)layers.size();   ci.ppEnabledLayerNames = layers.data();
    VK_CHECK(vkCreateInstance(&ci, nullptr, &G.inst));
    if (G.validation){
        auto f = (PFN_vkCreateDebugUtilsMessengerEXT)vkGetInstanceProcAddr(G.inst, "vkCreateDebugUtilsMessengerEXT");
        VkDebugUtilsMessengerCreateInfoEXT di{VK_STRUCTURE_TYPE_DEBUG_UTILS_MESSENGER_CREATE_INFO_EXT};
        di.messageSeverity = VK_DEBUG_UTILS_MESSAGE_SEVERITY_WARNING_BIT_EXT | VK_DEBUG_UTILS_MESSAGE_SEVERITY_ERROR_BIT_EXT;
        di.messageType = VK_DEBUG_UTILS_MESSAGE_TYPE_GENERAL_BIT_EXT | VK_DEBUG_UTILS_MESSAGE_TYPE_VALIDATION_BIT_EXT |
                         VK_DEBUG_UTILS_MESSAGE_TYPE_PERFORMANCE_BIT_EXT;
        di.pfnUserCallback = dbgCb;
        if (f) VK_CHECK(f(G.inst, &di, nullptr, &G.dbg));
    }
}

static void pickDeviceByLuid(){
    // CUDA side first
    CU_CHECK(cudaSetDevice(0));
    cudaDeviceProp p{}; CU_CHECK(cudaGetDeviceProperties(&p, 0));
    memcpy(G.cuLuid, p.luid, 8);
    // Vulkan side: match LUID (the iGPU also enumerates — D-034)
    uint32_t n = 0; VK_CHECK(vkEnumeratePhysicalDevices(G.inst, &n, nullptr));
    std::vector<VkPhysicalDevice> devs(n);
    VK_CHECK(vkEnumeratePhysicalDevices(G.inst, &n, devs.data()));
    for (auto d : devs){
        VkPhysicalDeviceIDProperties idp{VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_ID_PROPERTIES};
        VkPhysicalDeviceProperties2 p2{VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_PROPERTIES_2};
        p2.pNext = &idp;
        vkGetPhysicalDeviceProperties2(d, &p2);
        if (idp.deviceLUIDValid && memcmp(idp.deviceLUID, G.cuLuid, 8) == 0){
            G.phys = d;
            memcpy(G.vkLuid, idp.deviceLUID, 8);
            snprintf(G.devName, sizeof(G.devName), "%s", p2.properties.deviceName);
            G.luidMatch = true;
            return;
        }
    }
    fprintf(stderr, "G-DEVICE FAIL: no Vulkan physical device matches the CUDA device LUID\n");
    exit(2);
}

static void createDevice(){
    // queue family: graphics (implies transfer)
    uint32_t n=0; vkGetPhysicalDeviceQueueFamilyProperties(G.phys, &n, nullptr);
    std::vector<VkQueueFamilyProperties> qf(n);
    vkGetPhysicalDeviceQueueFamilyProperties(G.phys, &n, qf.data());
    bool found=false;
    for (uint32_t i=0;i<n;i++) if (qf[i].queueFlags & VK_QUEUE_GRAPHICS_BIT){ G.qfam=i; found=true; break; }
    if (!found){ fprintf(stderr,"no graphics queue family\n"); exit(2); }
    float prio=1.0f;
    VkDeviceQueueCreateInfo qi{VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO};
    qi.queueFamilyIndex=G.qfam; qi.queueCount=1; qi.pQueuePriorities=&prio;
    std::vector<const char*> exts = {
        VK_KHR_EXTERNAL_MEMORY_EXTENSION_NAME,
        VK_KHR_EXTERNAL_MEMORY_WIN32_EXTENSION_NAME,
        VK_KHR_EXTERNAL_SEMAPHORE_EXTENSION_NAME,
        VK_KHR_EXTERNAL_SEMAPHORE_WIN32_EXTENSION_NAME,
    };
    if (G.windowed) exts.push_back(VK_KHR_SWAPCHAIN_EXTENSION_NAME);
    VkPhysicalDeviceVulkan12Features f12{VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_VULKAN_1_2_FEATURES};
    f12.timelineSemaphore = VK_TRUE;
    VkDeviceCreateInfo ci{VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO};
    ci.pNext = &f12;
    ci.queueCreateInfoCount=1; ci.pQueueCreateInfos=&qi;
    ci.enabledExtensionCount=(uint32_t)exts.size(); ci.ppEnabledExtensionNames=exts.data();
    VK_CHECK(vkCreateDevice(G.phys, &ci, nullptr, &G.dev));
    vkGetDeviceQueue(G.dev, G.qfam, 0, &G.q);
    VkCommandPoolCreateInfo pi{VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO};
    pi.queueFamilyIndex=G.qfam; pi.flags=VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT;
    VK_CHECK(vkCreateCommandPool(G.dev, &pi, nullptr, &G.pool));
}

static uint32_t memType(uint32_t bits, VkMemoryPropertyFlags props){
    VkPhysicalDeviceMemoryProperties mp;
    vkGetPhysicalDeviceMemoryProperties(G.phys, &mp);
    for (uint32_t i=0;i<mp.memoryTypeCount;i++)
        if ((bits & (1u<<i)) && (mp.memoryTypes[i].propertyFlags & props) == props) return i;
    fprintf(stderr,"no suitable memory type\n"); exit(2);
}

static void createSharedBuffer(){
    VkExternalMemoryBufferCreateInfo ext{VK_STRUCTURE_TYPE_EXTERNAL_MEMORY_BUFFER_CREATE_INFO};
    ext.handleTypes = VK_EXTERNAL_MEMORY_HANDLE_TYPE_OPAQUE_WIN32_BIT;
    VkBufferCreateInfo bi{VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO};
    bi.pNext=&ext; bi.size=BYTES;
    bi.usage=VK_BUFFER_USAGE_TRANSFER_SRC_BIT|VK_BUFFER_USAGE_TRANSFER_DST_BIT;
    bi.sharingMode=VK_SHARING_MODE_EXCLUSIVE;
    VK_CHECK(vkCreateBuffer(G.dev, &bi, nullptr, &G.sharedBuf));
    VkMemoryRequirements mr; vkGetBufferMemoryRequirements(G.dev, G.sharedBuf, &mr);
    VkExportMemoryAllocateInfo exp{VK_STRUCTURE_TYPE_EXPORT_MEMORY_ALLOCATE_INFO};
    exp.handleTypes = VK_EXTERNAL_MEMORY_HANDLE_TYPE_OPAQUE_WIN32_BIT;
    VkMemoryAllocateInfo mai{VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO};
    mai.pNext=&exp; mai.allocationSize=mr.size;
    mai.memoryTypeIndex=memType(mr.memoryTypeBits, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
    VK_CHECK(vkAllocateMemory(G.dev, &mai, nullptr, &G.sharedMem));
    VK_CHECK(vkBindBufferMemory(G.dev, G.sharedBuf, G.sharedMem, 0));
    // export -> CUDA import
    auto getH = (PFN_vkGetMemoryWin32HandleKHR)vkGetDeviceProcAddr(G.dev, "vkGetMemoryWin32HandleKHR");
    if (!getH){ fprintf(stderr,"vkGetMemoryWin32HandleKHR missing\n"); exit(2); }
    VkMemoryGetWin32HandleInfoKHR gi{VK_STRUCTURE_TYPE_MEMORY_GET_WIN32_HANDLE_INFO_KHR};
    gi.memory=G.sharedMem; gi.handleType=VK_EXTERNAL_MEMORY_HANDLE_TYPE_OPAQUE_WIN32_BIT;
    HANDLE h=nullptr; VK_CHECK(getH(G.dev, &gi, &h));
    cudaExternalMemoryHandleDesc hd{};
    hd.type = cudaExternalMemoryHandleTypeOpaqueWin32;
    hd.handle.win32.handle = h;
    hd.size = mr.size;
    CU_CHECK(cudaImportExternalMemory(&G.extMem, &hd));
    CloseHandle(h);                                     // CUDA sample pattern: safe post-import
    cudaExternalMemoryBufferDesc bd{};
    bd.offset=0; bd.size=BYTES;
    void* dp=nullptr;
    CU_CHECK(cudaExternalMemoryGetMappedBuffer(&dp, G.extMem, &bd));
    G.dPix = (uchar4*)dp;
}

static VkSemaphore createTimeline(cudaExternalSemaphore_t* cuOut){
    VkExportSemaphoreCreateInfo exp{VK_STRUCTURE_TYPE_EXPORT_SEMAPHORE_CREATE_INFO};
    exp.handleTypes = VK_EXTERNAL_SEMAPHORE_HANDLE_TYPE_OPAQUE_WIN32_BIT;
    VkSemaphoreTypeCreateInfo ty{VK_STRUCTURE_TYPE_SEMAPHORE_TYPE_CREATE_INFO};
    ty.pNext=&exp; ty.semaphoreType=VK_SEMAPHORE_TYPE_TIMELINE; ty.initialValue=0;
    VkSemaphoreCreateInfo ci{VK_STRUCTURE_TYPE_SEMAPHORE_CREATE_INFO};
    ci.pNext=&ty;
    VkSemaphore s; VK_CHECK(vkCreateSemaphore(G.dev, &ci, nullptr, &s));
    auto getH = (PFN_vkGetSemaphoreWin32HandleKHR)vkGetDeviceProcAddr(G.dev, "vkGetSemaphoreWin32HandleKHR");
    if (!getH){ fprintf(stderr,"vkGetSemaphoreWin32HandleKHR missing\n"); exit(2); }
    VkSemaphoreGetWin32HandleInfoKHR gi{VK_STRUCTURE_TYPE_SEMAPHORE_GET_WIN32_HANDLE_INFO_KHR};
    gi.semaphore=s; gi.handleType=VK_EXTERNAL_SEMAPHORE_HANDLE_TYPE_OPAQUE_WIN32_BIT;
    HANDLE h=nullptr; VK_CHECK(getH(G.dev, &gi, &h));
    cudaExternalSemaphoreHandleDesc sd{};
    sd.type = cudaExternalSemaphoreHandleTypeTimelineSemaphoreWin32;
    sd.handle.win32.handle = h;
    CU_CHECK(cudaImportExternalSemaphore(cuOut, &sd));
    CloseHandle(h);
    return s;
}

static void cuWait(cudaExternalSemaphore_t s, uint64_t v){
    cudaExternalSemaphoreWaitParams p{}; p.params.fence.value=v;
    CU_CHECK(cudaWaitExternalSemaphoresAsync(&s, &p, 1, 0));
}
static void cuSignal(cudaExternalSemaphore_t s, uint64_t v){
    cudaExternalSemaphoreSignalParams p{}; p.params.fence.value=v;
    CU_CHECK(cudaSignalExternalSemaphoresAsync(&s, &p, 1, 0));
}

static void setupCommon(){
    createInstance();
    pickDeviceByLuid();
    createDevice();
    createSharedBuffer();
    G.tlCuda = createTimeline(&G.extCuda);
    G.tlVk   = createTimeline(&G.extVk);
}

// submit ONE prerecorded command buffer with timeline wait(cudaDone=v) / signal(vkDone=v)
// + optional binary wait/signal (windowed present path)
static void submitTimeline(VkCommandBuffer cb, uint64_t v,
                           VkSemaphore extraWaitBin = VK_NULL_HANDLE,
                           VkSemaphore extraSigBin = VK_NULL_HANDLE){
    uint64_t waitVals[2]  = { v, 0 };
    uint64_t sigVals[2]   = { v, 0 };
    VkSemaphore waits[2]  = { G.tlCuda, extraWaitBin };
    VkSemaphore sigs[2]   = { G.tlVk,   extraSigBin };
    VkPipelineStageFlags stages[2] = { VK_PIPELINE_STAGE_TRANSFER_BIT, VK_PIPELINE_STAGE_TRANSFER_BIT };
    uint32_t nw = extraWaitBin ? 2u : 1u, ns = extraSigBin ? 2u : 1u;
    VkTimelineSemaphoreSubmitInfo tsi{VK_STRUCTURE_TYPE_TIMELINE_SEMAPHORE_SUBMIT_INFO};
    tsi.waitSemaphoreValueCount = nw; tsi.pWaitSemaphoreValues = waitVals;
    tsi.signalSemaphoreValueCount = ns; tsi.pSignalSemaphoreValues = sigVals;
    VkSubmitInfo si{VK_STRUCTURE_TYPE_SUBMIT_INFO};
    si.pNext = &tsi;
    si.waitSemaphoreCount = nw; si.pWaitSemaphores = waits; si.pWaitDstStageMask = stages;
    si.commandBufferCount = 1; si.pCommandBuffers = &cb;
    si.signalSemaphoreCount = ns; si.pSignalSemaphores = sigs;
    VK_CHECK(vkQueueSubmit(G.q, 1, &si, VK_NULL_HANDLE));
}

static void hostWaitVkDone(uint64_t v){
    VkSemaphoreWaitInfo wi{VK_STRUCTURE_TYPE_SEMAPHORE_WAIT_INFO};
    wi.semaphoreCount=1; wi.pSemaphores=&G.tlVk; wi.pValues=&v;
    VK_CHECK(vkWaitSemaphores(G.dev, &wi, ~0ull));
}

// ----------------------------------------------------------------- headless (the declared face)
static int goldenFace(const std::string& declared, const char* path){
    std::string hash = blake2b::hash256_hex(declared);
    FILE* f = fopen(path, "rb");
    if (!f){ fprintf(stderr,"GOLDEN NOT FROZEN %s\n",hash.c_str()); printf("%s\n",hash.c_str()); return 2; }
    char want[128]={0};
    if (fscanf(f,"%127s",want)!=1){ fclose(f); fprintf(stderr,"GOLDEN FILE UNREADABLE\n"); return 2; }
    fclose(f);
    if (hash==std::string(want)){ fprintf(stderr,"GOLDEN OK %.8s\n",hash.c_str()); printf("%s\n",hash.c_str()); return 0; }
    fprintf(stderr,"GOLDEN MISMATCH have=%.8s want=%.8s\n",hash.c_str(),want);
    printf("%s\n",hash.c_str()); return 1;
}

static int headless(int frames, bool golden, bool json){
    G.validation = true;
    setupCommon();
    // readback buffer (host-visible)
    VkBuffer rb; VkDeviceMemory rbMem;
    {
        VkBufferCreateInfo bi{VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO};
        bi.size=BYTES; bi.usage=VK_BUFFER_USAGE_TRANSFER_DST_BIT; bi.sharingMode=VK_SHARING_MODE_EXCLUSIVE;
        VK_CHECK(vkCreateBuffer(G.dev,&bi,nullptr,&rb));
        VkMemoryRequirements mr; vkGetBufferMemoryRequirements(G.dev, rb, &mr);
        VkMemoryAllocateInfo mai{VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO};
        mai.allocationSize=mr.size;
        mai.memoryTypeIndex=memType(mr.memoryTypeBits,
            VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT|VK_MEMORY_PROPERTY_HOST_COHERENT_BIT);
        VK_CHECK(vkAllocateMemory(G.dev,&mai,nullptr,&rbMem));
        VK_CHECK(vkBindBufferMemory(G.dev,rb,rbMem,0));
    }
    void* rbMap=nullptr; VK_CHECK(vkMapMemory(G.dev,rbMem,0,BYTES,0,&rbMap));
    // one prerecorded command buffer: shared -> readback
    VkCommandBuffer cb;
    {
        VkCommandBufferAllocateInfo ai{VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO};
        ai.commandPool=G.pool; ai.level=VK_COMMAND_BUFFER_LEVEL_PRIMARY; ai.commandBufferCount=1;
        VK_CHECK(vkAllocateCommandBuffers(G.dev,&ai,&cb));
        VkCommandBufferBeginInfo bi{VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO};
        VK_CHECK(vkBeginCommandBuffer(cb,&bi));
        VkBufferCopy cp{0,0,BYTES};
        vkCmdCopyBuffer(cb, G.sharedBuf, rb, 1, &cp);
        VK_CHECK(vkEndCommandBuffer(cb));
    }
    dim3 bs(16,16), gs((W+15)/16,(H+15)/16);
    std::vector<uint8_t> host(BYTES);
    std::string h0, hMid, hLast;
    bool gSync = true, gRound = true;
    uint32_t prevCounter = 0;
    for (int n=1; n<=frames; n++){
        cuWait(G.extVk, (uint64_t)(n-1));               // buffer free (frame n-1 consumed)
        kPattern<<<gs,bs>>>(G.dPix, (uint32_t)n);
        CU_CHECK(cudaGetLastError());
        cuSignal(G.extCuda, (uint64_t)n);
        submitTimeline(cb, (uint64_t)n);
        hostWaitVkDone((uint64_t)n);
        const uint8_t* r = (const uint8_t*)rbMap;
        uint32_t counter = (uint32_t)r[0] | ((uint32_t)r[4])<<8 | ((uint32_t)r[8])<<16 | ((uint32_t)r[12])<<24;
        if (counter != (uint32_t)n || (n>1 && counter != prevCounter+1)) gSync = false;
        prevCounter = counter;
        if (n==1 || n==frames){
            // G-ROUNDTRIP: CUDA's direct view of the shared bytes vs the Vulkan-copied bytes
            CU_CHECK(cudaMemcpy(host.data(), G.dPix, BYTES, cudaMemcpyDeviceToHost));
            if (memcmp(host.data(), rbMap, BYTES) != 0) gRound = false;
        }
        if (n==1)        h0   = blake2b::hash256_hex(r, BYTES);
        if (n==frames/2) hMid = blake2b::hash256_hex(r, BYTES);
        if (n==frames)   hLast= blake2b::hash256_hex(r, BYTES);
    }
    VK_CHECK(vkQueueWaitIdle(G.q));                      // teardown only (allowed: not the steady loop)
    bool gVal = (G.vErrors==0 && G.vWarnings==0);
    bool verdict = gSync && gRound && G.luidMatch && gVal;
    char buf[1024];
    snprintf(buf,sizeof(buf),
        "{\"module\":\"interop\",\"ver\":\"0.1.0\",\"face\":\"headless\",\"W\":%u,\"H\":%u,"
        "\"frames\":%d,\"path\":\"buffer-copy\",\"h0\":\"%.16s\",\"hmid\":\"%.16s\",\"hlast\":\"%.16s\","
        "\"gates\":{\"roundtrip\":%d,\"sync\":%d,\"device_luid\":%d,\"validation_clean\":%d},"
        "\"verdict\":%d}",
        W,H,frames,h0.c_str(),hMid.c_str(),hLast.c_str(),
        gRound?1:0,gSync?1:0,G.luidMatch?1:0,gVal?1:0,verdict?1:0);
    std::string declared(buf);
    if (golden) return goldenFace(declared, "goldens/interop/golden.hash");
    if (json){ printf("%s\n", declared.c_str()); return verdict?0:1; }
    printf("interop v0.1.0 - R0 HEADLESS (Vulkan-CUDA external memory + timeline semaphores)\n");
    printf("-------------------------------------------------------\n");
    printf("  device                 %s (LUID match: %s)\n", G.devName, G.luidMatch?"PASS":"FAIL");
    printf("  frames                 %d @ %ux%u  (path: buffer-copy, declared)\n", frames, W, H);
    printf("  G-ROUNDTRIP            %s   (CUDA view == Vulkan-copied bytes, frames 1+%d)\n", gRound?"PASS":"FAIL", frames);
    printf("  G-SYNC                 %s   (strictly monotonic embedded counters)\n", gSync?"PASS":"FAIL");
    printf("  G-VALIDATION           %s   (%llu errors, %llu warnings)\n", gVal?"PASS":"FAIL",
           (unsigned long long)G.vErrors,(unsigned long long)G.vWarnings);
    printf("  frame hashes           h0=%.16s  mid=%.16s  last=%.16s\n", h0.c_str(),hMid.c_str(),hLast.c_str());
    printf("-------------------------------------------------------\n");
    printf("VERDICT: %s\n", verdict?"PASS":"FAIL");
    printf("declared hash: %s\n", blake2b::hash256_hex(declared).c_str());
    return verdict?0:1;
}

// ----------------------------------------------------------------- selftest
static int selftest(){
    G.validation = true;
    createInstance();
    pickDeviceByLuid();
    createDevice();
    printf("[selftest] device: %s\n", G.devName);
    printf("[selftest] LUID match (CUDA dev 0): PASS\n");
    printf("[selftest] timeline semaphores + external memory/semaphore win32: device created OK\n");
    bool gVal = (G.vErrors==0 && G.vWarnings==0);
    printf("[selftest] validation clean: %s\n", gVal?"PASS":"FAIL");
    printf("VERDICT: %s\n", gVal?"PASS":"FAIL");
    return gVal?0:1;
}

// ----------------------------------------------------------------- windowed ("the sim breathes")
static LRESULT CALLBACK wndProc(HWND h, UINT m, WPARAM w, LPARAM l){
    if (m==WM_DESTROY || (m==WM_KEYDOWN && w==VK_ESCAPE)){ PostQuitMessage(0); return 0; }
    return DefWindowProcA(h,m,w,l);
}

static int windowed(){
    G.windowed = true;
    setupCommon();
    // window
    HINSTANCE hi = GetModuleHandleA(nullptr);
    WNDCLASSA wc{}; wc.lpfnWndProc=wndProc; wc.hInstance=hi; wc.lpszClassName="TU_INTEROP";
    wc.hCursor=LoadCursorA(nullptr,IDC_ARROW);
    RegisterClassA(&wc);
    RECT rc{0,0,(LONG)W,(LONG)H};
    AdjustWindowRect(&rc, WS_OVERLAPPEDWINDOW, FALSE);
    HWND hwnd = CreateWindowA("TU_INTEROP","TINY UNIVERSE - R0 interop (Vulkan-CUDA)",
        WS_OVERLAPPEDWINDOW|WS_VISIBLE, CW_USEDEFAULT, CW_USEDEFAULT,
        rc.right-rc.left, rc.bottom-rc.top, nullptr, nullptr, hi, nullptr);
    // surface + swapchain
    VkSurfaceKHR surf;
    VkWin32SurfaceCreateInfoKHR sci{VK_STRUCTURE_TYPE_WIN32_SURFACE_CREATE_INFO_KHR};
    sci.hinstance=hi; sci.hwnd=hwnd;
    VK_CHECK(vkCreateWin32SurfaceKHR(G.inst,&sci,nullptr,&surf));
    VkBool32 sup=VK_FALSE;
    VK_CHECK(vkGetPhysicalDeviceSurfaceSupportKHR(G.phys, G.qfam, surf, &sup));
    if (!sup){ fprintf(stderr,"present not supported on the graphics queue\n"); return 2; }
    VkSurfaceCapabilitiesKHR caps;
    VK_CHECK(vkGetPhysicalDeviceSurfaceCapabilitiesKHR(G.phys,surf,&caps));
    VkSwapchainCreateInfoKHR sc{VK_STRUCTURE_TYPE_SWAPCHAIN_CREATE_INFO_KHR};
    sc.surface=surf;
    sc.minImageCount = caps.minImageCount+1;
    if (caps.maxImageCount && sc.minImageCount>caps.maxImageCount) sc.minImageCount=caps.maxImageCount;
    sc.imageFormat=VK_FORMAT_B8G8R8A8_UNORM;                 // Q-R0-1 resolution
    sc.imageColorSpace=VK_COLOR_SPACE_SRGB_NONLINEAR_KHR;
    sc.imageExtent={W,H};
    sc.imageArrayLayers=1;
    sc.imageUsage=VK_IMAGE_USAGE_TRANSFER_DST_BIT|VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT;
    sc.imageSharingMode=VK_SHARING_MODE_EXCLUSIVE;
    sc.preTransform=caps.currentTransform;
    sc.compositeAlpha=VK_COMPOSITE_ALPHA_OPAQUE_BIT_KHR;
    sc.presentMode=VK_PRESENT_MODE_FIFO_KHR;                 // vsync (pacing non-declared)
    sc.clipped=VK_TRUE;
    VkSwapchainKHR swap; VK_CHECK(vkCreateSwapchainKHR(G.dev,&sc,nullptr,&swap));
    uint32_t nImg=0; VK_CHECK(vkGetSwapchainImagesKHR(G.dev,swap,&nImg,nullptr));
    std::vector<VkImage> imgs(nImg);
    VK_CHECK(vkGetSwapchainImagesKHR(G.dev,swap,&nImg,imgs.data()));
    // prerecorded cmdbuf per swapchain image: UNDEFINED->DST, copy (RGBA source into BGRA
    // image: R/B swap is visually acceptable for the R0 pattern; exact swizzle is R1 polish),
    // ->PRESENT
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
            b1.image=imgs[i];
            b1.subresourceRange={VK_IMAGE_ASPECT_COLOR_BIT,0,1,0,1};
            b1.srcAccessMask=0; b1.dstAccessMask=VK_ACCESS_TRANSFER_WRITE_BIT;
            vkCmdPipelineBarrier(cbs[i],VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT,VK_PIPELINE_STAGE_TRANSFER_BIT,
                                 0,0,nullptr,0,nullptr,1,&b1);
            VkBufferImageCopy cp{};
            cp.bufferOffset=0; cp.bufferRowLength=0; cp.bufferImageHeight=0;
            cp.imageSubresource={VK_IMAGE_ASPECT_COLOR_BIT,0,0,1};
            cp.imageOffset={0,0,0}; cp.imageExtent={W,H,1};
            vkCmdCopyBufferToImage(cbs[i], G.sharedBuf, imgs[i], VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL, 1, &cp);
            VkImageMemoryBarrier b2=b1;
            b2.oldLayout=VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL; b2.newLayout=VK_IMAGE_LAYOUT_PRESENT_SRC_KHR;
            b2.srcAccessMask=VK_ACCESS_TRANSFER_WRITE_BIT; b2.dstAccessMask=0;
            vkCmdPipelineBarrier(cbs[i],VK_PIPELINE_STAGE_TRANSFER_BIT,VK_PIPELINE_STAGE_BOTTOM_OF_PIPE_BIT,
                                 0,0,nullptr,0,nullptr,1,&b2);
            VK_CHECK(vkEndCommandBuffer(cbs[i]));
        }
    }
    // binary semaphores for acquire/present
    VkSemaphore acq, ren;
    VkSemaphoreCreateInfo bi2{VK_STRUCTURE_TYPE_SEMAPHORE_CREATE_INFO};
    VK_CHECK(vkCreateSemaphore(G.dev,&bi2,nullptr,&acq));
    VK_CHECK(vkCreateSemaphore(G.dev,&bi2,nullptr,&ren));
    dim3 bs(16,16), gs((W+15)/16,(H+15)/16);
    uint64_t n=0; MSG msg{};
    ULONGLONG t0=GetTickCount64(); uint64_t f0=0;
    for (;;){
        while (PeekMessageA(&msg,nullptr,0,0,PM_REMOVE)){
            if (msg.message==WM_QUIT) goto done;
            TranslateMessage(&msg); DispatchMessageA(&msg);
        }
        n++;
        cuWait(G.extVk, n-1);
        kPattern<<<gs,bs>>>(G.dPix, (uint32_t)n);
        CU_CHECK(cudaGetLastError());
        cuSignal(G.extCuda, n);
        uint32_t idx=0;
        VkResult ar = vkAcquireNextImageKHR(G.dev, swap, ~0ull, acq, VK_NULL_HANDLE, &idx);
        if (ar!=VK_SUCCESS && ar!=VK_SUBOPTIMAL_KHR){ fprintf(stderr,"acquire failed %d\n",(int)ar); break; }
        submitTimeline(cbs[idx], n, acq, ren);
        VkPresentInfoKHR pi{VK_STRUCTURE_TYPE_PRESENT_INFO_KHR};
        pi.waitSemaphoreCount=1; pi.pWaitSemaphores=&ren;
        pi.swapchainCount=1; pi.pSwapchains=&swap; pi.pImageIndices=&idx;
        VkResult pr = vkQueuePresentKHR(G.q,&pi);
        if (pr!=VK_SUCCESS && pr!=VK_SUBOPTIMAL_KHR){ fprintf(stderr,"present failed %d\n",(int)pr); break; }
        ULONGLONG t=GetTickCount64();
        if (t-t0>=1000){
            char title[128];
            snprintf(title,sizeof(title),"TINY UNIVERSE - R0 interop  |  frame %llu  |  %.0f fps",
                     (unsigned long long)n, 1000.0*(double)(n-f0)/(double)(t-t0));
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
    for (int i=1;i<argc;i++){
        std::string a=argv[i];
        if      (a=="--json") json=true;
        else if (a=="--golden") golden=true;
        else if (a=="--selftest") self=true;
        else if (a=="--headless") head=true;
        else if (a=="--frames" && i+1<argc) frames=atoi(argv[++i]);
        else { fprintf(stderr,"usage: interop [--headless [--frames N] [--golden|--json]] [--selftest]\n"); return 2; }
    }
    if (self) return selftest();
    if (head) return headless(frames, golden, json);
    return windowed();
}
