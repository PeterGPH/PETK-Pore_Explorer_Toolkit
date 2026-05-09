/* Adapted from: https://developer.nvidia.com/blog/cuda-pro-tip-generate-custom-application-profile-timelines-nvtx/ */
#ifdef USE_NVTX
#include <nvToolsExt.h>

const uint32_t nvtx_colors[] = { 0xff00ff00, 0xff0000ff, 0xffffff00, 0xffff00ff, 0xff00ffff, 0xffff0000, 0xffffffff };
const int num_nvtx_colors = sizeof(nvtx_colors)/sizeof(uint32_t);

#define PUSH_NVTX(name,cid) { \
    int color_id = cid; \
    color_id = color_id%num_nvtx_colors;\
    nvtxEventAttributes_t eventAttrib = {0}; \
    eventAttrib.version = NVTX_VERSION; \
    eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE; \
    eventAttrib.colorType = NVTX_COLOR_ARGB; \
    eventAttrib.color = nvtx_colors[color_id]; \
    eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII; \
    eventAttrib.message.ascii = name; \
    nvtxRangePushEx(&eventAttrib); \
}
#define POP_NVTX nvtxRangePop();
#else
#define PUSH_NVTX(name,cid)
#define POP_NVTX
#endif
