#ifndef GPU_MANAGER_H
#define GPU_MANAGER_H

#include <cstdio>
#include <vector>
#include <cuda.h>
#include <cuda_runtime_api.h>

#include "useful.h"

#ifdef USE_NCCL
#include <nccl.h>
#define NCCLCHECK(cmd) do {					\
      ncclResult_t r = cmd;					\
      if (r!= ncclSuccess) {					\
	  printf("Failed, NCCL error %s:%d '%s'\n",             \
		 __FILE__,__LINE__,ncclGetErrorString(r));	\
	  exit(EXIT_FAILURE);					\
      }								\
  } while(0)
#endif

#ifndef gpuErrchk
#define delgpuErrchk
#define gpuErrchk(code) { if ((code) != cudaSuccess) {					       \
	fprintf(stderr,"CUDA Error: %s %s %d\n", cudaGetErrorString(code), __FILE__, __LINE__); \
    }}
#endif

#define NUMSTREAMS 8

// GPUs capable of Peer Access 
// (Currently unused)
struct GPUPeer {
	int gpu;
	std::vector<int> gpus;
	GPUPeer() : gpu(-1) {}
	GPUPeer(int gpu) : gpu(gpu) {}
};

class GPU {
    /* Class to represent individual GPUs on a node */
    friend class GPUManager;
private:
    bool may_timeout;
    unsigned int id;
    cudaStream_t streams[NUMSTREAMS];

    int last_stream;
    bool streams_created;
    void create_streams();
    void destroy_streams();

    cudaDeviceProp properties;

public:
    GPU(unsigned int id);
    ~GPU();

    inline const cudaStream_t& get_stream(unsigned int stream_id) {
	return streams[stream_id];
    }

    inline const cudaStream_t& get_next_stream() {
	if (last_stream == NUMSTREAMS-1) {
	    last_stream = 0;
	} else {
            last_stream +=1;
	}
	return streams[last_stream];
    };
};

class GPUManager {

private:
	static std::vector<GPU> allGpus, timeouts, notimeouts;
	static void init_devices();
	static int nGPUs;
	static bool is_safe;

	#ifdef USE_NCCL
	static void init_comms();
	static ncclComm_t* comms;
	#endif

public:	
	static size_t allGpuSize() { return allGpus.size(); }
	static std::vector<GPU> gpus;
	
	static bool safe() { return is_safe; }

	// init
	// Initializes gpus and properties vector
	// Bad things may happen if this is called more than once
	static void init();

	static void load_info();

	static void select_gpus(std::vector<unsigned int>& gpu_ids);
	// use
	// Use the GPU using local index 0..N (not cudaGetDevice index)
	static void use(int gpu_id);

	static void sync(int gpu_id);
	static void sync() {
	    if (gpus.size() > 1) {
		int curr;
		gpuErrchk( cudaGetDevice(&curr) );
		for (auto it = gpus.begin(); it != gpus.end(); ++it) {
		    gpuErrchk( cudaSetDevice(it->id) );
		    gpuErrchk( cudaDeviceSynchronize() );
		}
		gpuErrchk( cudaSetDevice(curr) );
	    } else gpuErrchk( cudaDeviceSynchronize() );
	}


	// current
	// @return the current GPU a thread is using
	static int current();
	
	// safe
	// @param whether gpus should contain GPUs that may timeout
	static void safe(bool make_safe);
	
	static int getInitialGPU();

        // 
    inline const cudaStream_t& get_next_stream() {
	return gpus[0].get_next_stream();
    };

    #ifdef USE_NCCL
    template<typename T>
    void nccl_broadcast(int root, std::vector<T*> send_d, std::vector<T*> recv_d, unsigned int size, int stream_id) {
	if (gpus.size() == 1) return;
	cudaStream_t stream = 0;
	NCCLCHECK(ncclGroupStart());
	for (size_t i = 0; i < gpus.size(); ++i) {
	    if (stream_id >= 0) stream = gpus[i].streams[stream_id];
	    NCCLCHECK( ncclBroadcast((const void*) send_d[i], (void*) recv_d[i],
				     size*sizeof(T)/sizeof(float), ncclFloat, root,
				     comms[i], stream) );
	}
	NCCLCHECK(ncclGroupEnd());
    }
    template<typename T>
	void nccl_broadcast(int root, std::vector<T*> send_d, std::vector<T*> recv_d, unsigned int size, cudaStream_t* streams) {
	if (gpus.size() == 1) return;
	NCCLCHECK(ncclGroupStart());
	for (size_t i = 0; i < gpus.size(); ++i) {
	    NCCLCHECK( ncclBroadcast((const void*) send_d[i], (void*) recv_d[i],
				     size*sizeof(T)/sizeof(float), ncclFloat, root,
				     comms[i], streams[i]) );
	}
	NCCLCHECK(ncclGroupEnd());
    }

    template<typename T>
    void nccl_reduce(int root, const std::vector<T*> send_d, const std::vector<T*> recv_d, const unsigned int size, const int stream_id) {
	if (gpus.size() == 1) return;
	cudaStream_t stream = 0;
	NCCLCHECK(ncclGroupStart());
	for (size_t i = 0; i < gpus.size(); ++i) {
	    if (stream_id >= 0) stream = gpus[i].streams[stream_id];
	    NCCLCHECK(ncclReduce((const void*) send_d[i], (void*) recv_d[i],
				 size*sizeof(T)/sizeof(float), ncclFloat, ncclSum, root, 
				 comms[i], stream));
	}
	NCCLCHECK(ncclGroupEnd());
    }
    #endif // USE_NCCL
    
};
#ifndef delgpuErrchk
#undef  delgpuErrchk
#undef  gpuErrchk
#endif

#endif
