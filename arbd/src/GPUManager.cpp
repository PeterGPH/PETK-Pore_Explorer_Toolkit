#include "GPUManager.h"

#ifndef gpuErrchk
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char* file, int line, bool abort=true) {
   if (code != cudaSuccess) {
      fprintf(stderr,"CUDA Error: %s %s:%d\n", cudaGetErrorString(code), __FILE__, line);
      if (abort) exit(code);
   }
}
#endif 

#define WITH_GPU(id,code) { int wg_curr; cudaGetDevice(&wg_curr); cudaSetDevice(id); code ; cudaSetDevice(wg_curr); }

int GPUManager::nGPUs = 0;
bool GPUManager::is_safe = true;
std::vector<GPU> GPUManager::allGpus, GPUManager::gpus, GPUManager::notimeouts;

GPU::GPU(unsigned int id) : id(id) {
    cudaSetDevice(id);
    cudaGetDeviceProperties(&properties, id);
    printf("[%d] %s ", id, properties.name);
    if (properties.kernelExecTimeoutEnabled) {
	printf("(may timeout) ");
	may_timeout = true;
    } else {
	may_timeout = false;
    }
    printf("| SM %d.%d, ", properties.major, properties.minor);
    printf("%.2fGHz, ", (float) properties.clockRate * 10E-7);
    printf("%.1fGB RAM\n", (float) properties.totalGlobalMem * 7.45058e-10);

    streams_created = false;
    // fflush(stdout);
    // gpuErrchk( cudaDeviceSynchronize() );
}
GPU::~GPU() {
    destroy_streams();
}

void GPU::create_streams() {
    int curr;
    gpuErrchk( cudaGetDevice(&curr) );
    gpuErrchk( cudaSetDevice(id) );

    if (streams_created) destroy_streams();
    last_stream = -1;
    for (int i = 0; i < NUMSTREAMS; i++) {
	// printf("  creating stream %d at %p\n", i, (void *) &streams[i]);
	gpuErrchk( cudaStreamCreate( &streams[i] ) );
	// gpuErrchk( cudaStreamCreateWithFlags( &(streams[i]) , cudaStreamNonBlocking ) );
    }
    streams_created = true;

    gpuErrchk( cudaSetDevice(id) );
    cudaSetDevice(curr);
}

void GPU::destroy_streams() {
    int curr;
    // printf("Destroying streams\n");
    if (cudaGetDevice(&curr) == cudaSuccess) { // Avoid errors when program is shutting down
	gpuErrchk( cudaSetDevice(id) );
	if (streams_created) {
	    for (int i = 0; i < NUMSTREAMS; i++) {
		// printf("  destroying stream %d at %p\n", i, (void *) &streams[i]);
		gpuErrchk( cudaStreamDestroy( streams[i] ) );
	    }
	}
	gpuErrchk( cudaSetDevice(curr) );
    }
    streams_created = false;
}


void GPUManager::init() {
    gpuErrchk(cudaGetDeviceCount(&nGPUs));
    printf("Found %d GPU(s)\n", nGPUs);
    for (int dev = 0; dev < nGPUs; dev++) {
	GPU g(dev);
	allGpus.push_back(g);
	if (!g.may_timeout) notimeouts.push_back(g);
    }
    is_safe = false;
    if (allGpus.size() == 0) {
	fprintf(stderr, "Error: Did not find a GPU\n");
	exit(1);
    }
}

void GPUManager::load_info() {
    init();
    gpus = allGpus;
    init_devices();
}

void GPUManager::init_devices() {
    printf("Initializing devices... ");
    for (unsigned int i = 0; i < gpus.size(); i++) {
    	if (i != gpus.size() - 1 && gpus.size() > 1)
    	    printf("%d, ", gpus[i].id);
    	else if (gpus.size() > 1)
    	    printf("and %d\n", gpus[i].id);
    	else
    	    printf("%d\n", gpus[i].id);

    	use(i);
    	cudaDeviceSetCacheConfig( cudaFuncCachePreferL1 );
    	gpus[i].create_streams();
    }
    use(0);
    gpuErrchk( cudaDeviceSynchronize() );
}

void GPUManager::select_gpus(std::vector<unsigned int>& gpu_ids) {
    gpus.clear();
    for (auto it = gpu_ids.begin(); it != gpu_ids.end(); ++it) {
	gpus.push_back( allGpus[*it] );
    }
    init_devices();
    #ifdef USE_NCCL
    init_comms();
    #endif
}

void GPUManager::use(int gpu_id) {
	gpu_id = gpu_id % (int) gpus.size();
	// printf("Setting device to %d\n",gpus[gpu_id].id);
	gpuErrchk( cudaSetDevice(gpus[gpu_id].id) );
	// printf("Done setting device\n");
}

void GPUManager::sync(int gpu_id) {
    WITH_GPU( gpus[gpu_id].id, 
    	      gpuErrchk( cudaDeviceSynchronize() ));
    // int wg_curr; 
    // gpuErrchk( cudaGetDevice(&wg_curr) );
    // gpuErrchk( cudaSetDevice(gpus[gpu_id].id) );
    // gpuErrchk( cudaSetDevice(wg_curr) );
}

int GPUManager::current() {
	int c;
	cudaGetDevice(&c);
	return c;
}

void GPUManager::safe(bool make_safe) {
	if (make_safe == is_safe) return;
	if (make_safe) {
		if (notimeouts.size() == 0) {
			printf("WARNING: No safe GPUs\n");
			return;
		}
		allGpus = notimeouts;
		is_safe = true;
	} else {
	    is_safe = false;
	}
}

int GPUManager::getInitialGPU() {
    // TODO: check the load on the gpus and select an unused one
    for (auto it = gpus.begin(); it != gpus.end(); ++it) {
	GPU& gpu = *it;
	if (!gpu.properties.kernelExecTimeoutEnabled)
	    return gpu.id; 
    }
    return 0;
}

#ifdef USE_NCCL
ncclComm_t* GPUManager::comms = NULL;
void GPUManager::init_comms() {
    if (gpus.size() == 1) return;
    int* gpu_ids = new int[gpus.size()];

    comms = new ncclComm_t[gpus.size()];
    int i = 0;
    for (auto &g: gpus) {
	gpu_ids[i] = g.id;
	++i;
    }
    NCCLCHECK(ncclCommInitAll(comms, gpus.size(), gpu_ids));
}
#endif
