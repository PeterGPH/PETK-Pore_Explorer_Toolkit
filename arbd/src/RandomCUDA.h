#ifndef RANDOM_CUDA_H
#define RANDOM_CUDA_H

// #include "/usr/include/linux/cuda.h"
// #include "/usr/local/encap/cuda-4.0/include/cuda_runtime.h"
// #include "/usr/local/encap/cuda-4.0/include/curand_kernel.h"
// #include "/usr/local/encap/cuda-4.0/include/curand.h"
#include <cuda.h>
#include <curand_kernel.h>
#include <curand.h>

#include "useful.h"
#include "ComputeForce.h"

#ifdef __CUDACC__
    #define HOST __host__
    #define DEVICE __device__
#else
    #define HOST 
    #define DEVICE 
#endif

#ifndef gpuErrchk
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true) {
	if (code != cudaSuccess) {
		fprintf(stderr,"CUDA Error: %s (%s:%d)\n", cudaGetErrorString(code), file, line);
		if (abort) exit(code);
	}
}
#endif

#define cuRandchk(ans) { cuRandAssert((ans), __FILE__, __LINE__); }
inline void cuRandAssert(curandStatus code, const char *file, int line, bool abort=true) {
	if (code != CURAND_STATUS_SUCCESS) {
		fprintf(stderr, "CURAND Error: %d (%s:%d)\n", code, file, line);
		if (abort) exit(code);
	}
}

class Random {
public:
	static const size_t RAND_N = 1024*4; // max random numbers stored

	curandState_t *states;
	curandGenerator_t generator;
	unsigned int *integer_h, *integer_d;
	float *uniform_h, *uniform_d;
	float *gaussian_h, *gaussian_d;
	size_t integer_n, uniform_n, gaussian_n;

public:

	Random() : states(NULL), integer_h(NULL), integer_d(NULL), uniform_h(NULL), uniform_d(NULL), gaussian_h(NULL), gaussian_d(NULL) { }
	Random(int num, unsigned long seed=0) : states(NULL), integer_h(NULL), integer_d(NULL), uniform_h(NULL), uniform_d(NULL), gaussian_h(NULL), gaussian_d(NULL) {		
		init(num, seed);
	}

	void init(int num, unsigned long seed);

	DEVICE inline float gaussian(int idx, int num) {
		// TODO do stuff
		if (idx < num)
			return curand_normal(&states[idx]);
		return 0.0f;
	}
	DEVICE inline float gaussian(curandState* state) {
		return curand_normal(state);
	}

	DEVICE inline Vector3 gaussian_vector(int idx, int num) {
		// TODO do stuff
		if (idx < num) {
			curandState localState = states[idx];
			Vector3 v = gaussian_vector(&localState);
			states[idx] = localState;
			return v;
		} else return Vector3(0.0f);			
	}
	DEVICE inline Vector3 gaussian_vector(curandState* state) {
		float x = gaussian(state);
		float y = gaussian(state);
		float z = gaussian(state);
		return Vector3(x, y, z);
	}

	unsigned int integer();
	unsigned int poisson(float lambda);
	float uniform();

	HOST inline float gaussian() {
	    if (gaussian_n < 1) {
		cuRandchk(curandGenerateNormal(generator, gaussian_d, RAND_N, 0, 1));
		gpuErrchk(cudaMemcpy(gaussian_h, gaussian_d, sizeof(float) * RAND_N, cudaMemcpyDeviceToHost));
	    }
	    return gaussian_h[--gaussian_n];
	}

	void reorder(int *a, int n);
};

#endif
