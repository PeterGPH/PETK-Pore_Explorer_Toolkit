#include "CudaUtil.cuh"
#include <cuda_runtime_api.h>

#ifndef CUDART_VERSION
#error CUDART_VERSION Undefined!
#elif (CUDART_VERSION < 9000)

#if __CUDA_ARCH__ < 300
volatile extern __shared__ int sh[];
__device__ int warp_bcast(int v, int leader) {
	// WARNING: might not be safe to call in divergent branches 
	const int tid = threadIdx.x;
	const int warpLane = tid % WARPSIZE;
	if (warpLane == leader)
		sh[tid/WARPSIZE] = v;
	return sh[tid/WARPSIZE];		
}	
#elif __CUDA_ARCH__ < 700
__device__ int warp_bcast(int v, int leader) {return __shfl(v, leader); }
#else
__device__ int warp_bcast(int v, int leader) {return __shfl_sync(v, leader); }
#endif

__device__ int atomicAggInc(int *ctr, int warpLane) {
	// https://devblogs.nvidia.com/parallelforall/cuda-pro-tip-optimized-filtering-warp-aggregated-atomics/
	int mask = __ballot(1);
	int leader = __ffs(mask)-1;

	int res;
	if ( warpLane == leader )
		res = atomicAdd(ctr, __popc(mask));
	res = warp_bcast(res,leader);
	return res + __popc( mask & ((1 << warpLane) - 1) );
}
#endif

__global__
void reduceVector(const int num, Vector3* __restrict__ vector, Vector3* netVector) {
	extern __shared__ Vector3 blockVector[];
	const int tid = threadIdx.x;

	// grid-stride loop over vector[]
	for (int i = threadIdx.x+blockIdx.x*blockDim.x; i < num; i+=blockDim.x*gridDim.x) {
		// assign vector to shared memory
		blockVector[tid] = vector[i];
		// blockVector[tid] = Vector3(0.0f);
		__syncthreads();
		
		
		// Reduce vectors in shared memory
		// http://www.cuvilib.com/Reduction.pdf
		for (int offset = blockDim.x/2; offset > 0; offset >>= 1) {
			if (tid < offset) {
				int oid = tid + offset;
				blockVector[tid] = blockVector[tid] + blockVector[oid];
			}
			__syncthreads();
		}

		if (tid == 0)
			atomicAdd( netVector, blockVector[0] );
	}
}
