#pragma once
#include "useful.h"
#define WARPSIZE 32

extern __device__ int warp_bcast(int v, int leader);

#ifndef CUDART_VERSION
#error CUDART_VERSION Undefined!
#elif (CUDART_VERSION < 9000)
extern __device__ int atomicAggInc(int *ctr, int warpLane);
#else
__device__ inline int atomicAggInc(int *ctr, int warpLane) {
    return atomicAdd(ctr, 1);
}
#endif

extern __global__
void reduceVector(const int num, Vector3* __restrict__ vector, Vector3* netVector);

__device__ inline void exclIntCumSum(int* in, const int n) {
	// 1) int* in must point to shared memory
	// 2) int n must be power of 2
	const int tid = threadIdx.x;
	// RBTODO: worry about possible bank conflicts http://www.eecs.umich.edu/courses/eecs570/hw/parprefix.pdf
	
	// build tree of sums
	int stride = 1;
	for (int d = n>>1; d > 0; d >>= 1) {
		__syncthreads();
		if (tid < d) {
			int id = 2*stride*(tid+1)-1;
			in[id] += in[id-stride];
		}
		stride *= 2;
	}
	if (tid == 0) in[n-1] = 0;		/* exclusive cumsum (starts at 0) */

	// traverse down tree and build 'scan'
	for (int d = 1; d < n; d*= 2) {
		stride >>= 1;
		__syncthreads();

		if (tid < d) { // RBTODO: this could be incorrect ==> test
			int id = 2*stride*(tid+1)-1;
			int t = in[id];
			in[id] += in[id-stride];
			in[id-stride] = t;
		}
	}
	__syncthreads();
}

__device__ inline void inclIntCumSum(int* in, const int n) {
	// 1) int* in must point to shared memory
	// 2) int n must be power of 2
	const int tid = threadIdx.x;
	
	// RBTODO: worry about possible bank conflicts http://www.eecs.umich.edu/courses/eecs570/hw/parprefix.pdf
	
	// build tree of sums
	int stride = 1;
	for (int d = n>>1; d > 0; d >>= 1) {
		__syncthreads();
		if (tid < d) {
			int id = 2*stride*(tid+1)-1;
			in[id] += in[id-stride];
		}
		stride *= 2;
	}
	// if (tid == 0) in[n-1] = 0;		/* exclusive cumsum (starts at 0) */

	// traverse down tree and build 'scan'
	for (int d = 1; d < n; d*= 2) {
		stride >>= 1;
		__syncthreads();

		if (tid < d) { // RBTODO: this could be incorrect ==> test
			int id = 2*stride*(tid+1)-1;
			in[id+stride] += in[id];
			/* int t = in[id]; */
			/* in[id] += in[id-stride]; */
			/* in[id-stride] = t; */
		}
	}
	__syncthreads();
}

__device__ inline void atomicAdd(Vector3* address, const Vector3 val) {
	atomicAdd( &(address->x), val.x);
	atomicAdd( &(address->y), val.y);
	atomicAdd( &(address->z), val.z);
}
__device__ inline void atomicAdd(ForceEnergy* address, const ForceEnergy val) {
    atomicAdd( &(address->f), val.f );
    atomicAdd( &(address->e), val.e);
}
