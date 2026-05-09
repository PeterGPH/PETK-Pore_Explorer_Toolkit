#ifndef KERNEL_H_
#define KERNEL_H_
#include "useful.h"

template<const int BlockSize>
static __global__ void BrownParticlesKineticEnergy(Vector3* P_n, int type[], BrownianParticleType* part[], 
                                                   float *vec_red, int num, int num_rb_attached_particles, int num_replicas)
{
    __shared__ __align__(4) float sdata[BlockSize];
    
    Vector3 p1, p2;
    float mass1, mass2;

    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*(BlockSize<<1) + tid;
    unsigned int gridSize = (BlockSize<<1)*gridDim.x;

    sdata[tid] = 0.f; 

    int n = (num*num_replicas);

    while (i < n) 
    {
	const int i1 = (i % num) +  (i/num)*(num+num_rb_attached_particles);
        const int t1 = type[i1];
        const BrownianParticleType& pt1 = *part[t1];

        p1    = P_n[i];
        mass1 = pt1.mass;
        
        if(i + BlockSize < n)
        {
	    const int i2 = ((i+BlockSize) % num) +  ((i+BlockSize)/num)*(num+num_rb_attached_particles);
            const int t2 = type[i2];
            const BrownianParticleType& pt2 = *part[t2];

            p2    = P_n[i+BlockSize];
            mass2 = pt2.mass;

            sdata[tid] += (p1.length2() / mass1 + p2.length2() / mass2); 
        }
        else
            sdata[tid] += p1.length2() / mass1;

        i += gridSize;
    }

    sdata[tid] *= 0.50f;

    __syncthreads();

    if (BlockSize == 512) 
    { 
        if (tid < 256) 
            sdata[tid] += sdata[tid + 256]; 
        __syncthreads();
       if (tid < 128)
            sdata[tid] += sdata[tid + 128];
        __syncthreads();
       if (tid < 64)
            sdata[tid] += sdata[tid + 64];
        __syncthreads();
       if (tid < 32)
            sdata[tid] += sdata[tid + 32];
        __syncthreads();
       if (tid < 16)
            sdata[tid] += sdata[tid + 16];
        __syncthreads();
       if (tid < 8)
            sdata[tid] += sdata[tid + 8];
        __syncthreads();
       if (tid < 4)
            sdata[tid] += sdata[tid + 4];
        __syncthreads();
       if (tid < 2)
            sdata[tid] += sdata[tid + 2];
        __syncthreads();
       if (tid < 1)
            sdata[tid] += sdata[tid + 1];
        __syncthreads();
    }
    else if (BlockSize == 256) 
    {
        if (tid < 128)
            sdata[tid] += sdata[tid + 128];
        __syncthreads();
       if (tid < 64)
            sdata[tid] += sdata[tid + 64];
        __syncthreads();
       if (tid < 32)
            sdata[tid] += sdata[tid + 32];
        __syncthreads();
       if (tid < 16)
            sdata[tid] += sdata[tid + 16];
        __syncthreads();
       if (tid < 8)
            sdata[tid] += sdata[tid + 8];
        __syncthreads();
       if (tid < 4)
            sdata[tid] += sdata[tid + 4];
        __syncthreads();
       if (tid < 2)
            sdata[tid] += sdata[tid + 2];
        __syncthreads();
       if (tid < 1)
            sdata[tid] += sdata[tid + 1];
        __syncthreads(); 
    }
    else if (BlockSize == 128) 
    {
       if (tid < 64)
            sdata[tid] += sdata[tid + 64];
        __syncthreads();
       if (tid < 32)
            sdata[tid] += sdata[tid + 32];
        __syncthreads();
       if (tid < 16)
            sdata[tid] += sdata[tid + 16];
        __syncthreads();
       if (tid < 8)
            sdata[tid] += sdata[tid + 8];
        __syncthreads();
       if (tid < 4)
            sdata[tid] += sdata[tid + 4];
        __syncthreads();
       if (tid < 2)
            sdata[tid] += sdata[tid + 2];
        __syncthreads();
       if (tid < 1)
            sdata[tid] += sdata[tid + 1];
        __syncthreads();
 
    }
    else if (BlockSize == 64)
    {
       if (tid < 32)
            sdata[tid] += sdata[tid + 32];
        __syncthreads();
       if (tid < 16)
            sdata[tid] += sdata[tid + 16];
        __syncthreads();
       if (tid < 8)
            sdata[tid] += sdata[tid + 8];
        __syncthreads();
       if (tid < 4)
            sdata[tid] += sdata[tid + 4];
        __syncthreads();
       if (tid < 2)
            sdata[tid] += sdata[tid + 2];
        __syncthreads();
       if (tid < 1)
            sdata[tid] += sdata[tid + 1];
        __syncthreads();

    }
    else if (BlockSize == 32)
    {
       if (tid < 16)
            sdata[tid] += sdata[tid + 16];
        __syncthreads();
       if (tid < 8)
            sdata[tid] += sdata[tid + 8];
        __syncthreads();
       if (tid < 4)
            sdata[tid] += sdata[tid + 4];
        __syncthreads();
       if (tid < 2)
            sdata[tid] += sdata[tid + 2];
        __syncthreads();
       if (tid < 1)
            sdata[tid] += sdata[tid + 1];
        __syncthreads();

    }
    __syncthreads();
    if (tid == 0) 
        vec_red[blockIdx.x] = sdata[0];
}
//The size must be power of 2, otherwise there is error
//This small kernel is to reduce the small vecotr obtained by
//the reduction kernel from the above.
//The grid size should be one.
//This small routine is to help do further reduction of a small vector

template<int BlockSize>
static __global__ void Reduction(float* dev_vec, float* result, int Size)
{
    __shared__ __align__(4) float data[BlockSize];
    const unsigned int tid = threadIdx.x;
    
    data[tid] = dev_vec[tid];
    size_t idx = tid + BlockSize;
    while(idx < Size)
    {
        data[tid] += dev_vec[idx];
        idx += BlockSize;
    }
    __syncthreads();

    int n = BlockSize;
    while(n > 1)
    {
        n = (n >> 1);
        if(tid < n)
            data[tid] += data[tid+n];
        __syncthreads();
    }
    if(tid == 0) 
        result[0] = data[0];
}
#endif
