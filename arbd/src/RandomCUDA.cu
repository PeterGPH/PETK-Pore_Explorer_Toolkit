#include "RandomCUDA.h"

__global__
void initKernel(unsigned long seed, curandState_t *state, int num);

void Random::init(int num, unsigned long seed) {
	if (states != NULL)
		gpuErrchk(cudaFree(states));
	gpuErrchk(cudaMalloc(&states, sizeof(curandState) * num));
	int nBlocks = num / NUM_THREADS + 1;
	initKernel<<< nBlocks, NUM_THREADS >>>(seed, states, num);
	gpuErrchk(cudaDeviceSynchronize());

	// Create RNG and set seed
	curandCreateGenerator(&generator, CURAND_RNG_PSEUDO_XORWOW);
	curandSetPseudoRandomGeneratorSeed(generator, seed);
	
	if (uniform_d != NULL) 
        {
            gpuErrchk(cudaFree(uniform_d));
            uniform_d = NULL;
        }
        if(integer_d!=NULL)
        {
            gpuErrchk(cudaFree(integer_d));
            integer_d = NULL;
        }
        if(gaussian_d!=NULL)
        {
            gpuErrchk(cudaFree(gaussian_d));
            gaussian_d = NULL;
        }
        if(integer_h!=NULL)
        {
	    delete[] integer_h;
            integer_h = NULL;
        }
        if(uniform_h!=NULL)
        {
	    delete[] uniform_h;
            uniform_h = NULL;
	}
        if(gaussian_h!=NULL)
        {
	    delete[] gaussian_h;
            gaussian_h = NULL;
        }
	gpuErrchk(cudaMalloc((void**)&uniform_d, sizeof(float) * RAND_N));
	gpuErrchk(cudaMalloc((void**)&integer_d, sizeof(unsigned int) * RAND_N));
	gpuErrchk(cudaMalloc((void**)&gaussian_d, sizeof(float) * RAND_N));
	integer_h = new unsigned int[RAND_N];
	uniform_h = new float[RAND_N];
	gaussian_h = new float[RAND_N];
	uniform_n = 0;
	integer_n = 0;
	gaussian_n = 0;
}

float Random::uniform() {
	if (uniform_n < 1) {
		cuRandchk(curandGenerateUniform(generator, uniform_d, RAND_N));
		gpuErrchk(cudaMemcpy(uniform_h, uniform_d, sizeof(float) * RAND_N, cudaMemcpyDeviceToHost));
		uniform_n = RAND_N;
	}
	return uniform_h[--uniform_n];
}

unsigned int Random::poisson(float lambda) {
	const float l = exp(-lambda);
	unsigned int k = 0;
	float p = uniform();
	while (p >= l) {
		p *= uniform();
		k = k + 1;
	}
	return k;
}

unsigned int Random::integer() {
	if (integer_n < 1) {
		curandGenerate(generator, integer_d, RAND_N);
		gpuErrchk(cudaMemcpy(integer_h, integer_d, sizeof(unsigned int) * RAND_N, cudaMemcpyDeviceToHost));
		integer_n = RAND_N;
	}
	return integer_h[--integer_n];
}

void Random::reorder(int a[], int n) {
	for (int i = 0; i < (n-1); ++i) {
		unsigned int j = i + (integer() % (n-i));
		if ( j == i )
			continue;
		std::swap<int>(a[i], a[j]);
		const int tmp = a[j];
		a[j] = a[i];
		a[i] = tmp;
	}
}

__global__ 
void initKernel(unsigned long seed, curandState_t *state, int num) {
       int idx = blockIdx.x * blockDim.x + threadIdx.x;
       int step = blockDim.x * gridDim.x;
       for(int i = idx; i < num; i=i+step)
       {
           curandState_t local;
           // curand_init(clock64()+seed,i,0,&local);
           //curand_init(clock64(),i,0,&state[i]);
	   curand_init(seed,i,0,&local);
           state[(size_t)i] = local;
       }

}
