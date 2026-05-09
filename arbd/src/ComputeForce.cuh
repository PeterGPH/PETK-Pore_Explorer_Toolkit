// ComputeForce.cuh
//
// Terrance Howard <heyterrance@gmail.com>
#pragma once
#include <cassert>

#include "CudaUtil.cuh"

#include "TabulatedMethods.cuh"

// From TabulatedMethods.cuh: constexpr float BD_PI = 3.1415927f; 
constexpr size_t MAX_CELLS_FOR_CELLNEIGHBORLIST = 1<<25;
constexpr size_t MAX_NLIST_PAIRS = 1<<27; // Reduce if ARBD crashes immediately with GPU memory allocation error

// texture<int,    1, cudaReadModeElementType>      NeighborsTex;
// texture<int,    1, cudaReadModeElementType> pairTabPotTypeTex;
//texture<int2,   1, cudaReadModeElementType>      pairListsTex;
// texture<float4, 1, cudaReadModeElementType>            PosTex;

__host__ __device__
EnergyForce ComputeForce::coulombForce(Vector3 r, float alpha,
																			 float start, float len) {
	float d = r.length();
	
	if (d >= start + len)
		return EnergyForce();
	if (d <= start) {
		float energy = alpha/d - alpha/start + 0.5f*alpha/(start*start)*len;
		Vector3 force = -alpha/(d*d*d)*r;

		return EnergyForce(energy, force);
	}

	// Switching.
	float c = alpha / (start * start);
	float energy = 0.5f/len * c * (start + len - d) * (start + len - d);
	Vector3 force = -c * (1.0f - (d - start) / len) / d * r;
	return EnergyForce(energy, force);
}

__host__ __device__
EnergyForce ComputeForce::coulombForceFull(Vector3 r, float alpha) {
	float d = r.length();
	return EnergyForce(alpha/d, -alpha/(d*d*d)*r);
}

__host__ __device__
EnergyForce ComputeForce::softcoreForce(Vector3 r, float eps, float rad6) {
	const float d2 = r.length2();
	const float d6 = d2*d2*d2;

	Vector3 force = -12*eps*(rad6*rad6/(d6*d6*d2) - rad6/(d6*d2))*r;
		if (isnan(force.x) or isnan(force.y) or isnan(force.z))
			printf(">>>> Damn.\n");

	if (d6 < rad6) {
		const float d6_2 = d6 * d6;
		const float rad6_2 = rad6 * rad6;
		float e = eps * ((rad6_2 / (d6_2)) - (2.0f * rad6 / d6)) + eps;
		Vector3 f = -12.0f * eps * (rad6_2 / (d6_2 * d2) - rad6 / (d6 * d2)) * r;
		return EnergyForce(e, f);
	}

	return EnergyForce();
}

__global__
void computeFullKernel(Vector3 force[], Vector3 pos[], int type[],
											 float tableAlpha[], float tableEps[], float tableRad6[],
											 int num, int numParts, BaseGrid* sys, float g_energies[],
											 int gridSize, int numReplicas, bool get_energy) {
	// Calculate the ID of each thread.
	// The combination of (X, Y) is unique among all threads
	// id_in_block is unique to a single thread block
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	float energy_local = 0.0f;

	// For all threads representing a valid pair of particles
	if (i < num * numReplicas) {
		const int repID = i / num;
		const int typei = type[i];
		const Vector3& posi = pos[i];
		float alpha, eps, rad6;
		Vector3 force_local(0.0f);
		int typej = -1;

		for (int j = repID * num; j < (repID + 1) * num; ++j) {
			if (i == j) continue;
			int newj = type[j];
			if (typej != newj) {
				// Get new values if type[j-1] != type[j]
				// Save time; avoid reading from memory
				typej = newj;
				alpha = tableAlpha[typei * numParts + typej];
				eps = tableEps[typei * numParts + typej];
				rad6 = tableRad6[typei * numParts + typej];
			}
			Vector3 dr = sys->wrapDiff(pos[j] - posi);
			EnergyForce fc = ComputeForce::coulombForceFull(dr, alpha);
			EnergyForce fh = ComputeForce::softcoreForce(dr, eps, rad6);
			// Only update the force in the X particle.
			// Another thread will handle the Y particle
			force_local += fc.f + fh.f;

			// Check if there is a bond between these two.
			// If there is, and the bond's flag is ADD, then add some more force to
			// the interaction.

			// Only update the energy once. The other thread handling x and y will
			// not update the energy.
			if (get_energy && j > i && i < num)
				energy_local += fc.e + fh.e;
		}
		force[i] = force_local;
		if (get_energy && i < num)
			g_energies[i] = energy_local;
	}
}


__global__
void computeSoftcoreFullKernel(Vector3 force[], Vector3 pos[], int type[],
															 float tableEps[], float tableRad6[],
															 int num, int numParts, BaseGrid* sys,
															 float g_energies[], int gridSize,
															 int numReplicas, bool get_energy) {
	// Calculate the ID of each thread.
	// The combination of (X, Y) is unique among all threads
	// id_in_block is unique to a single thread block
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	float energy_local = 0.0f;
	// For all threads representing a valid pair of particles
	if (i < num * numReplicas) {
		const int repID = i / num;
		const int typei = type[i];
		const Vector3& posi = pos[i];
		float eps, rad6;
		Vector3 force_local(0.0f);

		int typej = -1;
		for (int j = repID * num; j < (repID + 1) * num; ++j) {
			if (i == j) continue;
			int newj = type[j];
			if (typej != newj) {
				typej = newj;
				eps = tableEps[typei * numParts + typej];
				rad6 = tableRad6[typei * numParts + typej];
			}

			Vector3 dr = sys->wrapDiff(pos[j] - posi);
			EnergyForce fh = ComputeForce::softcoreForce(dr, eps, rad6);
			// Only update the force in the particle i
			// Another thread will handle the particle j
			force_local += fh.f;
			// Only update the energy once. The other thread handling x and y
			// will not update the energy
			if (get_energy && j > i)
				energy_local += fh.e;
		}
		force[i] = force_local;
		if (get_energy)
			g_energies[i] = energy_local;
	}
}

__global__
void computeElecFullKernel(Vector3 force[], Vector3 pos[], int type[],
													 float tableAlpha[], int num, int numParts,
													 BaseGrid* sys, float g_energies[],
													 int gridSize, int numReplicas,
													 bool get_energy) {
	// Calculate the ID of each thread.
	// The combination of (X, Y) is unique among all threads
	// id_in_block is unique to a single thread block
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	float energy_local = 0.0f;
	// For all threads representing a valid pair of particles
	if (i < num * numReplicas) {
		const int repID = i / num;
		const int typei = type[i];
		const Vector3 posi = pos[i];
		float alpha;

		Vector3 force_local(0.0f);
		int typej = -1;
		for (int j = repID * num; j < num * (repID-1); j++) {
			if (i == j) continue;
			int newj = type[j];
			if (typej != newj) {
				typej = newj;
				alpha = tableAlpha[typei * numParts + typej];
			}
			const Vector3 dr = sys->wrapDiff(pos[j] - posi);
			EnergyForce fc = ComputeForce::coulombForceFull(dr, alpha);
			// Only update the force in the X particle.
			// Another thread will handle the Y particle
			force_local += fc.f;
			// Only update the energy once. The other thread handling x and y
			// will not update the energy
			if (get_energy and j > i)
				energy_local += fc.e;
		}
		force[i] = force_local;
		if (get_energy)
			g_energies[i] = energy_local;
	}
}


/* const __device__ int maxPairs = 1 << 14; */

/* __global__ */
/* void pairlistTest(Vector3 pos[], int num, int numReplicas, */
/* 									BaseGrid* sys, CellDecomposition* decomp, */
/* 									const int nCells, const int blocksPerCell, */
/* 									int* g_numPairs, int* g_pairI, int* g_pairJ ) { */
/* 	const int gtid = threadIdx.x + blockIdx.x*blockDim.x; */
/* 	for (int i = gtid; i < gridDim.x*100; i+=blockDim.x) { */
/* 		assert( g_numPairs[i] == 0 ); */
/* 		assert( g_pairI[i] != NULL ); */
/* 		assert( g_pairJ[i] != NULL ); */
/* 	} */
/* } */

__device__ int* exSum;
void initExSum() {
    int tmp = 0;
    int* devPtr;
    cudaMalloc(&devPtr, sizeof(int));
    cudaMemcpyToSymbol(exSum, &devPtr, sizeof(int*));
    cudaMemcpy(devPtr, &tmp, sizeof(int), cudaMemcpyHostToDevice);

}
int getExSum() {
    int tmp;
    int* devPtr;
    cudaMemcpyFromSymbol(&devPtr, exSum, sizeof(int*));
    cudaMemcpy(&tmp, devPtr, sizeof(int), cudaMemcpyDeviceToHost);
    return tmp;
}
//
__device__
int computeCellNeighbor( const int3 cells, const int3 cell_idx, const int dx, const int dy, const int dz )
{
    int idx = cell_idx.x;
    int idy = cell_idx.y;
    int idz = cell_idx.z;

    int u = idx + dx;
    int v = idy + dy;
    int w = idz + dz;

    int nID;
    if (cells.x == 1 and u != 0) nID = -1;
    else if (cells.y == 1 and v != 0) nID =  -1;
    else if (cells.z == 1 and w != 0) nID = -1;
    else if (cells.x == 2 and (u < 0 || u > 1)) nID = -1;
    else if (cells.y == 2 and (v < 0 || v > 1)) nID = -1;
    else if (cells.z == 2 and (w < 0 || w > 1)) nID = -1;
    else
    {
	u = (u + cells.x) % cells.x;
	v = (v + cells.y) % cells.y;
	w = (w + cells.z) % cells.z;
	nID = w + cells.z * (v + cells.y * u);
    }

    return nID;
}

__global__ 
void createNeighborsList(const int3 *Cells,int* __restrict__ CellNeighborsList)
{
    const int tid = threadIdx.x + blockDim.x * blockIdx.x;
    const int3 cells = Cells[0]; 
    const int nCells = cells.x * cells.y * cells.z;
    const int Size   = blockDim.x * gridDim.x;
    int   nID;
    
    for (int cID = tid; cID < nCells; cID += Size) {

        int idz = cID %  cells.z;
        int idy = cID /  cells.z % cells.y;
        int idx = cID / (cells.z * cells.y);

        int count = 0;
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dz = -1; dz <= 1; ++dz) {

		    nID = computeCellNeighbor( cells, make_int3(idx,idy,idz), dx, dy, dz );
                    CellNeighborsList[size_t(count+27*cID)] = nID;
                    ++count;
                    //__syncthreads();
                }
            }
        }
    }
}
template<const int BlockSize,const int Size,const int N>
__global__ void createPairlists(Vector3* __restrict__ pos, const int num, const int numReplicas,
                                const BaseGrid* __restrict__ sys, const CellDecomposition* __restrict__ decomp,
                                const int nCells,int* g_numPairs, int2* g_pair, int numParts, const int* __restrict__ type,
                                int* __restrict__ g_pairTabPotType, const Exclude* __restrict__ excludes,
                                const int2* __restrict__ excludeMap, const int numExcludes, float pairlistdist2, cudaTextureObject_t PosTex, cudaTextureObject_t NeighborsTex)
{
    __shared__ float4 __align__(16) particle[N];
    __shared__ int     Index_i[N];

    const int TotalBlocks  = gridDim.x * gridDim.y;
    const int cells        = TotalBlocks / Size;
    const int cell_start   = ( blockIdx.x + gridDim.x * blockIdx.y) / Size;
    const int pid_start    = ((blockIdx.x + gridDim.x * blockIdx.y) % Size) * N;
    const int tid          =   threadIdx.x + blockDim.x * threadIdx.y
                                           + blockDim.x *  blockDim.y * threadIdx.z;
    const int warpLane     = tid % WARPSIZE;
    const int nReps        = gridDim.z;
    const int idx_  = tid % N;
    const int idx__ = tid / N;
    const int Step1 = Size * N;
    const int Step2 = Size / N; 

    const CellDecomposition::cell_t* __restrict__ cellInfo = decomp->getCells();

    for(int repID = blockIdx.z; repID < numReplicas; repID += nReps)
    {
        for(int cellid_i = cell_start; cellid_i < nCells; cellid_i += cells)
        {
            CellDecomposition::range_t rangeI = decomp->getRange(cellid_i,repID);
            int Ni = rangeI.last-rangeI.first;

            for(int pid_i = pid_start; pid_i < Ni; pid_i += Step1)
            {
                __syncthreads();
                if(tid + pid_i < Ni && tid < N)
                {
                    Index_i [tid] = cellInfo[rangeI.first+pid_i+tid].particle;
                    particle[tid] = tex1Dfetch<float4>(PosTex,Index_i[tid]);
                }
                __syncthreads();

                if(idx_ + pid_i < Ni)
                {
                    int ai = Index_i[idx_];
                    Vector3 A(particle[idx_]);

                    int2 ex_pair = make_int2(-1,-1);
                    if(numExcludes > 0 && excludeMap != NULL)
                    {
                        ex_pair = excludeMap[ai -repID * num];
                    }

                    //loop over neighbor directions
                    for(int idx = 0; idx < 27; ++idx)
                    {

			int currEx = ex_pair.x;
			int nextEx = (ex_pair.x >= 0) ? excludes[currEx].ind2 : -1;

			int neighbor_cell;
			if (nCells < MAX_CELLS_FOR_CELLNEIGHBORLIST) {
			    neighbor_cell = tex1Dfetch<int>(NeighborsTex,idx+27*cellid_i);
			} else {
			    int3 cells = decomp->nCells;
			    int3 cell_idx = make_int3(cellid_i %  cells.z,
						      cellid_i /  cells.z % cells.y,
						      cellid_i / (cells.z * cells.y));

			    int dz = (idx % 3) - 1;
			    int dy = ((idx/3) % 3) - 1;
			    int dx = ((idx/9) % 3) - 1;
			    neighbor_cell = computeCellNeighbor( decomp->nCells, cell_idx, dx, dy, dz );
			}

                        if(neighbor_cell < 0)
                        {
                            continue;
                        }

                        CellDecomposition::range_t rangeJ = decomp->getRange(neighbor_cell,repID);
                        int Nj = rangeJ.last-rangeJ.first;

                        // In each neighbor cell, loop over particles
                        for(int pid_j = idx__; pid_j < Nj; pid_j += Step2)
                        {
                            
                            int aj  = cellInfo[pid_j+rangeJ.first].particle;
                            if( aj <= ai)
                            {
                                continue;
                            }

                            while (nextEx >= 0 && nextEx < ( aj - repID * num))
                            {
                                nextEx = (currEx < ex_pair.y - 1) ? excludes[++currEx].ind2 : -1;
                            }

                            if (nextEx == (aj - repID * num))
                            {
                                #ifdef DEBUGEXCLUSIONS
                                atomicAggInc( exSum, warpLane );
                                #endif
                                nextEx = (currEx < ex_pair.y - 1) ? excludes[++currEx].ind2 : -1;
                                continue;
                            }

                            float4 b = tex1Dfetch<float4>(PosTex,aj);
                            Vector3 B(b.x,b.y,b.z);

                            float dr = (sys->wrapDiff(A-B)).length2();
                            if(dr <= pairlistdist2)
                            {
                                int gid = atomicAggInc( g_numPairs, warpLane );
				if (gid < MAX_NLIST_PAIRS) {
				    int pairType = type[ai] + type[aj] * numParts;

				    g_pair[gid] = make_int2(ai,aj);
				    g_pairTabPotType[gid] = pairType;
				} else {
				    if (gid == MAX_NLIST_PAIRS)
					printf("RAN OUT OF PAIRLIST SPACE\n");
				}
                            }
                        }
                    }
                }
            }
        }
    }
}

__global__
void createPairlists_debug(Vector3* __restrict__ pos, const int num, const int numReplicas,
                                const BaseGrid* __restrict__ sys, const CellDecomposition* __restrict__ decomp,
                                const int nCells,
                                int* g_numPairs, int2* g_pair,
                                int numParts, const int* __restrict__ type, int* __restrict__ g_pairTabPotType,
                                const Exclude* __restrict__ excludes, const int2* __restrict__ excludeMap, const int numExcludes,
                                float pairlistdist2)
{
    // TODO: loop over all cells with edges within pairlistdist2
    // Loop over threads searching for atom pairs
    //   Each thread has designated values in shared memory as a buffer
    //   A sync operation periodically moves data from shared to global
    const int tid = threadIdx.x;
    const int warpLane = tid % WARPSIZE; /* RBTODO: optimize */
    const int split = 32;                                   /* numblocks should be divisible by split */
    /* const int blocksPerCell = gridDim.x/split;  */
    const CellDecomposition::cell_t* __restrict__ cellInfo = decomp->getCells();
    for (int cID = 0 + (blockIdx.x % split); cID < nCells; cID += split)
    {
        for (int repID = 0; repID < numReplicas; repID++)
        {
            const CellDecomposition::range_t rangeI = decomp->getRange(cID, repID);
            for (int ci = rangeI.first + blockIdx.x/split; ci < rangeI.last; ci += gridDim.x/split)
            {
                const int ai = cellInfo[ci].particle;
                const CellDecomposition::cell_t celli = cellInfo[ci];
                const int ex_start = (numExcludes > 0 && excludeMap != NULL) ? excludeMap[ai -repID*num].x : -1;
                const int ex_end   = (numExcludes > 0 && excludeMap != NULL) ? excludeMap[ai -repID*num].y : -1;
                for(int x = -1; x <= 1; ++x) 
                {
                    for(int y = -1; y <= 1; ++y) 
                    {
                        for (int z = -1; z <= 1; ++z) 
                        {
                            const int nID = decomp->getNeighborID(celli, x, y, z);
                            //const int nID = CellNeighborsList[x+27*cID];//elli.id]; 
                            if (nID < 0) continue; // Initialize exclusions
                            // TODO: optimize exclusion code (and entire kernel)
                            int currEx = ex_start;
                            int nextEx = (ex_start >= 0) ? excludes[currEx].ind2 : -1;
                            //int ajLast = -1; // TODO: remove this sanity check
                            const CellDecomposition::range_t range = decomp->getRange(nID, repID);
                            for (int n = range.first + tid; n < range.last; n+=blockDim.x) 
                            {
                                const int aj = cellInfo[n].particle;
                                if (aj <= ai) continue;
                                // Skip excludes
                                // Implementation requires that aj increases monotonically
                                //assert( ajLast < aj ); ajLast = aj; // TODO: remove this sanity check
                                while (nextEx >= 0 && nextEx < (aj - repID * num)) // TODO get rid of this
                                    nextEx = (currEx < ex_end - 1) ? excludes[++currEx].ind2 : -1;
                                if (nextEx == (aj - repID * num))
                                {
                                    #ifdef DEBUGEXCLUSIONS
                                    atomicAggInc( exSum, warpLane );
                                    #endif
                                    nextEx = (currEx < ex_end - 1) ? excludes[++currEx].ind2 : -1;
                                    continue;
                                }
                                // TODO: Skip non-interacting types for efficiency
                                // Skip ones that are too far away
                                const float dr = (sys->wrapDiff(pos[aj] - pos[ai])).length2();
                                if (dr > pairlistdist2) continue;
                                // Add to pairlist
                                int gid = atomicAggInc( g_numPairs, warpLane );
                                int pairType = type[ai] + type[aj] * numParts;
                                g_pair[gid] = make_int2(ai,aj);
                                g_pairTabPotType[gid] = pairType;
                            }
                        }                      
                    }
                }
            }                              
        }
    }
}

// TODO: deprecate?
__global__
void computeKernel(Vector3 force[], Vector3 pos[], int type[],
									 float tableAlpha[], float tableEps[], float tableRad6[],
									 int num, int numParts, BaseGrid* sys,
									 CellDecomposition* decomp,
									 float g_energies[], float switchStart, float switchLen,
									 int gridSize, int numReplicas, bool get_energy) {
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	float energy_local = 0.0f;

	// i - index of the particle in the original, unsorted array
	if (i < num * numReplicas) {
		const int repID = i / num;
		const int typei = type[i];
		const Vector3 posi = pos[i];
		// TODO: Fix this: Find correct celli (add a new function to
		//       CellDecomposition, binary search over cells)
		CellDecomposition::cell_t celli = decomp->getCellForParticle(i);
		float alpha(0.0f), eps(0.0f), rad6(0.0f);
		Vector3 force_local(0.0f);

		const CellDecomposition::cell_t* pairs = decomp->getCells();
		for (int x = -1; x <= 1; ++x) {
			for (int y = -1; y <= 1; ++y) {
				for (int z = -1; z <= 1; ++z) {
					const int nID = decomp->getNeighborID(celli, x, y, z);
					// Skip if got wrong or duplicate cell.
					if (nID < 0) continue;

					const CellDecomposition::range_t range = decomp->getRange(nID, repID);

					int typej = -1;
					for (int n = range.first; n < range.last; ++n) {
						const int j = pairs[n].particle;
//						if (j < 0)
//							printf("%d -> pairs[%d].particle %d\n", i, n, j);
						if (j == i) continue;
						const int newj = type[j];
						// Update values.
						if (typej != newj) {
							typej = newj;
							alpha = tableAlpha[typei * numParts + typej];
							eps = tableEps[typei * numParts + typej];
							rad6 = tableRad6[typei * numParts + typej];
						}

						const Vector3 dr = sys->wrapDiff(pos[j] - posi);
						if (dr.length() < 1e-4) {
//							printf("dr = %g << 1: %d -> %d on [%d, %d)\n",
//										 dr.length(), i, j, range.first, range.last);
						}

						const EnergyForce fc =
								ComputeForce::coulombForce(dr, alpha, switchStart, switchLen);
						const EnergyForce fh = ComputeForce::softcoreForce(dr, eps, rad6);

						force_local += fc.f + fh.f;
						energy_local += 0.5f * (fc.e + fh.e);
					} 	// n
				} 		// z
			} 			// y
		} 				// x
		force[i] = force_local;
		if (isnan(force_local.x) or isnan(force_local.y) or isnan(force_local.z)) {
//			printf("Nan FORCE!\n");
			force[i] = Vector3(0.0f);
		}
		if (get_energy)
			g_energies[i] = energy_local;
	}
}

__device__ int pairForceCounter = 0;
__global__ void printPairForceCounter() {
	if (threadIdx.x + blockIdx.x == 0)
		printf("Computed the force for %d pairs\n", pairForceCounter);
}

template<const int BlockSize>
__device__ inline void _computeTabulatedKernel(Vector3* force, const BaseGrid* __restrict__ sys, 
					       float cutoff2, const int numPairs, const int2* __restrict__ g_pair, 
					       const int* __restrict__ g_pairTabPotType, TabulatedPotential** __restrict__ tablePot,
					       cudaTextureObject_t pairListsTex, cudaTextureObject_t PosTex, cudaTextureObject_t pairTabPotTypeTex
    )
{
    const int tid = threadIdx.x + blockDim.x * threadIdx.y
                                         + blockDim.x *  blockDim.y * threadIdx.z 
                                         + BlockSize  *( blockIdx.x + gridDim.x * blockIdx.y 
                                         + gridDim.x  * gridDim.y   * blockIdx.z );

    const int TotalThreads = BlockSize * gridDim.x * gridDim.y * gridDim.z;
    for (int i = tid; i < numPairs; i += TotalThreads) 
    {
        //int2 pair = g_pair[i];
        int2 pair = tex1Dfetch<int2>(pairListsTex,i);
        //int  ind  = tex1Dfetch(pairTabPotTypeTex,i); 

        int ai = pair.x;
        int aj = pair.y;
                        
        //int ind = g_pairTabPotType[i];

        Vector3 a(tex1Dfetch<float4>(PosTex, ai));
        Vector3 b(tex1Dfetch<float4>(PosTex, aj));
        Vector3 dr = sys->wrapDiff(b-a);
        
        float d2 = dr.length2();
        int  ind  = tex1Dfetch<int>(pairTabPotTypeTex,i);
        if (tablePot[ind] != NULL && d2 <= cutoff2) 
        {
            Vector3 f = tablePot[ind]->computef(dr,d2);
            atomicAdd( &force[ai],  f );
            atomicAdd( &force[aj], -f );
        }
    }
}

template<const int BlockSize>
__global__ void computeTabulatedKernel(Vector3* force, const BaseGrid* __restrict__ sys, 
                                       float cutoff2, const int* __restrict__ g_numPairs, const int2* __restrict__ g_pair, 
                                       const int* __restrict__ g_pairTabPotType, TabulatedPotential** __restrict__ tablePot,
				       cudaTextureObject_t pairListsTex, cudaTextureObject_t PosTex, cudaTextureObject_t pairTabPotTypeTex) {
    _computeTabulatedKernel<BlockSize>(force,sys,
				       cutoff2, *g_numPairs, g_pair,
				       g_pairTabPotType, tablePot,
				       pairListsTex, PosTex, pairTabPotTypeTex);
}

template<const int BlockSize>
__global__ void computeTabulatedKernel(Vector3* force, const BaseGrid* __restrict__ sys, 
                                       float cutoff2, const int2* __restrict__ g_pair, 
                                       const int* __restrict__ g_pairTabPotType, TabulatedPotential** __restrict__ tablePot,
				       cudaTextureObject_t pairListsTex, cudaTextureObject_t PosTex, cudaTextureObject_t pairTabPotTypeTex,
				       int start, int numPairs) {
    _computeTabulatedKernel<BlockSize>(force,sys,
				       cutoff2, numPairs, g_pair+start,
				       g_pairTabPotType+start, tablePot,
				       pairListsTex, PosTex, pairTabPotTypeTex);
} 

__global__ void clearEnergies(float* __restrict__  g_energies, int num) {
	for (int i = threadIdx.x+blockIdx.x*blockDim.x; i < num; i+=blockDim.x*gridDim.x) {
		g_energies[i] = 0.0f;
	}
}

__global__ void computeTabulatedEnergyKernel(Vector3* force, const Vector3* __restrict__ pos,
				const BaseGrid* __restrict__ sys, float cutoff2,
				const int* __restrict__ g_numPairs,	const int2* __restrict__ g_pair, const int* __restrict__ g_pairTabPotType, 	TabulatedPotential** __restrict__ tablePot, float* g_energies) {
	const int numPairs = *g_numPairs;
	for (int i = threadIdx.x+blockIdx.x*blockDim.x; i < numPairs; i+=blockDim.x*gridDim.x) {
		const int2 pair = g_pair[i];
		const int ai = pair.x;
		const int aj = pair.y;
		const int ind = g_pairTabPotType[i];

		// RBTODO: implement wrapDiff2, returns dr2 (???)
		Vector3 dr = pos[aj] - pos[ai];
		dr = sys->wrapDiff(dr);
		float d2 = dr.length2();
		// RBTODO: order pairs according to distance to reduce divergence // not actually faster
		
		if (tablePot[ind] != NULL && d2 <= cutoff2) { 
			EnergyForce fe = tablePot[ind]->compute(dr,d2);
			atomicAdd( &force[ai],  fe.f );
			atomicAdd( &force[aj], -fe.f );
			// RBTODO: reduce energies
			atomicAdd( &(g_energies[ai]), fe.e*0.5f );
			atomicAdd( &(g_energies[aj]), fe.e*0.5f );
		}
	}
}


// =============================================================================
// Kernel computes forces between Brownian particles (ions)
// NOT using cell decomposition
//
__global__
void computeTabulatedFullKernel(Vector3 force[], Vector3 pos[], int type[], TabulatedPotential* tablePot[], TabulatedPotential* tableBond[], int num, int numParts, BaseGrid *sys, Bond bonds[], int2 bondMap[], int numBonds, Exclude excludes[], int2 excludeMap[], int numExcludes, float g_energies[], int gridSize, int numReplicas, bool get_energy, Angle angles[]) 
{
	// Thread's unique ID.
	const int i = blockIdx.x * blockDim.x + threadIdx.x;

	// Initialize interaction energy (per particle)
	float energy_local = 0.0f;

	// Loop over ALL particles in ALL replicas
	if (i < num * numReplicas) {
		const int repID = i / num;

		// Each particle may have a varying number of bonds.
		// bondMap is an array with one element for each particle which keeps track
		// of where a particle's bonds are stored in the bonds array.
		// bondMap[i].x is the index in the bonds array where the ith particle's
		// bonds begin.
		// bondMap[i].y is the index in the bonds array where the ith particle's
		// bonds end.
		const int bond_start	= (bondMap != NULL) ? bondMap[i - repID * num].x : -1;
		const int bond_end 		= (bondMap != NULL) ? bondMap[i - repID * num].y : -1;

		// currBond is the index in the bonds array that we should look at next
		// currBond is initialized to bond_start because that is the first index of the
		// bonds array where this particle's bonds are stored
		int currBond = bond_start;

		// nextBond is the ID number of the next particle that this particle is bonded to
		// If this particle has at least one bond, then nextBond is initialized to be the
	 	// first particle that this particle is bonded to
		int nextBond = (bond_start >= 0) ? bonds[bond_start].ind2 : -1;

		// Same as for bonds, but for exclusions now
		const int ex_start 	= (excludeMap != NULL) ? excludeMap[i].x : -1;
		const int ex_end 		= (excludeMap != NULL) ? excludeMap[i].y : -1;
		int currEx = ex_start;
		int nextEx = (ex_start >= 0) ? excludes[ex_start].ind2 : -1;

		// Particle's type and position
		const int typei = type[i];
		const Vector3& posi = pos[i];

		// Initialize force_local - force on a particle (i)
		Vector3 force_local(0.0f);

		int typej = -1;
		int ind = -1;

		// Loop over ALL particles in a replica, where current particle belongs to
		const size_t first 	= repID * num;
		const size_t last	= first + num;
		for (int j = first; j < last; ++j) {
			if (i == j) continue;
			
			int newj = type[j];
			if (typej != newj) {
				typej = newj;
				ind = typei + typej * numParts;
			}

			Vector3 dr = sys->wrapDiff(pos[j] - posi);

			EnergyForce ft(0.0f, Vector3(0.0f));

			if (nextEx == (j - repID * num))
				nextEx = (currEx < ex_end - 1) ? excludes[++currEx].ind2 : -1;
			else if (tablePot[ind] != NULL)
				ft = tablePot[ind]->compute(dr);

			// If the next bond we want is the same as j, then there is a bond between
			// particles i and j.
			if (nextBond == (j - repID * num) and tableBond != NULL) {
				// If the user has specified the REPLACE option for this bond, then
				// overwrite the force we calculated from the regular tabulated
				// potential. If the user has specified the ADD option, then add the bond
				// force to the tabulated potential value.
				EnergyForce bond_ef = tableBond[bonds[currBond].tabFileIndex]->compute(dr);
				switch (bonds[currBond].flag) {
					case Bond::REPLACE:	ft  = bond_ef; break;
					case Bond::ADD:     ft += bond_ef; break;
				}

				// Increment currBond, so that we can find the index of the next particle
				// that this particle is bonded to
				if (currBond < bond_end - 1) nextBond = bonds[++currBond].ind2;
				else nextBond = -1;
			}

			force_local += ft.f;
			if (get_energy and j > i)
				energy_local += ft.e;

		} // Loop over all particles (i != j) in replica repID

		force[i] = force_local;
		if (get_energy and i < num)
			g_energies[i] = energy_local;
	}
}

__global__
void computeAngles(Vector3 force[], Vector3 pos[],
									 Angle angles[], TabulatedAnglePotential* tableAngle[],
									 int numAngles, int num, BaseGrid* sys,
									 float g_energies[], bool get_energy) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	float energy_local = 0.0f;
	Vector3 force_local(0.0f);
	if (idx < num) {
		for (int i = 0; i < numAngles; ++i) {
			Angle& a = angles[i];
			const int ind = a.getIndex(idx);
			if (ind >= 0) {
				EnergyForce ef = tableAngle[a.tabFileIndex]->computeOLD(&a, pos, sys, ind);
				force_local += ef.f;
				if (ind == 1 and get_energy)
					energy_local += ef.e;
			}
		}
		force[idx] += force_local;
		if (get_energy)
			g_energies[idx] += energy_local;
	}
}

// TODO: add kernels for energy calculations
//__global__ void computeTabulatedBonds(Vector3* force,
//				Vector3* __restrict__ pos,
//				BaseGrid* __restrict__ sys,
//				int numBonds, int3* __restrict__ bondList_d, TabulatedPotential** tableBond) {
__global__
void computeTabulatedBonds(Vector3* force, Vector3* __restrict__ pos, BaseGrid* __restrict__ sys, 
int numBonds, int3* __restrict__ bondList_d, TabulatedPotential** tableBond, float* energy, bool get_energy)
{
	// Loop over ALL bonds in ALL replicas
	for (int bid = threadIdx.x+blockIdx.x*blockDim.x; bid<numBonds; bid+=blockDim.x*gridDim.x) {
		// Initialize interaction energy (per particle)
		// float energy_local = 0.0f;
		
		int i = bondList_d[bid].x;
		int j = bondList_d[bid].y;

		// Find the distance between particles i and j,
		// wrapping this value if necessary
		const Vector3 dr = sys->wrapDiff(pos[j] - pos[i]);

		//Vector3 force_local = tableBond[ bondList_d[bid].z ]->computef(dr,dr.length2());
	        EnergyForce fe_local = tableBond[ bondList_d[bid].z ]->compute(dr,dr.length2());	
		//atomicAdd( &force[i], force_local );
		//atomicAdd( &force[j], -force_local );
		atomicAdd( &force[i], fe_local.f );
                atomicAdd( &force[j], -fe_local.f );

		if (get_energy)
		{
		 	//TODO: clarification on energy computation needed, consider changing.
		 	atomicAdd( &energy[i], fe_local.e*0.5f);
		        atomicAdd( &energy[j], fe_local.e*0.5f);
		}
	}
}

__global__
void computeTabulatedAngles(Vector3* force,
				Vector3* __restrict__ pos,
				BaseGrid* __restrict__ sys,
				int numAngles, int4* __restrict__ angleList_d, TabulatedAnglePotential** tableAngle, float* energy, bool get_energy) {
	// Loop over ALL angles in ALL replicas
	for (int i = threadIdx.x+blockIdx.x*blockDim.x; i<numAngles; i+=blockDim.x*gridDim.x) {
		int4& ids = angleList_d[i];
		computeAngle(tableAngle[ ids.w ], sys, force, pos, ids.x, ids.y, ids.z, energy, get_energy);
	// if (get_energy)
	// {
	//     //TODO: clarification on energy computation needed, consider changing.
	//     atomicAdd( &g_energies[i], energy_local);
	//     //atomicAdd( &g_energies[j], energy_local);
	// }
	}
}

__global__
void computeTabulatedBondAngles(Vector3* force,
				Vector3* __restrict__ pos,
				BaseGrid* __restrict__ sys,
				int numBondAngles, int4* __restrict__ bondAngleList_d, TabulatedAnglePotential** tableAngle,
				TabulatedPotential** tableBond,
				float* energy, bool get_energy) {
	// Loop over ALL angles in ALL replicas
	for (int i = threadIdx.x+blockIdx.x*blockDim.x; i<numBondAngles; i+=blockDim.x*gridDim.x) {
		int atom1 = bondAngleList_d[2*i].x;
		int atom2 = bondAngleList_d[2*i].y;
		int atom3 = bondAngleList_d[2*i].z;
		int atom4 = bondAngleList_d[2*i].w;

		int angleInd1 = bondAngleList_d[2*i+1].x;
		int bondInd   = bondAngleList_d[2*i+1].y;
		int angleInd2 = bondAngleList_d[2*i+1].z;

		computeBondAngle(tableAngle[ angleInd1 ], tableBond[ bondInd ], tableAngle[ angleInd2 ], sys, force, pos, atom1, atom2, atom3, atom4, energy, get_energy);
	}
}

__global__
void computeProductPotentials(Vector3* force,
			      Vector3* __restrict__ pos,
			      BaseGrid* __restrict__ sys,
			      int numProductPotentials,
			      int* __restrict__ productPotentialParticles,
			      SimplePotential* __restrict__ potentialList,
			      uint2* __restrict__ productPotential_list,
			      unsigned short* __restrict__ productCount,
			      float* energy, bool get_energy) {
    /*
      productPotential_list[i].x : index of first potential in potentialList for i_th productPotential
      productPotential_list[i].y : index of first atom in productPotentialParticles for i_th productPotential

      for three potentials, angle, bond, angle, we would have the following atomic indices in productPotentialParticles:
        pot1 : productPotential_list[i].y, productPotential_list[i].y + 1 , productPotential_list[i].y + 2
        pot2 : productPotential_list[i].y + 3, productPotential_list[i].y + 4
	and
        pot3 : productPotential_list[i].y + 5, productPotential_list[i].y + 6, productPotential_list[i].y + 7

      productCount[i] : number of potentials in the i_th productPotential
    */

    // CRAPPY NAIVE IMPLEMENTATION
    constexpr int MAX_XPOTS = 4;
    float2 energy_and_deriv[MAX_XPOTS];
    float tmp_force;

    for (int i = threadIdx.x+blockIdx.x*blockDim.x; i<numProductPotentials; i+=blockDim.x*gridDim.x) {
	unsigned short num_pots = productCount[i];

	unsigned int part_idx = productPotential_list[i].y;
#pragma unroll
	for (unsigned short int j = 0; j < MAX_XPOTS; ++j) {
	    if (j == num_pots) break;
	    SimplePotential& p = potentialList[ productPotential_list[i].x + j ];

	    // Hidden branch divergence in compute_value => sort potentials by type before running kernel
	    float tmp = p.compute_value(pos,sys, &productPotentialParticles[part_idx]);
	    energy_and_deriv[j] = p.compute_energy_and_deriv(tmp);
	    part_idx += p.type==BOND? 2: p.type==ANGLE? 3: 4;
	}

	part_idx = productPotential_list[i].y;
#pragma unroll
	for (unsigned short int j = 0; j < MAX_XPOTS; ++j) {
	    if (j == num_pots) break;
	    tmp_force = energy_and_deriv[j].y;
#pragma unroll
	    for (unsigned short int k = 0; k < MAX_XPOTS; ++k) {
		if (k == num_pots) break;
		if (j == k) continue;
		tmp_force *= energy_and_deriv[k].x;
	    }
	    SimplePotential& p = potentialList[ productPotential_list[i].x + j ];
	    if (tmp_force != 0) {
		// TODO add energy
		p.apply_force(pos,sys, force, &productPotentialParticles[part_idx], tmp_force);
	    }
	    part_idx += p.type==BOND? 2: p.type==ANGLE? 3: 4;
	}
    }
}


__global__
void computeDihedrals(Vector3 force[], Vector3 pos[],
											Dihedral dihedrals[],
											TabulatedDihedralPotential* tableDihedral[],
											int numDihedrals, int num, BaseGrid* sys, float g_energies[],
											bool get_energy) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	// float energy_local = 0.0f;
	Vector3 force_local(0.0f);

	if (i < numDihedrals) {
		// RBTODO: optimize
		Dihedral& d = dihedrals[i];

		const Vector3 ab = sys->wrapDiff( pos[d.ind1] - pos[d.ind2] );
		const Vector3 bc = sys->wrapDiff( pos[d.ind2] - pos[d.ind3] );
		const Vector3 cd = sys->wrapDiff( pos[d.ind3] - pos[d.ind4] );
		
		//const float distab = ab.length();
		const float distbc = bc.length();
		//const float distcd = cd.length();
	
		Vector3 crossABC = ab.cross(bc);
		Vector3 crossBCD = bc.cross(cd);
		Vector3 crossX = bc.cross(crossABC);

		const float cos_phi = crossABC.dot(crossBCD) / (crossABC.length() * crossBCD.length());
		const float sin_phi = crossX.dot(crossBCD) / (crossX.length() * crossBCD.length());
		
		const float angle = -atan2(sin_phi, cos_phi);

	
		Vector3 f1, f2, f3; // forces
		f1 = -distbc * crossABC.rLength2() * crossABC;
		f3 = -distbc * crossBCD.rLength2() * crossBCD;
		f2 = -(ab.dot(bc) * bc.rLength2()) * f1 - (bc.dot(cd) * bc.rLength2()) * f3;
	
		// Shift "angle" by "PI" since    -PI < dihedral < PI
		// And our tabulated potential data: 0 < angle < 2 PI
		float& dangleInv = tableDihedral[d.tabFileIndex]->angle_step_inv;
		float t = (angle + BD_PI) * dangleInv;
		int home = (int) floorf(t);
		t = t - home;

		int size = tableDihedral[d.tabFileIndex]->size;
		home = home % size;
		int home1 = (home + 1) % size;

		//================================================
		// Linear interpolation
		float * pot = tableDihedral[d.tabFileIndex]->pot;
		float U0 = pot[home];       // Potential
		float dU = pot[home1] - U0; // Change in potential
		
		float energy = dU * t + U0;
		float f = -dU * dangleInv;
		//================================================
		// TODO: add an option for cubic interpolation [Probably not]

		if (crossABC.rLength() > 1.0f || crossBCD.rLength() > 1.0f)
			// avoid singularity when one angle is straight 
			f = 0.0f;

		f1 *= f;
		f2 *= f;
		f3 *= f;

		atomicAdd( &force[d.ind1], f1 );
		atomicAdd( &force[d.ind2], f2-f1 );
		atomicAdd( &force[d.ind3], f3-f2 );
		atomicAdd( &force[d.ind4], -f3 );

		if (get_energy) {
			atomicAdd( &g_energies[d.ind1], energy );
			atomicAdd( &g_energies[d.ind2], energy );
			atomicAdd( &g_energies[d.ind3], energy );
			atomicAdd( &g_energies[d.ind4], energy );
		}
	}
}


    // void computeTabulatedDihedrals(Vector3* __restrict__ force, Vector3* __restrict__ pos, int num,
    // 			    int numParts, BaseGrid* __restrict__ sys, int4* __restrict__ dihedralList_d,
    // 			    int* __restrict__ dihedralPotList_d,
    // 			    int numDihedrals, int numReplicas, float* __restrict g_energies,
    // 			    bool get_energy, TabulatedDihedralPotential** __restrict__ tableDihedral) {

__global__
void computeTabulatedDihedrals(Vector3* force, const Vector3* __restrict__ pos,
			       const BaseGrid* __restrict__ sys,
			       int numDihedrals, const int4* const __restrict__ dihedralList_d,
			       const int* __restrict__ dihedralPotList_d, TabulatedDihedralPotential** tableDihedral, float* energy, bool get_energy) {

	// int currDihedral = blockIdx.x * blockDim.x + threadIdx.x; // first particle ID

    // Loop over ALL dihedrals in ALL replicas
	for (int i = threadIdx.x+blockIdx.x*blockDim.x; i < numDihedrals; i+=blockDim.x*gridDim.x) {
		const int4& ids = dihedralList_d[i];
		const int& id = dihedralPotList_d[i];
		computeDihedral(tableDihedral[ id ], sys, force, pos, ids.x, ids.y, ids.z, ids.w, energy, get_energy);

	// if (get_energy)
	// {
	//     //TODO: clarification on energy computation needed, consider changing.
	//     atomicAdd( &g_energies[i], energy_local);
	//     //atomicAdd( &g_energies[j], energy_local);
	// }
    }
}

__global__
void computeHarmonicRestraints(Vector3* force, const Vector3* __restrict__ pos,
			       const BaseGrid* __restrict__ sys,
			       int numRestraints, const int2* const __restrict__ restraintList,
			       const Vector3* __restrict__ r0, const float* __restrict__ k) {

    // Loop over ALL dihedrals in ALL replicas
    for (int i = threadIdx.x+blockIdx.x*blockDim.x; i < numRestraints; i+=blockDim.x*gridDim.x) {
	const int& id = restraintList[i].x;
	const int& rid = restraintList[i].y;

	const Vector3 dr = sys->wrapDiff(pos[id]-r0[rid]);
	Vector3 f = -k[rid]*dr;
	atomicAdd( &force[ id ], f );
    }
}
