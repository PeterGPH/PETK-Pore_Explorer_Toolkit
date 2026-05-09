// CellDecomposition.cu
//
// Terrance Howard <heyterrance@gmail.com>

#include "CellDecomposition.h"

// *****************************************************************************
// Error Check

inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true) {
	if (code != cudaSuccess) {
		fprintf(stderr,"CUDA Error: %s %s %d\n",
						cudaGetErrorString(code), file, line);
		if (abort) exit(code);
	}
}

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

#define gpuKernelCheck() {kernelCheck( __FILE__, __LINE__); }
inline void kernelCheck(const char* file, int line)
{
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        std::fprintf(stderr,"Error: %s in %s %d\n", cudaGetErrorString(err),file, line);
        assert(1==2);
    }
    //gpuErrchk(cudaDeviceSynchronize());
}

// *****************************************************************************
// CUDA Kernel Definitions

__global__ void decomposeKernel(Vector3 pos[],
																CellDecomposition::cell_t cells[],
																Vector3 origin, float cutoff,
																int3 nCells, size_t num, int numReplicas);

__global__
void make_rangesKernel(CellDecomposition::cell_t cells[], int tmp[],
											 size_t num, int numCells, int numReplicas);

__global__
void bind_rangesKernel(CellDecomposition::range_t ranges[], int tmp[],
											 int numCells, int numReplicas);

// *****************************************************************************
// CellDecomposition Implementations

CellDecomposition::CellDecomposition(Matrix3 box, Vector3 origin,
																		 float cutoff, int numReplicas) :
		BaseGrid(box, origin, cutoff), cutoff(cutoff), numReplicas(numReplicas),
		cells(NULL), cells_d(NULL), unsorted_cells(NULL), unsorted_cells_d(NULL),
		ranges(NULL), ranges_d(NULL) {
	const Vector3 dim = getExtent();
	nCells.x = int((dim.x - 1) / cutoff) + 1;
	nCells.y = int((dim.y - 1) / cutoff) + 1;
	nCells.z = int((dim.z - 1) / cutoff) + 1;
	numCells = nCells.x * nCells.y * nCells.z;
	printf("Created Cell Decomposition (%lu, %lu, %lu)\n",
			nCells.x, nCells.y, nCells.z);
}


CellDecomposition* CellDecomposition::copyToCUDA() {
	cell_t* tmp_cells = this->cells;
	cell_t* tmp_unsorted = this->unsorted_cells;

	this->cells = this->cells_d;
	this->unsorted_cells = this->unsorted_cells_d;

	const size_t sz = sizeof(CellDecomposition);
	CellDecomposition *c_d = NULL;
	gpuErrchk(cudaMalloc(&c_d, sz));
	gpuErrchk(cudaMemcpy(c_d, this, sz, cudaMemcpyHostToDevice));

	this->cells = tmp_cells;
	this->unsorted_cells = tmp_unsorted;

	return c_d;
}

void CellDecomposition::decompose_d(Vector3 pos_d[], size_t num) {
	const size_t cells_sz = sizeof(cell_t) * num * numReplicas;
	const size_t numCellRep = numCells * numReplicas;

	if (cells_d == NULL) {
		gpuErrchk(cudaMalloc(&cells_d, cells_sz));
		gpuErrchk(cudaMalloc(&unsorted_cells_d, cells_sz));
		gpuErrchk(cudaMalloc(&ranges_d, sizeof(range_t) * numCellRep));
		unsorted_cells = new cell_t[num * numReplicas];
		cells = new cell_t[num * numReplicas];
		ranges = new range_t[numCellRep];
	}

	// Pair particles with cells.
	size_t nBlocks = (num * numReplicas) / NUM_THREADS + 1;
	decomposeKernel<<< nBlocks, NUM_THREADS >>>(pos_d, cells_d, origin, cutoff,
																							nCells, num, numReplicas);
	gpuErrchk(cudaDeviceSynchronize());
	gpuErrchk(cudaMemcpy(unsorted_cells_d, cells_d, cells_sz,
											 cudaMemcpyDeviceToDevice));
	gpuErrchk(cudaMemcpyAsync(unsorted_cells, unsorted_cells_d, cells_sz,
														cudaMemcpyDeviceToHost));

	// Sort cells.
	thrust::device_ptr<cell_t> c_d(cells_d);
	thrust::sort(c_d, c_d + num * numReplicas);
	gpuErrchk(cudaMemcpyAsync(cells, cells_d, cells_sz, cudaMemcpyDeviceToHost));
	//Han-Yi Chou
        //gpuErrchk(cudaMemcpy(cells, cells_d, cells_sz, cudaMemcpyDeviceToHost));
	const size_t nMax = std::max(2lu * numCells, num);
	nBlocks = (nMax * numReplicas) / NUM_THREADS + 1;

	// Create ranges for cells.
	int* temp_ranges = NULL;
	gpuErrchk(cudaMalloc(&temp_ranges, 2 * sizeof(int) * numCellRep));
	gpuErrchk(cudaMemset(temp_ranges, -1, 2 * sizeof(int) * numCellRep));
	make_rangesKernel<<< nBlocks, NUM_THREADS >>>(cells_d, temp_ranges,
																								num, numCells, numReplicas);
	gpuErrchk( cudaDeviceSynchronize() );

	// Copy temp_ranges to ranges_d
	bind_rangesKernel<<< nBlocks, NUM_THREADS >>>(ranges_d, temp_ranges,
																								numCells, numReplicas);
	gpuErrchk(cudaMemcpy(ranges, ranges_d, numCellRep, cudaMemcpyDeviceToHost));
	gpuErrchk( cudaFree(temp_ranges) );


}

// *****************************************************************************
// CUDA Kernels

__global__
void decomposeKernel(Vector3 *pos, CellDecomposition::cell_t *cells,
										 Vector3 origin, float cutoff, int3 nCells,
										 size_t num, int numReplicas) {
	const size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < num * numReplicas) {
		const int repID =  idx / num;
		const Vector3& p = pos[idx];
		const int id = CellDecomposition::getCellID(p, origin, cutoff, nCells);
		const int3 r = CellDecomposition::getCellPos(id, nCells);
		cells[idx] = CellDecomposition::cell_t(idx, id, r, repID);
	}
}

__global__
void make_rangesKernel(CellDecomposition::cell_t cells[], int tmp[],
											 size_t num, int numCells, int numReplicas) {
	const size_t idx = blockIdx.x * blockDim.x + threadIdx.x; 

	if (idx < num * numReplicas) {
		const int repID = cells[idx].repID;
		assert(repID == idx/num);
		const int cellID = cells[idx].id + repID * numCells; // cellID in tmp array

		// Get positions in tmp array.
		const int first = cellID * 2;
		const int last = first + 1;
		const int particle = idx % num;

		if (particle == 0)
			tmp[first] = idx;

		if (particle == num - 1)
			tmp[last] = idx + 1;

		const int prev_id = idx - 1;
		if (prev_id >= 0
		    and cells[prev_id].repID == repID
		    and cellID != cells[prev_id].id + repID * numCells)
			tmp[first] = idx;
		
		const int next_id = idx + 1;
		if (next_id < num * numReplicas
		    and cells[next_id].repID == repID
		    and cellID != cells[next_id].id + repID * numCells)
			tmp[last] = idx + 1;
	}
}

__global__
void bind_rangesKernel(CellDecomposition::range_t ranges[], int tmp[],
											 int numCells, int numReplicas) {
        const size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < numCells * numReplicas)
		ranges[idx] = CellDecomposition::range_t(tmp[2*idx], tmp[2*idx+1]);
	/* Print range of each cell. Skip over empty cells
	__syncthreads();
	if (idx == 0) {
		for (int i = 0; i < numCells * numReplicas; ++i) {
			if (ranges[i].first == -1 and ranges[i].last == -1) continue;
			printf("cell %d : [%d, %d)\n", i, ranges[i].first, ranges[i].last);
		}
	}
	// */
}
