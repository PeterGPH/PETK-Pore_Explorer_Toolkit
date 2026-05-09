// CellDecomposition.h (2013)
// contains CellDecomposition class and two related structs
//
// Authors: Terrance Howard <howard33@illinois.edu>,
//          Justin Dufresne <jdufres1@friars.providence.edu>
//
// "When I wrote this, only God and myself understood what I was thinking.
//  Now, only God knows."

#ifndef CELL_DECOMPOSITION_H
#define CELL_DECOMPOSITION_H

#include <vector>
#include <algorithm> // std::sort

#include <thrust/sort.h>
#include <thrust/device_ptr.h>
#include <thrust/partition.h>

#include <cuda.h>	 	// cudaMalloc, cudaMemcpy
#include <vector_types.h>	// int3

#include "useful.h" // Vector3, Matrix3
#include "BaseGrid.h"

#if defined(__CUDACC__) // NVCC
   #define MY_ALIGN(n) __align__(n)
#elif defined(__GNUC__) // GCC
  #define MY_ALIGN(n) __attribute__((aligned(n)))
#elif defined(_MSC_VER) // MSVC
  #define MY_ALIGN(n) __declspec(align(n))
#else
  #error "Please provide a definition for MY_ALIGN macro for your host compiler!"
#endif

class CellDecomposition : public BaseGrid {
public:
	// range_t
	// Contains first and last exclusive indices in cells array.
	struct range_t {
	public:
		HOST DEVICE
		inline range_t() : first(-1), last(-1) { }

		HOST DEVICE
		inline range_t(int first, int last) : first(first), last(last) { }

	public:
		int first, last; // [first, last)
	};

	// cell_t
	// Contains replica id, particle id, cell id, and position of cell
	struct cell_t {
	public:
		HOST DEVICE
		inline cell_t() : particle(-1), id(-1) { }

		HOST DEVICE
		inline cell_t(int particle, int id, const int3& r, int repID) :
				particle(particle), repID(repID), id(id), pos(r) { }

		HOST DEVICE
		inline bool operator<(const cell_t& p) const {
				if (repID != p.repID) return repID < p.repID;
				if (id != p.id) return id < p.id;
				return particle < p.particle;
		}
	
	public:
		int particle; // id of particle
		int repID;
		int id; // location in CellDecomposition's cells array
		int3 pos; // position of cell in grid
	};

public:
	CellDecomposition(Matrix3 box, Vector3 origin, float cutoff, int numReplicas);

	// Place particles in cells and create a range for each cell.
	// Decompose on the GPU.
	void decompose_d(Vector3 *pos_d, size_t num);

	// copyToCUDA
	// Return a copy of CellDecomposition to GPU.
	CellDecomposition* copyToCUDA();

	HOST DEVICE
	inline float getCutoff() const { return cutoff; }

	HOST DEVICE
	inline size_t size() const { return numCells; }

	HOST DEVICE
	inline const cell_t& getCell(int ind) const { return cells[ind]; }

	HOST DEVICE
	inline const cell_t& getCellForParticle(int particle) const {
		return unsorted_cells[particle];
	}

	// Return cell array
	HOST DEVICE
	inline const cell_t* getCells() const {
		return cells;
	}
        //Han-Yi Chou
        HOST DEVICE
        inline const cell_t* getCells_d() const {
                return cells_d;
        }

	/*
	HOST DEVICE
	inline const range_t& getRange(const cell_t& c) const {
		return ranges_d[c.id + c.repID * numCells];
	}
	*/

	HOST DEVICE
	inline const range_t& getRange(int ind, int repID) const {
		return ranges_d[ind + repID * numCells];
	}

	HOST DEVICE
	inline int getCellID(const Vector3 &r0) const {
		const Vector3 r = r0 - origin;
		//const int x = int(r.x / cutoff);
		//const int y = int(r.y / cutoff);
		//const int z = int(r.z / cutoff);
		const int x = floorf(r.x / cutoff);
                const int y = floorf(r.y / cutoff);
                const int z = floorf(r.z / cutoff);
		return getCellID(x, y, z, nCells);
	}

	HOST DEVICE
	static inline int getCellID(const Vector3& r0, const Vector3& origin,
															float cutoff, int3 nCells) {
		const Vector3 r = r0 - origin;
		//const int x = int(r.x / cutoff);
		//const int y = int(r.y / cutoff);
		//const int z = int(r.z / cutoff);
		const int x = floorf(r.x / cutoff);
                const int y = floorf(r.y / cutoff);
                const int z = floorf(r.z / cutoff);
		return getCellID(x, y, z, nCells);
	}

	// Return position of cell in grid.
	HOST DEVICE
	inline int3 getCellPos(int id) const {
		return getCellPos(id, nCells);
	}

	HOST DEVICE
	static inline int3 getCellPos(int id, int3 nCells) {
		int3 p;
		p.z = id % nCells.z;
		p.y = (id / nCells.z) % nCells.y;
		p.x = id / (nCells.z * nCells.y);
		return p;
	}

	// Return ID of cell in position (i, j, k) relative to c.
	// Return -1 if wrapping to an adjacent cell.
	HOST DEVICE
	inline int getNeighborID(const cell_t& c, int i, int j, int k) const {
		if (i == 0 and j == 0 and k == 0)
			return c.id;
		int u = i + c.pos.x;
		int v = j + c.pos.y;
		int w = k + c.pos.z;
		if (nCells.x == 1 and u != 0) return -1;
		if (nCells.y == 1 and v != 0) return -1;
		if (nCells.z == 1 and w != 0) return -1;
		if (nCells.x == 2 and (u < 0 || u > 1)) return -1;
		if (nCells.y == 2 and (v < 0 || v > 1)) return -1;
		if (nCells.z == 2 and (w < 0 || w > 1)) return -1;
		return getCellID(u, v, w, nCells);
	}
/*
        HOST DEVICE
inline int getNeighborID(int idx, int dx, int dy, int dz) const
{
    if(dx == 0 and dy == 0 and dz == 0)
        return idx;
    int idx_z = idx % nCells.z;
    int idx_y = idx / nCells.z % nCells.y;
    int idx_x = idx / (nCells.z * nCells.y);

    int u = (dx + idx_x + nCells.x) % nCells.x;
    int v = (dy + idx_y + nCells.y) % nCells.y;
    int w = (dz + idx_z + nCells.z) % nCells.z;
    if (nCells.x == 1 and u != 0) return -1;
    if (nCells.y == 1 and v != 0) return -1;
    if (nCells.z == 1 and w != 0) return -1;
    if (nCells.x == 2 and (u < 0 || u > 1)) return -1;
    if (nCells.y == 2 and (v < 0 || v > 1)) return -1;
    if (nCells.z == 2 and (w < 0 || w > 1)) return -1;
    return getCellID(u, v, w, nCells);
}
*/
public:
	int3 nCells;

private:
	// Wrap an integer with inclusive lower and upper bounds.
	HOST DEVICE
	static inline int wrapInt(int k, int lower, int upper) {
		int range = upper - lower + 1;
		if (k < lower)
			k += range * ((lower - k) / range + 1);
		return lower + (k - lower) % range;
	}
	
	// Calculate a cell's id given a position in the grid.
	HOST DEVICE
	static inline int getCellID(int i, int j, int k, int3 nCells) {
		i = wrapInt(i, 0, nCells.x - 1);
		j = wrapInt(j, 0, nCells.y - 1);
		k = wrapInt(k, 0, nCells.z - 1);
		return k + nCells.z * (j + (nCells.y * i));
	}


private:
	static const unsigned int NUM_THREADS = 256;

	int numCells;
	int numReplicas;

	cell_t* cells;
	cell_t* cells_d;
	cell_t* unsorted_cells;
	cell_t* unsorted_cells_d;
	range_t* ranges;
	range_t* ranges_d;

	float cutoff;

	// build_ranges
	// @param	number of particles
	// used by decompose()
	void build_ranges(size_t num);

};

#endif
