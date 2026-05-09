// BrownianParticleType.h (2002)
// Contains BrownianParticleType and TypeDecomposition classes
//
// Author: Jeff Comer <jcomer2@illinois.edu>
// Edited (2013) by Terrance Howard <howard33@illinois.edu>,
//                  Justin Dufresne <jdufres1@friars.providence.edu>

#ifndef BROWNIANPARTICLETYPE_H
#define BROWNIANPARTICLETYPE_H

#include <vector>

#include "Reservoir.h"
#include "BaseGrid.h"
#include "CellDecomposition.h"

// Stores particle type's potential grid and other information
class BrownianParticleType {
	private:
		// clear
		// Deletes all members
		void clear();

		// copy
		// Copies all members
		// @param BrownianParticleType to copy
		void copy(const BrownianParticleType& src);

	public:
		BrownianParticleType(const String& name = "") :
				name(name), num(0),
				diffusion(0.0f), radius(1.0f), charge(0.0f), eps(0.0f), meanPmf(NULL),
				numPartGridFiles(-1), reservoir(NULL), pmf(NULL), pmf_scale(NULL), pmf_boundary_conditions(NULL),
				diffusionGrid(NULL),
				forceXGrid(NULL), forceYGrid(NULL), forceZGrid(NULL){ }

		BrownianParticleType(const BrownianParticleType& src) { copy(src); }

		~BrownianParticleType() { clear(); }

		BrownianParticleType& operator=(const BrownianParticleType& src);

    void set_boundary_conditions( int num, BoundaryCondition* bcs ) {
	if (num <= 0) return;

	if (pmf_boundary_conditions != NULL) {
	    delete[] pmf_boundary_conditions;
	}

	pmf_boundary_conditions = new BoundaryCondition[num];
	for (int i=0; i < num; ++i) {
	    pmf_boundary_conditions[i] = bcs[i];
	}
    }

		// crop
		// Crops all BaseGrid members
		// @param  boundries to crop to (x0, y0, z0) -> (x1, y1, z1);
		//         whether to change the origin
		// @return success of function (if false nothing was done)
		//bool crop(int x0, int y0, int z0, int x1, int y1, int z1, bool keep_origin);

public:
		String name;
		int num; // number of particles of this type
                float mass; // mass of brownian particles Han-Yi Chou
                Vector3 transDamping; // translational damping coefficient Han-Yi Chou
		float diffusion;
		float radius;
		float charge;
		float eps;
		//float meanPmf;
		float *meanPmf;
                int   numPartGridFiles;
                float mu; //for Nose-Hoover Langevin dynamics

		Reservoir* reservoir;
		BaseGrid** pmf;
    float* pmf_scale;
    BoundaryCondition* pmf_boundary_conditions;
		BaseGrid* diffusionGrid;
		BaseGrid* forceXGrid;
		BaseGrid* forceYGrid;
		BaseGrid* forceZGrid;
};

/*
// Spatially decomposes BrownianParticleTypes
class TypeDecomposition {
	private:
		size_t num_cells_;
		size_t num_parts_;
		std::vector<BrownianParticleType*> parts_; // 2D array; parts_[cell][particle_type]

		TypeDecomposition() {}

	public:
		TypeDecomposition(const CellDecomposition &decomp,
				const BrownianParticleType *parts, size_t num_parts);

		~TypeDecomposition();

		// Getters
		const BrownianParticleType* at(size_t i) const;
		const BrownianParticleType* operator[](size_t i) const { return at(i); }

		const std::vector<BrownianParticleType*>& parts() const { return parts_; }

		int num_cells() const { return num_cells_; }
		int num_parts() const { return num_parts_; }
};
*/
#endif
