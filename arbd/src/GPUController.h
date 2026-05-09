// provides interface between main CPU loop and various GPUs
//   -- holds data for each GPU

#pragma once
#include "useful.h"

class GPUcontroller {
public:
	GPUcontroller(const Configuration& c, const long int randomSeed,
			bool debug, int numReplicas = 0);
	~GPUcontroller();

	static bool DEBUG;

private:  

	void copyToCUDA();

	
private:
	const Configuration& conf;
	int numReplicas;

	// Integrator variables
	BaseGrid* sys;
	ComputeForce* internal;
	Vector3* forceInternal;

	// CUDA device variables
	Vector3 *pos_d, *forceInternal_d, *force_d;
	int *type_d;
	BrownianParticleType **part_d;
	BaseGrid *sys_d, *kTGrid_d;
	Random *randoGen_d;
	Bond* bonds_d;
	int2* bondMap_d;
	Exclude* excludes_d;
	int2* excludeMap_d;
	Angle* angles_d;
	Dihedral* dihedrals_d;
	
}
	
	
