// GrandBrownTown.cuh
//
// Terrance Howard <heyterrance@gmail.com>
#pragma once
//#define MDSTEP
#define Unit1 4.18679994e4
#define Unit2 2.046167337
//#define Debug
#include "CudaUtil.cuh"
#include "RigidBodyType.h"
#include "RigidBodyGrid.h"

__device__
Vector3 step(Vector3& r0, float kTlocal, Vector3 force, float diffusion, Vector3 diffGrad,
						 float timestep, BaseGrid *sys, Random *randoGen, int num);

inline __device__
ForceEnergy compute_position_dependent_force(
    const Vector3* __restrict__ pos, Vector3* __restrict__ forceInternal,
    const int* __restrict__ type, BrownianParticleType** part,
    const float electricField, const int scheme, const int idx)
{
    int t = type[idx];
    Vector3 r0 = pos[idx];
    const BrownianParticleType& pt = *part[t];
    Vector3 forceExternal = Vector3(0.0f, 0.0f, pt.charge * electricField);

    ForceEnergy fe(0.f, 0.f);
    for(int i = 0; i < pt.numPartGridFiles; ++i)
    {
	ForceEnergy tmp(0.f, 0.f);
	if(!scheme) {
	    BoundaryCondition bc = pt.pmf_boundary_conditions[i];
	    INTERPOLATE_FORCE(tmp, pt.pmf[i]->interpolateForceDLinearly, bc, r0)
		} else
	    tmp = pt.pmf[i]->interpolateForceD(r0);
	fe.f += tmp.f * pt.pmf_scale[i];
	fe.e += tmp.e * pt.pmf_scale[i];
    }
    // if(get_energy)
    // 	energy[idx] += fe.e;

#ifndef FORCEGRIDOFF
    // Add a force defined via 3D FORCE maps (not 3D potential maps)
    if(!scheme)
    {
	if (pt.forceXGrid != NULL) fe.f.x += pt.forceXGrid->interpolatePotentialLinearly(r0);
	if (pt.forceYGrid != NULL) fe.f.y += pt.forceYGrid->interpolatePotentialLinearly(r0);
	if (pt.forceZGrid != NULL) fe.f.z += pt.forceZGrid->interpolatePotentialLinearly(r0);
    }
    else
    {
	if (pt.forceXGrid != NULL) fe.f.x += pt.forceXGrid->interpolatePotential(r0);
	if (pt.forceYGrid != NULL) fe.f.y += pt.forceYGrid->interpolatePotential(r0);
	if (pt.forceZGrid != NULL) fe.f.z += pt.forceZGrid->interpolatePotential(r0);
    }
#endif
    fe.f = fe.f + forceExternal;
    return fe;
}


////The kernel is for Nose-Hoover Langevin dynamics
template <bool need_position_dep_force = false>
__global__ void 
updateKernelNoseHooverLangevin(Vector3* __restrict__ pos, Vector3* __restrict__ momentum, float* random, 
                               Vector3* __restrict__ forceInternal, int type[], BrownianParticleType* part[], 
                               float kT, BaseGrid* kTGrid, float electricField, int tGridLength, float timestep, 
                               int num, int num_rb_attached_particles, BaseGrid* sys, Random* randoGen, int numReplicas, int scheme)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < num * numReplicas)
    {
	idx = (idx % num) + (idx/num) * (num+num_rb_attached_particles);

	ForceEnergy fe;
	if (need_position_dep_force) {
	    fe = compute_position_dependent_force(
		pos, forceInternal, type, part, electricField, scheme, idx );
	}
	
        int t = type[idx];
        Vector3 r0  = pos[idx];
        Vector3 p0  = momentum[idx];
        float   ran = random[idx];

        const BrownianParticleType& pt = *part[t];

	if (need_position_dep_force) {
	    forceInternal[idx] = forceInternal[idx] + fe.f;
	}
        Vector3 force = forceInternal[idx];
        #ifdef Debug
        forceInternal[idx] = -force;
        #endif
        // Get local kT value
        float kTlocal;
        if(!scheme)
            kTlocal = (tGridLength == 0) ? kT : kTGrid->interpolatePotentialLinearly(r0); /* periodic */
        else
            kTlocal = (tGridLength == 0) ? kT : kTGrid->interpolatePotential(r0); /* periodic */

        // Update the particle's position using the calculated values for time, force, etc.
        float mass  = pt.mass;
        float mu    = pt.mu;
        Vector3 gamma   = pt.transDamping;
        float rando = (*randoGen).gaussian(idx, num * numReplicas);

        float tmp   = sqrtf(kTlocal * (1.f - expf(-2.f * gamma.x * timestep)) / mu);

        if (pt.diffusionGrid != NULL)
        {

            Vector3 gridCenter = pt.diffusionGrid->origin +
            pt.diffusionGrid->basis.transform( Vector3(0.5*pt.diffusionGrid->nx, 0.5*pt.diffusionGrid->ny,
                                                       0.5*pt.diffusionGrid->nz));
            Vector3 p2 = r0 - gridCenter;
            p2 = sys->wrapDiff( p2 ) + gridCenter;
            ForceEnergy diff;
            if(!scheme)
                diff = pt.diffusionGrid->interpolateForceDLinearly<periodic>(p2);
            else
                diff = pt.diffusionGrid->interpolateForceD(p2);
            gamma = Vector3(kTlocal / (mass * diff.e));
        }

        #ifdef MDSTEP
        force = Vector3(-r0.x, -r0.y, -r0.z);
        #endif

        p0  = p0  + 0.5f * timestep * force * Unit1;

        r0  = r0  + 0.5f * timestep * p0 * 1e4 / mass;
        //r0 = sys->wrap(r0);

        ran = ran + 0.5f * (p0.length2() / mass * 0.238845899f - 3.f * kTlocal) / mu;

        p0  = expf(-0.5f * ran) * p0;

        ran = expf(-gamma.x * timestep) * ran + tmp * rando;

        p0  = expf(-0.5f * ran) * p0;

        ran = ran + 0.5f * (p0.length2() / mass * 0.238845899 - 3.f * kTlocal) / mu;
 
        r0  = r0  + 0.5f * timestep * p0 * 1e4 / mass;
        r0 = sys->wrap(r0);

        pos[idx] = r0;
        momentum[idx] = p0;
        random[idx]= ran;     
    }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//This is the kernel for BAOAB (velocity verlet algorithm) which is symplectic. For more information,
//please infer http://www.MolecularDynamics.info
//The original BBK kernel is no longer used since the random numbers should be reused 
//which is not possible in GPU code.
//Han-Yi Chou
template <bool need_position_dep_force = false>
__global__ void updateKernelBAOAB(Vector3* pos, Vector3* momentum, Vector3* __restrict__ forceInternal,
                                  int type[], BrownianParticleType* part[], float kT, BaseGrid* kTGrid, 
                                  float electricField,int tGridLength, float timestep,
				  int num, int num_rb_attached_particles, BaseGrid* sys,
                                  Random* randoGen, int numReplicas, int scheme)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;


    if (idx < num * numReplicas)
    {
	idx = (idx % num) + (idx/num) * (num+num_rb_attached_particles);

	ForceEnergy fe;
	if (need_position_dep_force) {
	    fe = compute_position_dependent_force(
		pos, forceInternal, type, part, electricField, scheme, idx );
	}
	// if (get_energy) energy[idx] += fe.e;

        int t = type[idx];
        Vector3 r0 = pos[idx];
        Vector3 p0 = momentum[idx];
        const BrownianParticleType& pt = *part[t];

	if (need_position_dep_force) { // only needed at ts == 1
	    forceInternal[idx] = forceInternal[idx] + fe.f;
	}
	Vector3 force = forceInternal[idx];
	
#ifdef Debug
        forceInternal[idx] = -force;
#endif

        // Get local kT value
        float kTlocal;
        if(!scheme)
            kTlocal = (tGridLength == 0) ? kT : kTGrid->interpolatePotentialLinearly(r0); /* periodic */
        else
            kTlocal = (tGridLength == 0) ? kT : kTGrid->interpolatePotential(r0); /* periodic */

        // Update the particle's position using the calculated values for time, force, etc.
        float mass      = pt.mass;
        Vector3 gamma   = pt.transDamping;
        Vector3 rando = (*randoGen).gaussian_vector(idx, num * numReplicas);
        //printf("%f %f %f\n", rando.x, rando.y, rando.z);
        if (pt.diffusionGrid != NULL) 
        {

            Vector3 gridCenter = pt.diffusionGrid->origin +
            pt.diffusionGrid->basis.transform( Vector3(0.5*pt.diffusionGrid->nx, 0.5*pt.diffusionGrid->ny,
                                                       0.5*pt.diffusionGrid->nz));
            Vector3 p2 = r0 - gridCenter;
            p2 = sys->wrapDiff( p2 ) + gridCenter;
            ForceEnergy diff;
            if(!scheme)
                diff = pt.diffusionGrid->interpolateForceDLinearly<periodic>(p2);
            else
                diff = pt.diffusionGrid->interpolateForceD(p2);
            gamma = Vector3(kTlocal / (mass * diff.e));
        }

        #ifdef MDSTEP
        force = Vector3(-r0.x, -r0.y, -r0.z);
        #endif

        p0 = p0 + 0.5f * timestep * force * Unit1;
        r0 = r0 + 0.5f * timestep / mass * p0 * 1e4;

        p0.x = expf(-timestep * gamma.x) * p0.x + sqrtf(mass * kTlocal * (1.f-expf(-2.f*timestep*gamma.x))) * rando.x * Unit2;
        p0.y = expf(-timestep * gamma.y) * p0.y + sqrtf(mass * kTlocal * (1.f-expf(-2.f*timestep*gamma.y))) * rando.y * Unit2;
        p0.z = expf(-timestep * gamma.z) * p0.z + sqrtf(mass * kTlocal * (1.f-expf(-2.f*timestep*gamma.z))) * rando.z * Unit2;

        r0 = r0 + 0.5f * timestep * p0 * 1e4 / mass;
        r0 = sys->wrap(r0);
        
        pos[idx]      = r0;
        momentum[idx] = p0;

        //if(idx == 0)
          //  printf("%f %f %f\n", pos[idx].x,pos[idx].y,pos[idx].z);
    }
}

//update momentum in the last step of BAOAB integrator for the Langevin dynamics. Han-Yi Chou
__global__ void LastUpdateKernelBAOAB(Vector3* pos,Vector3* momentum, Vector3* __restrict__ forceInternal,
                                      int type[], BrownianParticleType* part[], float kT, BaseGrid* kTGrid, 
                                      float electricField, int tGridLength, float timestep, int num, int num_rb_attached_particles,
                                      BaseGrid* sys, Random* randoGen, int numReplicas, float* __restrict__ energy, bool get_energy,int scheme)
{
    int idx  = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);

    if (idx < num * numReplicas)
    {
	idx = (idx % num) + (idx/num) * (num+num_rb_attached_particles);

	ForceEnergy fe = compute_position_dependent_force(
	    pos, forceInternal, type, part, electricField, scheme, idx );
	if (get_energy) energy[idx] += fe.e;

        Vector3 r0 = pos[idx];
        Vector3 p0 = momentum[idx];

        Vector3 force = forceInternal[idx] + fe.f;
        forceInternal[idx] = forceInternal[idx] + fe.f; // Adding position-dep force so it doesn't need to be calculated in the "plain" kernel 
#ifdef Debug
        forceInternal[idx] = -force;
#endif

        #ifdef MDSTEP
        force = Vector3(-r0.x, -r0.y, -r0.z);
        #endif

        p0 = p0 + 0.5f * timestep * force * Unit1;
        momentum[idx] = p0;
    }
}

//Update kernel for Brownian dynamics
__global__
void updateKernel(Vector3* pos, Vector3* __restrict__ forceInternal, int type[], 
                  BrownianParticleType* part[],float kT, BaseGrid* kTGrid, float electricField, 
                  int tGridLength, float timestep, int num, int num_rb_attached_particles, BaseGrid* sys,
		  Random* randoGen, int numReplicas, float* energy, bool get_energy, int scheme) 
{
	// Calculate this thread's ID
	int idx = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
        
	// TODO: Make this a grid-stride loop to make efficient reuse of RNG states 
	// Loop over ALL particles in ALL replicas
	if (idx < num * numReplicas) 
        {
	    idx = (idx % num) + (idx/num) * (num+num_rb_attached_particles);
		const int t = type[idx];
		Vector3   p = pos[idx];

		const BrownianParticleType& pt = *part[t];
                
	 	/* printf("atom %d: forceInternal: %f %f %f\n", idx, forceInternal[idx].x, forceInternal[idx].y, forceInternal[idx].z);  */

		ForceEnergy fe = compute_position_dependent_force(
		    pos, forceInternal, type, part, electricField, scheme, idx );

		// Compute total force:
		//	  Internal:  interaction between particles
		//	  External:  electric field (now this is basically a constant vector)
		//	  forceGrid: ADD force due to PMF or other potentials defined in 3D space
		Vector3 force = forceInternal[idx] + fe.f;
                #ifdef Debug
                forceInternal[idx] = -force;
                #endif


		// Get local kT value
		float kTlocal;
                if(!scheme)
                    kTlocal = (tGridLength == 0) ? kT : kTGrid->interpolatePotentialLinearly(p); /* periodic */
                else
                    kTlocal = (tGridLength == 0) ? kT : kTGrid->interpolatePotential(p); /* periodic */

		// Update the particle's position using the calculated values for time, force, etc.
		float diffusion = pt.diffusion;
		Vector3 diffGrad = Vector3(0.0f);
		// printf("force: %f %f %f %f %f %f\n", p.x, p.y, p.z,
		//        fe.f.x, fe.f.y, fe.f.z);
		//        // force.x, force.y, force.z);

		if (pt.diffusionGrid != NULL) 
                {
			// printf("atom %d: pos: %f %f %f\n", idx, p.x, p.y, p.z);
			// p = pt.diffusionGrid->wrap(p); // illegal mem access; no origin/basis?

			Vector3 gridCenter = pt.diffusionGrid->origin +
				pt.diffusionGrid->basis.transform( Vector3(0.5*pt.diffusionGrid->nx,
									   0.5*pt.diffusionGrid->ny,
									   0.5*pt.diffusionGrid->nz)); 
			Vector3 p2 = p - gridCenter;
			p2 = sys->wrapDiff( p2 ) + gridCenter;			
			/* p2 = sys->wrap( p2 ); */
			/* p2 = p2 - gridCenter; */
			/* printf("atom %d: ps2: %f %f %f\n", idx, p2.x, p2.y, p2.z); */
                        ForceEnergy diff;
                        if(!scheme)	
			    diff = pt.diffusionGrid->interpolateForceDLinearly<periodic>(p2);
                        else
                            diff = pt.diffusionGrid->interpolateForceD(p2);
			diffusion = diff.e;
			diffGrad = diff.f;
		}

		// if (idx == 0) {
		// 	printf("force: "); force.print();
		// }
		
		Vector3 tmp = step(p, kTlocal, force, diffusion, -diffGrad, timestep, sys, randoGen, 
                                   num * numReplicas);
		// assert( tmp.length() < 10000.0f );
		pos[idx] = tmp;

                if(get_energy)
                {
                    float en_local = 0.f;
                    for(int i = 0; i < pt.numPartGridFiles; ++i)
                    {
			float en_tmp = 0.0f;
                        if(!scheme)
                            en_tmp = pt.pmf[i]->interpolatePotentialLinearly(tmp);
                        else
                            en_tmp = pt.pmf[i]->interpolatePotential(tmp);
			en_tmp *= pt.pmf_scale[i];
                    }
                    energy[idx] += en_local;
                }		
	}
}
/*
//This is the BBK Langevin integrator for Langevin dynamics Han-Yi Chou
__device__ inline void step(Vector3& r0, Vector3& p0, float kTlocal, Vector3 force, float diffusion, Vector3 diffGrad,
                            float mass, Vector3& gamma, float timestep, BaseGrid *sys, Random *randoGen, int num)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    Vector3 rando = randoGen->gaussian_vector(idx, num);
    float tmp = sqrtf(diffusion * timestep);

    p0.x = (1.0f - 0.50f * timestep * gamma.x) * p0.x + (0.50f * timestep * force.x * Unit1 + (tmp * sqrtf(gamma.x) * rando.x) * mass * Unit2);
    p0.y = (1.0f - 0.50f * timestep * gamma.y) * p0.y + (0.50f * timestep * force.y * Unit1 + (tmp * sqrtf(gamma.y) * rando.y) * mass * Unit2);
    p0.z = (1.0f - 0.50f * timestep * gamma.z) * p0.z + (0.50f * timestep * force.z * Unit1 + (tmp * sqrtf(gamma.z) * rando.z) * mass * Unit2);

    r0 = r0 + (timestep * p0) / mass * 1e4;
    Vector3 r = r0;
    r0 = sys->wrap(r);
}

//This is the BBK integrator for updating momentum in Langevin dynamics Han-Yi Chou
__device__ inline void step(Vector3& p0, float kTlocal, Vector3 force, float diffusion, Vector3 diffGrad,
                            float mass, Vector3& gamma, float timestep, BaseGrid *sys, Random *randoGen, int num)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    Vector3 rando = randoGen->gaussian_vector(idx, num);
    float tmp = sqrtf(diffusion * timestep);
    
    p0.x = (p0.x + (0.50 * timestep * force.x * Unit1 + (tmp * sqrtf(gamma.x) * rando.x) * mass * Unit2)) / (1.0f+0.5f*timestep*gamma.x);
    p0.y = (p0.y + (0.50 * timestep * force.y * Unit1 + (tmp * sqrtf(gamma.y) * rando.y) * mass * Unit2)) / (1.0f+0.5f*timestep*gamma.y);
    p0.z = (p0.z + (0.50 * timestep * force.z * Unit1 + (tmp * sqrtf(gamma.z) * rando.z) * mass * Unit2)) / (1.0f+0.5f*timestep*gamma.z);
}
*/
//For Brownian dynamics
__device__
inline Vector3 step(Vector3& r0, float kTlocal, Vector3 force, float diffusion,
						Vector3 diffGrad, float timestep, BaseGrid *sys,
						Random *randoGen, int num) {
	const int idx = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
	// TODO: improve performance by storing state locally, then sending it back to GPU
	Vector3 rando = (*randoGen).gaussian_vector(idx, num);

	diffusion *= timestep;
	Vector3 r = r0 + (diffusion / kTlocal) * force
							+ timestep * diffGrad
							+ sqrtf(2.0f * diffusion) * rando;
	// Wrap about periodic boundaries
	return sys->wrap(r);
}

template<bool print=false>
__global__
void updateGroupSites(Vector3 pos[], int* groupSiteData, int num, int numGroupSites, int numReplicas) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    // TODO: improve naive implementation so that each thread loads memory single pos elemment

    // For all threads representing a valid pair of particles
    if (i < numGroupSites*numReplicas) {
	pos[num*numReplicas + i] = Vector3(0.0f); 
    }

    // For all threads representing a valid pair of particles
    if (i < numGroupSites*numReplicas) {
	const int imod = i % numGroupSites;
	const int rep = i/numGroupSites;
	const int start  = groupSiteData[imod];
	const int finish = groupSiteData[imod+1];
	float weight = 1.0 / (finish-start);

	Vector3 tmp = Vector3(0.0f);

	for (int j = start; j < finish; j++) {
	    const int aj = groupSiteData[j] + num*rep;
	    tmp += weight * pos[aj];
	}
	if (print) {
	    printf("GroupSite %d (rep %d/%d) COM (start,finish, x,y,z): (%d,%d, %f,%f,%f)\n",i, rep, numReplicas, start, finish, tmp.x, tmp.y, tmp.z);
	}
	pos[num*numReplicas + i] = tmp;
    }
}

template<bool print=false>
__global__
void distributeGroupSiteForces(Vector3 force[], int* groupSiteData, int num, int numGroupSites, int numReplicas) {
    // TODO, handle groupsite energies
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    // For all threads representing a valid pair of particles
    if (i < numGroupSites*numReplicas) {
	const int imod = i % numGroupSites;
	const int rep = i/numGroupSites;
	const int start  = groupSiteData[imod];
	const int finish = groupSiteData[imod+1];
	float weight = 1.0 / (finish-start);

	const Vector3 tmp = weight*force[num*numReplicas+i];
	if (print) {
	    printf("GroupSite %d Force rep %d/%d: %f %f %f\n",i, rep, numReplicas, tmp.x, tmp.y, tmp.z);
	}

	for (int j = start; j < finish; j++) {
	    const int aj = groupSiteData[j] + num*rep;
	    atomicAdd( force+aj, tmp );
	}
    }
}

__global__ void devicePrint(RigidBodyType* rb[]) {
	// printf("Device printing\n");
	int i = 0;
	printf("RigidBodyType %d: numGrids = %d\n", i, rb[i]->numPotGrids);
	// printf("  RigidBodyType %d: potGrid: %p\n", i, rb[i]->rawPotentialGrids);
	// int j = 0;
	// printf("  RigidBodyType %d: potGrid[%d]: %p\n", i, j, &(rb[i]->rawPotentialGrids[j]));
	// printf("  RigidBodyType %d: potGrid[%d] size: %d\n", i, j, rb[i]->rawPotentialGrids[j].getSize());
	// BaseGrid g = rb[i]->rawPotentialGrids[j];
	// for (int k = 0; k < rb[i]->rawPotentialGrids[j].size(); k++)
	// for (int k = 0; k < rb[i]->rawPotentialGrids[j].getSize(); k++)
	// 	printf("    rbType_d[%d]->potGrid[%d].val[%d]: %g\n",
	// 				 i, j, k, rb[i]->rawPotentialGrids[j].val[k]);
	// i, j, k, rb[i]->rawPotentialGrids[j]).val[k];
	
}

// __global__ void devicePrint(RigidBodyType* rb[]) {
// 	// printf("Device printing\n");
// 	int i = 0;
// 	printf("RigidBodyType %d: numGrids = %d\n", i, rb[i]->numPotGrids);
// 	printf("RigidBodyType %d: potGrid: %p\n", i, rb[i]->rawPotentialGrids);
// 	int j = 0;
// 	printf("RigidBodyType %d: potGrid[%d]: %p\n", i, &(rb[i]->rawPotentialGrids[j]));
// 	BaseGrid g = rb[i]->rawPotentialGrids[j];
// 	// for (int k = 0; k < rb[i]->rawPotentialGrids[j].size(); k++)
// 	for (int k = 0; k < g->getSize(); k++)
// 		printf("rbType_d[%d]->potGrid[%d].val[%d]: %g\n",
// 					 i, j, k, g.val[k]);
// 	// i, j, k, rb[i]->rawPotentialGrids[j]).val[k];
	
// }


// __device__ Vector3* totalForce;
// Vector3* totalForce_h;
// void initTotalForce() {
//     cudaMalloc( &totalForce_h, sizeof(Vector3) );
//     cudaMemcpyToSymbol(totalForce, &totalForce_h, sizeof(totalForce_h));
// }

__global__ void compute_position_dependent_force_for_rb_attached_particles(
    const Vector3* __restrict__ pos,
    Vector3* __restrict__ forceInternal, float* __restrict__ energy,
    const int* __restrict__ type, BrownianParticleType** __restrict__ part,
    const float electricField, const int num, const int num_rb_attached_particles,
    const int numReplicas, const int scheme)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < num_rb_attached_particles * numReplicas)
    {
	idx = num + (idx % num_rb_attached_particles) + (idx/num_rb_attached_particles) * (num+num_rb_attached_particles);
	ForceEnergy fe = compute_position_dependent_force(
	    pos, forceInternal, type, part, electricField, scheme, idx );
	atomicAdd( &forceInternal[idx], fe.f );
	atomicAdd( &energy[idx], fe.e );
    }
}
__global__ void compute_position_dependent_force_for_rb_attached_particles(
    const Vector3* __restrict__ pos, Vector3* __restrict__ forceInternal,
    const int* __restrict__ type, BrownianParticleType** __restrict__ part,
    const float electricField, const int num, const int num_rb_attached_particles,
    const int numReplicas, const int scheme)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < num_rb_attached_particles * numReplicas)
    {
	idx = num + (idx % num_rb_attached_particles) + (idx/num_rb_attached_particles) * (num+num_rb_attached_particles);
	ForceEnergy fe = compute_position_dependent_force(
	    pos, forceInternal, type, part, electricField, scheme, idx );
	atomicAdd( &forceInternal[idx], fe.f );
    }
}
