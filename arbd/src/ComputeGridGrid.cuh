// Included in RigidBodyController.cu
#pragma once
#include "useful.h"
#define NUMTHREADS 128
#define WARPSIZE 32

class RigidBodyGrid;
class BaseGrid;

extern __global__
void computeGridGridForce(const RigidBodyGrid* rho, const RigidBodyGrid* u,
				const Matrix3 basis_rho, const Matrix3 basis_u_inv,
				const Vector3 origin_rho_minus_center_u, const Vector3 center_u_minus_origin_u,
				ForceEnergy* retForce, Vector3 * retTorque, int scheme, BaseGrid* sys_d);

extern __global__
void computePmfGridForce(const RigidBodyGrid* rho, const RigidBodyGrid* u,
			 const Matrix3 basis_rho, const Matrix3 basis_u_inv,
			 const Vector3 origin_rho_minus_origin_u,
			 ForceEnergy* retForce, Vector3 * retTorque, int scheme);

extern __global__
void computePartGridForce(const Vector3* __restrict__ pos, Vector3* particleForce,
				const int num, const int* __restrict__ particleIds,
				const RigidBodyGrid* __restrict__ u,
				const Matrix3 basis_u_inv, const Vector3 center_u, const Vector3 origin_u,
				ForceEnergy* __restrict__ retForceTorque, float* energy, bool get_energy, int scheme, BaseGrid* sys_d);

extern __global__
void createPartlist(const Vector3* __restrict__ pos,
				const int numTypeParticles, const int* __restrict__ typeParticles_d,
		    const int attached_particle_start, const int attached_particle_end,
				int* numParticles_d, int* particles_d,
				const Vector3 gridCenter, const float radius2, BaseGrid* sys_d);
	
extern __global__
void printRigidBodyGrid(const RigidBodyGrid* rho);
