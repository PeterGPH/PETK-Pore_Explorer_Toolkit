// Included in RigidBodyController.cu
#include "ComputeGridGrid.cuh"
#include "RigidBodyGrid.h"
#include "CudaUtil.cuh"
//RBTODO: add __restrict__, benchmark (Q: how to restrict member data?)

class GridPositionTransformer {
public:
    __device__ GridPositionTransformer(const Vector3 o, const Vector3 c, BaseGrid* s) :
	o(o), c(c), s(s) { }
    __device__ inline Vector3 operator() (Vector3 pos) const {
	return s->wrapDiff(pos + o) + c;
    }
private:
    const Vector3 o;
    const Vector3 c;
    const BaseGrid* s;
};

//class PmfPositionTransformer : public BasePositionTransformer {
class PmfPositionTransformer {
public:
    __device__ PmfPositionTransformer(const Vector3 o) : o(o) { }
    __device__ inline Vector3 operator() (Vector3 pos) const {
	return pos + o;
    }
private:
    const Vector3 o;
};

template <typename T>
__device__
inline void common_computeGridGridForce(const RigidBodyGrid* rho, const RigidBodyGrid* u, const Matrix3 basis_rho, const Matrix3 basis_u_inv, const T& transformer,
					ForceEnergy* retForce, Vector3 * retTorque, int scheme)
{

	extern __shared__ ForceEnergy s[];
	ForceEnergy *force = s;
	//Vector3 *torque = &s[NUMTHREADS];
        ForceEnergy *torque = &s[NUMTHREADS];

  // RBTODO: http://devblogs.nvidia.com/parallelforall/cuda-pro-tip-write-flexible-kernels-grid-stride-loops
	const int tid = threadIdx.x;
	const int r_id = blockIdx.x * blockDim.x + threadIdx.x;

	force[tid] = ForceEnergy(0.f,0.f);
	torque[tid] = ForceEnergy(0.f,0.f);
	if (r_id < rho->getSize()) { // skip threads with nothing to do
		// RBTODO: reduce registers used;
		//   commenting out interpolateForceD still uses ~40 registers
		//   -- the innocuous-looking fn below is responsible; consumes ~17 registers!
	    Vector3 r_pos= rho->getPosition(r_id); /* i,j,k value of voxel */

	    r_pos = basis_rho.transform( r_pos );
	    const Vector3 u_ijk_float = basis_u_inv.transform( transformer( r_pos ) );
		// RBTODO: Test for non-unit delta
		/* Vector3 tmpf  = Vector3(0.0f); */
		/* float tmpe = 0.0f; */
		/* const ForceEnergy fe = ForceEnergy( tmpf, tmpe); */

                ForceEnergy fe;
                if(!scheme)
		    fe = u->interpolateForceDLinearly( u_ijk_float ); /* in coord frame of u */
                else
                    fe = u->interpolateForceD( u_ijk_float );

		force[tid] = fe;
                //force[tid].e = fe.e;

		const float r_val = rho->val[r_id]; /* maybe move to beginning of function?  */
		force[tid].f = basis_u_inv.transpose().transform( r_val*(force[tid].f) ); /* transform to lab frame, with correct scaling factor */
                force[tid].e = r_val;

		// Calculate torque about origin_rho in the lab frame
		torque[tid].f = r_pos.cross(force[tid].f);
	}

	// Reduce force and torques
	// http://www.cuvilib.com/Reduction.pdf
	// RBTODO optimize further, perhaps
	// assert( NUMTHREADS==32 || NUMTHREADS==64 || NUMTHREADS==128 || NUMTHREADS==256 || NUMTHREADS==512 );
	__syncthreads();
	for (int offset = blockDim.x/2; offset > 0; offset >>= 1) {
		if (tid < offset) {
			int oid = tid + offset;
                        //if(get_energy)
                            //force[tid].e = force[tid].e + force[oid].e;
			force[tid] = force[tid] + force[oid];
			torque[tid] = torque[tid] + torque[oid];
		}
		__syncthreads();
	}

	if (tid == 0) {
		atomicAdd( retForce, force[0] ); // apply force to particle
		atomicAdd( retTorque, torque[0].f ); // apply force to particle
	}
}

__global__
void computeGridGridForce(const RigidBodyGrid* rho, const RigidBodyGrid* u, const Matrix3 basis_rho, const Matrix3 basis_u_inv, const Vector3 origin_rho_minus_center_u, const Vector3 center_u_minus_origin_u,
			ForceEnergy* retForce, Vector3 * retTorque, int scheme, BaseGrid* sys_d)
{
    GridPositionTransformer transformer = GridPositionTransformer(origin_rho_minus_center_u, center_u_minus_origin_u, sys_d);
    common_computeGridGridForce<GridPositionTransformer>(rho, u, basis_rho, basis_u_inv, transformer, retForce, retTorque, scheme);
}

__global__
void computePmfGridForce(const RigidBodyGrid* rho, const RigidBodyGrid* u, const Matrix3 basis_rho, const Matrix3 basis_u_inv, const Vector3 origin_rho_minus_origin_u,
			 ForceEnergy* retForce, Vector3 * retTorque, int scheme)
{
    PmfPositionTransformer transformer = PmfPositionTransformer(origin_rho_minus_origin_u);
    common_computeGridGridForce<PmfPositionTransformer>(rho, u, basis_rho, basis_u_inv, transformer, retForce, retTorque, scheme);
}

__global__
void computePartGridForce(const Vector3* __restrict__ pos, Vector3* particleForce,
				const int num, const int* __restrict__ particleIds, 
				const RigidBodyGrid* __restrict__ u,
				const Matrix3 basis_u_inv, const Vector3 center_u, const Vector3 origin_u,
				ForceEnergy* __restrict__ retForceTorque, float* __restrict__ energy, bool get_energy, int scheme, BaseGrid* sys_d) {

	extern __shared__ ForceEnergy s[];
	ForceEnergy *force  = s;
	ForceEnergy *torque = &s[NUMTHREADS];
  	
	const int tid = threadIdx.x;
	const int i = blockIdx.x * blockDim.x + threadIdx.x;

	force[tid]  = ForceEnergy(0.f, 0.f);
	torque[tid] = ForceEnergy(0.f,0.f);
	if (i < num) {
		const int id = particleIds[i];
		Vector3 p = sys_d->wrapDiff(pos[id]-center_u) + center_u - origin_u;
		const Vector3 u_ijk_float = basis_u_inv.transform( p );

                ForceEnergy fe;
                if(!scheme)                       
		    fe = u->interpolateForceDLinearly( u_ijk_float ); /* in coord frame of u */
                else
                    fe = u->interpolateForceD( u_ijk_float );
                
		force[tid] = fe;
                //force[tid].e = fe.e;
                if(get_energy)
                    atomicAdd(&energy[id], fe.e);
		force[tid].f = basis_u_inv.transpose().transform( force[tid].f ); /* transform to lab frame */
		atomicAdd( &particleForce[id], force[tid].f ); // apply force to particle
		
		// Calculate torque about origin_u in the lab frame
		torque[tid].f = p.cross(force[tid].f);				// RBTODO: test sign
	}

	// Reduce force and torques
	// assert( NUMTHREADS==32 || NUMTHREADS==64 || NUMTHREADS==128 || NUMTHREADS==256 || NUMTHREADS==512 );
	__syncthreads();
	for (int offset = blockDim.x/2; offset > 0; offset >>= 1) {
		if (tid < offset) {
			int oid = tid + offset;
                        //if(get_energy)
                            //force[tid].e = force[tid].e + force[oid].e;
			force[tid] = force[tid] + force[oid];
			torque[tid] = torque[tid] + torque[oid];
		}
		__syncthreads();
	}
	
	if (tid == 0) {
		retForceTorque[2*blockIdx.x] = force[0];
		retForceTorque[2*blockIdx.x+1] = torque[0];
	}
}


__global__
void createPartlist(const Vector3* __restrict__ pos,
				const int numTypeParticles, const int* __restrict__ typeParticles_d,
		    const int attached_particle_start, const int attached_particle_end,
				int* numParticles_d, int* particles_d,
				const Vector3 gridCenter, const float radius2, BaseGrid* sys_d) {
	const int tid = threadIdx.x;
	const int warpLane = tid % WARPSIZE; /* RBTODO: optimize */
	
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < numTypeParticles) {
		int aid = typeParticles_d[i];
		if (aid < attached_particle_start || aid >= attached_particle_end) { 
		    float dist = (sys_d->wrapDiff(pos[aid] - gridCenter)).length2();
		
		    if (dist <= radius2) {
			int tmp = atomicAggInc(numParticles_d, warpLane);
			particles_d[tmp] = aid;
		    }
		}
	}
}		

__global__
void printRigidBodyGrid(const RigidBodyGrid* rho) {
  printf("Printing an RB of size %d\n",rho->size);
  for (int i=0; i < rho->size; i++)
	printf("  val[%d] = %f\n", i, rho->val[i]);
}
