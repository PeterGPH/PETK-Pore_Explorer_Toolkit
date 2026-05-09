/*===========================\
| RigidBody Class for device |
\===========================*/
#pragma once

#include "useful.h"
#include "RandomCPU.h"		/* for BD integration; RBTODO: fix this */
#include "GPUManager.h"

#ifdef __CUDACC__
    #define HOST __host__
    #define DEVICE __device__
#else
    #define HOST 
    #define DEVICE 
#endif

#include "RigidBodyType.h"
#include "RigidBodyController.h"

class Configuration;
class BaseGrid;
typedef float BigReal;					/* strip this out later */
typedef Vector3 Force;


class RigidBody { // host side representation of rigid bodies
	friend class RigidBodyController;
	/*––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––.
	| See Appendix A of: Dullweber, Leimkuhler and McLaclan. "Symplectic        |
	| splitting methods for rigid body molecular dynamics". J Chem Phys. (1997) |
	`––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––*/
		public:
    RigidBody(String name, const Configuration& c, const RigidBodyType& t, RigidBodyController* RBC,
	      int attached_particle_start, int attached_particle_end);

    RigidBody(const RigidBody& rb);
    // RigidBody(const RigidBody& rb) : RigidBody(rb.name, *rb.c, *rb.t) {};
	void init();
	/* HOST DEVICE RigidBody(RigidBodyType t); */
	~RigidBody();

	int appendNumParticleBlocks( std::vector<int>* blocks );

    void update_particle_positions(Vector3* pos_d, Vector3* force_d, float* energy_d);

	HOST DEVICE void addForce(Force f); 
	HOST DEVICE void addTorque(Force t);
        HOST DEVICE void addEnergy(float e);
	HOST DEVICE void addLangevin(Vector3 w1, Vector3 w2);
        HOST inline void setKinetic(float e) { kinetic = e; };	
	HOST DEVICE inline void clearForce() { force = Force(0.0f); energy = 0.f;}
	//HOST DEVICE inline void clearForce() { force = ForceEnergy(0.f, 0.f); }
	HOST DEVICE inline void clearTorque() { torque = Force(0.0f); }

	// HOST DEVICE void integrate(Vector3& old_trans, Matrix3& old_rot, int startFinishAll);
	// HOST DEVICE void integrate(Vector3& old_trans, Matrix3& old_rot, int startFinishAll);
	void integrateDLM(BaseGrid* sys, int startFinishAll);
	void integrate(BaseGrid* sys, int startFinishAll);	

	// HOST DEVICE inline String getKey() const { return key; }
	// HOST DEVICE inline String getKey() const { return t->name; }
	HOST DEVICE inline String getKey() const { return name; }
	
	HOST DEVICE inline Vector3 transformBodyToLab(Vector3 v) const { return orientation*v + position; }
	HOST DEVICE inline Vector3 getPosition() const { return position; }
	HOST DEVICE inline Matrix3 getOrientation() const { return orientation; }
	// HOST DEVICE inline Matrix3 getBasis() const { return orientation; }
	HOST DEVICE inline BigReal getMass() const { return t->mass; }
	//HOST DEVICE inline Vector3 getVelocity() const { return momentum/t->mass; }
	HOST DEVICE inline Vector3 getVelocity() const { return momentum; }
        HOST float getEnergy() { return energy; }
        HOST float getKinetic(){ return kinetic; }
	//HOST DEVICE inline Vector3 getAngularVelocity() const { 
	//	return Vector3( angularMomentum.x / t->inertia.x,
	//								 angularMomentum.y / t->inertia.y,
									 //angularMomentum.z / t->inertia.z );
	//}
	HOST DEVICE inline Vector3 getAngularVelocity() const { 
              return Vector3( angularMomentum.x, angularMomentum.y, angularMomentum.z);
        }

        void initializeParticleLists();
	void updateParticleList(Vector3* pos_d, BaseGrid* sys_d);
	void callGridParticleForceKernel(Vector3* pos_d, Vector3* force_d, int s, float* energy, bool get_energy, int scheme, BaseGrid* sys, BaseGrid* sys_d, ForceEnergy* forcestorques_d, const std::vector<int>& forcestorques_offset, int& fto_idx);
	void apply_attached_particle_forces(const Vector3* force);
	void applyGridParticleForces(BaseGrid* sys, ForceEnergy* forcestorques, const std::vector<int>& forcestorques_offset, int& fto_idx);
	
	bool langevin;
	Vector3 torque; // lab frame (except in integrate())
        
private:
	static GPUManager gpuman;

	RigidBodyController* RBC;
	inline Vector3 getRandomGaussVector() { 
	    return RBC->getRandomGaussVector();
	}

	// String key;
	String name;
	/* static const SimParameters * simParams; */
	Vector3 position;		  /* position of center of mass */
	// Q = orientation.transpose(); in Dullweber et al
	Matrix3 orientation;					/* rotation that brings RB coordinates into the lab frame */

	Vector3 momentum;		 /* in lab frame */
	Vector3 angularMomentum; // angular momentum along corresponding principal axes
        Vector3 W1,W2;
 
	// Langevin
	Vector3 langevinTransFriction; /* RBTODO: make this work with a grid */
	Vector3 langevinRotFriction;
	BigReal Temp;

	/* Vector3 transDampingCoeff; */
	/* Vector3 transForceCoeff; */
	/* Vector3 rotDampingCoeff; */
	/* Vector3 rotTorqueCoeff;     */

	// integration
	const Configuration* c;
	const RigidBodyType* t;
	float timestep;					
	Vector3 force;  // lab frame
        float energy; //potential energy
        float kinetic; 
	bool isFirstStep; 
	
	int* numParticles;		  /* particles affected by potential grids */
	int** possible_particles_d;		 	
	int** particles_d;		 	
	const cudaStream_t** particleForceStreams;

    int attached_particle_start, attached_particle_end;
    
	/*–––––––––––––––––––––––––––––––––––––––––.
	| units "kcal_mol/AA * ns" "(AA/ns) * amu" |
	`–––––––––––––––––––––––––––––––––––––––––*/
	BigReal impulse_to_momentum; /* should be const, but copy constructor failed */

	HOST DEVICE inline void applyRotation(const Matrix3& R);
	HOST DEVICE inline Matrix3 Rx(BigReal t);
	HOST DEVICE inline Matrix3 Ry(BigReal t);
	HOST DEVICE inline Matrix3 Rz(BigReal t);
	HOST DEVICE inline Matrix3 eulerToMatrix(const Vector3 e);
        float Temperature();
        void  Boltzmann(unsigned long int);
};

