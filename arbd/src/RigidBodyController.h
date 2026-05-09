#pragma once

#include <vector>
#include <fstream>
/* #include "RandomCUDA.h" */
#include <cuda.h>
#include <cuda_runtime.h>
#include "useful.h"
#include "BaseGrid.h"
#include "RigidBodyGrid.h"
#include "GPUManager.h"

#define NUMSTREAMS 8

// #include "RigidBody.h"

class RigidBodyType;
class RigidBody;
class Configuration;
class ForceEnergy;
// class RandomCPU;
#include "RandomCPU.h"

// TODO: performance: create RigidBodyGridPair so pairlistdist check is done per grid pair, not per RB pair
class RigidBodyForcePair  {
	friend class RigidBodyController;

public:
	RigidBodyForcePair(RigidBodyType* t1, RigidBodyType* t2,
										 RigidBody* rb1, RigidBody* rb2,
					std::vector<int> gridKeyId1, std::vector<int> gridKeyId2, bool isPmf, int updatePeriod) :
	updatePeriod(updatePeriod), type1(t1), type2(t2), rb1(rb1), rb2(rb2),
		gridKeyId1(gridKeyId1), gridKeyId2(gridKeyId2), isPmf(isPmf)
		{
			printf("    Constructing RB force pair...\n");
			/* initialize(); */
			// printf("    done constructing RB force pair\n");
		}
	RigidBodyForcePair(const RigidBodyForcePair& o) :
		updatePeriod(o.updatePeriod), type1(o.type1), type2(o.type2), rb1(o.rb1), rb2(o.rb2),
		gridKeyId1(o.gridKeyId1), gridKeyId2(o.gridKeyId2), isPmf(o.isPmf) {
		printf("    Copying RB force pair...\n");
		/* initialize(); */
	}
	RigidBodyForcePair& operator=(RigidBodyForcePair& o) {
		printf("    Copying assigning RB force pair...\n");
		swap(*this,o);
		return *this;
	}	
	~RigidBodyForcePair();

	bool isOverlapping(BaseGrid* sys) const;

private:
	int initialize();
	void swap(RigidBodyForcePair& a, RigidBodyForcePair& b);

	int updatePeriod;
	
	RigidBodyType* type1;
	RigidBodyType* type2;
	RigidBody* rb1;
	RigidBody* rb2;
	
	std::vector<int> gridKeyId1;
	std::vector<int> gridKeyId2;
	std::vector<int> numBlocks;

	bool isPmf;
	
	//std::vector<Vector3*> forces;
	//std::vector<Vector3*> forces_d;
	std::vector<ForceEnergy*> forces;
        std::vector<ForceEnergy*> forces_d;
	std::vector<Vector3*> torques;
	std::vector<Vector3*> torques_d;

	static int nextStreamID; 
	std::vector<int> streamID;
	static cudaStream_t* stream;
	static void createStreams();

	static int lastStreamID;
	static RigidBodyForcePair* lastRbForcePair;
	static int lastRbGridID;
	
	void callGridForceKernel(int pairId, int s,int scheme, BaseGrid* sys_d);
	void retrieveForcesForGrid(const int i);
	void processGPUForces(BaseGrid*);
	Matrix3 getBasis1(const int i);
	Matrix3 getBasis2(const int i);
	Vector3 getOrigin1(const int i);
	Vector3 getOrigin2(const int i);
	Vector3 getCenter2(const int i);


	static GPUManager gpuman;
};

class RigidBodyController {
public:
	/* DEVICE RigidBodyController(const NamdState *s, int reductionTag, SimParameters *sp); */
	RigidBodyController();
        ~RigidBodyController();
	RigidBodyController(const Configuration& c, const char* outArg, unsigned long int seed, int repID);

        void AddLangevin();
        void SetRandomTorques();
	void integrate(BaseGrid* sys, int step);
        void integrateDLM(BaseGrid* sys, int step);
	void updateForces(Vector3* pos_d, Vector3* force_d, int s, float* energy, bool get_energy, int scheme, BaseGrid* sys, BaseGrid* sys_d, int num, int num_rb_attached_particles);
	void updateParticleLists(Vector3* pos_d, BaseGrid* sys_d);
    void update_attached_particle_positions(Vector3* pos_d, Vector3* force_d, float* energy_d, BaseGrid* sys_d, int num, int num_rb_attached_particles, int numReplicas);
        void clearForceAndTorque(); 
        void KineticEnergy();
        void print(int step);
        //void printEnergyData(std::fstream &file);
        float getEnergy(float (RigidBody::*get)());
private:
	bool loadRBCoordinates(const char* fileName);
    bool load_restart_coordinates(const char* filename);
	void initializeForcePairs();

	//void print(int step);
	void printLegend(std::ofstream &file);
	void printData(int step, std::ofstream &file);
public:
	RigidBodyType** rbType_d;

	inline Vector3 getRandomGaussVector() {
	    return random->gaussian_vector();
	}
	/* RequireReduction *gridReduction; */

	BaseGrid* grids;
	RigidBodyGrid* grids_d;
	
private:
	std::ofstream trajFile;
	
	const Configuration& conf;
	char outArg[128];
	
	RandomCPU* random;
	/* RequireReduction *gridReduction; */
	
	Vector3* trans; // would have made these static, but
	Matrix3* rot;  	// there are errors on rigidBody->integrate
	std::vector< std::vector<RigidBody> > rigidBodyByType;
	std::vector< RigidBodyForcePair > forcePairs;

	void construct_grids();
	void destruct_grids();

        //float* rb_energy;
	Vector3* attached_particle_forces;
	ForceEnergy* particleForces;
	ForceEnergy* particleForces_d;
	std::vector<int> particleForceNumBlocks;
	std::vector<int> particleForce_offset;
	int totalParticleForceNumBlocks;
};
