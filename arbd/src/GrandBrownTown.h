// Author: Jeff Comer <jcomer2@illinois.edu>
#ifndef GRANDBROWNTOWN_H
#define GRANDBROWNTOWN_H

#ifdef __CUDACC__
    #define HOST __host__
    #define DEVICE __device__
#else
    #define HOST
    #define DEVICE
#endif

#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <locale.h> // setlocale
#include <sstream> // std::stringstream
#include <string> // std::string
#include <vector> // std::vector

#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>

#include "GPUManager.h"
#include "useful.h"
#include "BaseGrid.h"
#include "OverlordGrid.h"
#include "Reader.h"
#include "RandomCUDA.h"
#include "ComputeForce.h"
#include "BrownianParticleType.h"
#include "TrajectoryWriter.h"
#include "JamesBond.h"
#include "Exclude.h"
#include "Angle.h"
#include "Configuration.h"
#include "Dihedral.h"
/* #include "RigidBody.h" */
/* #include "RigidBodyType.h" */
/* #include "RigidBodyGrid.h" */
#include "RigidBodyController.h"
#include "WKFUtils.h"

// IMD
#include "vmdsock.h"
#include "imd.h"

//#include "analyticForce.h"

// using namespace std;

// #define FORCEGRIDOFF

class GrandBrownTown {
public:
	GrandBrownTown(const char* configFile, const char* outArg,
			bool debug, bool imd_on, unsigned int imd_port, int numReplicas = 0);
	GrandBrownTown(const Configuration& c, const char* outArg,
			bool debug, bool imd_on, unsigned int imd_port, int numReplicas = 0);
	~GrandBrownTown();

	void run();
	static bool DEBUG;

private:  

	// Given the numbers of each particle, populate the type list.
	void populate();

	// Count the number of atoms in the restart file.
	int countRestart(const char* fileName);

	void writeRestart(int repID) const;
        void writeMomentumRestart(int repID) const;

	void initialCondCen();
	void initialCond();

	// A couple old routines for getting particle positions.
	Vector3 findPos(int typ);
	Vector3 findPos(int typ, float minZ);
	
	bool readTableFile(const String& value, int currTab);
	bool readBondFile(const String& value, int currBond);
	bool readAngleFile(const String& value, int currAngle);

	void newCurrent(int repID) const;
	void writeCurrent(int repID, float t) const;
	void writeCurrentSegment(int repID, float t, float segZ) const;
	void getDebugForce();
	
	void copyRandToCUDA();
	void copyToCUDA();

        //Compute the kinetic energy in general. Han-Yi Chou
        float KineticEnergy();
        //float RotKineticEnergy();

        //Initialize the Nose-Hoover auxilliary variables
        void InitNoseHooverBath(int N);
        //curandState_t *randoDevice;

	void init_cuda_group_sites();

public:
	// Compute the current in nanoamperes.
	float current(float t) const;

	// Compute the current in nanoamperes for a restricted segment (-segZ < z < segZ).
	float currentSegment(float t, float segZ, int carrier) const;

	int getReservoirCount(int partInd, int resInd) const;

	IndexList getReservoirList(int partInd, int resInd) const;
		
	// Find an open position to place a particle.
	Vector3 freePosition(Vector3 r0, Vector3 r1, float minDist);

private:
	static GPUManager gpuman;
	const Configuration& conf;
	int numReplicas;
	
	// IMD variables
	bool imd_on;
	unsigned int imd_port;
	Vector3* imdForces;
	
	// Output variables
	std::vector<std::string> outCurrFiles;
	std::vector<std::string> restartFiles;
	std::vector<std::string> outFilePrefixes;

        //Hna-Yi Chou Langevin Dynamics
        std::vector<std::string> restartMomentumFiles;
        std::vector<std::string> outMomentumFilePrefixes;//, outForceFilePrefixes;

	std::vector<TrajectoryWriter*> writers;

        //For momentum, i.e. Langevin dynamic Han-Yi Chou
        std::vector<TrajectoryWriter*> momentum_writers;
        //std::vector<TrajectoryWriter*> force_writers;

	Vector3 sysDim;

	// Integrator variables
	BaseGrid* sys;
	Random *randoGen;
	ComputeForce* internal;
	Vector3* forceInternal;

	// Particle variables
	String* partsFromFile;
	int* indices;
	int numPartsFromFile;
	Bond* bonds;
	int numCap; 		// max number of particles
	int num; 			// number of particles
	Vector3* pos; 		// particle positions
        Vector3* momentum;      // particle momentum Han-Yi Chou
        float *random;
        //Vector3* force;
	int* type; 			// particle types: 0, 1, ... -> num * numReplicas
	String* name; 		// particle types: POT, CLA, ... -> num * numReplicas
	int* serial; 		// particle serial numbers
	int currSerial; 	// the serial number of the next new particle
	Vector3* posLast; 	// previous positions of particles
        Vector3* momLast;
	float timeLast; 	// used with posLast
	float minimumSep; 	// minimum separation allowed with placing new particles

	std::vector<RigidBodyController*> RBC;
	Vector3* rbPos; 		// rigid body positions
	
	// CUDA device variables
	//Vector3 *pos_d, *forceInternal_d, *force_d;
	//int *type_d;
	BrownianParticleType **part_d;
	BaseGrid *sys_d, *kTGrid_d;
	Random* randoGen_d;
	//Bond* bonds_d;
	//int2* bondMap_d;
	//Exclude* excludes_d;
	//int2* excludeMap_d;
	//Angle* angles_d;
	//Dihedral* dihedrals_d;

	// System parameters
	String outputName;
	float timestep;
	long int steps;
	unsigned long int seed;
	String temperatureGridFile;
	String inputCoordinates;
	String restartCoordinates;
	int numberFluct;
	int interparticleForce;
	int tabulatedPotential;
	int fullLongRange;
	float kT;
	float temperature;
	float coulombConst;
	float electricField;
	float cutoff;
	float switchLen;
	int outputPeriod;
	int outputEnergyPeriod;
	int outputFormat;
	float currentSegmentZ;
	int numberFluctPeriod;
	int decompPeriod;
	int numCapFactor;
	BaseGrid* kTGrid;
	BaseGrid* tGrid;
	BaseGrid* sigmaT;

	// Other parameters.
	float switchStart;
	float maxInitialPot;
	float initialZ;

	// Particle parameters.
	BrownianParticleType* part;
	int numParts;
	int numBonds;
	int numExcludes;
	int numAngles;
	int numDihedrals;

    int num_rb_attached_particles;

	int numGroupSites;
	int* groupSiteData_d;

	String partFile;
	String bondFile;
	String excludeFile;
	String angleFile;
	String dihedralFile;
	bool readPartsFromFile;
	bool readBondsFromFile;
	bool readExcludesFromFile;
	bool readAnglesFromFile;
	bool readDihedralsFromFile;
	String* partGridFile;
	String* partDiffusionGridFile;
	String* partForceXGridFile;
	String* partForceYGridFile;
	String* partForceZGridFile;
	String* partTableFile;
	String* partReservoirFile;
	int* partTableIndex0;
	int* partTableIndex1;

	String* bondTableFile;
	int numTabBondFiles;
	int2* bondMap;
	int3 *bondList;

	Exclude* excludes;
	int2* excludeMap;
	String excludeRule;
	int excludeCapacity;

	Angle* angles;
	String* angleTableFile;
	int numTabAngleFiles;
	int4 *angleList;

	Dihedral* dihedrals;
	String* dihedralTableFile;
	int numTabDihedralFiles;
	int4 *dihedralList;
	int  *dihedralPotList;

	int numBondAngles;
	BondAngle* bondAngles;
	int4 *bondAngleList;

        //Han-Yi Chou
        String particle_dynamic;
        String rigidbody_dynamic;
        String particle_langevin_integrator;
        int ParticleInterpolationType;
        int RigidBodyInterpolationType;
	void updateNameList();

	void remember(float t);

	void deleteParticles(IndexList& p);

	void addParticles(int n, int typ);

	// Add particles randomly within the region defined by r0 and r1.
	void addParticles(int n, int typ, Vector3 r0, Vector3 r1);

	// Add particles randomly within the region defined by r0 and r1.
	// Maintains a minimum distance of minDist between particles.
	void addParticles(int n, int typ, Vector3 r0, Vector3 r1, float minDist);

	// Add or delete particles in the reservoirs.
	// Reservoirs are not wrapped.
	void updateReservoirs();

};

#endif
