// Configuration.h (2013)
// Loads .brown config file that can be shared between simulations
// To be used by GrandBrownTown to initialize its members
//
// Authors: Terrance Howward <howard33@illinois.edu>
//          Justin Dufresne <jdufres1@friars.providence.edu>
//

#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include <algorithm> // sort
#include <vector>
#include <map>

#include <cuda.h>
#include <cuda_runtime.h>

#include "useful.h" // String, Vector3
#include "BrownianParticleType.h"
#include "BaseGrid.h"
#include "OverlordGrid.h"
#include "ComputeForce.h"
#include "Reader.h"
#include "TrajectoryWriter.h"
#include "TabulatedPotential.h"
#include "TabulatedAngle.h"
#include "ProductPotential.h"
#include "GPUManager.h"
#include "RigidBodyType.h"
#include "RigidBody.h"

// Units:
//    Energy: kcal/mol (6.947694e-24 kJ)
//    Temperature: Kelvin
//    Time: nanoseconds
//    Length: angstroms
//    Momentum: Da * \mu m / ns

// Forward declerations
class Angle;
class Dihedral;
struct Restraint;

class Configuration {
	struct compare {
		bool operator()(const String& lhs, const String& rhs);
		bool operator()(const Bond& lhs, const Bond& rhs);
		bool operator()(const Exclude& lhs, const Exclude& rhs);
		bool operator()(const Angle& lhs, const Angle& rhs);
		bool operator()(const Dihedral& lhs, const Dihedral& rhs);
		bool operator()(const BondAngle& lhs, const BondAngle& rhs);
		bool operator()(const ProductPotentialConf& lhs, const ProductPotentialConf& rhs);
	};

	void setDefaults();
	Vector3 stringToVector3(String s);
	Matrix3 stringToMatrix3(String s);

	int readParameters(const char* config_file);
	void readAngles();
	void readAtoms();
	void readGroups();
	void readBonds();
	void readExcludes();
	void addExclusion(int ind1, int ind2);
	void buildExcludeMap();
	void readDihedrals();
	void readRestraints();
	void readBondAngles();


	bool readTableFile(const String& value, int currTab);
	bool readBondFile(const String& value, int currBond);
	bool readAngleFile(const String& value, int currAngle);
	bool readDihedralFile(const String& value, int currDihedral);

	bool readBondAngleFile(const String& value, const String& bondfile1, const String& bondfile2, int currBondAngle);

	// Given the numbers of each particle, populate the type list.
	void populate();

	void loadRestart(const char* file_name);
	bool loadCoordinates(const char* file_name);
	int countRestart(const char* file_name);

	void getDebugForce();

        //Han-Yi Chou
        bool Boltzmann(const Vector3& com_v,int N);
        bool loadMomentum(const char* file_name);
        void loadRestartMomentum(const char* file_name);
        void Print();
        void PrintMomentum();
public:
	Configuration(const char * config_file, int simNum = 0, bool debug=false);
	~Configuration();

    int find_particle_type(const char* s) const {
	for (int j = 0; j < numParts; j++) {
	    // printf("Searching particle %d (%s) =? %s\n", j, part[j].name.val(), s);
	    if (strcmp(s,part[j].name.val()) == 0)
		return j;
	}
	return -1;
    }

	void copyToCUDA();

	// Output variables
	Vector3 sysDim;
	BaseGrid* sys;
	// temporary variables
	Vector3 origin, size, basis1, basis2, basis3;


	bool loadedCoordinates;
        bool loadedMomentum;

	// Device Variables
	//int *type_d;
	BrownianParticleType **part_d;
	BaseGrid *sys_d, *kTGrid_d;
	//Bond* bonds_d;
	//int2* bondMap_d;
	//Exclude* excludes_d;
	//int2* excludeMap_d;
	//Angle* angles_d;
	//Dihedral* dihedrals_d;

	// number of simulations
	int simNum;

	// Particle variables
	String* partsFromFile;
	int* indices;
	int numPartsFromFile;
	Bond* bonds;
	int numCap; // max number of particles
	int num; // current number of particles
    int num_rb_attached_particles;
        Vector3* pos; //  position of each particle
        Vector3* momentum; //momentum of each brownian particles Han-Yi Chou
        Vector3  COM_Velocity; //center of mass velocity Han-Yi Chou
	int* type; // type of each particle
	int* serial; // serial number of each particle
	int currSerial; // the serial number of the next new particle
	String* name; // name of each particle
	Vector3* posLast; // used for current computation
        Vector3* momLast; //used for Lagevin dynamics
	float timeLast; // used with posLast
	float minimumSep; // minimum separation allowed with placing new particles


	// RigidBody variables
	/* int numRB; */
	/* std::vector< std::vector<RigidBody> > rbs; */
	
	// System parameters
	String outputName;
	float timestep;
	long int steps;
	long int seed;
	// String kTGridFile;
	String temperatureGridFile;
	String inputCoordinates;
        String inputMomentum; //Han-Yi Chou
	String inputRBCoordinates;
	String restartRBCoordinates;
	int copyReplicaCoordinates;
	String restartCoordinates;
        String restartMomentum; //Han-Yi Chou
	int numberFluct;
	int interparticleForce;
	int tabulatedPotential;
	int fullLongRange;
	float kT;
	float temperature;
	float coulombConst;
	float electricField;
	float cutoff;
	float pairlistDistance;
	float switchLen;
	float imdForceScale;
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
	unsigned long randoSeed;

	// Other parameters.
	int rigidBodyGridGridPeriod;
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
	int numBondAngles;
	int numRestraints;
	int* numPartsOfType;
	String partFile;
	String bondFile;
	String excludeFile;
	String angleFile;
	String dihedralFile;
	String restraintFile;
	String bondAngleFile;
	bool readPartsFromFile;
	bool readGroupSitesFromFile;
	bool readBondsFromFile;
	bool readExcludesFromFile;
	bool readAnglesFromFile;
	bool readDihedralsFromFile;
	bool readBondAnglesFromFile;
	bool readRestraintsFromFile;
	//String* partGridFile;
	String **partGridFile;
	//float* partGridFileScale;
	float **partGridFileScale;
        //int *numPartGridFiles;
    std::map<std::string,BaseGrid> part_grid_dictionary;
    std::map<std::string,BaseGrid*> part_grid_dictionary_d;
	std::vector< std::vector<String> > partRigidBodyGrid;
	String* partDiffusionGridFile;
	String* partForceXGridFile;
	String* partForceYGridFile;
	String* partForceZGridFile;
	String* partTableFile;
	String* partReservoirFile;
	int* partTableIndex0;
	int* partTableIndex1;

	String groupSiteFile;
	int numGroupSites;
	std::vector<std::vector<int>> groupSiteData;

	String* bondTableFile;
	int numTabBondFiles;
	int2* bondMap;
	
	Exclude* excludes;
	int2* excludeMap;
	String excludeRule;
	int excludeCapacity;

	Angle* angles;
	String* angleTableFile;
	int numTabAngleFiles;

	Dihedral* dihedrals;
	String* dihedralTableFile;
	int numTabDihedralFiles;

	BondAngle* bondAngles;

	Restraint* restraints;

	void readProductPotentials();
	String productPotentialFile;
	int numProductPotentials;
	bool readProductPotentialsFromFile;
        ProductPotentialConf* productPotentials;
	XpotMap simple_potential_ids;
        std::vector<SimplePotential> simple_potentials;

        //Han-Yi Chou
        String ParticleDynamicType;
        String RigidBodyDynamicType;
        String ParticleLangevinIntegrator;
	// RigidBody parameters.
	RigidBodyType* rigidBody;
	int numRigidTypes;
        int ParticleInterpolationType;
        int RigidBodyInterpolationType;
};

#endif
