// RigidBodyType.h (2015)
// Author: Chris Maffeo <cmaffeo2@illinois.edu>

#pragma once
#include <vector>
/* #include <thrust/host_vector.h> */
/* #include <thrust/device_vector.h> */
/* #include "Reservoir.h" */
/* #include "BaseGrid.h" */
/* #include "RigidBodyGrid.h" */
#include "useful.h"

#include <cstdio>

class Reservoir;
class BaseGrid;
class RigidBodyGrid;
class Configuration;
class RigidBodyController;
class RigidBody;

class RigidBodyType {
    friend class RigidBody;
public:
RigidBodyType(const String& name = "", const Configuration* conf = NULL ) :
	name(name), conf(conf), num(0),
		reservoir(NULL), mass(1.0f), inertia(), transDamping(),
		rotDamping(), initPos(), initRot(Matrix3(1.0f)), initMomentum(Vector3(0.f)), initAngularMomentum(Vector3(0.f)),
		numPotGrids(0), numDenGrids(0), numPmfs(0), numParticles(NULL) { }
	~RigidBodyType() { clear(); }
private:
	// Deletes all members
	void clear();
	// void copy(const RigidBodyType& src);

	void addGrid(String s, std::vector<String> &keys, std::vector<String> &files);
	void addScaleFactor(String s, std::vector<String> &keys, std::vector<float> &vals);
	
public:
	/* RigidBodyType& operator=(const RigidBodyType& src); */
	void copyGridsToDevice();
	
    void append_attached_particle_file(String s) { attached_particle_files.push_back(s); }
    void attach_particles();
    size_t num_attached_particles() const { return attached_particle_types.size() ;}
    const std::vector<int>& get_attached_particle_types() const { return attached_particle_types; }

	void addPotentialGrid(String s);
	void addDensityGrid(String s);
	void addPMF(String s);
	void scalePotentialGrid(String s);
	void scaleDensityGrid(String s);
	void scalePMF(String s);

	void setDampingCoeffs(float timestep);

	void initializeParticleLists();
	// TODO: privatize
public:
	String name;
private:
	const Configuration* conf;
	std::vector<String> attached_particle_files;
	std::vector<int>attached_particle_types;
private:
    std::vector<Vector3>attached_particle_positions;
    Vector3* attached_particle_positions_d;

public:
	int num; // number of particles of this type

	Reservoir* reservoir;

	float mass;
	Vector3 inertia;
	Vector3 transDamping;
	Vector3 rotDamping;
	Vector3 transForceCoeff;
	Vector3 rotTorqueCoeff;

	Vector3 initPos;	
	Matrix3 initRot;
        Vector3 initMomentum;
        Vector3 initAngularMomentum;
	
	std::vector<String> potentialGridKeys;
	std::vector<String> densityGridKeys;
	std::vector<String> pmfKeys;

	std::vector<String> potentialGridFiles;
	std::vector<String> densityGridFiles;
	std::vector<String> pmfFiles;

	std::vector<String> potentialGridScaleKeys;
	std::vector<String> densityGridScaleKeys;
	std::vector<String> pmfScaleKeys;

	std::vector<float> potentialGridScale;
	std::vector<float> densityGridScale;
	std::vector<float> pmfScale;
	
	// RBTODO: clear std::vectors after initialization, (but keep offsets)
	// duplicates of std::vector grids for device
public:
	int numPotGrids;
	int numDenGrids;
	int numPmfs;

	int* numParticles;		  /* particles affected by potential grids */
	int** particles;		 	
	int** particles_d;		 	


	/* RigidBodyGrid* rawPotentialGrids; */
	/* RigidBodyGrid* rawDensityGrids; */
	/* BaseGrid* rawPmfs; */
	/* Matrix3* rawPotentialBases; */
	/* Matrix3* rawDensityBases; */
	/* Vector3* rawPotentialOrigins; */
	/* Vector3* rawDensityOrigins;		 */

	
	// device pointers
	/* RigidBodyGrid** rawPotentialGrids_d; */
	/* RigidBodyGrid** rawDensityGrids_d; */
	/* RigidBodyGrid** rawPmfs_d; */

	size_t* potential_grid_idx;
	size_t* density_grid_idx;
	size_t* pmf_grid_idx;

	size_t* potential_grid_idx_d;
	size_t* density_grid_idx_d;
	size_t* pmf_grid_idx_d;
	
	RigidBodyController* RBC;
};
