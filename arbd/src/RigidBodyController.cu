/* #include "RigidBody.h" */
#include <iomanip>
#include "RigidBodyController.h"
#include "Configuration.h"
#include "RigidBodyType.h"
#include "RigidBodyGrid.h"
#include "ComputeGridGrid.cuh"

// #include "GPUManager.h"
// GPUManager RigidBodyController::gpuman;

#include <cuda_profiler_api.h>

// #include <vector>
#include "Debug.h"

#include "RandomCPU.h"							/* RBTODO: fix this? */

#ifndef gpuErrchk
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, String file, int line, bool abort=true) {
   if (code != cudaSuccess) {
       fprintf(stderr,"CUDA Error: %s %s %d\n", cudaGetErrorString(code), __FILE__, line);
      if (abort) exit(code);
   }
}
#endif

// allocate and initialize an array of stream handles
cudaStream_t *RigidBodyForcePair::stream = (cudaStream_t *) malloc(NUMSTREAMS * sizeof(cudaStream_t));
int RigidBodyForcePair::nextStreamID = 0;        /* used during stream init */
int RigidBodyForcePair::lastRbGridID = -1; /* used to schedule kernel interaction */
RigidBodyForcePair* RigidBodyForcePair::lastRbForcePair = NULL;
/* #include <cuda.h> */
/* #include <cuda_runtime.h> */
/* #include <curand_kernel.h> */

RigidBodyController::RigidBodyController(const Configuration& c, const char* prefix, unsigned long int seed, int repID) : conf(c)
{
        char str[8];
        sprintf(str, "%d", repID);
        strcpy(outArg, prefix);
        strcat(outArg, ".");
        strcat(outArg, str);

	gpuErrchk(cudaDeviceSynchronize()); /* RBTODO: this should be extraneous */
	construct_grids();
	for (int i = 0; i < conf.numRigidTypes; i++)
	    conf.rigidBody[i].initializeParticleLists();

	int numRB = 0;
	// grow list of rbs

	int attached_particle_offset = 0;
	for (int i = 0; i < conf.numRigidTypes; i++) {			
		numRB += conf.rigidBody[i].num;
		int attached_particle_in_type = conf.rigidBody[i].num_attached_particles();
		std::vector<RigidBody> tmp;
		// RBTODO: change conf.rigidBody to conf.rigidBodyType
		const int jmax = conf.rigidBody[i].num;
		for (int j = 0; j < jmax; j++) {
			String name = conf.rigidBody[i].name;
			if (jmax > 1) {
			    char stmp[128];
			    snprintf(stmp, 128, "#%d", j);
			    name.add( stmp );
			}
			RigidBody r(name, conf, conf.rigidBody[i], this,
				    attached_particle_offset, attached_particle_offset+attached_particle_in_type);
			attached_particle_offset += attached_particle_in_type;

			int nb = r.appendNumParticleBlocks( &particleForceNumBlocks );
			tmp.push_back( r );
		}
		rigidBodyByType.push_back(tmp);
	}
	attached_particle_forces = new Vector3[attached_particle_offset];

	totalParticleForceNumBlocks = 0;
	for (int i=0; i < particleForceNumBlocks.size(); ++i) {
	    particleForce_offset.push_back(2*totalParticleForceNumBlocks);
	    totalParticleForceNumBlocks += particleForceNumBlocks[i];
	}

	gpuErrchk(cudaMallocHost(&(particleForces), sizeof(ForceEnergy) * 2*totalParticleForceNumBlocks))
	gpuErrchk(cudaMalloc(&(particleForces_d), sizeof(ForceEnergy) * 2*totalParticleForceNumBlocks))

	if (conf.restartRBCoordinates.length() > 0)
	    load_restart_coordinates(conf.restartRBCoordinates.val());
	else if (conf.inputRBCoordinates.length() > 0)
	    loadRBCoordinates(conf.inputRBCoordinates.val());
	
	random = new RandomCPU(conf.seed + repID + 1); /* +1 to avoid using same seed as RandomCUDA */
	
	
	initializeForcePairs();	// Must run after construct_grids()
	gpuErrchk(cudaDeviceSynchronize()); /* RBTODO: this should be extraneous */
}

RigidBodyController::~RigidBodyController() {
	for (int i = 0; i < rigidBodyByType.size(); i++)
		rigidBodyByType[i].clear();
	rigidBodyByType.clear();

	delete [] attached_particle_forces;
	delete random;
}

struct GridKey {
    String name;
    float scale;
    GridKey(const String& name, const float& scale) :
	name(name), scale(scale) { }

    bool operator==(const GridKey& o) const { return name == o.name && scale == o.scale; }
};

void RigidBodyController::construct_grids() {
    // typedef std::tuple<String, float> GridKey;
    
    // Build dictionary to reuse grids across all types, first finding scale factors
    std::vector<GridKey> all_files;
    std::vector<GridKey>::iterator itr;

    for (int t_idx = 0; t_idx < conf.numRigidTypes; ++t_idx)
    {
	// TODO: don't duplicate the code below three times
	RigidBodyType& t = conf.rigidBody[t_idx];
	t.RBC = this;

	t.numPotGrids = t.potentialGridFiles.size();
	t.numDenGrids = t.densityGridFiles.size();
	t.numPmfs = t.pmfFiles.size();

	t.potential_grid_idx = new size_t[t.numPotGrids]; // TODO; don't allocate here
	t.density_grid_idx = new size_t[t.numDenGrids]; // TODO; don't allocate here
	t.pmf_grid_idx = new size_t[t.numPmfs]; // TODO; don't allocate here
	for (size_t i = 0; i < t.potentialGridFiles.size(); ++i)
	{

	    String& filename = t.potentialGridFiles[i];
	    String& name = t.potentialGridKeys[i];
	    float scale = 1.0f;
	    for (size_t j = 0; j < t.potentialGridScaleKeys.size(); ++j)
	    {
		if (name == t.potentialGridScaleKeys[j])
		    scale = t.potentialGridScale[j];
	    }

	    GridKey key = GridKey(filename, scale);
	    size_t key_idx;
	    // Find key if it exists
	    itr = std::find(all_files.begin(), all_files.end(), key);
	    if (itr == all_files.end())
	    {
		key_idx = all_files.size();
		all_files.push_back( key );
	    }
	    else 
	    {
		key_idx = std::distance(all_files.begin(), itr);
	    }

	    // Assign index into all_files to RigidBodyType
	    t.potential_grid_idx[i] = key_idx;

	}

	// Density
	for (size_t i = 0; i < t.densityGridFiles.size(); ++i)
	{

	    String& filename = t.densityGridFiles[i];
	    String& name = t.densityGridKeys[i];
	    float scale = 1.0f;
	    for (size_t j = 0; j < t.densityGridScaleKeys.size(); ++j)
	    {
		if (name == t.densityGridScaleKeys[j])
		    scale = t.densityGridScale[j];
	    }

	    GridKey key = GridKey(filename, scale);
	    size_t key_idx;
	    // Find key if it exists
	    itr = std::find(all_files.begin(), all_files.end(), key);
	    if (itr == all_files.end())
	    {
		key_idx = all_files.size();
		all_files.push_back( key );
	    }
	    else 
	    {
		key_idx = std::distance(all_files.begin(), itr);
	    }

	    // Assign index into all_files to RigidBodyType
	    t.density_grid_idx[i] = key_idx;
	}

	//PMF	
	for (size_t i = 0; i < t.pmfFiles.size(); ++i)
	{

	    String& filename = t.pmfFiles[i];
	    String& name = t.pmfKeys[i];
	    float scale = 1.0f;
	    for (size_t j = 0; j < t.pmfScaleKeys.size(); ++j)
	    {
		if (name == t.pmfScaleKeys[j])
		    scale = t.pmfScale[j];
	    }

	    GridKey key = GridKey(filename, scale);
	    size_t key_idx;
	    // Find key if it exists
	    itr = std::find(all_files.begin(), all_files.end(), key);
	    if (itr == all_files.end())
	    {
		key_idx = all_files.size();
		all_files.push_back( key );
	    }
	    else 
	    {
		key_idx = std::distance(all_files.begin(), itr);
	    }

	    // Assign index into all_files to RigidBodyType
	    t.pmf_grid_idx[i] = key_idx;
	}
	
	// TODO: have RBType manage this allocation
	gpuErrchk(cudaMalloc(&t.potential_grid_idx_d, sizeof(size_t)*t.numPotGrids ));
	gpuErrchk(cudaMalloc(&t.density_grid_idx_d, sizeof(size_t)*t.numDenGrids ));
	gpuErrchk(cudaMalloc(&t.pmf_grid_idx_d, sizeof(size_t)*t.numPmfs ));

	gpuErrchk(cudaMemcpy(t.potential_grid_idx_d, t.potential_grid_idx, sizeof(size_t)*t.numPotGrids, cudaMemcpyHostToDevice ));
	gpuErrchk(cudaMemcpy(t.density_grid_idx_d, t.density_grid_idx, sizeof(size_t)*t.numDenGrids, cudaMemcpyHostToDevice ));
	gpuErrchk(cudaMemcpy(t.pmf_grid_idx_d, t.pmf_grid_idx, sizeof(size_t)*t.numPmfs, cudaMemcpyHostToDevice ));

    }
    
    // Store grids 
    grids = new BaseGrid[all_files.size()];
    gpuErrchk(cudaMalloc( &grids_d, sizeof(RigidBodyGrid)*all_files.size() ));
    
    // Read and scale grids, then copy to GPU
    for (size_t i = 0; i < all_files.size(); ++i)
    {
	GridKey& key = all_files[i];
	BaseGrid& g0 = grids[i];
	g0 = BaseGrid(key.name);
	g0.scale(key.scale);

	RigidBodyGrid g = RigidBodyGrid();
	g.nx = g0.nx;
	g.ny = g0.ny;
	g.nz = g0.nz;
	g.size = g0.size;
	g.val = g0.val;
	   
	// Copy to GPU, starting with grid data
	float* tmp;
	size_t sz = sizeof(float) * g.getSize();
	gpuErrchk(cudaMalloc( &tmp, sz)); 
	gpuErrchk(cudaMemcpy( tmp, g.val, sz, cudaMemcpyHostToDevice));

	// Set grid pointer to device 
	g.val = tmp;

	// Copy grid
	sz = sizeof(RigidBodyGrid);
	// gpuErrchk(cudaMalloc(&ptr_d, sz));
	gpuErrchk(cudaMemcpy(&grids_d[i], &g, sz, cudaMemcpyHostToDevice));

	// Restore pointer
	g.val = NULL;
	tmp = NULL;
    }
}	

void RigidBodyController::destruct_grids() {
    // TODO

}

bool RigidBodyController::loadRBCoordinates(const char* fileName) {
	char line[STRLEN];
	FILE* inp = fopen(fileName, "r");

	if (inp == NULL) {
		printf("GrandBrownTown: load RB coordinates: File '%s' does not exist\n", fileName);
		exit(-1);	   
	}

	int imax = rigidBodyByType.size();
	int i = 0;
	int jmax = rigidBodyByType[i].size();
	int j = 0;

	while (fgets(line, STRLEN, inp) != NULL) {
		// Ignore comments.
		int len = strlen(line);
		if (line[0] == '#') continue;
		if (len < 2) continue;

		String s(line);
		int numTokens = s.tokenCount();
		if (numTokens < 3+9) {
			printf("GrandBrownTown: load RB coordinates: Invalid coordinate file line: %s\n", line);
			fclose(inp);	
			exit(-1);
		}
                if(conf.RigidBodyDynamicType == String("Langevin") && numTokens < 18)
                {
                    std::cout << "Warning the initial momentum set by random number" << std::endl;
                }

		String* tokenList = new String[numTokens];
		s.tokenize(tokenList);
		if (tokenList == NULL) {
			printf("GrandBrownTown: load RB coordinates: Invalid coordinate file line: %s\n", line);
			fclose(inp);
			exit(-1);
		}

		RigidBody& rb = rigidBodyByType[i][j];
		rb.position = Vector3(
			(float) strtod(tokenList[0],NULL), (float) strtod(tokenList[1],NULL), (float) strtod(tokenList[2],NULL));
		rb.orientation = Matrix3(
			(float) strtod(tokenList[3],NULL), (float) strtod(tokenList[4],NULL), (float) strtod(tokenList[5],NULL),
			(float) strtod(tokenList[6],NULL), (float) strtod(tokenList[7],NULL), (float) strtod(tokenList[8],NULL),
			(float) strtod(tokenList[9],NULL), (float) strtod(tokenList[10],NULL), (float) strtod(tokenList[11],NULL));

	        if(conf.RigidBodyDynamicType == String("Langevin") && numTokens >= 18)
                {
                    rb.momentum = Vector3((float)strtod(tokenList[12],NULL), (float) strtod(tokenList[13],NULL), (float) strtod(tokenList[14],NULL));
                    rb.angularMomentum = Vector3((float)strtod(tokenList[15],NULL), (float) strtod(tokenList[16],NULL), (float) strtod(tokenList[17],NULL));
                }

	        if(conf.RigidBodyDynamicType == String("Langevin") && numTokens < 18)
		    rb.Boltzmann(conf.seed);
               
		delete[] tokenList;

		j++;
		if (j == jmax) {
			i++;
			if (i == imax)
				break;
			j=0;
			jmax = rigidBodyByType[i].size();
		}
	}
	fclose(inp);
	return true;
}
bool RigidBodyController::load_restart_coordinates(const char* filename) {
	char line[STRLEN];
	FILE* inp = fopen(filename, "r");

	if (inp == NULL) {
		printf("ARBD: load RB coordinates: File '%s' does not exist\n", filename);
		exit(-1);	   
	}

	int imax = rigidBodyByType.size();
	int i = 0;
	int jmax = rigidBodyByType[i].size();
	int j = 0;

	while (fgets(line, STRLEN, inp) != NULL) {
		// Ignore comments.
		int len = strlen(line);
		if (line[0] == '#') continue;
		if (len < 2) continue;

		String s(line);
		int numTokens = s.tokenCount();
		if (numTokens < 2+3+9) {
			printf("ARBD: invalid RB restart coordinate line: %s\n", line);
			fclose(inp);	
			exit(-1);
		}
                if(conf.RigidBodyDynamicType == String("Langevin") && numTokens < 20)
                {
                    std::cout << "Warning the initial momentum set by random number" << std::endl;
                }

		String* tokenList = new String[numTokens];
		s.tokenize(tokenList);
		if (tokenList == NULL) {
			printf("ARBD: invalid RB restart coordinate line: %s\n", line);
			fclose(inp);
			exit(-1);
		}

		RigidBody& rb = rigidBodyByType[i][j];
		printf("Assinging positions %d %d %f\n",i,j, (float) strtod(tokenList[2],NULL));
		rb.position = Vector3(
			(float) strtod(tokenList[2],NULL), (float) strtod(tokenList[3],NULL), (float) strtod(tokenList[4],NULL));
		rb.orientation = Matrix3(
			(float) strtod(tokenList[5],NULL), (float) strtod(tokenList[6],NULL), (float) strtod(tokenList[7],NULL),
			(float) strtod(tokenList[8],NULL), (float) strtod(tokenList[9],NULL), (float) strtod(tokenList[10],NULL),
			(float) strtod(tokenList[11],NULL), (float) strtod(tokenList[12],NULL), (float) strtod(tokenList[13],NULL));

	        if(conf.RigidBodyDynamicType == String("Langevin") && numTokens >= 20)
                {
                    rb.momentum = Vector3((float)strtod(tokenList[14],NULL), (float) strtod(tokenList[15],NULL), (float) strtod(tokenList[16],NULL));
                    rb.angularMomentum = Vector3((float)strtod(tokenList[17],NULL), (float) strtod(tokenList[18],NULL), (float) strtod(tokenList[19],NULL));
                }

	        if(conf.RigidBodyDynamicType == String("Langevin") && numTokens < 20)
		    rb.Boltzmann(conf.seed);

		delete[] tokenList;

		j++;
		if (j == jmax) {
			i++;
			if (i == imax)
				break;
			j=0;
			jmax = rigidBodyByType[i].size();
		}
	}
	fclose(inp);
	if (i < imax || j < jmax) {
	    printf("ARBD: RB restart file did not contain the correct number of lines: %s\n", filename);
	    exit(-1);
	}
	return true;
}

		

void RigidBodyController::initializeForcePairs() {
	// Loop over all pairs of rigid body types
	//   the references here make the code more readable, but they may incur a performance loss
	RigidBodyForcePair::createStreams();
	// printf("Initializing force pairs\n");
	for (int ti = 0; ti < conf.numRigidTypes; ti++) {
		RigidBodyType& t1 = conf.rigidBody[ti];
		for (int tj = ti; tj < conf.numRigidTypes; tj++) {
			RigidBodyType& t2 = conf.rigidBody[tj];


			const std::vector<String>& keys1 = t1.densityGridKeys; 
			const std::vector<String>& keys2 = t2.potentialGridKeys;

			// printf("  Working on type pair ");
			// t1.name.printInline(); printf(":"); t2.name.print();
			
			// Loop over all pairs of grid keys (e.g. "Elec")
			std::vector<int> gridKeyId1;
			std::vector<int> gridKeyId2;
			
			// printf("  Grid keys %d:%d\n",keys1.size(),keys2.size());

			bool paired = false;
			for(int k1 = 0; k1 < keys1.size(); k1++) {
				for(int k2 = 0; k2 < keys2.size(); k2++) {
				    // printf("    checking grid keys ");
				    //	keys1[k1].printInline(); printf(":"); keys2[k2].print();
					
					if ( keys1[k1] == keys2[k2] ) {
						gridKeyId1.push_back( t1.density_grid_idx[k1] );
						gridKeyId2.push_back( t2.potential_grid_idx[k2] );
						paired = true;
					}
				}
			}
			
			if (paired) {
				// found matching keys => calculate force between all grid pairs
				std::vector<RigidBody>& rbs1 = rigidBodyByType[ti];
				std::vector<RigidBody>& rbs2 = rigidBodyByType[tj];

				// Loop over rigid bodies of these types
				for (int i = 0; i < rbs1.size(); i++) {
					for (int j = (ti==tj ? i+1 : 0); j < rbs2.size(); j++) {
						RigidBody* rb1 = &(rbs1[i]);
						RigidBody* rb2 = &(rbs2[j]);

						// printf("    pushing RB force pair for %d:%d\n",i,j);
						RigidBodyForcePair fp = RigidBodyForcePair(&(t1),&(t2),rb1,rb2,gridKeyId1,gridKeyId2, false, conf.rigidBodyGridGridPeriod );
						gpuErrchk(cudaDeviceSynchronize()); /* RBTODO: this should be extraneous */
						forcePairs.push_back( fp ); 
						// printf("    done pushing RB force pair for %d:%d\n",i,j);
					}
				}
			}
		}
	}

	// add Pmfs (not a true pairwise RB interaction; hacky implementation)
	for (int ti = 0; ti < conf.numRigidTypes; ti++) {
		RigidBodyType& t1 = conf.rigidBody[ti];

		const std::vector<String>& keys1 = t1.densityGridKeys; 
		const std::vector<String>& keys2 = t1.pmfKeys;
		std::vector<int> gridKeyId1;
		std::vector<int> gridKeyId2;
		
		// Loop over all pairs of grid keys (e.g. "Elec")
		bool paired = false;
		for(int k1 = 0; k1 < keys1.size(); k1++) {
			for(int k2 = 0; k2 < keys2.size(); k2++) {
				if ( keys1[k1] == keys2[k2] ) {
				    gridKeyId1.push_back( t1.density_grid_idx[k1] );
				    gridKeyId2.push_back( t1.pmf_grid_idx[k2] );
				    paired = true;
				}
			}
		}	
		if (paired) {
			// found matching keys => calculate force between all grid pairs
			std::vector<RigidBody>& rbs1 = rigidBodyByType[ti];
			
			// Loop over rigid bodies of these types
			for (int i = 0; i < rbs1.size(); i++) {
					RigidBody* rb1 = &(rbs1[i]);
					RigidBodyForcePair fp = RigidBodyForcePair(&(t1),&(t1),rb1,rb1,gridKeyId1,gridKeyId2, true, conf.rigidBodyGridGridPeriod);
					gpuErrchk(cudaDeviceSynchronize()); /* RBTODO: this should be extraneous */
					forcePairs.push_back( fp ); 
			}
		}
	}

	// Initialize device data for RB force pairs after std::vector is done growing
	for (int i = 0; i < forcePairs.size(); i++)
		forcePairs[i].initialize();
			
}

void RigidBodyController::update_attached_particle_positions(Vector3* pos_d, Vector3* force_d, float* energy_d, BaseGrid* sys_d, int num, int num_rb_attached_particles, int numReplicas) {
    for (int i = 0; i < rigidBodyByType.size(); i++) {
	for (int j = 0; j < rigidBodyByType[i].size(); j++) {
	    rigidBodyByType[i][j].update_particle_positions(pos_d, force_d, energy_d);
	}
    }
    gpuErrchk(cudaMemset((void*) force_d, 0, num_rb_attached_particles*sizeof(Vector3)));
    gpuErrchk(cudaMemset((void*) energy_d, 0, num_rb_attached_particles*sizeof(float)));
}


void RigidBodyController::updateParticleLists(Vector3* pos_d, BaseGrid* sys_d) {
	for (int i = 0; i < rigidBodyByType.size(); i++) {
		for (int j = 0; j < rigidBodyByType[i].size(); j++) {
			rigidBodyByType[i][j].updateParticleList(pos_d, sys_d);
		}
	}
}

void RigidBodyController::clearForceAndTorque()
{
    // clear old forces
    for (int i = 0; i < rigidBodyByType.size(); i++) 
    {
        for (int j = 0; j < rigidBodyByType[i].size(); j++) 
        {
            RigidBody& rb = rigidBodyByType[i][j];
            rb.clearForce();
            rb.clearTorque();
        }
    }
}

void RigidBodyController::updateForces(Vector3* pos_d, Vector3* force_d, int s, float* energy, bool get_energy, int scheme, BaseGrid* sys, BaseGrid* sys_d, int num, int num_rb_attached_particles)
{
	//if (s <= 1)
		//gpuErrchk( cudaProfilerStart() );
	
	// Grid–particle forces	
	int pfo_idx = 0;
	for (int i = 0; i < rigidBodyByType.size(); i++) {
		for (int j = 0; j < rigidBodyByType[i].size(); j++) {
			RigidBody& rb = rigidBodyByType[i][j];
			rb.callGridParticleForceKernel( pos_d, force_d, s, energy, get_energy, scheme, sys, sys_d, particleForces_d, particleForce_offset, pfo_idx );
		}
	}

	// RBTODO: launch kernels ahead of time and sync using event and memcpyAsync 
	gpuErrchk( cudaDeviceSynchronize() );
	cudaMemcpy(particleForces, particleForces_d, sizeof(ForceEnergy)*2*totalParticleForceNumBlocks, cudaMemcpyDeviceToHost);
	cudaMemcpyAsync(attached_particle_forces, &force_d[num], sizeof(Vector3)*(num_rb_attached_particles), cudaMemcpyDeviceToHost);

	pfo_idx=0;
	for (int i = 0; i < rigidBodyByType.size(); i++) {
		for (int j = 0; j < rigidBodyByType[i].size(); j++) {
			RigidBody& rb = rigidBodyByType[i][j];
			rb.applyGridParticleForces(sys, particleForces, particleForce_offset, pfo_idx);
;
		}
	}

	{
	    gpuErrchk( cudaDeviceSynchronize() );
	    for (auto &rbv: rigidBodyByType) {
		for (auto &rb: rbv) {
		    rb.apply_attached_particle_forces( attached_particle_forces );
		}
	    }
	}


	// Grid–Grid forces
	if ( ((s % conf.rigidBodyGridGridPeriod) == 0 || s == 1 ) && forcePairs.size() > 0) {
		for (int i=0; i < forcePairs.size(); i++) {
			// TODO: performance: make this check occur less frequently
		    if (forcePairs[i].isOverlapping(sys)) {
				forcePairs[i].callGridForceKernel(i, s, scheme, sys_d);
		    }
		}
		
		// each kernel call is followed by async memcpy for previous; now get last
		RigidBodyForcePair* fp = RigidBodyForcePair::lastRbForcePair;
                if(RigidBodyForcePair::lastRbGridID >= 0)
                {
		    fp->retrieveForcesForGrid( fp->lastRbGridID );
		    fp->lastRbGridID = -1;
                }
		// stream sync was slower than device sync
		/* for (int i = 0; i < NUMSTREAMS; i++) { */
		/* 	const cudaStream_t &s = RigidBodyForcePair::stream[i]; */
		/* 	gpuErrchk(cudaStreamSynchronize( s ));  */
		/* } */
		gpuErrchk(cudaDeviceSynchronize());
		for (int i=0; i < forcePairs.size(); i++)
			if (forcePairs[i].isOverlapping(sys))
				forcePairs[i].processGPUForces(sys);
	}
}

void RigidBodyController::SetRandomTorques()
{
    for (int i = 0; i < rigidBodyByType.size(); i++)
    {
        for (int j = 0; j < rigidBodyByType[i].size(); j++)
        {
            RigidBody& rb = rigidBodyByType[i][j];
            rb.W1 = random->gaussian_vector();
            rb.W2 = random->gaussian_vector();
        }
    }           
}

void RigidBodyController::AddLangevin()
{
    for (int i = 0; i < rigidBodyByType.size(); i++)
    {
        for (int j = 0; j < rigidBodyByType[i].size(); j++)
        {
            RigidBody& rb = rigidBodyByType[i][j];

            //printf("%f %f %f\n",rb.W1.x,rb.W1.y,rb.W1.z);
            //printf("%f %f %f\n",rb.W2.x,rb.W2.y,rb.W2.z);

            rb.addLangevin(rb.W1,rb.W2);
        }
    }
}

void RigidBodyController::integrateDLM(BaseGrid* sys, int step) 
{
    // tell RBs to integrate
    for (int i = 0; i < rigidBodyByType.size(); i++) 
    {
        for (int j = 0; j < rigidBodyByType[i].size(); j++) 
        {
            RigidBody& rb = rigidBodyByType[i][j];
            rb.integrateDLM(sys, step);
        }
    }
}


//Chris original part for Brownian motion
void RigidBodyController::integrate(BaseGrid* sys, int step) 
{
 	// tell RBs to integrate
	if ( step % conf.outputPeriod == 0 ) 
        {       /* PRINT & INTEGRATE */
		if (step == 0) 
                {	// first step so only start this cycle
			print(step);
			for (int i = 0; i < rigidBodyByType.size(); i++)
                        {
				for (int j = 0; j < rigidBodyByType[i].size(); j++)
                                {
					RigidBody& rb = rigidBodyByType[i][j];
					rb.integrate(sys, 0);	
				}
			}
		} 
                else 
                {       // finish last cycle
			for (int i = 0; i < rigidBodyByType.size(); i++)
                        {
				for (int j = 0; j < rigidBodyByType[i].size(); j++)
                                {
					RigidBody& rb = rigidBodyByType[i][j];
					rb.integrate(sys, 1);	
				}
			}
			//print(step);

			// start this cycle
			/*for (int i = 0; i < rigidBodyByType.size(); i++) {
				for (int j = 0; j < rigidBodyByType[i].size(); j++) {
					RigidBody& rb = rigidBodyByType[i][j];
					rb.integrate(sys, 0);	
				}
			}*/
		}
	} 
        else 
        {	/* INTEGRATE ONLY */
		if (step == 0) 
                {		// first step so only start this cycle
			print(step);
			for (int i = 0; i < rigidBodyByType.size(); i++)
                        {
				for (int j = 0; j < rigidBodyByType[i].size(); j++)
                                {
					RigidBody& rb = rigidBodyByType[i][j];
					rb.integrate(sys, 0);	
				}
			}
		} 
                else 
                {       // integrate end of last step and start of this one
			for (int i = 0; i < rigidBodyByType.size(); i++) 
                        {
				for (int j = 0; j < rigidBodyByType[i].size(); j++)
                                {
					RigidBody& rb = rigidBodyByType[i][j];
					rb.integrate(sys, 2);	
				}
			}
		}
	}
}

void RigidBodyController::KineticEnergy()
{
    //float e = 0.;
    //int num = 0;
    for (int i = 0; i < rigidBodyByType.size(); i++) 
    {
        for (int j = 0; j < rigidBodyByType[i].size(); j++) 
        {
            RigidBody& rb = rigidBodyByType[i][j];
            rb.setKinetic(rb.Temperature());
            //rb.kinetic=tmp;
            //e += tmp;
            //num += 1;
        }
    }
    //return e;
    /*if(num > 0)
        return e / num;
    else
        return 0.;*/
}

void RigidBodyForcePair::createStreams() {
	for (int i = 0; i < NUMSTREAMS; i++)
		gpuErrchk( cudaStreamCreate( &(stream[i]) ) );
		// gpuErrchk( cudaStreamCreateWithFlags( &(stream[i]) , cudaStreamNonBlocking ) );
}
bool RigidBodyForcePair::isOverlapping(BaseGrid* sys) const {
	if (isPmf) return true;
	// float pairlistDist = 2.0f; /* TODO: get from conf */
	float rbDist = sys->wrapDiff((rb1->getPosition() - rb2->getPosition())).length();
	for (int i = 0; i < gridKeyId1.size(); ++i) {
		const int k1 = gridKeyId1[i];
		const int k2 = gridKeyId2[i];
		float d1 = type1->RBC->grids[k1].getRadius() + type1->RBC->grids[k1].getCenter().length();
		float d2 = type2->RBC->grids[k2].getRadius() + type2->RBC->grids[k2].getCenter().length();
		if (rbDist < d1+d2)
			return true;
	}
	return false;
}
Vector3 RigidBodyForcePair::getOrigin1(const int i) {
	const int k1 = gridKeyId1[i];
	return rb1->transformBodyToLab( type1->RBC->grids[k1].getOrigin() );
}
Vector3 RigidBodyForcePair::getOrigin2(const int i) {
	const int k2 = gridKeyId2[i];
	Vector3 o = type2->RBC->grids[k2].getOrigin();
	if (!isPmf)
	    return rb2->transformBodyToLab( o );
	else
	    return o;
}		
Vector3 RigidBodyForcePair::getCenter2(const int i) {
    Vector3 c;
    if (!isPmf)
	c = rb2->getPosition();
    else {
	const int k2 = gridKeyId2[i];
	Vector3 o = type2->RBC->grids[k2].getCenter();
    }
    return c;
}
Matrix3 RigidBodyForcePair::getBasis1(const int i) {
	const int k1 = gridKeyId1[i];
	return rb1->getOrientation()*type1->RBC->grids[k1].getBasis();
}
Matrix3 RigidBodyForcePair::getBasis2(const int i) {
	const int k2 = gridKeyId2[i];
	Matrix3 b = type2->RBC->grids[k2].getBasis();
	if (!isPmf)
	    return rb2->getOrientation()*b;
	else
	    return b;
}

// RBTODO: bundle several rigidbodypair evaluations in single kernel call
void RigidBodyForcePair::callGridForceKernel(int pairId, int s, int scheme, BaseGrid* sys_d) 
{
	// get the force/torque between a pair of rigid bodies
	/* printf("  Updating rbPair forces\n"); */
	const int numGrids = gridKeyId1.size();

	/* if (s%10 != 0) */
	/* 	pairId = -1000; */

	// RBTODO: precompute certain common transformations and pass in kernel call
	for (int i = 0; i < numGrids; i++) {
		const int nb = numBlocks[i];
		const int k1 = gridKeyId1[i];
		const int k2 = gridKeyId2[i];
		const cudaStream_t &s = gpuman.gpus[0].get_stream(streamID[i]);

		/*
			ijk: index of grid value
			r: postion of point ijk in real space
			B: grid Basis
			o: grid origin
			R: rigid body orientation
			c: rigid body center

			B': R.B 
			c': R.o + c

  		/.––––––––––––––––––.
	  	| r = R.(B.ijk+o)+c |
	  	| r = B'.ijk + c'   |
	  	`––––––––––––––––––./
		*/
		Matrix3 B1 = getBasis1(i);
		// Vector3 c = getOrigin1(i) - getOrigin2(i);
		Vector3 center_u = getCenter2(i);
		Matrix3 B2 = getBasis2(i).inverse();
                
		// RBTODO: get energy
		if (!isPmf) {								/* pair of RBs */
			computeGridGridForce<<< nb, NUMTHREADS, 2*sizeof(ForceEnergy)*NUMTHREADS, s>>>
				(&type1->RBC->grids_d[k1], &type2->RBC->grids_d[k2],
				 B1, B2, getOrigin1(i) - center_u, center_u - getOrigin2(i),
				 forces_d[i], torques_d[i], scheme, sys_d);
		} else {										/* RB with a PMF */
			computeGridGridForce<<< nb, NUMTHREADS, 2*sizeof(ForceEnergy)*NUMTHREADS, s>>>
				(&type1->RBC->grids_d[k1], &type2->RBC->grids_d[k2],
				 B1, B2, getOrigin1(i) - center_u, center_u-getOrigin2(i),
				 forces_d[i], torques_d[i], scheme, sys_d);
		}
		// retrieveForcesForGrid(i); // this is slower than approach below, unsure why
		
		if (lastRbGridID >= 0)
			lastRbForcePair->retrieveForcesForGrid(lastRbGridID);
		lastRbForcePair = this;
		lastRbGridID = i;
	}
}

void RigidBodyForcePair::retrieveForcesForGrid(const int i) {
	// i: grid ID (less than numGrids)
        const cudaStream_t &s = gpuman.gpus[0].get_stream(streamID[i]);
	// const int nb = numBlocks[i];
	const int nb = 1;

	gpuErrchk(cudaMemcpyAsync(forces[i], forces_d[i], sizeof(ForceEnergy)*nb, cudaMemcpyDeviceToHost, s));
	gpuErrchk(cudaMemcpyAsync(torques[i], torques_d[i], sizeof(Vector3)*nb, cudaMemcpyDeviceToHost, s));
}
void RigidBodyForcePair::processGPUForces(BaseGrid* sys) {
	
	const int numGrids = gridKeyId1.size();
	Vector3 f = Vector3(0.f);
	Vector3 t = Vector3(0.f);
        float energy = 0.f;
	for (int i = 0; i < numGrids; i++) {
	    // const int nb = numBlocks[i];
	    const int nb = 1;

		//Vector3 tmpF = Vector3(0.0f);
		ForceEnergy tmpF = ForceEnergy(0.f, 0.f);
		Vector3 tmpT = Vector3(0.f);
			
		for (int j = 0; j < nb; j++) {
			tmpF = tmpF + forces[i][j];
			tmpT = tmpT + torques[i][j];
		}
		
		// tmpT is the torque calculated about the origin of density grid
		//   so here we transform torque to be about rb1
		Vector3 o1 = getOrigin1(i);
		tmpT = tmpT - (rb1->getPosition() - o1).cross( tmpF.f );

		// clear forces on GPU
		gpuErrchk(cudaMemset((void*)(forces_d[i]),0,nb*sizeof(ForceEnergy)));
		gpuErrchk(cudaMemset((void*)(torques_d[i]),0,nb*sizeof(Vector3)));

		// sum energies,forces and torques
                energy += tmpF.e;
		f = f + tmpF.f;
		t = t + tmpT;
	}

	f *= updatePeriod;
	t *= updatePeriod;
	
	rb1->addForce( f );
	rb1->addTorque( t );
        if(isPmf)
            rb1->addEnergy( energy );
	//if (!isPmf) {
	else 
        {
		const Vector3 t2 = -t + sys->wrapDiff(rb2->getPosition()-rb1->getPosition()).cross( f );
		rb2->addForce( -f );
		rb2->addTorque( t2 );
                rb1->addEnergy(energy*.5);
                rb2->addEnergy(energy*.5);
	}
        
	// printf("force: %s\n", f.toString().val());
	// printf("torque: %s\n", t.toString().val());
}

void RigidBodyController::print(int step) {
	// modeled after outputExtendedData() in Controller.C
    if (conf.numRigidTypes <= 0) return;
	if ( step >= 0 ) {
		// Write RIGID BODY trajectory file
		if ( step % conf.outputPeriod == 0 ) {
			if ( ! trajFile.rdbuf()->is_open() ) {
	      // open file
			    // printf("OPENING RIGID BODY TRAJECTORY FILE\n");
				// RBTODO: backup_file(simParams->rigidBodyTrajectoryFile);

				char fname[140];
				strcpy(fname,outArg);
				strcat(fname, ".rb-traj");
	      trajFile.open(fname);
				
	      while (!trajFile) {
					/* if ( errno == EINTR ) {
						printf("Warning: Interrupted system call opening RIGIDBODY trajectory file, retrying.\n");
						trajFile.clear();
						trajFile.open(simParams->rigidBodyTrajectoryFile);
						continue;
					}
					*/ 
					//char err_msg[257];
					printf("Error opening RigidBody trajectory file %s",fname);
					exit(1);
	      }
	      trajFile << "# RigidBody trajectory file" << std::endl;
	      printLegend(trajFile);
			}
			// printf("WRITING RIGID BODY COORDINATES AT STEP %d\n",step);
			printData(step,trajFile);
			trajFile.flush();    
		}
                if(step % conf.outputEnergyPeriod == 0)
                {
                
                }
    
		// Write restart File
		/* if ( simParams->restartFrequency && */
		/* 		 ((step % simParams->restartFrequency) == 0) && */
		/* 		 (step != simParams->firstTimestep) )	{ */
		if ( step % conf.outputPeriod == 0 && step != 0 ){
		    // printf("RIGID BODY: WRITING RESTART FILE AT STEP %d\n", step);
			char fname[140];
			strcpy(fname,outArg);
			strcat(fname, ".rigid");
			// RBTODO: NAMD_backup_file(fname,".old"); /*  */
			std::ofstream restartFile(fname);
			while (!restartFile) {
				/* RBTODO 
	      if ( errno == EINTR ) {
					printf("Warning: Interrupted system call opening rigid body restart file, retrying.\n");
					restartFile.clear();
					restartFile.open(fname);
					continue;
	      }
				*/
	      printf("Error opening rigid body restart file %s",fname);
	      exit(1); // NAMD_err(err_msg);
			}
			restartFile << "# RigidBody restart file" << std::endl;
			printLegend(restartFile);
			printData(step,restartFile);
			if (!restartFile) {
	      printf("Error writing rigid body restart file %s",fname);
	      exit(-1); // NAMD_err(err_msg);
			} 
		}
	}
}

void RigidBodyController::printLegend(std::ofstream &file) {
        file << "#$LABELS step RigidBodyKey"
		 << " posX  posY  posZ"
		 << " rotXX rotXY rotXZ"
		 << " rotYX rotYY rotYZ"
		 << " rotZX rotZY rotZZ"
		 << " velX  velY  velZ"
		 << " angVelX angVelY angVelZ" << std::endl;
}
void RigidBodyController::printData(int step,std::ofstream &file) {
	// tell RBs to integrate
	for (int i = 0; i < rigidBodyByType.size(); i++) {
		for (int j = 0; j < rigidBodyByType[i].size(); j++) {
			const RigidBody& rb = rigidBodyByType[i][j];
			
			Vector3 v =  rb.getPosition();
			Matrix3 t =  rb.getOrientation();
			file << step <<" "<< rb.getKey()
					 <<" "<< v.x <<" "<< v.y <<" "<< v.z;
			file << std::setprecision(10) <<" "<< t.exx <<" "<< t.exy <<" "<< t.exz
					 <<" "<< t.eyx <<" "<< t.eyy <<" "<< t.eyz
					 <<" "<< t.ezx <<" "<< t.ezy <<" "<< t.ezz;
			v = rb.getVelocity();
			file << std::setprecision(10) <<" "<< v.x <<" "<< v.y <<" "<< v.z;
			v = rb.getAngularVelocity();
			file << std::setprecision(10) <<" "<< v.x <<" "<< v.y <<" "<< v.z
					 << std::endl;
		}
	}
}

float RigidBodyController::getEnergy(float (RigidBody::*Get)())
{
    float e = 0.f;
    for (int i = 0; i < rigidBodyByType.size(); i++)
    {
        for(int j = 0; j < rigidBodyByType[i].size(); j++) 
        { 
            RigidBody& rb = rigidBodyByType[i][j];
            //e += rb.getKinetic();
            e += (rb.*Get)();
        }
    }
    return e;
}

#if 0
void RigidBodyController::printEnergyData(std::fstream &file)
{
    if(file.is_open())
    {

        for (int i = 0; i < rigidBodyByType.size(); i++) 
        {
            for(int j = 0; j < rigidBodyByType[i].size(); j++)
            {
                const RigidBody& rb = rigidBodyByType[i][j];
                file << "Kinetic Energy " << rb.getKey() << ": " << rb.getKinetic() << " (kT)" << std::endl;
                file << " Potential Energy " << rb.getKey() << ": " << rb.getEnergy() << " (kcal/mol)" << std::endl;
            }
       }
    }
    else
    {
        std::cout << " Error in opening files\n"; 
    }      
}
#endif
int RigidBodyForcePair::initialize() {
    // printf("    Initializing (streams for) RB force pair...\n");

	const int numGrids = gridKeyId1.size();
	// RBTODO assert gridKeysIds are same size 

	// allocate memory for forces/torques
	for (int i = 0; i < numGrids; i++) {
		const int k1 = gridKeyId1[i];
		const int sz = type1->RBC->grids[k1].getSize();
		int nb = sz / NUMTHREADS + ((sz % NUMTHREADS == 0) ? 0:1 );
		streamID.push_back( nextStreamID % NUMSTREAMS );
		nextStreamID++;

		numBlocks.push_back(nb);

		nb = 1;
		//forces.push_back( new Vector3[nb] );
		forces.push_back( new ForceEnergy[nb]);
		torques.push_back( new Vector3[nb] );

		//forces_d.push_back( new Vector3[nb] ); // RBTODO: correct?
		forces_d.push_back( new ForceEnergy[nb]);
		torques_d.push_back( new Vector3[nb] );

		// allocate device memory for numBlocks of torque, etc.
    // printf("      Allocating device memory for forces/torques\n");
		gpuErrchk(cudaMalloc(&(forces_d[i]), sizeof(ForceEnergy) * nb));
		gpuErrchk(cudaMalloc(&(torques_d[i]), sizeof(Vector3) * nb));
	}
	gpuErrchk(cudaDeviceSynchronize());
	// printf("    Done initializing RB force pair\n");
	return nextStreamID;
}

void RigidBodyForcePair::swap(RigidBodyForcePair& a, RigidBodyForcePair& b) {
	using std::swap;
	swap(a.type1, b.type1);
	swap(a.type2, b.type2);
	swap(a.rb1, b.rb1);
	swap(a.rb2, b.rb2);

	swap(a.gridKeyId1, b.gridKeyId1);
	swap(a.gridKeyId2, b.gridKeyId2);

	swap(a.numBlocks, b.numBlocks);

	swap(a.forces,    b.forces);
	swap(a.forces_d,  b.forces_d);
	swap(a.torques,   b.torques);
	swap(a.torques_d, b.torques_d);
}


RigidBodyForcePair::~RigidBodyForcePair() {
    //printf("    Destructing RB force pair\n");
	const int numGrids = gridKeyId1.size();

	// printf("      numGrids = %d\n",numGrids);

	// RBTODO assert gridKeysIds are same size 

	// allocate memory for forces/torques
	if (streamID.size() > 0) {
		for (int i = 0; i < numGrids; i++) {
			const int k1 = gridKeyId1[i];
			const int nb = numBlocks[i];

			// free device memory for numBlocks of torque, etc.
			// printf("      Freeing device memory for forces/torques\n");
			gpuErrchk(cudaFree( forces_d[i] ));	
			gpuErrchk(cudaFree( torques_d[i] ));
		}
	}
	streamID.clear();
	numBlocks.clear();
	forces.clear();
	forces_d.clear();
	torques.clear();
	torques_d.clear();
}



