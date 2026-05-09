#include <cassert>
#include "Configuration.h"
#include "RigidBodyType.h"
#include "Reservoir.h"
#include "BaseGrid.h"
#include "RigidBodyGrid.h"

#ifndef gpuErrchk
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true) {
   if (code != cudaSuccess) {
      fprintf(stderr,"CUDA Error: %s %s %d\n", cudaGetErrorString(code), __FILE__, line);
      if (abort) exit(code);
   }
}
#endif

void RigidBodyType::clear() {
	num = 0;											// RBTODO: not 100% sure about this
	if (reservoir != NULL) delete reservoir;
	reservoir = NULL;
	// pmf = NULL;
	mass = 1.0;

	// TODO: make sure that this actually removes grid data
	potentialGridKeys.clear();
	densityGridKeys.clear();
	pmfKeys.clear();

	
	if (numParticles != NULL) {
		for (int i=0; i < numPotGrids; ++i) {
			printf("CLEARING\n");
			if (numParticles[i] > 0) {
				delete[] particles[i];
				gpuErrchk(cudaFree( particles_d[i] ));
			}
		}
		delete[] numParticles;
		delete[] particles;
		delete[] particles_d;
		numParticles = NULL;
	}
}


// void RigidBodyType::setDampingCoeffs(float timestep, float tmp_mass, Vector3 tmp_inertia, float tmp_transDamping, float tmp_rotDamping) {
void RigidBodyType::setDampingCoeffs(float timestep) { /* MUST ONLY BE CALLED ONCE!!! */
	/*–––––––––––––––––––––––––––––––––––––––––––––––––––––––.
	| DiffCoeff = kT / dampingCoeff mass                     |
	|                                                        |
	| type->DampingCoeff has units of (1/ps)                 |
	|                                                        |
	| f[kcal/mol AA] = - dampingCoeff * momentum[amu AA/ns]  |
	|                                                        |
	| units "(1/ns) * (amu AA/ns)" "kcal_mol/AA" * 2.390e-09 |
	`–––––––––––––––––––––––––––––––––––––––––––––––––––––––*/

	/*––––––––––––––––––––––––––––––––––––––––––––––––––––.
	| < f(t) f(t') > = 2 kT dampingCoeff mass delta(t-t') |
	|                                                     |
	|  units "sqrt( k K (1/ns) amu / ns )" "kcal_mol/AA"  |
	|    * 2.1793421e-06                                  |
	`––––––––––––––––––––––––––––––––––––––––––––––––––––*/
	// RBTODO: make units consistent with rest of RB code 
	float Temp = 295; /* RBTODO: temperature should be read from grid? Or set in uniformly in config file */
	transForceCoeff = 2.1793421e-06 * Vector3::element_sqrt( 2*Temp*mass*transDamping/timestep );

	// setup for langevin
	// langevin = rbParams->langevin;
	// if (langevin) {
	// T = - dampingCoeff * angularMomentum

	// < f(t) f(t') > = 2 kT dampingCoeff inertia delta(t-t')
	rotTorqueCoeff = 2.1793421e-06 *
		Vector3::element_sqrt( 2*Temp* Vector3::element_mult(inertia,rotDamping) / timestep );


	transDamping = 2.3900574e-9 * transDamping;
	rotDamping = 2.3900574e-9 * rotDamping;

}
	
void RigidBodyType::attach_particles() {
    for (const auto& filename: attached_particle_files) {
	const size_t line_char = 256;
	FILE* inp = fopen(filename.val(), "r");
	char line[line_char];

	// If the particle file cannot be found, exit the program
	if (inp == NULL) {
	    printf("ERROR: Could not open `%s'.\n", filename.val());
	    fclose(inp);
	    exit(1);
	}

	// Get and process all lines of input
	while (fgets(line, line_char, inp) != NULL) {
		// Lines in the particle file that begin with # are comments
		if (line[0] == '#') continue;

		String s(line);
		int numTokens = s.tokenCount();

		// Break the line down into pieces (tokens) so we can process them individually
		String* tokenList = new String[numTokens];
		s.tokenize(tokenList);

		// Legitimate GROUP input lines have at least 3 tokens:
		// Particle_type | x | y | z
		// A line without exactly six tokens should be discarded.
		if (numTokens != 4) {
		    printf("Error: Invalid attached particle file line: %s\n", tokenList[0].val());
		    fclose(inp);
		    exit(-1);
		}

		// Make sure the index of this particle is unique.
		// NOTE: The particle list is sorted by index.
		int type_idx = conf->find_particle_type( tokenList[0].val() );
		if (type_idx < 0) {
		    printf("Error: Unrecognized particle type: %s\n", line);
		    fclose(inp);
		    exit(-1);
		}
		attached_particle_types.push_back( type_idx );
		attached_particle_positions.push_back( Vector3(atof(tokenList[1].val()), atof(tokenList[2].val()), atof(tokenList[3].val())) );
	}
	fclose(inp);
    }
    size_t sz = sizeof(Vector3)*attached_particle_positions.size();
    gpuErrchk(cudaMalloc( &(attached_particle_positions_d), sz ));
    gpuErrchk(cudaMemcpyAsync( attached_particle_positions_d, &attached_particle_positions[0], sz, cudaMemcpyHostToDevice));
}

void RigidBodyType::addGrid(String s, std::vector<String> &keys, std::vector<String> &files) {
	// tokenize and return
	int numTokens = s.tokenCount();
	if (numTokens != 2) {
		printf("ERROR: could not add Grid.\n"); // TODO improve this message
		exit(1);
	}
	String* token = new String[numTokens];
	s.tokenize(token);
	keys.push_back( String(token[0]) );
	files.push_back( String(token[1]) );
	delete[] token;
}
void RigidBodyType::addPotentialGrid(String s) {
    addGrid(s, potentialGridKeys, potentialGridFiles);
}
void RigidBodyType::addDensityGrid(String s) {
    addGrid(s, densityGridKeys, densityGridFiles);
}
void RigidBodyType::addPMF(String s) {
    addGrid(s, pmfKeys, pmfFiles);
}

void RigidBodyType::addScaleFactor(String s, std::vector<String> &keys, std::vector<float> &vals) {
	// tokenize and return
	int numTokens = s.tokenCount();
	if (numTokens != 2) {
		printf("ERROR: could not add Grid.\n"); // TODO improve this message
		exit(1);
	}
	String* token = new String[numTokens];
	s.tokenize(token);
	String key = token[0];
	float v = (float) strtod(token[1], NULL);
	keys.push_back( key );
	vals.push_back( v );
	delete[] token;
}
void RigidBodyType::scalePotentialGrid(String s) {
    addScaleFactor(s, potentialGridScaleKeys, potentialGridScale);
}
void RigidBodyType::scaleDensityGrid(String s) {
    addScaleFactor(s, densityGridScaleKeys, densityGridScale);
}
void RigidBodyType::scalePMF(String s) {
    addScaleFactor(s, pmfScaleKeys, pmfScale);
}

void RigidBodyType::initializeParticleLists() {
	if (numPotGrids < 1) return;

	numParticles = new int[numPotGrids];
	particles = new int*[numPotGrids];
	particles_d = new int*[numPotGrids];

	// Loop over potential grids
	for (int i = 0; i < numPotGrids; ++i) {
		String& gridName = potentialGridKeys[i];
		numParticles[i] = 0;

		// Count the particles interacting with potential grid i
		// Loop over particle types
		for (int j = 0; j < conf->numParts; ++j) {
			// Loop over rigid body grid names associated with particle type
			const std::vector<String>& gridNames = conf->partRigidBodyGrid[j];
			for (int k = 0; k < gridNames.size(); ++k) {
				if (gridNames[k] == gridName) {
					numParticles[i] += conf->numPartsOfType[j];
				}
			}
		}

		if (numParticles[i] > 0) {

		    // allocate array of particle ids for the potential grid 
		    particles[i] = new int[numParticles[i]];
		    int pid = 0;
		
		    // Loop over particle types to count the number of particles
		    for (int j = 0; j < conf->numParts; ++j) {

			// Build temporary id array of type j particles
			int tmp[conf->numPartsOfType[j]];
			int currId = 0;
			for (int aid = 0; aid < conf->num + conf->num_rb_attached_particles ; ++aid) {
			    if (conf->type[aid] == j)
				tmp[currId++] = aid;
			}
			if (currId == 0) continue;

			// Loop over rigid body grid names associated with particle type
			const std::vector<String>& gridNames = conf->partRigidBodyGrid[j];
			for (int k = 0; k < gridNames.size(); ++k) {
			    if (gridNames[k] == gridName) {
				// Copy type j particles to particles[i]
				memcpy( &(particles[i][pid]), tmp, sizeof(int)*currId );
				assert(currId == conf->numPartsOfType[j]);
				pid += conf->numPartsOfType[j];
			    }
			}
		    }

		    // Initialize device data
		    size_t sz = sizeof(int) * numParticles[i];
		    gpuErrchk(cudaMalloc( &(particles_d[i]), sz ));
		    gpuErrchk(cudaMemcpyAsync( particles_d[i], particles[i], sz, cudaMemcpyHostToDevice));
		}
	}
}
