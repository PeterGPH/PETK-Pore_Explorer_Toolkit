#include "Configuration.h"
#include "Angle.h"
#include "Dihedral.h"
#include "Restraint.h"
#include "ProductPotential.h"
#include <cmath>
#include <cassert>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <string>
#include <iostream>
using namespace std;

#ifndef gpuErrchk
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true) {
   if (code != cudaSuccess) {
      fprintf(stderr,"CUDA Error: %s %s %d\n", cudaGetErrorString(code), __FILE__, line);
      if (abort) exit(code);
   }
}
#endif

namespace
{
    template<class T> 
    void convertString(const String& token, void* data)
    {
        exit(1);
    }

    template<> 
    void convertString<float>(const String& token, void* data)
    {
        float* tmp = (float*)data;
        *tmp = atof(token);
    }

    template<>
    void convertString<String>(const String& token, void* data)
    {
        String* tmp = (String*)data;
        *tmp = token;
    }

    template<class T>
    void stringToArray(String* str, int& size, T** array)
    {
        register int num;
        String *token;
        num =  str->tokenCount();
        size = num;
        *array = new T[num];
        token  = new String[num];
        str->tokenize(token);

        for(int i = 0; i < num; ++i)
            convertString<T>(token[i], (*array)+i);
        delete [] token;
    }
}
Configuration::Configuration(const char* config_file, int simNum, bool debug) :
		simNum(simNum) {
	// Read the parameters.
	//type_d = NULL;
	kTGrid_d = NULL;
	//bonds_d = NULL;
	//bondMap_d = NULL;
	//excludes_d = NULL;
	//excludeMap_d = NULL;
	//angles_d = NULL;
	//dihedrals_d = NULL;
	setDefaults();
	readParameters(config_file);

	// Get the number of particles
	// printf("\nCounting particles specified in the ");
	if (restartCoordinates.length() > 0) {
    // Read them from the restart file.
	    // printf("restart file.\n");
		num = countRestart(restartCoordinates.val());
		if (copyReplicaCoordinates <= 0) {
		    num /= simNum;
		}
  } else {
    if (readPartsFromFile) readAtoms();
    if (numPartsFromFile > 0) {
      // Determine number of particles from input file (PDB-style)
	// printf("input file.\n");
      num = numPartsFromFile;
    } else {
      // Sum up all particles in config file
	// printf("configuration file.\n");
      //int num0 = 0;
      num = 0;
      for (int i = 0; i < numParts; i++) num += part[i].num;
      //num = num0;
    }
  } // end result: variable "num" is set

	// Count particles associated with rigid bodies
	num_rb_attached_particles = 0;
	if (numRigidTypes > 0) {
	    // grow list of rbs
	    for (int i = 0; i < numRigidTypes; i++) {
		RigidBodyType &rbt = rigidBody[i];
		rbt.attach_particles();
		num_rb_attached_particles += rbt.num * rbt.num_attached_particles();
	    }
	}
	assert( num_rb_attached_particles == 0 || simNum == 1 ); // replicas not yet implemented
	// num = num+num_rb_attached_particles;


	// Set the number capacity
	printf("\n%d particles\n", num);
	printf("%d particles attached to RBs\n", num_rb_attached_particles);

	if (numCap <= 0) numCap = numCapFactor*num; // max number of particles
	if (numCap <= 0) numCap = 20;

	if (readGroupSitesFromFile) readGroups();
	printf("%d groups\n", numGroupSites);

	// Allocate particle variables.
	// Each replica works with num+num_rb_attached_particles in array
	pos = new Vector3[ (num+num_rb_attached_particles) * simNum];

        //Han-Yi Chou
        if (ParticleDynamicType == String("Langevin") || ParticleDynamicType == String("NoseHooverLangevin"))
            momentum = new Vector3[(num+num_rb_attached_particles) * simNum];

	type   = new int[(num+num_rb_attached_particles) * simNum];
	serial = new int[(num+num_rb_attached_particles) * simNum];
	posLast = new Vector3[(num+num_rb_attached_particles) * simNum];

	{
	    int pidx = 0;
	    for (int i = 0; i < numRigidTypes; i++) { // Loop over RB types
		RigidBodyType &rbt = rigidBody[i];
		for (int j = 0; j < rbt.num; ++j) { // Loop over RBs
		    for (const int& t: rbt.get_attached_particle_types()) {
			type[num+pidx] = t;
			serial[num+pidx] = num+pidx;
			pidx++;
		    }
		}
	    }
	}	
	
        //Han-Yi Chou
        if(ParticleDynamicType == String("Langevin") || ParticleDynamicType == String("NoseHooverLangevin"))
           momLast = new Vector3[(num+num_rb_attached_particles) * simNum];
	name = new String[(num+num_rb_attached_particles) * simNum];
	currSerial = 0;


  // Now, load the coordinates
	loadedCoordinates = false;
        loadedMomentum    = false; //Han-Yi Chou

  //I need kT here Han-Yi Chou
  kT = temperature * 0.0019872065f; // `units "k K" "kcal_mol"`
  //kT = temperature * 0.593f;
 // If we have a restart file - use it
	if (restartCoordinates.length() > 0) {
		loadRestart(restartCoordinates.val()); 
		printf("Loaded %d restart coordinates from `%s'.\n", num, restartCoordinates.val());
		printf("Particle numbers specified in the configuration file will be ignored.\n");
		loadedCoordinates = true;
                //Han-Yi Chou Langevin dynamic
                if (ParticleDynamicType == String("Langevin") || ParticleDynamicType == String("NoseHooverLangevin"))
                {
                    if (restartMomentum.length() > 0)
                    {
                        loadRestartMomentum(restartMomentum.val());
                        printf("Loaded %d restart momentum from `%s'.\n", num, restartMomentum.val());
                        printf("Particle numbers specified in the configuration file will be ignored.\n");
                        loadedMomentum = true;
                    }
                    else
                    {
                        printf("Warning: There is no restart momentum file when using restart coordinates in Langevin Dynamics\n");
                        printf("Initialize with Boltzmann distribution\n");
                        loadedMomentum = Boltzmann(COM_Velocity, num * simNum);
                    }
               }
	} 
        else 
        {
		// Load coordinates from a file?
		if (numPartsFromFile > 0) {
			loadedCoordinates = true;
			for (int i = 0; i < num; i++) {
				int numTokens = partsFromFile[i].tokenCount();

				// Break the line down into pieces (tokens) so we can process them individually
				String* tokenList = new String[numTokens];
				partsFromFile[i].tokenize(tokenList);

				int currType = find_particle_type(tokenList[2]);
				if (currType == -1) {
				    printf("Error: Unable to find particle type %s\n", tokenList[2].val());
				    exit(1);

				}
				for (int j = 0; j < numParts; j++)
					if (tokenList[2] == part[j].name)
						currType = j;

				for (int s = 0; s < simNum; ++s)
				    type[i + s*num] = currType;

				serial[i] = currSerial++;

				pos[i] = Vector3(atof(tokenList[3].val()), atof(tokenList[4].val()), atof(tokenList[5].val()));
                                //Han-Yi Chou
                                if (ParticleDynamicType == String("Langevin") || ParticleDynamicType == String("NoseHooverLangevin"))
                                {
                                    loadedMomentum = true;
                                    if(numTokens == 9)
                                        momentum[i] = Vector3(atof(tokenList[6].val()), atof(tokenList[7].val()), atof(tokenList[8].val()));
                                    else
                                    {
                                        printf("Error occurs in %s at line %d. Please specify momentum\n", __FILE__, __LINE__);
                                        assert(1==2);
                                    }
                                }
			}
			delete[] partsFromFile;
			partsFromFile = NULL;
                        //Han-Yi Chou
                        for(int i = 1; i < simNum; ++i)
                            for(int j = 0; j < num; ++j)
                                serial[j + num * i] = currSerial++;
                }
                else 
                {
	            // Not loading coordinates from a file
	            populate();
	            if (inputCoordinates.length() > 0) 
                    {
		        printf("Loading coordinates from %s ... ", inputCoordinates.val());
			loadedCoordinates = loadCoordinates(inputCoordinates.val());
			if (loadedCoordinates)
			    printf("done!\n");
	            }
                    if(ParticleDynamicType == String("Langevin") || ParticleDynamicType == String("NoseHooverLangevin"))
                    {
                        if (inputMomentum.length() > 0) 
                        {
                            printf("Loading momentum from %s ... ", inputMomentum.val());
                            loadedMomentum = loadMomentum(inputMomentum.val());
                            if (loadedMomentum)
                                printf("done!\n");
                        }
                        else
                            loadedMomentum = Boltzmann(COM_Velocity, (num * simNum));
                    }
                }
            }
        //Check initialize momentum
        //if(ParticleDynamicType == String("Langevin"))
            //PrintMomentum();
	/* Initialize exclusions */
	excludeCapacity = 256;
	numExcludes = 0;
	excludes = new Exclude[excludeCapacity];

	if (readExcludesFromFile) readExcludes();
	if (readBondsFromFile) readBonds();
	if (readAnglesFromFile) readAngles();
	if (readDihedralsFromFile) readDihedrals();
	if (readRestraintsFromFile) readRestraints();
	if (readBondAnglesFromFile) readBondAngles();
	if (readProductPotentialsFromFile) readProductPotentials();


	if (temperatureGridFile.length() != 0) {
		printf("\nFound temperature grid file: %s\n", temperatureGridFile.val());
		tGrid = new BaseGrid(temperatureGridFile.val());
		printf("Loaded `%s'.\n", temperatureGridFile.val());
		printf("Grid size %s.\n", tGrid->getExtent().toString().val());

		// TODO: ask Max Belkin what this is about and how to remove hard-coded temps
		float ToSo = 1.0f / (295.0f * 4.634248239f); // 1 / (To * sigma(To))
		sigmaT = new BaseGrid(*tGrid);
		sigmaT->shift(-122.8305f);
		sigmaT->scale(0.0269167f);
		sigmaT->mult(*tGrid);
		sigmaT->scale(ToSo);

		kTGrid = new BaseGrid(*tGrid);
		float factor = 0.0019872065f; // `units "k K" "kcal_mol"`
		kTGrid->scale(factor);
		// char outFile[256];
		// char comment[256]; sprintf(comment,"KTGrid");
		// sprintf(outFile,"kTGrid.dx");
		// kTGrid->write(outFile, comment);
	}

	printf("\nFound %d particle types.\n", numParts);

	printf("Loading the potential grids...\n");
	// First load a single copy of each grid
	for (int i = 0; i < numParts; i++) 
        {
	    for(int j = 0; j < part[i].numPartGridFiles; ++j)
	    {
		std::string fname(partGridFile[i][j].val(), partGridFile[i][j].length());
		if (part_grid_dictionary.count( fname ) == 0)
		{
		    int len = fname.length();
		    if (len >= 3 && fname[len-3]=='.' && fname[len-2]=='d' && fname[len-1]=='x')
		    {
			part_grid_dictionary.insert({fname, BaseGrid(fname.c_str())});
		    }
		    else if  (len >= 4 && fname[len-4]=='.' && fname[len-3]=='d' && fname[len-2]=='e' && fname[len-1]=='f')
		    {
			assert(1==2); // Throw exception because this implementation needs to be revisited
/*                OverlordGrid* over = new OverlordGrid[part[i].numPartGridFiles];
		  part[i].meanPmf = new float[part[i].numPartGridFiles];
		  for(int j = 0; j < part[i].numPartGridFiles; ++j)
		  {
		  map = partGridFile[i][j];
		  len = map.length();
		  if (!(len >= 4 && map[len-4]=='.' && map[len-3]=='d' && map[len-2]=='e' && map[len-1]=='f'))
		  {
		  cout << "currently do not support different format " << endl;
		  exit(1);
		  }

		  String rootGrid = OverlordGrid::readDefFirst(map);
		  over[j] = OverlordGrid(rootGrid.val());
		  int count = over->readDef(map);
		  printf("Loaded system def file `%s'.\n", map.val());
		  printf("Found %d unique grids.\n", over->getUniqueGridNum());
		  printf("Linked %d subgrids.\n", count);
		  part[i].meanPmf[j] = part[i].pmf[j].mean();
		  }
		  part[i].pmf = static_cast<BaseGrid*>(over);
*/
		    } else {
			printf("WARNING: Unrecognized gridFile extension. Must be *.def or *.dx.\n");
			exit(-1);
		    }
		}
	    }
	}

	std::map<std::string,float> grid_mean_dict;
	for (const auto& pair : part_grid_dictionary)
	{
	    grid_mean_dict.insert({pair.first, pair.second.mean()});
	}

	// Then assign grid addresses to particles
	for (int i = 0; i < numParts; i++)
        {
	    part[i].pmf     = new BaseGrid*[part[i].numPartGridFiles];
	    part[i].pmf_scale = new float[part[i].numPartGridFiles];
	    part[i].meanPmf = new float[part[i].numPartGridFiles];
	    for(int j = 0; j < part[i].numPartGridFiles; ++j)
	    {
		part[i].pmf[j] = &(part_grid_dictionary.find( std::string(partGridFile[i][j]) )->second);
		part[i].pmf_scale[j] = partGridFileScale[i][j];
		part[i].meanPmf[j] = grid_mean_dict.find( std::string(partGridFile[i][j]) )->second * part[i].pmf_scale[j];
	    }
		if (partForceXGridFile[i].length() != 0) {
			part[i].forceXGrid = new BaseGrid(partForceXGridFile[i].val());
			printf("Loaded `%s'.\n", partForceXGridFile[i].val());
			printf("Grid size %s.\n", part[i].forceXGrid->getExtent().toString().val());
		}

		if (partForceYGridFile[i].length() != 0) {
			part[i].forceYGrid = new BaseGrid(partForceYGridFile[i].val());
			printf("Loaded `%s'.\n", partForceYGridFile[i].val());
			printf("Grid size %s.\n", part[i].forceYGrid->getExtent().toString().val());
		}

		if (partForceZGridFile[i].length() != 0) {
			part[i].forceZGrid = new BaseGrid(partForceZGridFile[i].val());
			printf("Loaded `%s'.\n", partForceZGridFile[i].val());
			printf("Grid size %s.\n", part[i].forceZGrid->getExtent().toString().val());
		}

		if (partDiffusionGridFile[i].length() != 0) {
			part[i].diffusionGrid = new BaseGrid(partDiffusionGridFile[i].val());
			printf("Loaded `%s'.\n", partDiffusionGridFile[i].val());
			printf("Grid size %s.\n", part[i].diffusionGrid->getExtent().toString().val());
		}

		if (temperatureGridFile.length() != 0) {
			if (partDiffusionGridFile[i].length() != 0) {
				part[i].diffusionGrid->mult(*sigmaT);
			} else {
				part[i].diffusionGrid = new BaseGrid(*sigmaT);
				part[i].diffusionGrid->scale(part[i].diffusion);
				// char outFile[256];
				// char comment[256]; sprintf(comment,"Diffusion for particle type %d", i);
				// sprintf(outFile,"diffusion%d.dx",i);
				// part[i].diffusionGrid->write(outFile, comment);
			}
		}
           
	}

    // Load reservoir files if any
    for (int i = 0; i < numParts; i++) {
        if (partReservoirFile[i].length() != 0) {
            printf("\nLoading the reservoirs for %s... \n", part[i].name.val());
            part[i].reservoir = new Reservoir(partReservoirFile[i].val());
            int nRes = part[i].reservoir->length();
            printf("\t -> %d reservoir(s) found in `%s'.\n", nRes, partReservoirFile[i].val());
        }
    }

    // Get the system dimensions
    // from the dimensions of supplied 3D potential maps
    if (size.length2() > 0) {	// use size if it's defined
	if (basis1.length2() > 0 || basis2.length2() > 0 || basis3.length2() > 0)
	    printf("WARNING: both 'size' and 'basis' were specified... using 'size'\n"); 
	basis1 = Vector3(size.x,0,0);
	basis2 = Vector3(0,size.y,0);
	basis3 = Vector3(0,0,size.z);
    }
    if (basis1.length2() > 0 && basis2.length2() > 0 && basis3.length2() > 0) {
	sys = new BaseGrid( Matrix3(basis1,basis2,basis3), origin, 1, 1, 1 );
    } else {
	// TODO: use largest system in x,y,z
	sys = *part[0].pmf;
    }
    sysDim = sys->getExtent();

// RBTODO: clean this mess up
	/* // RigidBodies... */
	/* if (numRigidTypes > 0) { */
	/* 	printf("\nCounting rigid bodies specified in the configuration file.\n"); */
	/* 	numRB = 0; */

	/* 	// grow list of rbs */
	/* 	for (int i = 0; i < numRigidTypes; i++) {			 */
	/* 		numRB += rigidBody[i].num; */

	/* 		std::vector<RigidBody> tmp; */
	/* 		for (int j = 0; j < rigidBody[i].num; j++) { */
	/* 			tmp.push_back( new RigidBody( this, rigidBody[i] ) ); */
	/* 		} */

	/* 		rbs.push_back(tmp); */
	/* 	} */
		// // state data
		// rbPos = new Vector3[numRB * simNum];
		// type = new int[numRB * simNum];

	/* } */
	/* printf("Initial RigidBodies: %d\n", numRB); */


	// Create exclusions from the exclude rule, if it was specified in the config file
	if (excludeRule != String("")) {
		int oldNumExcludes = numExcludes;
		Exclude* newExcludes = makeExcludes(bonds, bondMap, num, numBonds, excludeRule, numExcludes);
		if (excludes == NULL) {
			excludes = new Exclude[numExcludes];
		} else if (numExcludes >= excludeCapacity) {
			Exclude* tempExcludes = excludes;
			excludes = new Exclude[numExcludes];
			for (int i = 0; i < oldNumExcludes; i++)
				excludes[i] = tempExcludes[i];
			delete [] tempExcludes;
		}
		for (int i = oldNumExcludes; i < numExcludes; i++)
			excludes[i] = newExcludes[i - oldNumExcludes];
		printf("Built %d exclusions.\n",numExcludes);
	}

	{ // Add exclusions for RB attached particles
	    std::vector<Exclude> ex;
	    int start = num;
	    for (int i = 0; i < numRigidTypes; i++) { // Loop over RB types
		RigidBodyType &rbt = rigidBody[i];
		const int nap = rbt.num_attached_particles();
		for (int j = 0; j < rbt.num; ++j) { // Loop over RBs
		    for (int ai = 0; ai < nap-1; ++ai) {
			for (int aj = ai+1; aj < nap; ++aj) {
			    ex.push_back( Exclude( ai+start, aj+start ) );
			}
		    }
		    start += nap;
		}
	    }
	    // copy
	    int oldNumExcludes = numExcludes;
	    numExcludes = numExcludes + ex.size();
	    if (excludes == NULL) {
		excludes = new Exclude[numExcludes];
	    } else if (numExcludes >= excludeCapacity) {
		Exclude* tempExcludes = excludes;
		excludes = new Exclude[numExcludes];
		for (int i = 0; i < oldNumExcludes; i++)
		    excludes[i] = tempExcludes[i];
		delete [] tempExcludes;
	    }
	    for (int i = oldNumExcludes; i < numExcludes; i++)
		excludes[i] = ex[i - oldNumExcludes];
	}

	printf("Built %d exclusions.\n",numExcludes);		
	buildExcludeMap();

	// Count number of particles of each type
	numPartsOfType = new int[numParts];
	for (int i = 0; i < numParts; ++i) {
		numPartsOfType[i] = 0;
	}
	for (int i = 0; i < num+num_rb_attached_particles; ++i) {
		++numPartsOfType[type[i]];
	}

	// Some geometric stuff that should be gotten rid of.
	Vector3 buffer = (sys->getCenter() + 2.0f*sys->getOrigin())/3.0f;
	initialZ = buffer.z;

	// Set the initial conditions.
	// Do the initial conditions come from restart coordinates?
	// inputCoordinates are ignored if restartCoordinates exist.
	/*
	if (restartCoordinates.length() > 0) {
		loadRestart(restartCoordinates.val());
		printf("Loaded %d restart coordinates from `%s'.\n", num, restartCoordinates.val());
		printf("Particle numbers specified in the configuration file will be ignored.\n");
	} else {
		// Set the particle types.

		// Load coordinates from a file?
		if (numPartsFromFile > 0) {
			for (int i = 0; i < num; i++) {
				int numTokens = partsFromFile[i].tokenCount();

				// Break the line down into pieces (tokens) so we can process them individually
				String* tokenList = new String[numTokens];
				partsFromFile[i].tokenize(tokenList);
				int currType = 0;
				for (int j = 0; j < numParts; j++)
					if (tokenList[2] == part[j].name)
						currType = j;
				type[i] = currType;
				serial[i] = currSerial;
				currSerial++;

				pos[i] = Vector3(atof(tokenList[3].val()), atof(tokenList[4].val()), atof(tokenList[5].val()));
			}
			if (partsFromFile != NULL) {
				delete[] partsFromFile;
				partsFromFile = NULL;
			}
		} else if (inputCoordinates.length() > 0) {
			populate();
			printf("Loading coordinates from %s.\n", inputCoordinates.val());
			bool loaded = loadCoordinates(inputCoordinates.val());
			if (loaded) 
				printf("Loaded initial coordinates from %s.\n", inputCoordinates.val());
		}
	}
	*/
	

	// Get the maximum particle radius.
	minimumSep = 0.0f;
	for (int i = 0; i < numParts; ++i)
		minimumSep = std::max(minimumSep, part[i].radius);
	minimumSep *= 2.5f; // Make it a little bigger.

	// Default outputEnergyPeriod
	if (outputEnergyPeriod < 0)
		outputEnergyPeriod = 10 * outputPeriod;
	
	// If we are running with debug ON, ask the user which force computation to use
	if (debug)
		getDebugForce();

	printf("\n");
	switchStart = cutoff - switchLen;
	if (fullLongRange == 0)
		printf("Cutting off the potential from %.10g to %.10g.\n", switchStart, switchStart+switchLen);
	
	if (fullLongRange != 0)
		printf("No cell decomposition created.\n");

}

Configuration::~Configuration() {
	// System state
	delete[] pos;
        //Han-Yi Chou
        if (ParticleDynamicType == String("Langevin") || ParticleDynamicType == String("NoseHooverLangevin"))
            delete[] momentum;

	delete[] posLast;
        //Han-Yi Chou
        if (ParticleDynamicType == String("Langevin") || ParticleDynamicType == String("NoseHooverLangevin"))
            delete[] momLast;

	delete[] type;
	delete[] name;
	
	// Particle parameters
	delete[] part;
	//delete[] partGridFile;
	//delete[] partGridFileScale;
	for(int i = 0; i < numParts; ++i)
        {
            if(partGridFile[i] != NULL) 
            {
                delete[] partGridFile[i];
                partGridFile[i] = NULL;
            }
            if(partGridFileScale[i] != NULL)
            {
                delete[] partGridFileScale[i];
                partGridFileScale[i] = NULL;
            }
        }
        delete [] partGridFile;
        delete [] partGridFileScale;
        //delete numPartGridFiles;
	delete[] partForceXGridFile;
	delete[] partForceYGridFile;
	delete[] partForceZGridFile;
	delete[] partDiffusionGridFile;
	delete[] partReservoirFile;
	partRigidBodyGrid.clear();
	
	// TODO: plug memory leaks
	if (partsFromFile != NULL) delete[] partsFromFile;
	if (bonds != NULL) delete[] bonds;
	if (bondMap != NULL) delete[] bondMap;
	if (excludes != NULL) delete[] excludes;
	if (excludeMap != NULL) delete[] excludeMap;
	if (angles != NULL) delete[] angles;
	if (dihedrals != NULL) delete[] dihedrals;
	if (bondAngles != NULL) delete[] bondAngles;
	if (productPotentials != NULL) delete[] productPotentials;

	delete[] numPartsOfType;
	  
	// Table parameters
	delete[] partTableFile;
	delete[] partTableIndex0;
	delete[] partTableIndex1;

	delete[] bondTableFile;

	delete[] angleTableFile;

	delete[] dihedralTableFile;

	//if (type_d != NULL) {
		//gpuErrchk(cudaFree(type_d));
		gpuErrchk(cudaFree(sys_d));
		gpuErrchk(cudaFree(kTGrid_d));
		gpuErrchk(cudaFree(part_d));
		//gpuErrchk(cudaFree(bonds_d));
		//gpuErrchk(cudaFree(bondMap_d));
		//gpuErrchk(cudaFree(excludes_d));
		//gpuErrchk(cudaFree(excludeMap_d));
		//gpuErrchk(cudaFree(angles_d));
		//gpuErrchk(cudaFree(dihedrals_d));
	//}
}

void Configuration::copyToCUDA() {
    printf("Copying particle grids to GPU %d\n", GPUManager::current());
    for (const auto& pair : part_grid_dictionary)
    {
	// Copy PMF
	const BaseGrid& g = pair.second;
	BaseGrid *g_d = g.copy_to_cuda();
	part_grid_dictionary_d.insert({pair.first, g_d});
	// printf("Assigning grid for %s to %p (originally %p)\n", pair.first.c_str(), (void *) part_grid_dictionary_d[pair.first], (void *) g_d);
    }

    printf("Copying particle data to GPU %d\n", GPUManager::current());

	BrownianParticleType **part_addr = new BrownianParticleType*[numParts];

	// Copy the BaseGrid objects and their member variables/objects
	gpuErrchk(cudaMalloc(&part_d, sizeof(BrownianParticleType*) * numParts));
	// TODO: The above line fails when there is not enough memory. If it fails, stop.
	
	for (int i = 0; i < numParts; i++) 
        {
		BrownianParticleType *b = new BrownianParticleType(part[i]);
		// Copy PMF pointers
		if (part[i].pmf != NULL) 
                {
		    {
			BaseGrid** tmp_d = new BaseGrid*[part[i].numPartGridFiles];
			BaseGrid** tmp   = new BaseGrid*[part[i].numPartGridFiles];
			for(int j = 0; j < part[i].numPartGridFiles; ++j) {
			    // printf("Retrieving grid for %s (at %p)\n", partGridFile[i][j].val(), (void *) part_grid_dictionary_d[std::string(partGridFile[i][j])]);
			    tmp[j] = part_grid_dictionary_d[std::string(partGridFile[i][j])];
			}
			gpuErrchk(cudaMalloc(&tmp_d, sizeof(BaseGrid*)*part[i].numPartGridFiles));
			gpuErrchk(cudaMemcpy(tmp_d, tmp, sizeof(BaseGrid*)*part[i].numPartGridFiles,
					     cudaMemcpyHostToDevice));
			b->pmf = tmp_d;
		    }

		    {
			float *tmp;
			gpuErrchk(cudaMalloc(&tmp, sizeof(float)*part[i].numPartGridFiles));
			gpuErrchk(cudaMemcpy(tmp, part[i].pmf_scale, sizeof(float)*part[i].numPartGridFiles,
					     cudaMemcpyHostToDevice));
			b->pmf_scale = tmp;
		    }

		    {
			float *tmp;
			gpuErrchk(cudaMalloc(&tmp, sizeof(float)*part[i].numPartGridFiles));
			gpuErrchk(cudaMemcpy(tmp, part[i].meanPmf, sizeof(float)*part[i].numPartGridFiles, 
					     cudaMemcpyHostToDevice));
			b->meanPmf = tmp;
		    }

		    {
			BoundaryCondition *tmp;
			size_t s = sizeof(BoundaryCondition)*part[i].numPartGridFiles;
			gpuErrchk(cudaMalloc(&tmp, s));
			gpuErrchk(cudaMemcpy(tmp, part[i].pmf_boundary_conditions, s, cudaMemcpyHostToDevice));
			b->pmf_boundary_conditions = tmp;
		    }
                    
		}

		// Copy the diffusion grid
		if (part[i].diffusionGrid != NULL) {
		    b->diffusionGrid = part[i].diffusionGrid->copy_to_cuda();
		} else {
		    b->diffusionGrid = NULL;
		}
		
		//b->pmf = pmf;
		gpuErrchk(cudaMalloc(&part_addr[i], sizeof(BrownianParticleType)));
		gpuErrchk(cudaMemcpyAsync(part_addr[i], b, sizeof(BrownianParticleType),
				cudaMemcpyHostToDevice));
	}
	// RBTODO: moved this out of preceding loop; was that correct?
	gpuErrchk(cudaMemcpyAsync(part_d, part_addr, sizeof(BrownianParticleType*) * numParts,
				cudaMemcpyHostToDevice));

	// kTGrid_d
	kTGrid_d = NULL;
	if (temperatureGridFile.length() > 0) {
		gpuErrchk(cudaMalloc(&kTGrid_d, sizeof(BaseGrid)));
		gpuErrchk(cudaMemcpyAsync(kTGrid_d, kTGrid, sizeof(BaseGrid), cudaMemcpyHostToDevice));
	}

	// type_d and sys_d
	gpuErrchk(cudaMalloc(&sys_d, sizeof(BaseGrid)));
	gpuErrchk(cudaMemcpyAsync(sys_d, sys, sizeof(BaseGrid), cudaMemcpyHostToDevice));
	/*gpuErrchk(cudaMalloc(&type_d, sizeof(int) * num * simNum));
	gpuErrchk(cudaMemcpyAsync(type_d, type, sizeof(int+num_rb_attached_particles) * num * simNum, cudaMemcpyHostToDevice));
	
	if (numBonds > 0) {
		// bonds_d
		gpuErrchk(cudaMalloc(&bonds_d, sizeof(Bond) * numBonds));
		gpuErrchk(cudaMemcpyAsync(bonds_d, bonds, sizeof(Bond) * numBonds, cudaMemcpyHostToDevice));
		
		// bondMap_d
		gpuErrchk(cudaMalloc(&bondMap_d, sizeof(int2) * num));
		gpuErrchk(cudaMemcpyAsync(bondMap_d, bondMap, sizeof(int2) * num, cudaMemcpyHostToDevice));
	}

	if (numExcludes > 0) {
		// excludes_d
		gpuErrchk(cudaMalloc(&excludes_d, sizeof(Exclude) * numExcludes));
		gpuErrchk(cudaMemcpyAsync(excludes_d, excludes, sizeof(Exclude) * numExcludes,
				cudaMemcpyHostToDevice));
		
		// excludeMap_d
		gpuErrchk(cudaMalloc(&excludeMap_d, sizeof(int2) * (num));
		gpuErrchk(cudaMemcpyAsync(excludeMap_d, excludeMap, sizeof(int2) * num,
				cudaMemcpyHostToDevice));
	}

	if (numAngles > 0) {
		// angles_d
		gpuErrchk(cudaMalloc(&angles_d, sizeof(Angle) * numAngles));
		gpuErrchk(cudaMemcpyAsync(angles_d, angles, sizeof(Angle) * numAngles,
				cudaMemcpyHostToDevice));
	}

	if (numDihedrals > 0) {
		// dihedrals_d
		gpuErrchk(cudaMalloc(&dihedrals_d, sizeof(Dihedral) * numDihedrals));
		gpuErrchk(cudaMemcpyAsync(dihedrals_d, dihedrals,
												 		  sizeof(Dihedral) * numDihedrals,
														 	cudaMemcpyHostToDevice));
	}*/
	gpuErrchk(cudaDeviceSynchronize());
}

void Configuration::setDefaults() {
    // System parameters
	outputName = "out";
	timestep = 1e-5f;
	rigidBodyGridGridPeriod = 1;
	steps = 100;

	unsigned long int r0 = clock();
	for (int i = 0; i < 4; i++)
	    r0 *= r0 + 1;
	seed = time(NULL) + r0;

	origin = Vector3(0,0,0);
	size = Vector3(0,0,0);
	basis1 = Vector3(0,0,0);
	basis2 = Vector3(0,0,0);
	basis3 = Vector3(0,0,0);
	
	inputCoordinates = "";
	restartCoordinates = "";
        //Han-Yi Chou
        inputMomentum = "";
        restartMomentum = "";
	copyReplicaCoordinates = 1;
	numberFluct = 0;
	numberFluctPeriod = 200;
	interparticleForce = 1;
	tabulatedPotential = 0;
	fullLongRange = 0;
	//	kTGridFile = ""; // Commented out for an unknown reason
	temperature = 295.0f;
	temperatureGridFile = "";
	coulombConst = 566.440698f/92.0f;
	electricField = 0.0f;
	cutoff = 10.0f;
	switchLen = 2.0f;
	pairlistDistance = 2.0f;
	imdForceScale = 1.0f;
	outputPeriod = 200;
	outputEnergyPeriod = -1;
	outputFormat = TrajectoryWriter::formatDcd;
	currentSegmentZ = -1.0f;
	numCap = 0;
	decompPeriod = 10;
	readPartsFromFile = 0;
	numPartsFromFile = 0;
	partsFromFile = NULL;
	readBondsFromFile = false;
	numGroupSites = 0;
	readGroupSitesFromFile = false;
	

	numBonds = 0;
	bonds = NULL;
	bondMap = NULL;
	numTabBondFiles = 0;
	readExcludesFromFile = false;
	numExcludes = 0;
	excludeCapacity = 256;
	excludes = NULL;
	excludeMap = NULL;
	excludeRule = "";
	readAnglesFromFile = false;
	numAngles = 0;
	angles = NULL;
	numTabAngleFiles = 0;
	readDihedralsFromFile = false;
	numDihedrals = 0;
	dihedrals = NULL;
	numTabDihedralFiles = 0;

	readBondAnglesFromFile = false;
	numBondAngles = 0;
	bondAngles = NULL;

	readProductPotentialsFromFile = false;
	numProductPotentials = 0;
	productPotentials = NULL;
	simple_potential_ids = XpotMap();
	simple_potentials = std::vector<SimplePotential>();

	readRestraintsFromFile = false;
	numRestraints = 0;
	restraints = NULL;

        //Han-Yi Chou default values
        ParticleDynamicType  = String("Brown");
        RigidBodyDynamicType = String("Brown");
        COM_Velocity = Vector3(0.f,0.f,0.f);
        ParticleLangevinIntegrator = String("BAOAB"); //The default is BAOAB

	// Hidden parameters
	// Might be parameters later
	numCapFactor = 5;

        ParticleInterpolationType = 0;
        RigidBodyInterpolationType = 0;
}

int Configuration::readParameters(const char * config_file) {
	Reader config(config_file);
	printf("Read config file %s\n", config_file);

	// Get the number of particles.
	const int numParams = config.length();
	numParts = config.countParameter("particle");
	numRigidTypes = config.countParameter("rigidBody");

	// Allocate the particle variables.
	part = new BrownianParticleType[numParts];
	//partGridFile = new String[numParts];
	//partGridFileScale = new float[numParts];
	partGridFile       = new String*[numParts];
        //partGridFileScale = new float[numParts];
        partGridFileScale  = new float*[numParts];
        //int numPartGridFiles = new int[numParts];

	partForceXGridFile = new String[numParts];
	partForceYGridFile = new String[numParts];
	partForceZGridFile = new String[numParts];
	partDiffusionGridFile = new String[numParts];
	partReservoirFile = new String[numParts];
	partRigidBodyGrid.resize(numParts);
	
	// Allocate the table variables.
	partTableFile = new String[numParts*numParts];
	partTableIndex0 = new int[numParts*numParts];
	partTableIndex1 = new int[numParts*numParts];

	// Allocate rigid body types
	rigidBody = new RigidBodyType[numRigidTypes];
	
	// Set a default
	/*
	for (int i = 0; i < numParts; ++i) {
	    partGridFileScale[i] = 1.0f;
	}*/

        for(int i = 0; i < numParts; ++i)
        {
            partGridFile[i] = NULL;
            partGridFileScale[i] = NULL;
            //part[i].numPartGridFiles = -1;
        }
        //for(int i = 0; i < numParts; ++i)
          //  cout << part[i].numPartGridFiles << endl;

	int btfcap = 10;
	bondTableFile = new String[btfcap];

	int atfcap = 10;
	angleTableFile = new String[atfcap];

	int dtfcap = 10;
	dihedralTableFile = new String[dtfcap];

	int currPart = -1;
	int currTab = -1;
	int currBond = -1;
	int currAngle = -1;
	int currDihedral = -1;
	int currRB = -1;

	int partClassPart =  0;
	int partClassRB   =  1;
	int currPartClass = -1;				// 0 => particle, 1 => rigidBody



	for (int i = 0; i < numParams; i++) {
		String param = config.getParameter(i);
		String value = config.getValue(i);
		// printf("Parsing %s: %s\n", param.val(), value.val());
		if (param == String("outputName"))
			outputName = value;
		else if (param == String("timestep"))
			timestep = (float) strtod(value.val(), NULL);
		else if (param == String("rigidBodyGridGridPeriod"))
			rigidBodyGridGridPeriod = atoi(value.val());
		else if (param == String("steps"))
			steps = atol(value.val());
		else if (param == String("seed"))
			seed = atoi(value.val());
		else if (param == String("origin"))
		    origin = stringToVector3( value );
		else if (param == String("systemSize"))
		    size = stringToVector3( value );
		else if (param == String("basis1"))
		    basis1 = stringToVector3( value );
		else if (param == String("basis2"))
		    basis2 = stringToVector3( value );
		else if (param == String("basis3"))
		    basis3 = stringToVector3( value );
		else if (param == String("inputCoordinates"))
			inputCoordinates = value;
		else if (param == String("restartCoordinates"))
			restartCoordinates = value;
                //Han-Yi Chou
                else if (param == String("inputMomentum"))
                        inputMomentum = value;
                else if (param == String("restartMomentum"))
                        restartMomentum = value;
		else if (param == String("copyReplicaCoordinates"))
		        copyReplicaCoordinates = atoi(value.val());
		else if (param == String("temperature"))
			temperature =  (float) strtod(value.val(),NULL);
		else if (param == String("temperatureGrid"))
			temperatureGridFile = value;
		else if (param == String("numberFluct"))
			numberFluct = atoi(value.val());
		else if (param == String("numberFluctPeriod"))
			numberFluctPeriod = atoi(value.val());
		else if (param == String("interparticleForce"))
			interparticleForce = atoi(value.val());
		else if (param == String("fullLongRange") || param == String("fullElect") )
			fullLongRange = atoi(value.val());
		else if (param == String("coulombConst"))
			coulombConst = (float) strtod(value.val(), NULL);
		else if (param == String("electricField"))
			electricField = (float) strtod(value.val(), NULL);
		else if (param == String("cutoff"))
			cutoff = (float) strtod(value.val(), NULL);
		else if (param == String("switchLen"))
			switchLen = (float) strtod(value.val(), NULL);
		else if (param == String("pairlistDistance"))
			pairlistDistance = (float) strtod(value.val(), NULL);
		else if (param == String("scaleIMDForce"))
			imdForceScale = (float) strtod(value.val(), NULL);		
		else if (param == String("outputPeriod"))
			outputPeriod = atoi(value.val());
		else if (param == String("outputEnergyPeriod"))
			outputEnergyPeriod = atoi(value.val());
		else if (param == String("outputFormat"))
			outputFormat = TrajectoryWriter::getFormatCode(value);
		else if (param == String("currentSegmentZ"))
			currentSegmentZ = (float) strtod(value.val(), NULL);
		else if (param == String("numCap"))
			numCap = atoi(value.val());
		else if (param == String("decompPeriod"))
			decompPeriod = atoi(value.val());

                //Han-Yi Chou
                else if (param == String("ParticleDynamicType"))
                    ParticleDynamicType = value;
                else if (param == String("RigidBodyDynamicType"))
                    RigidBodyDynamicType = value;
                else if (param == String("ParticleLangevinIntegrator"))
                    ParticleLangevinIntegrator = value;
                else if (param == String("ParticleInterpolationType"))
                    ParticleInterpolationType = atoi(value.val());
                else if (param == String("RigidBodyInterpolationType"))
                    RigidBodyInterpolationType = atoi(value.val());
		// PARTICLES
		else if (param == String("particle")) {
		    part[++currPart] = BrownianParticleType(value);
		    currPartClass = partClassPart;
		}
                else if (param == String("mu")) { // for Nose-Hoover Langevin
		    if (currPart < 0) exit(1);
		    part[currPart].mu = (float) strtod(value.val(), NULL);
		} else if (param == String("forceXGridFile")) {
		    if (currPart < 0) exit(1);
		    partForceXGridFile[currPart] = value;
		} else if (param == String("forceYGridFile")) {
		    if (currPart < 0) exit(1);
		    partForceYGridFile[currPart] = value;
		} else if (param == String("forceZGridFile")) {
		    if (currPart < 0) exit(1);
		    partForceZGridFile[currPart] = value;
		} else if (param == String("diffusionGridFile")) {
		    if (currPart < 0) exit(1);
		    partDiffusionGridFile[currPart] = value;
		} else if (param == String("diffusion")) {
		    if (currPart < 0) exit(1);
		    part[currPart].diffusion = (float) strtod(value.val(), NULL);
		} else if (param == String("charge")) {
		    if (currPart < 0) exit(1);
		    part[currPart].charge = (float) strtod(value.val(), NULL);
		} else if (param == String("radius")) {
		    if (currPart < 0) exit(1);
		    part[currPart].radius = (float) strtod(value.val(), NULL);
		} else if (param == String("eps")) {
		    if (currPart < 0) exit(1);
		    part[currPart].eps = (float) strtod(value.val(), NULL);
		} else if (param == String("reservoirFile")) {
		    if (currPart < 0) exit(1);
		    partReservoirFile[currPart] = value;
		}
		else if (param == String("tabulatedPotential"))
			tabulatedPotential = atoi(value.val());
		else if (param == String("tabulatedFile"))
			readTableFile(value, ++currTab);
		else if (param == String("tabulatedBondFile")) {
			if (numTabBondFiles >= btfcap) {
				String* temp = bondTableFile;
				btfcap *= 2;	
				bondTableFile = new String[btfcap];
				for (int j = 0; j < numTabBondFiles; j++)
					bondTableFile[j] = temp[j];
				delete[] temp;
			}
			if (readBondFile(value, ++currBond))
				numTabBondFiles++;
		} else if (param == String("inputParticles")) {
			if (readPartsFromFile) {
				printf("WARNING: More than one particle file specified. Ignoring new file.\n");
			} else {
				partFile = value;
				readPartsFromFile = true;
				loadedCoordinates = true;
			}
		} else if (param == String("inputGroups")) {
			if (readGroupSitesFromFile) {
				printf("WARNING: More than one group file specified. Ignoring new file.\n");
			} else {
				groupSiteFile = value;
				readGroupSitesFromFile = true;
			}
		} else if (param == String("inputBonds")) {
			if (readBondsFromFile) {
				printf("WARNING: More than one bond file specified. Ignoring new bond file.\n");
			} else {
				bondFile = value;				
				readBondsFromFile = true;
			}
		} else if (param == String("inputExcludes")) {
			if (readExcludesFromFile) {
				printf("WARNING: More than one exclude file specified. Ignoring new exclude file.\n");
			} else {
			    printf("inputExclude %s\n", value.val());
				excludeFile = value;				
				readExcludesFromFile = true;
			}
		} else if (param == String("exclude") or param == String("exclusion")) {
			excludeRule = value; 
		} else if (param == String("inputAngles")) {
			if (readAnglesFromFile) {
				printf("WARNING: More than one angle file specified. Ignoring new angle file.\n");
			} else {
				angleFile = value;
				readAnglesFromFile = true;
			}
		} else if (param == String("inputBondAngles")) {
			if (readBondAnglesFromFile) {
				printf("WARNING: More than one bondangle file specified. Ignoring new bondangle file.\n");
			} else {
			        bondAngleFile = value;
				readBondAnglesFromFile = true;
			}
		} else if (param == String("inputProductPotentials")) {
			if (readBondAnglesFromFile) {
				printf("WARNING: More than one product potential file specified. Ignoring new file.\n");
			} else {
			        productPotentialFile = value;
				readProductPotentialsFromFile = true;
			}
		} else if (param == String("tabulatedAngleFile")) {
			if (numTabAngleFiles >= atfcap) {
				String* temp = angleTableFile;
				atfcap *= 2;	
				angleTableFile = new String[atfcap];
				for (int j = 0; j < numTabAngleFiles; j++)
					angleTableFile[j] = temp[j];
				delete[] temp;
			}
			if (readAngleFile(value, ++currAngle))
				numTabAngleFiles++;
		} else if (param == String("inputDihedrals")) {
			if (readDihedralsFromFile) {
				printf("WARNING: More than one dihedral file specified. Ignoring new dihedral file.\n");
			} else {
				dihedralFile = value;
				readDihedralsFromFile = true;
			}
		} else if (param == String("tabulatedDihedralFile")) {
			if (numTabDihedralFiles >= dtfcap) {
				String * temp = dihedralTableFile;
				dtfcap *= 2;
				dihedralTableFile = new String[dtfcap];
				for (int j = 0; j < numTabDihedralFiles; j++)
					dihedralTableFile[j] = temp[j];
				delete[] temp;
			}
			if (readDihedralFile(value, ++currDihedral))
				numTabDihedralFiles++;
		} else if (param == String("inputRestraints")) {
			if (readRestraintsFromFile) {
				printf("WARNING: More than one restraint file specified. Ignoring new restraint file.\n");
			} else {
				restraintFile = value;
				readRestraintsFromFile = true;
			}
		} else if (param == String("gridFileScale")) {
		    if (currPart < 0) exit(1);
			//partGridFileScale[currPart] = (float) strtod(value.val(), NULL);
			  stringToArray<float>(&value, part[currPart].numPartGridFiles, 
                                                      &partGridFileScale[currPart]);
		} else if (param == String("gridFileBoundaryConditions")) {
		    if (currPart < 0) exit(1);
		    register size_t num = value.tokenCount();
		    if (num > 0) {
			String *tokens  = new String[num];
			BoundaryCondition *data = new BoundaryCondition[num];
			value.tokenize(tokens);
			for(size_t i = 0; i < num; ++i) {
			    tokens[i].lower();
			    if (tokens[i] == "dirichlet")
				data[i] = dirichlet;
			    else if (tokens[i] == "neumann")
				data[i] = neumann;
			    else if (tokens[i] == "periodic")
				data[i] = periodic;
			    else {
				fprintf(stderr,"WARNING: Unrecognized gridFile boundary condition \"%s\". Using Dirichlet.\n", tokens[i].val() );
				data[i] = dirichlet;
			    }
			}
			delete[] tokens;
			part[currPart].set_boundary_conditions(num, data);
			delete[] data;
		    }
		} else if (param == String("rigidBodyPotential")) {
		    if (currPart < 0) exit(1);
		    partRigidBodyGrid[currPart].push_back(value);
		}
                //Han-Yi Chou initial COM velocity for total particles
                else if (param == String("COM_Velocity"))
                    COM_Velocity = stringToVector3(value);

		// RIGID BODY
		else if (param == String("rigidBody")) {
			// part[++currPart] = BrownianParticleType(value);
			rigidBody[++currRB] = RigidBodyType(value, this);
			currPartClass = partClassRB;
		}
		else if (param == String("inertia")) {
		    if (currRB < 0) exit(1);
			rigidBody[currRB].inertia = stringToVector3( value );
		} else if (param == String("rotDamping")) {
		    if (currRB < 0) exit(1);
			rigidBody[currRB].rotDamping = stringToVector3( value );
		} else if (param == String("attachedParticles")) {
			rigidBody[currRB].append_attached_particle_file(value);
		} else if (param == String("densityGrid")) {
		    if (currRB < 0) exit(1);
			rigidBody[currRB].addDensityGrid(value);
		} else if (param == String("potentialGrid")) {
		    if (currRB < 0) exit(1);
			rigidBody[currRB].addPotentialGrid(value);
		} else if (param == String("densityGridScale")) {
		    if (currRB < 0) exit(1);
			rigidBody[currRB].scaleDensityGrid(value);
		} else if (param == String("potentialGridScale")) {
		    if (currRB < 0) exit(1);
			rigidBody[currRB].scalePotentialGrid(value);
		} else if (param == String("pmfScale")) {
		    if (currRB < 0) exit(1);
			rigidBody[currRB].scalePMF(value);
		} else if (param == String("position")) {
		    if (currRB < 0) exit(1);
			rigidBody[currRB].initPos = stringToVector3( value );
		} else if (param == String("orientation")) {
		    if (currRB < 0) exit(1);
			rigidBody[currRB].initRot = stringToMatrix3( value );
                } else if (param == String("momentum")) {
                        rigidBody[currRB].initMomentum = stringToVector3(value);
                } else if (param == String("angularMomentum")) {
		    if (currRB < 0) exit(1);
                        rigidBody[currRB].initAngularMomentum = stringToVector3(value);
		}
		else if (param == String("inputRBCoordinates"))
			inputRBCoordinates = value;
		else if (param == String("restartRBCoordinates"))
		        restartRBCoordinates = value;
		
		// COMMON
		else if (param == String("num")) {
		    if (currPartClass == partClassPart) {
			if (currPart < 0) exit(1);
			part[currPart].num = atoi(value.val());
		    } else if (currPartClass == partClassRB) {
			if (currRB < 0) exit(1);
			rigidBody[currRB].num = atoi(value.val());
		    }
		}
                //set mass here Han-Yi Chou
                else if (param == String("mass"))
                {
                    if (currPartClass == partClassPart) {
			if (currPart < 0) exit(1);
                        part[currPart].mass    = (float) strtod(value.val(),NULL);
                    } else if (currPartClass == partClassRB) {
			if (currRB < 0) exit(1);
                        rigidBody[currRB].mass = (float) strtod(value.val(),NULL);
		    }
                }
                //set damping here, using anisotropic damping, i.e. data type Vector3 Han-Yi Chou
                else if (param == String("transDamping"))
                {
                    if (currPartClass == partClassPart) {
			if (currPart < 0) exit(1);
                        part[currPart].transDamping    = stringToVector3(value);
		    } else if (currPartClass == partClassRB) {
			if (currRB < 0) exit(1);
                        rigidBody[currRB].transDamping = stringToVector3(value);
		    }
                }
		else if (param == String("gridFile")) {
			if (currPartClass == partClassPart)
                        {
			    if (currPart < 0) exit(1);
                                printf("Applying grid file '%s'\n", value.val());
				stringToArray<String>(&value, part[currPart].numPartGridFiles, 
                                                             &partGridFile[currPart]);
				const int& num = part[currPart].numPartGridFiles;
				partGridFileScale[currPart] = new float[num];
                                for(int i = 0; i < num; ++i) {
                                    // printf("%s ", partGridFile[currPart]->val());
				    partGridFileScale[currPart][i] = 1.0f;
				}

				// Set default boundary conditions for grids
				BoundaryCondition *bc = part[currPart].pmf_boundary_conditions;
				if (bc == NULL) {
				    bc = new BoundaryCondition[num];
				    for(int i = 0; i < num; ++i) {
					bc[i] = dirichlet;
				    }
				    part[currPart].pmf_boundary_conditions = bc;
				}
                        }
			else if (currPartClass == partClassRB) {
			    if (currRB < 0) exit(1);
				rigidBody[currRB].addPMF(value);
			}
		}
		// UNKNOWN
		else {
			printf("ERROR: Unrecognized keyword `%s'.\n", param.val());
			exit(1);
		}
	}

	// extra configuration for RB types
	for (int i = 0; i < numRigidTypes; i++)
		rigidBody[i].setDampingCoeffs(timestep);

        //For debugging purpose Han-Yi Chou
        //Print();
	return numParams;
}
//Han-Yi Chou
void Configuration::Print()
{
    printf("The dynamic type for particle is %s \n", ParticleDynamicType.val());
    for(int i = 0; i < numParts; ++i)
    {
        printf("The type %d has mass %f \n", i,part[i].mass);
        printf("The diffusion coefficient is %f \n", part[i].diffusion);
        printf("The translational damping is %f %f %f \n", part[i].transDamping.x, part[i].transDamping.y, part[i].transDamping.z);
    }
    printf("Done with check for Langevin");
    //assert(1==2);
}

void Configuration::PrintMomentum()
{
    for(int i = 0; i < num; ++i)
    {
        printf("%f %f %f\n", momentum[i].x, momentum[i].y, momentum[i].z);
    }
    //assert(1==2);
}
Vector3 Configuration::stringToVector3(String s) {
	// tokenize and return
	int numTokens = s.tokenCount();
	if (numTokens != 3) {
		printf("ERROR: could not convert input to Vector3.\n"); // TODO improve this message
		exit(1);
	}
	String* token = new String[numTokens];
	s.tokenize(token);
	Vector3 v( (float) strtod(token[0], NULL),
						 (float) strtod(token[1], NULL),
						 (float) strtod(token[2], NULL) );
	return v;
}
Matrix3 Configuration::stringToMatrix3(String s) {
	// tokenize and return
	int numTokens = s.tokenCount();
	if (numTokens != 9) {
		printf("ERROR: could not convert input to Matrix3.\n"); // TODO improve this message
		exit(1);
	}
	String* token = new String[numTokens];
	s.tokenize(token);
	Matrix3 m( (float) strtod(token[0], NULL),
						 (float) strtod(token[1], NULL),
						 (float) strtod(token[2], NULL),
						 (float) strtod(token[3], NULL),
						 (float) strtod(token[4], NULL),
						 (float) strtod(token[5], NULL),
						 (float) strtod(token[6], NULL),
						 (float) strtod(token[7], NULL),
						 (float) strtod(token[8], NULL) );
	return m;
}

void Configuration::readAtoms() {
	// Open the file
	FILE* inp = fopen(partFile.val(), "r");
	char line[256];

	// If the particle file cannot be found, exit the program
	if (inp == NULL) {
		printf("ERROR: Could not open `%s'.\n", partFile.val());
		bool found = true;
		for (int i = 0; i < numParts; i++)
			if (part[i].num == 0)
				found = false;
		// assert(false); // TODO probably relax constraint that particle must be found; could just be in RB
		if (!found) {
			printf("ERROR: Number of particles not specified in config file.\n");
			exit(1);
		}
		printf("Using default coordinates file\n");
		return;
	}

	// Our particle array has a starting capacity of 256
	// We will expand this later if we need to.
	int capacity = 256;
	numPartsFromFile = 0;
	partsFromFile = new String[capacity];
	indices = new int[capacity];
	indices[0] = 0;

	// Get and process all lines of input
	while (fgets(line, 256, inp) != NULL) {
		// Lines in the particle file that begin with # are comments
		if (line[0] == '#') continue;
		      
		String s(line);
		int numTokens = s.tokenCount();
		      
		// Break the line down into pieces (tokens) so we can process them individually
		String* tokenList = new String[numTokens];
		s.tokenize(tokenList);

		// Legitimate ATOM input lines have 6 tokens: 
		// ATOM | Index | Name | X-coord | Y-coord | Z-coord
		// A line without exactly six tokens should be discarded.
                if (ParticleDynamicType == String("Langevin") || ParticleDynamicType == String("NoseHooverLangevin")) {
		    if (numTokens != 9) {
			printf("Error: Invalid particle file line: %s\n", line);
			exit(-1);
		    }
		} else {
		    if (numTokens != 6) {
			printf("Error: Invalid particle file line: %s\n", line);
			exit(-1);
		    }
		}

		// Ensure that this particle's type was defined in the config file.
		// If not, discard this line.
		bool found;
		for (int j = 0; j < numParts; j++) {
			// If this particle type exists, add a new one to the list
			if (part[j].name == tokenList[2]) {
				found = true;
				part[j].num++;
			}
		}

		// If the particle's type does not exist according to the config file, discard it.
		if (!found) {
			printf("WARNING Unknown particle type %s found and discarded.\n", tokenList[2].val());
			continue;
		}

		// If we don't have enough room in our particle array, we need to expand it.
		if (numPartsFromFile >= capacity) {
			// Temporary pointers to the old arrays
			String* temp = partsFromFile;	
			int* temp2 = indices;

			// Double the capacity
			capacity *= 2;

			// Create pointers to new arrays which are twice the size of the old ones
			partsFromFile = new String[capacity];
			indices = new int[capacity];
		
			// Copy the old values into the new arrays
			for (int j = 0; j < numPartsFromFile; j++) {
				partsFromFile[j] = temp[j];
				indices[j] = temp2[j];
			}

			// delete the old arrays
			delete[] temp;
			delete[] temp2;
		}
		// Make sure the index of this particle is unique.
		// NOTE: The particle list is sorted by index. 
		bool uniqueID = true;		
		int key = atoi(tokenList[1].val());
		int mid = 0;

		// If the index is greater than the last index in the list, 
		// this particle belongs at the end of the list. Since the 
		// list is kept sorted, we know this is okay.
		if (numPartsFromFile == 0 || key > indices[numPartsFromFile - 1]) {
			indices[numPartsFromFile] = key;
			partsFromFile[numPartsFromFile++] = line;
		}
		// We need to do a binary search to figure out if
		// the index already exists in the list. 
		// The assumption is that input files SHOULD have their indices sorted in 
		// ascending order, so we shouldn't actually use the binary search 
		// or the sort (which is pretty time consuming) very often.
		else {
			int low = 0, high = numPartsFromFile - 1;
			
			while (low <= high) {
				mid = (int)((high - low) / 2 + low);
				int curr = indices[mid];
				if (curr < key) {
					low = mid + 1;
				} else if (curr > key) {
					high = mid - 1;
				} else {
					// For now, particles with non-unique IDs are simply not added to the array
					// Other possible approaches which are not yet implemented:
					// 1: Keep track of these particles and assign them new IDs after you have
					//    already added all of the other particles. 	
					// 2: Get rid of ALL particles with that ID, even the ones that have already 
					//    been added.
					printf("WARNING: Non-unique ID found: %s\n", line);
					uniqueID = false;
					break;
				}
			}
			if (uniqueID) {
				// Add the particle to the end of the array, then sort it. 
				indices[numPartsFromFile] = key;
				partsFromFile[numPartsFromFile++] = line;
				std::sort(indices, indices + numPartsFromFile);
				std::sort(partsFromFile, partsFromFile + numPartsFromFile, compare());		
			}
		}
	}
}
void Configuration::readGroups() {
	// Open the file
    const size_t line_char = 16384;
	FILE* inp = fopen(groupSiteFile.val(), "r");
	char line[line_char];

	// If the particle file cannot be found, exit the program
	if (inp == NULL) {
		printf("ERROR: Could not open `%s'.\n", partFile.val());
		exit(1);
	}

	// Our particle array has a starting capacity of 256
	// We will expand this later if we need to.
	// int capacity = 256;
	numGroupSites = 0;

	// partsFromFile = new String[capacity];
	// indices = new int[capacity];
	// indices[0] = 0;

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
		// GROUP | Atom_1_idx | Atom_2_idx | ...
		// A line without exactly six tokens should be discarded.
		if (numTokens < 3) {
		    printf("Error: Invalid group file line: %s\n", line);
		    exit(-1);
		}

		// Make sure the index of this particle is unique.
		// NOTE: The particle list is sorted by index. 
		std::vector<int> tmp;
		for (int i=1; i < numTokens; ++i) {
		    const int ai = atoi(tokenList[i].val());
		    if (ai >= num+num_rb_attached_particles) {
			printf("Error: Attempted to include invalid particle in group: %s\n", line);
			exit(-1);
		    } else if (ai >= num) {
			printf("WARNING: including RB particles in group with line: %s\n", line);
		    }
		    tmp.push_back( ai );
		}
		groupSiteData.push_back(tmp);
		numGroupSites++;
	}
}

void Configuration::readBonds() {
	// Open the file
	FILE* inp = fopen(bondFile.val(), "r");
	char line[256];

	// If the particle file cannot be found, exit the program
	if (inp == NULL) {
		printf("WARNING: Could not open `%s'.\n", bondFile.val());
		printf("         This simulation will not use particle bonds.\n");
		return;
	}

	// Our particle array has a starting capacity of 256
	// We will expand this later if we need to.
	int capacity = 256;
	numBonds = 0;
	bonds = new Bond[capacity];

	// Get and process all lines of input
	while (fgets(line, 256, inp) != NULL) {
		
		// Lines in the particle file that begin with # are comments
		if (line[0] == '#') continue;
		      
		String s(line);
		int numTokens = s.tokenCount();

		// Break the line down into pieces (tokens) so we can process them individually
		String* tokenList = new String[numTokens];
		s.tokenize(tokenList);

		// Legitimate BOND input lines have 4 tokens: 
		// BOND | OPERATION_FLAG | INDEX1 | INDEX2 | FILENAME 
		// A line without exactly five tokens should be discarded.
		if (numTokens != 5) {
			printf("WARNING: Invalid bond file line: %s\n", line);
			continue;
		}

		String op = tokenList[1];
		int ind1 = atoi(tokenList[2].val());
		int ind2 = atoi(tokenList[3].val());
		String file_name = tokenList[4];

		if (ind1 == ind2) {
			printf("WARNING: Invalid bond file line: %s\n", line);
			continue;
		}

		if (ind1 < 0 || ind1 >= num+num_rb_attached_particles+numGroupSites ||
		    ind2 < 0 || ind2 >= num+num_rb_attached_particles+numGroupSites) {
			printf("ERROR: Bond file line '%s' includes invalid index\n", line);
			exit(1);
		}

		
		// If we don't have enough room in our bond array, we need to expand it.
		if (numBonds+1 >= capacity) { // "numBonds+1" because we are adding two bonds to array
			// Temporary pointer to the old array
			Bond* temp = bonds;	

			// Double the capacity
			capacity *= 2;

			// Create pointer to new array which is twice the size of the old one
			bonds = new Bond[capacity];
		
			// Copy the old values into the new array
			for (int j = 0; j < numBonds; j++)
				bonds[j] = temp[j];

			// delete the old array
			delete[] temp;
		}
		// Add the bond to the bond array
		// We must add it twice: Once for (ind1, ind2) and once for (ind2, ind1)
		
		// RBTODO: add ind1/2 to exclusion list here iff op == REPLACE

		if (op == "REPLACE")
		    addExclusion(ind1, ind2);

		Bond* b = new Bond(op, ind1, ind2, file_name);
		bonds[numBonds++] = *b;
		b = new Bond(op, ind2, ind1, file_name);
		bonds[numBonds++] = *b;
		delete[] tokenList;
	}	
	// Call compareBondIndex with qsort to sort the bonds by BOTH ind1 AND ind2
	std::sort(bonds, bonds + numBonds, compare());

	/* Each particle may have a varying number of bonds
	 * bondMap is an array with one element for each particle
	 * which keeps track of where a particle's bonds are stored
	 * in the bonds array.
	 * bondMap[i].x is the index in the bonds array where the ith particle's bonds begin
	 * bondMap[i].y is the index in the bonds array where the ith particle's bonds end
	 */
	bondMap = new int2[num+num_rb_attached_particles+numGroupSites];
	for (int i = 0; i < num+num_rb_attached_particles+numGroupSites; i++) {
		bondMap[i].x = -1;
		bondMap[i].y = -1;
	}
	int currPart = -1;
	int lastPart = -1;
	for (int i = 0; i < numBonds; i++) {
		if (bonds[i].ind1 != currPart) {
			currPart = bonds[i].ind1;
			bondMap[currPart].x = i;
			if (lastPart >= 0) bondMap[lastPart].y = i;
			lastPart = currPart;
		}
	}
	if (bondMap[lastPart].x > 0)
		bondMap[lastPart].y = numBonds;
}

void Configuration::readExcludes()
{
	// Open the file
	FILE* inp = fopen(excludeFile.val(), "r");
	char line[256];

	// If the exclusion file cannot be found, exit the program
	if (inp == NULL) {
		printf("WARNING: Could not open `%s'.\n", excludeFile.val());
		printf("This simulation will not use exclusions.\n");
		return;
	}


	// Get and process all lines of input
	while (fgets(line, 256, inp) != NULL) {
		// Lines in the particle file that begin with # are comments
		if (line[0] == '#') continue;
		      
		String s(line);
		int numTokens = s.tokenCount();
		      
		// Break the line down into pieces (tokens) so we can process them individually
		String* tokenList = new String[numTokens];
		s.tokenize(tokenList);

		// Legitimate EXCLUDE input lines have 3 tokens: 
		// BOND | INDEX1 | INDEX2
		// A line without exactly three tokens should be discarded.
		if (numTokens != 3) {
			printf("WARNING: Invalid exclude file line: %s\n", line);
			continue;
		}
		int ind1 = atoi(tokenList[1].val());
		int ind2 = atoi(tokenList[2].val());
		addExclusion(ind1, ind2);
		delete[] tokenList;
	}
}
void Configuration::addExclusion(int ind1, int ind2) {
    if (ind1 >= num+num_rb_attached_particles || ind2 >= num+num_rb_attached_particles) {
	printf("WARNING: Attempted to add an exclusion for an out-of-range particle index (%d or %d >= %d).\n", ind1, ind2, num+num_rb_attached_particles);
	return;
    }
		
    // If we don't have enough room in our bond array, we need to expand it.
    if (numExcludes >= excludeCapacity) {
	// Temporary pointer to the old array
	Exclude* temp = excludes;	

	// Double the capacity
	excludeCapacity *= 2;

	// Create pointer to new array which is twice the size of the old one
	excludes = new Exclude[excludeCapacity];
		
	// Copy the old values into the new array
	for (int j = 0; j < numExcludes; j++)
	    excludes[j] = temp[j];

	// delete the old array
	delete[] temp;
    }

    // Add the bond to the exclude array
    // We must add it twice: Once for (ind1, ind2) and once for (ind2, ind1)
    Exclude ex1(ind1, ind2);
    excludes[numExcludes++] = ex1;
    Exclude ex2(ind2, ind1);
    excludes[numExcludes++] = ex2;
    
}    

void Configuration::buildExcludeMap() {
    // Call compareExcludeIndex with qsort to sort the excludes by BOTH ind1 AND ind2
    std::sort(excludes, excludes + numExcludes, compare());

    /* Each particle may have a varying number of excludes
     * excludeMap is an array with one element for each particle
     * which keeps track of where a particle's excludes are stored
     * in the excludes array.
     * excludeMap[i].x is the index in the excludes array where the ith particle's excludes begin
     * excludeMap[i].y is the index in the excludes array where the ith particle's excludes end
     */
    excludeMap = new int2[num+num_rb_attached_particles];
    for (int i = 0; i < num+num_rb_attached_particles; i++) {
	excludeMap[i].x = -1;
	excludeMap[i].y = -1;
    }
    int currPart = -1;
    int lastPart = -1;
    for (int i = 0; i < numExcludes; i++) {
	if (excludes[i].ind1 != currPart) {
	    currPart = excludes[i].ind1;
	    assert(currPart < num+num_rb_attached_particles);
	    excludeMap[currPart].x = i;
	    if (lastPart >= 0)
		excludeMap[lastPart].y = i;
	    lastPart = currPart;
	}
    }
    if (excludeMap[lastPart].x > 0)
	excludeMap[lastPart].y = numExcludes;
}

void Configuration::readAngles() {
	FILE* inp = fopen(angleFile.val(), "r");
	char line[256];
	int capacity = 256;
	numAngles = 0;
	angles = new Angle[capacity];

	// If the angle file cannot be found, exit the program
	if (inp == NULL) {
		printf("WARNING: Could not open `%s'.\n", angleFile.val());
		printf("This simulation will not use angles.\n");
		return;
	}

	while(fgets(line, 256, inp)) {
		if (line[0] == '#') continue;
		String s(line);
		int numTokens = s.tokenCount();
		String* tokenList = new String[numTokens];
		s.tokenize(tokenList);
		
		// Legitimate ANGLE inputs have 5 tokens
		// ANGLE | INDEX1 | INDEX2 | INDEX3 | FILENAME
		// Any angle input line without exactly 5 tokens should be discarded
		if (numTokens != 5) {
			printf("WARNING: Invalid angle input line: %s\n", line);
			continue;
		}		
		
		// Discard any empty line
		if (tokenList == NULL) 
			continue;
		
		int ind1 = atoi(tokenList[1].val());
		int ind2 = atoi(tokenList[2].val());
		int ind3 = atoi(tokenList[3].val());
		String file_name = tokenList[4];
		//printf("file_name %s\n", file_name.val());
		if (ind1 >= num+num_rb_attached_particles+numGroupSites or ind2 >= num+num_rb_attached_particles+numGroupSites or ind3 >= num+num_rb_attached_particles+numGroupSites)
			continue;

		if (numAngles >= capacity) {
			Angle* temp = angles;
			capacity *= 2;
			angles = new Angle[capacity];
			for (int i = 0; i < numAngles; i++)
				angles[i] = temp[i];
			delete[] temp;
		}

		Angle a(ind1, ind2, ind3, file_name);
		angles[numAngles++] = a;
		delete[] tokenList;
	}
	std::sort(angles, angles + numAngles, compare());	

	// for(int i = 0; i < numAngles; i++)
	// 	angles[i].print();
}

void Configuration::readDihedrals() {
	FILE* inp = fopen(dihedralFile.val(), "r");
	char line[256];
	int capacity = 256;
	numDihedrals = 0;
	dihedrals = new Dihedral[capacity];

	// If the dihedral file cannot be found, exit the program
	if (inp == NULL) {
		printf("WARNING: Could not open `%s'.\n", dihedralFile.val());
		printf("This simulation will not use dihedrals.\n");
		return;
	}

	while(fgets(line, 256, inp)) {
		if (line[0] == '#') continue;
		String s(line);
		int numTokens = s.tokenCount();
		String* tokenList = new String[numTokens];
		s.tokenize(tokenList);
		
		// Legitimate DIHEDRAL inputs have 6 tokens
		// DIHEDRAL | INDEX1 | INDEX2 | INDEX3 | INDEX4 | FILENAME
		// Any angle input line without exactly 6 tokens should be discarded
		if (numTokens != 6) {
			printf("WARNING: Invalid dihedral input line: %s\n", line);
			continue;
		}		
		
		// Discard any empty line
		if (tokenList == NULL) 
			continue;
		
		int ind1 = atoi(tokenList[1].val());
		int ind2 = atoi(tokenList[2].val());
		int ind3 = atoi(tokenList[3].val());
		int ind4 = atoi(tokenList[4].val());
		String file_name = tokenList[5];
		//printf("file_name %s\n", file_name.val());
		if (ind1 >= num+num_rb_attached_particles+numGroupSites or
		    ind2 >= num+num_rb_attached_particles+numGroupSites or
		    ind3 >= num+num_rb_attached_particles+numGroupSites or
		    ind4 >= num+num_rb_attached_particles+numGroupSites)
			continue;

		if (numDihedrals >= capacity) {
			Dihedral* temp = dihedrals;
			capacity *= 2;
			dihedrals = new Dihedral[capacity];
			for (int i = 0; i < numDihedrals; ++i)
				dihedrals[i] = temp[i];
			delete[] temp;
		}

		Dihedral d(ind1, ind2, ind3, ind4, file_name);
		dihedrals[numDihedrals++] = d;
		delete[] tokenList;
	}
	std::sort(dihedrals, dihedrals + numDihedrals, compare());	

	// for(int i = 0; i < numDihedrals; i++)
	// 	dihedrals[i].print();
}

void Configuration::readBondAngles() {
	FILE* inp = fopen(bondAngleFile.val(), "r");
	char line[256];
	int capacity = 256;
	numBondAngles = 0;
	bondAngles = new BondAngle[capacity];

	// If the angle file cannot be found, exit the program
	if (inp == NULL) {
		printf("WARNING: Could not open `%s'.\n", bondAngleFile.val());
		printf("This simulation will not use angles.\n");
		return;
	}

	while(fgets(line, 256, inp)) {
		if (line[0] == '#') continue;
		String s(line);
		int numTokens = s.tokenCount();
		String* tokenList = new String[numTokens];
		s.tokenize(tokenList);

		// Legitimate BONDANGLE inputs have 8 tokens
		// BONDANGLE | INDEX1 | INDEX2 | INDEX3 | INDEX4 | ANGLE_FILENAME | BOND_FILENAME1 | BOND_FILENAME2
		if (numTokens != 8) {
			printf("WARNING: Invalid bond_angle input line: %s\n", line);
			continue;
		}

		// Discard any empty line
		if (tokenList == NULL)
			continue;

		int ind1 = atoi(tokenList[1].val());
		int ind2 = atoi(tokenList[2].val());
		int ind3 = atoi(tokenList[3].val());
		int ind4 = atoi(tokenList[4].val());
		String file_name1 = tokenList[5];
		String file_name2 = tokenList[6];
		String file_name3 = tokenList[7];
		//printf("file_name %s\n", file_name.val());
		if (ind1 >= num or ind2 >= num or ind3 >= num or ind4 >= num)
			continue;

		if (numBondAngles >= capacity) {
			BondAngle* temp = bondAngles;
			capacity *= 2;
			bondAngles = new BondAngle[capacity];
			for (int i = 0; i < numBondAngles; i++)
				bondAngles[i] = temp[i];
			delete[] temp;
		}

		BondAngle a(ind1, ind2, ind3, ind4, file_name1, file_name2, file_name3);
		bondAngles[numBondAngles++] = a;
		delete[] tokenList;
	}
	std::sort(bondAngles, bondAngles + numBondAngles, compare());

	// for(int i = 0; i < numAngles; i++)
	// 	angles[i].print();
}

void Configuration::readProductPotentials() {
	FILE* inp = fopen(productPotentialFile.val(), "r");
	char line[256];
	int capacity = 256;
	numProductPotentials = 0;
	productPotentials = new ProductPotentialConf[capacity];

	// If the angle file cannot be found, exit the program
	if (inp == NULL) {
		printf("WARNING: Could not open `%s'.\n", productPotentialFile.val());
		printf("This simulation will not use product potentials.\n");
		return;
	}
	printf("DEBUG: READING PRODUCT POTENTAL FILE\n");
	std::vector<std::vector<int>> indices;
	std::vector<int> tmp;
	std::vector<String> pot_names;

	while(fgets(line, 256, inp)) {
		if (line[0] == '#') continue;
		String s(line);
		int numTokens = s.tokenCount();
		String* tokenList = new String[numTokens];
		s.tokenize(tokenList);

		indices.clear();
		tmp.clear();
		pot_names.clear();		    

		printf("\rDEBUG: reading line %d",numProductPotentials+1);

		// Legitimate ProductPotential inputs have at least 7 tokens
		// BONDANGLE | INDEX1 | INDEX2 | INDEX3 | [TYPE1] | POT_FILENAME1 | INDEX4 | INDEX5 | [TYPE2] POT_FILENAME2 ...
		if (numTokens < 5) {
		    printf("WARNING: Invalid product potential input line (too few tokens %d): %s\n", numTokens, line);
			continue;
		}

		// Discard any empty line
		if (tokenList == NULL)
			continue;

		SimplePotentialType type = BOND; // initialize to suppress warning
		bool type_specified = false;
		for (int i = 1; i < numTokens; ++i) {
		    char *end;
		    // printf("DEBUG: Working on token %d '%s'\n", i, tokenList[i].val());

		    // Try to convert token to integer
		    int index = (int) strtol(tokenList[i].val(), &end, 10);
		    if (tokenList[i].val() == end || *end != '\0' || errno == ERANGE) {
			// Failed to convert token to integer; therefore it must be a potential name or type

			// Try to match a type
			String n = tokenList[i];
			n.lower();
			if (n == "bond") { type = DIHEDRAL; type_specified = true; }
			else if (n == "angle")  { type = DIHEDRAL; type_specified = true; }
			else if (n == "dihedral")  { type = DIHEDRAL; type_specified = true; }
			else if (n == "vecangle") { type = VECANGLE; type_specified = true; }
			else { // Not a type, therefore a path to a potential
			    n = tokenList[i];
			    indices.push_back(tmp);
			    pot_names.push_back( n );
			    // TODO: Key should be tuple of (type,n)
			    std::string n_str = std::string(n.val());
			    if ( simple_potential_ids.find(n_str) == simple_potential_ids.end() ) {
				// Could not find fileName in dictionary, so read and add it
				unsigned int s = tmp.size();
				if (s < 2 || s > 4) {
				    printf("WARNING: Invalid product potential input line (indices of potential %d == %d): %s\n", i, s, line);
				    continue;
				}
				simple_potential_ids[ n_str ] = simple_potentials.size();
				if (not type_specified) type = s==2? BOND: s==3? ANGLE: DIHEDRAL;
				simple_potentials.push_back( SimplePotential(n.val(), type) );
			    }
			    tmp.clear();
			    type_specified = false;

			}
		    } else {
			if (index >= num) {
			    continue;
			}
			tmp.push_back(index);
		    }
		}

		if (numProductPotentials >= capacity) {
			ProductPotentialConf* temp = productPotentials;
			capacity *= 2;
			productPotentials = new ProductPotentialConf[capacity];
			for (int i = 0; i < numProductPotentials; i++)
				productPotentials[i] = temp[i];
			delete[] temp;
		}

		ProductPotentialConf a(indices, pot_names);
		productPotentials[numProductPotentials++] = a;
		delete[] tokenList;
	}
	printf("\nDEBUG: Sorting\n");
	std::sort(productPotentials, productPotentials + numProductPotentials, compare());

	// for(int i = 0; i < numAngles; i++)
	// 	angles[i].print();
}


void Configuration::readRestraints() {
	FILE* inp = fopen(restraintFile.val(), "r");
	char line[256];
	int capacity = 16;
	numRestraints = 0;
	restraints = new Restraint[capacity];

	// If the restraint file cannot be found, exit the program
	if (inp == NULL) {
		printf("WARNING: Could not open `%s'.\n", restraintFile.val());
		printf("  This simulation will not use restraints.\n");
		return;
	}

	while(fgets(line, 256, inp)) {
		if (line[0] == '#') continue;
		String s(line);
		int numTokens = s.tokenCount();
		String* tokenList = new String[numTokens];
		s.tokenize(tokenList);

		// inputs have 6 tokens
		// RESTRAINT | INDEX1 | k | x0 | y0 | z0
		if (numTokens != 6) {
			printf("WARNING: Invalid restraint input line: %s\n", line);
			continue;
		}

		// Discard any empty line
		if (tokenList == NULL) continue;

		int   id = atoi(tokenList[1].val());
		float k  = (float) strtod(tokenList[2].val(), NULL);
		float x0 = (float) strtod(tokenList[3].val(), NULL);
		float y0 = (float) strtod(tokenList[4].val(), NULL);
		float z0 = (float) strtod(tokenList[5].val(), NULL);

		if (id >= num + num_rb_attached_particles + numGroupSites) continue;

		if (numRestraints >= capacity) {
			Restraint* temp = restraints;
			capacity *= 2;
			restraints = new Restraint[capacity];
			for (int i = 0; i < numRestraints; ++i)
				restraints[i] = temp[i];
			delete[] temp;
		}

		Restraint tmp(id, Vector3(x0,y0,z0), k);
		restraints[numRestraints++] = tmp;
		delete[] tokenList;
	}
	// std::sort(restraints, restraints + numRestraints, compare());
}

//populate the type list and serial list
void Configuration::populate() {
    for (int repID = 0; repID < simNum; ++repID) {
                const int offset = repID * num;
                int pn = 0;
                int p = 0;
                for (int i = 0; i < num; ++i) {
                        type[i + offset] = p;
                        serial[i + offset] = currSerial++;

                        if (++pn >= part[p].num) {
                                p++;
                                pn = 0;
                        }
                }
        }
}

bool Configuration::readBondFile(const String& value, int currBond) {
	int numTokens = value.tokenCount();
	if (numTokens != 1) {
		printf("ERROR: Invalid tabulatedBondFile: %s, numTokens = %d\n", value.val(), numTokens);
		return false;
	}

	String* tokenList = new String[numTokens];
	value.tokenize(tokenList);
	if (tokenList == NULL) {
		printf("ERROR: Invalid tabulatedBondFile: %s; tokenList is NULL\n", value.val());
		return false;
	}

	bondTableFile[currBond] = tokenList[0];

	// printf("Tabulated Bond Potential: %s\n", bondTableFile[currBond].val() );

	return true;
}

bool Configuration::readAngleFile(const String& value, int currAngle) {
	int numTokens = value.tokenCount();
	if (numTokens != 1) {
		printf("ERROR: Invalid tabulatedAngleFile: %s, numTokens = %d\n", value.val(), numTokens);
		return false;
	}

	String* tokenList = new String[numTokens];
	value.tokenize(tokenList);
	if (tokenList == NULL) {
		printf("ERROR: Invalid tabulatedAngleFile: %s; tokenList is NULL\n", value.val());
		return false;
	}

	angleTableFile[currAngle] = tokenList[0];

	// printf("Tabulated Angle Potential: %s\n", angleTableFile[currAngle].val() );

	return true;
}

bool Configuration::readDihedralFile(const String& value, int currDihedral) {
	int numTokens = value.tokenCount();
	if (numTokens != 1) {
		printf("ERROR: Invalid tabulatedDihedralFile: %s, numTokens = %d\n", value.val(), numTokens);
		return false;
	}

	String* tokenList = new String[numTokens];
	value.tokenize(tokenList);
	if (tokenList == NULL) {
		printf("ERROR: Invalid tabulatedDihedralFile: %s; tokenList is NULL\n", value.val());
		return false;
	}

	dihedralTableFile[currDihedral] = tokenList[0];

	// printf("Tabulated Dihedral Potential: %s\n", dihedralTableFile[currDihedral].val() );

	return true;
}
//Load the restart coordiantes only
void Configuration::loadRestart(const char* file_name) {
	char line[STRLEN];
	FILE* inp = fopen(file_name, "r");

	if (inp == NULL) {
		printf("GrandBrownTown:loadRestart File `%s' does not exist\n", file_name);
		exit(-1);
	}

	int count = 0;
	while (fgets(line, STRLEN, inp) != NULL) {
		// Ignore comments.
		int len = strlen(line);
		if (line[0] == '#') continue;
		if (len < 2) continue;

		String s(line);
		int numTokens = s.tokenCount();
		if (numTokens != 4) {
			printf("GrandBrownTown:loadRestart Invalid coordinate file line: %s\n", line);
			fclose(inp);	
			exit(-1);
		}

		String* tokenList = new String[numTokens];
		s.tokenize(tokenList);
		if (tokenList == NULL) {
			printf("GrandBrownTown:loadRestart Invalid coordinate file line: %s\n", line);
			fclose(inp);
			exit(-1);
		}

		int typ = atoi(tokenList[0]);
		float x = (float) strtod(tokenList[1],NULL);
		float y = (float) strtod(tokenList[2],NULL);
		float z = (float) strtod(tokenList[3],NULL);

		pos[count] = Vector3(x,y,z);
		type[count] = typ;
		serial[count] = currSerial;
		currSerial++;
		if (typ < 0 || typ >= numParts) {
			printf("GrandBrownTown:countRestart Invalid particle type: %d\n", typ);
			fclose(inp);
			exit(-1);
		}

		count++;
		delete[] tokenList;
	}

	fclose(inp);    
}
//Han-Yi Chou
//First the resart coordinates should be loaded
void Configuration::loadRestartMomentum(const char* file_name) 
{
    char line[STRLEN];
    FILE* inp = fopen(file_name, "r");

    if (inp == NULL) 
    {
        printf("GrandBrownTown:loadRestart File `%s' does not exist\n", file_name);
        exit(-1);
    }
    if(!loadedCoordinates)
    {
        printf("First load the restart coordinates\n");
        assert(1==2);
    }
    int count = 0;
    while (fgets(line, STRLEN, inp) != NULL) 
    {
        // Ignore comments.
        int len = strlen(line);
        if (line[0] == '#') continue;
        if (len < 2) continue;

        String s(line);
        int numTokens = s.tokenCount();
        if (numTokens != 4) 
        {
            printf("GrandBrownTown:loadRestart Invalid momentum file line: %s\n", line);
            fclose(inp);
            exit(-1);
        }

        String* tokenList = new String[numTokens];
        s.tokenize(tokenList);
        if (tokenList == NULL) 
        {
            printf("GrandBrownTown:loadRestart Invalid momentum file line: %s\n", line);
            fclose(inp);
            exit(-1);
        }

        int typ = atoi(tokenList[0]);
        float x = (float) strtod(tokenList[1],NULL);
        float y = (float) strtod(tokenList[2],NULL);
        float z = (float) strtod(tokenList[3],NULL);

        if (typ < 0 || typ >= numParts) 
        {
            printf("GrandBrownTown:countRestart Invalid particle type : %d\n", typ);
            fclose(inp);
            exit(-1);
        }

        if(typ != type[count])
        {
            printf("Inconsistent in momentum file with the position file\n");
            fclose(inp);
            exit(-1);
        }
        momentum[count] = Vector3(x,y,z);
        ++count;
        delete[] tokenList;
    }
    fclose(inp);
}

bool Configuration::loadCoordinates(const char* file_name) {
	char line[STRLEN];
	FILE* inp = fopen(file_name, "r");

	if (inp == NULL) {
	    printf("ERROR: Could not open file for reading: %s\n", file_name);
	    exit(-1);
	    return false;
	}

	int count = 0;
	while (fgets(line, STRLEN, inp) != NULL) {
		// Ignore comments.
		int len = strlen(line);
		if (line[0] == '#') continue;
		if (len < 2) continue;

		String s(line);
		int numTokens = s.tokenCount();
		if (numTokens != 3) {
			printf("ERROR: Invalid coordinate file line: %s\n", line);
			fclose(inp);	
			return false;
		}

		String* tokenList = new String[numTokens];
		s.tokenize(tokenList);
		if (tokenList == NULL) {
			printf("ERROR: Invalid coordinate file line: %s\n", line);
			fclose(inp);
			return false;
		}

		if (count >= num*simNum) {
			printf("WARNING: Too many coordinates in coordinate file %s.\n", file_name);
			fclose(inp);
			return true;
		}

		float x = (float) strtod(tokenList[0],NULL);
		float y = (float) strtod(tokenList[1],NULL);
		float z = (float) strtod(tokenList[2],NULL);
		pos[count] = Vector3(x,y,z);
		count++;

		delete[] tokenList;
	}
	fclose(inp);

	if (count < num) {
		printf("ERROR: Too few coordinates in coordinate file.\n");
		return false;
	}
	return true;
}
//Han-Yi Chou The function populate should be called before entering this function
bool Configuration::loadMomentum(const char* file_name) 
{
    char line[STRLEN];
    FILE* inp = fopen(file_name, "r");

    if (inp == NULL) 
        return false;

    int count = 0;
    while (fgets(line, STRLEN, inp) != NULL) 
    {
        // Ignore comments.
        int len = strlen(line);
        if (line[0] == '#') 
            continue;
        if (len < 2) 
            continue;

        String s(line);
        int numTokens = s.tokenCount();
        if (numTokens != 3) 
        {
            printf("ERROR: Invalid momentum file line: %s\n", line);
            fclose(inp);
            return false;
        }

        String* tokenList = new String[numTokens];
        s.tokenize(tokenList);
        if (tokenList == NULL) 
        {
            printf("ERROR: Invalid momentum file line: %s\n", line);
            fclose(inp);
            return false;
        }

        if (count >= num) 
        {
            printf("WARNING: Too many momentum in momentum file %s.\n", file_name);
            fclose(inp);
            return false;
        }

        float x = (float) strtod(tokenList[0],NULL);
        float y = (float) strtod(tokenList[1],NULL);
        float z = (float) strtod(tokenList[2],NULL);
        momentum[count] = Vector3(x,y,z);
        ++count;
        delete[] tokenList;
    }
    fclose(inp);

    if (count < num) 
    {
        printf("ERROR: Too few momentum in momentum file.\n");
        return false;
    }
    return true;
}

// Count the number of atoms in the restart file.
int Configuration::countRestart(const char* file_name) {
	char line[STRLEN];
	FILE* inp = fopen(file_name, "r");

	if (inp == NULL) {
		printf("ERROR: countRestart File `%s' does not exist\n", file_name);
		exit(-1);
	}

	int count = 0;
	while (fgets(line, STRLEN, inp) != NULL) {
		int len = strlen(line);
		// Ignore comments.
		if (line[0] == '#') continue;
		if (len < 2) continue;

		String s(line);
		int numTokens = s.tokenCount();
		if (numTokens != 4) {
			printf("ERROR: countRestart Invalid coordinate file line: %s\n", line);
			fclose(inp);	
			exit(-1);
		}

		String* tokenList = new String[numTokens];
		s.tokenize(tokenList);
		if (tokenList == NULL) {
			printf("ERROR: countRestart Invalid coordinate file line: %s\n", line);
			fclose(inp);
			exit(-1);
		}

		int typ = atoi(tokenList[0]);
		// float x = strtod(tokenList[1],NULL);
		// float y = strtod(tokenList[2],NULL);
		// float z = strtod(tokenList[3],NULL);
		if (typ < 0 || typ >= numParts) {
			printf("ERROR: countRestart Invalid particle type: %d\n", typ);
			fclose(inp);
			exit(-1);
		}

		count++;
		delete[] tokenList;
	}

	fclose(inp);    
	return count;
}

bool Configuration::readTableFile(const String& value, int currTab) {
	int numTokens = value.tokenCount('@');
	if (numTokens != 3) {
		printf("ERROR: Invalid tabulatedFile: %s\n", value.val());
		return false;
	}

	String* tokenList = new String[numTokens];
	value.tokenize(tokenList, '@');
	if (tokenList == NULL) {
		printf("ERROR: Invalid tabulatedFile: %s\n", value.val());
		return false;
	}

	if (currTab >= numParts*numParts) {
	    printf("ERROR: Number of tabulatedFile entries exceeded %d*%d particle types.\n", numParts,numParts);
	    exit(1);
	}

	partTableIndex0[currTab] = atoi(tokenList[0]);
	partTableIndex1[currTab] = atoi(tokenList[1]);
	partTableFile[currTab] = tokenList[2];

	// printf("Tabulated Potential: %d %d %s\n", partTableIndex0[currTab],
	// 		partTableIndex1[currTab], partTableFile[currTab].val() );
	delete[] tokenList;
	return true;
}

void Configuration::getDebugForce() {
	// Allow the user to choose which force computation to use
	printf("\n");
	printf("(1) ComputeFull [Default]          (2) ComputeSoftcoreFull\n");
	printf("(3) ComputeElecFull                (4) Compute (Decomposed)\n");
	printf("(5) ComputeTabulated (Decomposed)  (6) ComputeTabulatedFull\n");

	printf("WARNING: ");
	if (tabulatedPotential) {
		if (fullLongRange) printf("(6) was specified by config file\n");
		else printf("(5) was specified by config file\n");
	} else {
		if (fullLongRange != 0) printf("(%d) was specified by config file\n", fullLongRange);
		else printf("(4) was specified by config file\n");
	}

	char buffer[256];
	int choice;
	while (true) {
		printf("Choose a force computation (1 - 6): ");
		fgets(buffer, 256, stdin);
		bool good = sscanf(buffer, "%d", &choice) && (choice >= 1 && choice <= 6);
		if (good)
			break;
	}
	switch(choice) {
		case 1:
			tabulatedPotential = 0;
			fullLongRange = 1;
			break;
		case 2:
			tabulatedPotential = 0;
			fullLongRange = 2;
			break;
		case 3:
			tabulatedPotential = 0;
			fullLongRange = 3;
			break;
		case 4:
			tabulatedPotential = 0;
			fullLongRange = 0;
			break;
		case 5:
			tabulatedPotential = 1;
			fullLongRange = 0;
			break;
		case 6:
			tabulatedPotential = 1;
			fullLongRange = 1;
			break;
		default:
			tabulatedPotential = 0;
			fullLongRange = 1;
			break;
	}
	printf("\n");
}
//Han-Yi Chou setting boltzman distribution of momentum with a given center of mass velocity
//Before using this code, make sure the array type list and serial list are both already initialized
bool Configuration::Boltzmann(const Vector3& v_com, int N)
{
    int count = 0;
    Vector3 total_momentum = Vector3(0.);

    RandomCPU random = RandomCPU(seed + 2); /* +2 to avoid using same seed elsewhere */

    for(int i = 0; i < N; ++i)
    {
        int typ = type[i];
        double M = part[typ].mass;
        double sigma = sqrt(kT * M) * 2.046167337e4;
   
        Vector3 tmp = random.gaussian_vector() * sigma;
        tmp = tmp * 1e-4;
        total_momentum += tmp;
        momentum[(size_t)count] = tmp;
        ++count;
    }
    if(N > 1)
    {
        total_momentum = total_momentum / (double)N;

        for(int i = 0; i < N; ++i)
        {
            int typ = type[i];
            double M = part[typ].mass;
        
            momentum[i] = momentum[i] - total_momentum + M * v_com;
        }
    }

    
    return true;
} 

//////////////////////////
// Comparison operators //
//////////////////////////
bool Configuration::compare::operator()(const String& lhs, const String& rhs) {
	String* list_lhs = new String[lhs.tokenCount()];
	String* list_rhs = new String[rhs.tokenCount()];
	lhs.tokenize(list_lhs);
	rhs.tokenize(list_rhs);
	int key_lhs = atoi(list_lhs[1].val());
	int key_rhs = atoi(list_rhs[1].val());
	delete[] list_lhs;
	delete[] list_rhs;
	return key_lhs < key_rhs;
}

bool Configuration::compare::operator()(const Bond& lhs, const Bond& rhs) {
	int diff = lhs.ind1 - rhs.ind1;
	if (diff != 0)
		return lhs.ind1 < rhs.ind1;
	return lhs.ind2 < rhs.ind2;
}

bool Configuration::compare::operator()(const Exclude& lhs, const Exclude& rhs) {
	int diff = lhs.ind1 - rhs.ind1;
	if (diff != 0)
		return lhs.ind1 < rhs.ind1;
	return lhs.ind2 < rhs.ind2;
}

bool Configuration::compare::operator()(const Angle& lhs, const Angle& rhs) {
	int diff = lhs.ind1 - rhs.ind1;
	if (diff != 0)
		return lhs.ind1 < rhs.ind1;
	diff = lhs.ind2 - rhs.ind2;
	if (diff != 0)
		return lhs.ind2 < rhs.ind2;
	return lhs.ind3 < rhs.ind3;
}

bool Configuration::compare::operator()(const Dihedral& lhs, const Dihedral& rhs) {
	int diff = lhs.ind1 - rhs.ind1;
	if (diff != 0) 
		return lhs.ind1 < rhs.ind1;
	diff = lhs.ind2 - rhs.ind2;
	if (diff != 0) 
		return lhs.ind2 < rhs.ind2;
	diff = lhs.ind3 - rhs.ind3;
	if (diff != 0) 
		return lhs.ind3 < rhs.ind3;
	return lhs.ind4 < rhs.ind4;
}

bool Configuration::compare::operator()(const BondAngle& lhs, const BondAngle& rhs) {
	int diff = lhs.ind1 - rhs.ind1;
	if (diff != 0)
		return lhs.ind1 < rhs.ind1;
	diff = lhs.ind2 - rhs.ind2;
	if (diff != 0)
		return lhs.ind2 < rhs.ind2;
	diff = lhs.ind3 - rhs.ind3;
	if (diff != 0) 
		return lhs.ind3 < rhs.ind3;
	return lhs.ind4 < rhs.ind4;
}

bool Configuration::compare::operator()(const ProductPotentialConf& lhs, const ProductPotentialConf& rhs) {
    int diff = rhs.indices.size() - lhs.indices.size();
    if (diff != 0) return diff > 0;

    for (unsigned int i = 0; i < lhs.indices.size(); ++i) {
	diff = rhs.indices[i].size() - lhs.indices[i].size();
	if (diff != 0) return diff > 0;
    }

    for (unsigned int i = 0; i < lhs.indices.size(); ++i) {
	for (unsigned int j = 0; j < lhs.indices[i].size(); ++j) {
	    diff = rhs.indices[i][j] - lhs.indices[i][j];
	    if (diff != 0) return diff > 0;
	}
    }
    return true;
}
