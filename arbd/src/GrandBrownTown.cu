#include "GrandBrownTown.h"
#include "GrandBrownTown.cuh"
/* #include "ComputeGridGrid.cuh" */
#include "WKFUtils.h"
#include "BrownParticlesKernel.h"
#include "nvtx_defs.h"
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <thrust/device_ptr.h>
#include <fstream>
#include <cuda_profiler_api.h>

#ifdef _OPENMP
#include <omp.h>
//#else
//typedef int omp_int_t;
//inline omp_int_t omp_get_thread_num() { return 0; }
//inline omp_int_t omp_get_max_threads() { return 1; }
#endif

#ifndef gpuErrchk
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true) {
	if (code != cudaSuccess) {
		fprintf(stderr,"CUDA Error: %s %s %d\n", cudaGetErrorString(code), file, line);
	if (abort) exit(code);
	}
}
#endif

bool GrandBrownTown::DEBUG;

cudaEvent_t START, STOP;

GPUManager GrandBrownTown::gpuman = GPUManager();

GrandBrownTown::GrandBrownTown(const Configuration& c, const char* outArg,
		bool debug, bool imd_on, unsigned int imd_port, int numReplicas) :
	imd_on(imd_on), imd_port(imd_port), numReplicas(numReplicas),
	//conf(c), RBC(RigidBodyController(c,outArg)) {
	conf(c) {

        RBC.resize(numReplicas);      
        for(int i = 0; i < numReplicas; ++i)
        {
            RigidBodyController* rb = new RigidBodyController(c, outArg, seed, i);
            RBC[i] = rb;
        }

        //printf("%d\n",__LINE__);
        //Determine which dynamic. Han-Yi Chou
        particle_dynamic  = c.ParticleDynamicType;
        rigidbody_dynamic = c.RigidBodyDynamicType;
        ParticleInterpolationType = c.ParticleInterpolationType;
        RigidBodyInterpolationType = c.RigidBodyInterpolationType;
        //particle_langevin_integrator = c.ParticleLangevinIntegrator;
        printf("%d\n",__LINE__);
	for (int i = 0; i < numReplicas; ++i) 
        {
		std::stringstream curr_file, restart_file, out_prefix;

		if (numReplicas > 1) {
		    curr_file << outArg << '.' << i << ".curr";
		    restart_file   << outArg << '.' << i << ".restart";
		    out_prefix << outArg << '.' << i;
		} else {
		    curr_file << outArg << ".curr";
		    restart_file   << outArg << ".restart";
		    out_prefix << outArg;
		}

                outCurrFiles.push_back(curr_file.str());
                restartFiles.push_back(restart_file.str());
                outFilePrefixes.push_back(out_prefix.str());

                //Han-Yi Chou for flush out the momentum
                if(particle_dynamic == String("Langevin") || particle_dynamic == String("NoseHooverLangevin"))
                {
                    std::stringstream restart_file_p, out_momentum_prefix, out_force_prefix;
                    restart_file_p << outArg << '.' << i << ".momentum.restart";
                    out_momentum_prefix << outArg << '.' << i << ".momentum";
                    //out_force_prefix << outArg << ".force." << i;

                    restartMomentumFiles.push_back(restart_file_p.str());//Han-Yi Chou
                    outMomentumFilePrefixes.push_back(out_momentum_prefix.str());
                    //outForceFilePrefixes.push_back(out_force_prefix.str());
                }           
	}

	GrandBrownTown::DEBUG = debug;
	sysDim = c.sysDim;
	sys = c.sys;

	// Particle variables
	partsFromFile = c.partsFromFile;
	indices = c.indices;
	numPartsFromFile = c.numPartsFromFile;  // number of particle types
	bonds = c.bonds;
	numCap = c.numCap;                      // max number of particles
	num = c.num;                            // current number of particles
	num_rb_attached_particles = c.num_rb_attached_particles;
	numGroupSites = c.numGroupSites;

	// Allocate arrays of positions, types and serial numbers
	pos    = new Vector3[(num+num_rb_attached_particles+numGroupSites) * numReplicas];  // [HOST] array of particles' positions.
        // Allocate arrays of momentum Han-Yi Chou
        if(particle_dynamic == String("Langevin") || particle_dynamic == String("NoseHooverLangevin"))
        {
            momentum = new Vector3[num * numReplicas]; //[HOST] array of particles' momentum
            if(particle_dynamic == String("NoseHooverLangevin"))
                random = new float[num * numReplicas];
        }
        //printf("%d\n",__LINE__);
        //for debug
        //force = new Vector3[num * numReplicas];

	type   = new     int[(num+num_rb_attached_particles) * numReplicas];  // [HOST] array of particles' types.
	serial = new     int[(num+num_rb_attached_particles) * numReplicas];  // [HOST] array of particles' serial numbers.

	// Allocate things for rigid body
	// RBC = RigidBodyController(c);
  // printf("About to devicePrint\n");
	// devicePrint<<<1,1>>>(&(c.rigidBody[0]));
	// devicePrint<<<1,1>>>(RBC.rbType_d);
	cudaDeviceSynchronize();
	// printf("Done with devicePrint\n");


	
	// Replicate identical initial conditions across all replicas
	for (int r = 0; r < numReplicas; ++r) {
	    std::copy(c.type, c.type + num+num_rb_attached_particles,
		      type + r*(num+num_rb_attached_particles));
	    std::copy(c.serial, c.serial + num + num_rb_attached_particles,
		      serial + r*(num+num_rb_attached_particles));
	  if (c.copyReplicaCoordinates > 0)
	    std::copy(c.pos, c.pos + num, pos + r*num);
	}
        if (c.copyReplicaCoordinates <= 0)
          std::copy(c.pos, c.pos + numReplicas*num, pos);

        //printf("%d\n",__LINE__); 
        //Han-Yi Chou
        if(particle_dynamic == String("Langevin") || particle_dynamic == String("NoseHooverLangevin"))
            std::copy(c.momentum,c.momentum+num*numReplicas,momentum);

        //printf("%d\n",__LINE__);

	currSerial = c.currSerial;  // serial number of the next new particle
	name = c.name;              // list of particle types! useful when 'numFluct == 1'
	posLast = c.posLast;        // previous positions of particles  (used for computing ionic current)
        momLast = c.momLast;
	timeLast = c.timeLast;      // previous time (used with posLast)
	minimumSep = c.minimumSep;  // minimum separation allowed when placing new particles

	// System parameters
	outputName = c.outputName;
	timestep = c.timestep;
	steps = c.steps;
	seed = c.seed;
	temperatureGridFile = c.temperatureGridFile;
	inputCoordinates = c.inputCoordinates;
	restartCoordinates = c.restartCoordinates;
	numberFluct = c.numberFluct;
	interparticleForce = c.interparticleForce;
	tabulatedPotential = c.tabulatedPotential;
	readBondsFromFile = c.readBondsFromFile;
	fullLongRange = c.fullLongRange;
	kT = c.kT;
	temperature = c.temperature;
	coulombConst = c.coulombConst;
	electricField = c.electricField;
	cutoff = c.cutoff;
	switchLen = c.switchLen;
	outputPeriod = c.outputPeriod;
	outputEnergyPeriod = c.outputEnergyPeriod;
	outputFormat = c.outputFormat;
	currentSegmentZ = c.currentSegmentZ;
	numberFluctPeriod = c.numberFluctPeriod;
	decompPeriod = c.decompPeriod;
	numCapFactor = c.numCapFactor;
	kTGrid = c.kTGrid;
	tGrid = c.tGrid;
	sigmaT = c.sigmaT;

	// Parameter files
	partTableFile = c.partTableFile;
	bondTableFile = c.bondTableFile;
	angleTableFile = c.angleTableFile;
	dihedralTableFile = c.dihedralTableFile;

	// Other parameters.
	switchStart = c.switchStart;
	maxInitialPot = c.maxInitialPot;
	initialZ = c.initialZ;

	// Particle parameters.
	part = c.part;
	numParts = c.numParts;
	numBonds = c.numBonds;
	numExcludes = c.numExcludes;
	numAngles = c.numAngles;
	numDihedrals = c.numDihedrals;
	partTableIndex0 = c.partTableIndex0;
	partTableIndex1 = c.partTableIndex1;

	numBondAngles = c.numBondAngles;

	numTabBondFiles = c.numTabBondFiles;
	bondMap = c.bondMap;
	// TODO: bondList = c.bondList;

	excludes = c.excludes;
        part = c.part;
	excludeMap = c.excludeMap;
	excludeRule = c.excludeRule;
	excludeCapacity = c.excludeCapacity;

	angles = c.angles;
	numTabAngleFiles = c.numTabAngleFiles;

	dihedrals = c.dihedrals;
	numTabDihedralFiles = c.numTabDihedralFiles;

	bondAngles = c.bondAngles;

	// Device parameters
	//type_d = c.type_d;
	part_d = c.part_d;
	sys_d = c.sys_d;
	kTGrid_d = c.kTGrid_d;
	//bonds_d = c.bonds_d;
	//bondMap_d = c.bondMap_d;
	//excludes_d = c.excludes_d;
	//excludeMap_d = c.excludeMap_d;
	//angles_d = c.angles_d;
	//dihedrals_d = c.dihedrals_d;

	printf("Setting up random number generator with seed %lu\n", seed);
	randoGen = new Random(num * numReplicas, seed);
	copyRandToCUDA();

        if(particle_dynamic == String("NoseHooverLangevin"))
            InitNoseHooverBath(num * numReplicas);

	// "Some geometric stuff that should be gotten rid of." -- Jeff Comer
	Vector3 buffer = (sys->getCenter() + 2.0f * sys->getOrigin())/3.0f;
	initialZ = buffer.z;

	// Load random coordinates if necessary Han-Yi Chou
	if (!c.loadedCoordinates) {
		//printf("Populating\n");
		//populate(); Han-Yi Chou, Actually the list is already populated 
	    initialCond();
		printf("Setting random initial conditions.\n");
	}

	// Prepare internal force computation
	 //internal = new ComputeForce(num, part, numParts, sys, switchStart, switchLen, coulombConst,
	 //			    fullLongRange, numBonds, numTabBondFiles, numExcludes, numAngles, numTabAngleFiles,
	 //			    numDihedrals, numTabDihedralFiles, c.pairlistDistance, numReplicas);
	internal = new ComputeForce(c, numReplicas);

	//MLog: I did the other halve of the copyToCUDA function from the Configuration class here, keep an eye on any mistakes that may occur due to the location.
	internal -> copyToCUDA(c.simNum, c.type, c.bonds, c.bondMap, c.excludes, c.excludeMap, c.angles, c.dihedrals, c.restraints, c.bondAngles, c.simple_potential_ids, c.simple_potentials, c.productPotentials );
	if (numGroupSites > 0) init_cuda_group_sites();

	// TODO: check for duplicate potentials 
	if (c.tabulatedPotential) {
		printf("Loading %d tabulated non-bonded potentials...\n", numParts*numParts);
		for (int p = 0; p < numParts*numParts; p++) {
			if (partTableFile[p].length() > 0) {
				int type0 = partTableIndex0[p];
				int type1 = partTableIndex1[p];

				internal->addTabulatedPotential(partTableFile[p].val(), type0, type1);
				// printf("  Loaded %s for types %s and %s.\n", partTableFile[p].val(),
				// 		part[type0].name.val(), part[type1].name.val());
			}
		}
	}
	printf("Using %d non-bonded exclusions\n",c.numExcludes/2);

	if (c.readBondsFromFile) {
		printf("Loading %d tabulated bond potentials...\n", numTabBondFiles);
		for (int p = 0; p < numTabBondFiles; p++)
			if (bondTableFile[p].length() > 0) {
				//MLog: make sure to add to all GPUs
			    // printf("...loading %s\n",bondTableFile[p].val());
			    internal->addBondPotential(bondTableFile[p].val(), p, bonds, bondAngles);
				// printf("%s\n",bondTableFile[p].val());
			} else {
			    printf("...skipping %s (\n",bondTableFile[p].val());
			    internal->addBondPotential(bondTableFile[p].val(), p, bonds, bondAngles);
			}
			    
	}

	if (c.readAnglesFromFile) {
		printf("Loading %d tabulated angle potentials...\n", numTabAngleFiles);
		for (int p = 0; p < numTabAngleFiles; p++)
			if (angleTableFile[p].length() > 0)
			{
				//MLog: make sure to do this for every GPU
			    internal->addAnglePotential(angleTableFile[p].val(), p, angles, bondAngles);
			}
	}

	if (c.readDihedralsFromFile) {
		printf("Loading %d tabulated dihedral potentials...\n", numTabDihedralFiles);
		for (int p = 0; p < numTabDihedralFiles; p++)
			if (dihedralTableFile[p].length() > 0)
				internal->addDihedralPotential(dihedralTableFile[p].val(), p, dihedrals);
	}

	auto _get_index = [this](int idx, int replica) {
	    // Convenient lambda function to deal with increasingly complicated indexing
	    auto num = this->num;
	    auto numReplicas = this->numReplicas;
	    auto num_rb_attached_particles = this->num_rb_attached_particles;
	    auto numGroupSites = this->numGroupSites;
	    idx = (idx < num+num_rb_attached_particles) ? idx + replica*(num+num_rb_attached_particles)
		: (idx-num-num_rb_attached_particles) + numReplicas*(num+num_rb_attached_particles) + replica * numGroupSites;
	    return idx;
	};

	//Mlog: this is where we create the bondList.
	if (numBonds > 0) {
		bondList = new int3[ (numBonds / 2) * numReplicas ];
		int j = 0;

		for(int k = 0 ; k < numReplicas; k++)
		{
			for(int i = 0; i < numBonds; ++i)
			{
				if(bonds[i].ind1 < bonds[i].ind2)
				{
					if (bonds[i].tabFileIndex == -1) {
						fprintf(stderr,"Error: bondfile '%s' was not read with tabulatedBondFile command.\n", bonds[i].fileName.val());
						exit(1);
					}
					bondList[j] = make_int3( _get_index(bonds[i].ind1, k), _get_index(bonds[i].ind2, k), bonds[i].tabFileIndex );
					// cout << "Displaying: bondList["<< j <<"].x = " << bondList[j].x << ".\n"
					// << "Displaying: bondList["<< j <<"].y = " << bondList[j].y << ".\n"
					// << "Displaying: bondList["<< j <<"].z = " << bondList[j].z << ".\n";
					++j;
				}
			}
		}
	}
	// internal->createBondList(bondList);

	if (numAngles > 0) {
	angleList = new int4[ (numAngles) * numReplicas ];
	for(int k = 0 ; k < numReplicas; k++) {
	    for(int i = 0; i < numAngles; ++i) {
			if (angles[i].tabFileIndex == -1) {
				fprintf(stderr,"Error: anglefile '%s' was not read with tabulatedAngleFile command.\n", angles[i].fileName.val());
				exit(1);
			}
			angleList[i+k*numAngles] = make_int4( _get_index(angles[i].ind1,k), _get_index(angles[i].ind2,k), _get_index(angles[i].ind3,k), angles[i].tabFileIndex );
	    }
	}
	}
	
	if (numDihedrals > 0) {
	dihedralList = new int4[ (numDihedrals) * numReplicas ];
	dihedralPotList = new  int[ (numDihedrals) * numReplicas ];
	for(int k = 0 ; k < numReplicas; k++) {
	    for(int i = 0; i < numDihedrals; ++i) {
			if (dihedrals[i].tabFileIndex == -1) {
				fprintf(stderr,"Error: dihedralfile '%s' was not read with tabulatedDihedralFile command.\n", dihedrals[i].fileName.val());
				exit(1);
			}
			dihedralList[i+k*numDihedrals] = make_int4( _get_index(dihedrals[i].ind1,k), _get_index(dihedrals[i].ind2,k), _get_index(dihedrals[i].ind3,k), _get_index(dihedrals[i].ind4,k) );
		dihedralPotList[i+k*numDihedrals] = dihedrals[i].tabFileIndex;
	    }
	}
	}

	if (numBondAngles > 0) {
	bondAngleList = new int4[ (numBondAngles*2) * numReplicas ];
	for(int k = 0 ; k < numReplicas; k++) {
	    for(int i = 0; i < numBondAngles; ++i) {
			if (bondAngles[i].tabFileIndex1 == -1) {
				fprintf(stderr,"Error: bondanglefile '%s' was not read with tabulatedAngleFile command.\n", bondAngles[i].angleFileName1.val());
				exit(1);
			}
			if (bondAngles[i].tabFileIndex2 == -1) {
				fprintf(stderr,"Error: bondanglefile1 '%s' was not read with tabulatedBondFile command.\n", bondAngles[i].bondFileName.val());
				exit(1);
			}
			if (bondAngles[i].tabFileIndex3 == -1) {
				fprintf(stderr,"Error: bondanglefile2 '%s' was not read with tabulatedBondFile command.\n", bondAngles[i].angleFileName2.val());
				exit(1);
			}
			int idx = i+k*numBondAngles;
			bondAngleList[idx*2]   = make_int4( _get_index(bondAngles[i].ind1,k),
							    _get_index(bondAngles[i].ind2,k),
							    _get_index(bondAngles[i].ind3,k),
							    _get_index(bondAngles[i].ind4,k) );
			bondAngleList[idx*2+1] = make_int4( bondAngles[i].tabFileIndex1, bondAngles[i].tabFileIndex2, bondAngles[i].tabFileIndex3, -1 );
	    }
	}
	}

	// TODO: remove older copy-to-cuda
	int2* restraintList = NULL;
	if (c.numRestraints > 0) {
	    restraintList = new int2[c.numRestraints * numReplicas];
	    int j = 0;
	    for(int k = 0 ; k < numReplicas; k++)
	    {
		for(int i = 0; i < c.numRestraints; ++i)
		{
		    restraintList[j] = make_int2( _get_index(c.restraints[i].id, k), i );
		    // cout << "Displaying: bondList["<< j <<"].x = " << bondList[j].x << ".\n"
		    // << "Displaying: bondList["<< j <<"].y = " << bondList[j].y << ".\n"
		    // << "Displaying: bondList["<< j <<"].z = " << bondList[j].z << ".\n";
		    ++j;
		}
	    }
	}
	
	internal->copyBondedListsToGPU(bondList,angleList,dihedralList,dihedralPotList,bondAngleList,restraintList);
	if (c.numRestraints > 0) delete[] restraintList;
	
	forceInternal = new Vector3[(num+num_rb_attached_particles+numGroupSites)*numReplicas];
	if (fullLongRange != 0)
	    printf("No cell decomposition created.\n");

	// Prepare the trajectory output writer.
	for (int repID = 0; repID < numReplicas; ++repID) {
		TrajectoryWriter *w = new TrajectoryWriter(outFilePrefixes[repID].c_str(), TrajectoryWriter::getFormatName(outputFormat),
							   sys->getBox(), num, timestep, outputPeriod);
                
		writers.push_back(w);
	}

        //Preparing the writers for momentum if necessary Han-Yi Chou
        if(particle_dynamic == String("Langevin") || particle_dynamic == String("NoseHooverLangevin"))
        {
            for (int repID = 0; repID < numReplicas; ++repID) 
            {

                TrajectoryWriter *w = new TrajectoryWriter(outMomentumFilePrefixes[repID].c_str(), TrajectoryWriter::getFormatName(outputFormat),
                                                           sys->getBox(), num, timestep, outputPeriod);
                momentum_writers.push_back(w);
            }
        }
	updateNameList();
}

GrandBrownTown::~GrandBrownTown() {
	delete[] forceInternal;
        forceInternal = NULL;
	delete[] pos;
        pos = NULL;
        //Han-Yi Chou
        if(particle_dynamic == String("Langevin") || particle_dynamic == String("NoseHooverLangevin"))
        {
            delete[] momentum;
            momentum = NULL;
            if(particle_dynamic == String("NoseHooverLangevin"))
            {
                delete[] random;
                random = NULL;
            }
        }
        //for debug
        //delete[] force;

	delete[] type;
	delete[] serial;
	//delete randoGen;

	if (numBonds > 0)
		delete[] bondList;
	if (numAngles > 0)
		delete[] angleList;
	if (numDihedrals > 0) {
		delete[] dihedralList;
		delete[] dihedralPotList;
	}
        if(randoGen->states != NULL)
            {
                gpuErrchk(cudaFree(randoGen->states));
                randoGen->states = NULL;
            }
            if(randoGen->integer_h != NULL)
            {
                delete[] randoGen->integer_h;
                randoGen->integer_h = NULL;
            }
            if(randoGen->integer_d != NULL)
            {
                gpuErrchk(cudaFree(randoGen->integer_d));
                randoGen->integer_d = NULL;
            }
            if(randoGen->uniform_h != NULL)
            {
                delete[] randoGen->uniform_h;
                randoGen->uniform_h = NULL;
            }
            if(randoGen->uniform_d != NULL)
            {
                gpuErrchk(cudaFree(randoGen->uniform_d));
                randoGen->uniform_d = NULL;
            }
            //curandDestroyGenerator(randoGen->generator);
            delete randoGen;
            gpuErrchk(cudaFree(randoGen_d));
            for(std::vector<RigidBodyController*>::iterator iter = RBC.begin(); iter != RBC.end(); ++iter)
            {
                //(*iter)->~RigidBodyController();
                delete *iter;
            }
            RBC.clear();
	// Auxillary objects
	delete internal;
        internal = NULL;
	for (int i = 0; i < numReplicas; ++i)
        {
		delete writers[i];
                writers[i] = NULL;
        }

        if(particle_dynamic == String("Langevin") || particle_dynamic == String("NoseHooverLangevin"))
        {
            for (int i = 0; i < numReplicas; ++i)
            {
                delete momentum_writers[i];
                momentum_writers[i]=NULL;
             }
            //for (int i = 0; i < numReplicas; ++i)
                //delete force_writers[i];

        }
	//gpuErrchk(cudaFree(pos_d));
	//gpuErrchk(cudaFree(forceInternal_d));
	//gpuErrchk(cudaFree(randoGen_d));
	//gpuErrchk( cudaFree(bondList_d) );

	if (imd_on)
		delete[] imdForces;
	
		
}

//Nose Hoover is now implement for particles.
void GrandBrownTown::run()
{

    // Open the files for recording ionic currents
    for (int repID = 0; repID < numReplicas; ++repID) 
    {
        writers[repID]->newFile(pos + repID*(num+num_rb_attached_particles), name, 0.0f, num); // 'pos + (repID*num)' == array-to-pointer decay
        if(particle_dynamic == String("Langevin") || particle_dynamic == String("NoseHooverLangevin"))
            momentum_writers[repID]->newFile((momentum + repID * num), name, 0.0f, num); // 'pos + (repID*num)' == array-to-pointer decay
        //random_writers[repID]->newFile(random + (repID * num), name, 0.0f, num);
    }

    // Initialize timers (util.*)
    wkf_timerhandle timer0, timerS;
    timer0 = wkf_timer_create();
    timerS = wkf_timer_create();

    #ifdef USE_NCCL
    cudaStream_t* nccl_broadcast_streams = new cudaStream_t[gpuman.gpus.size()];
    for (int i=0; i< gpuman.gpus.size(); ++i) nccl_broadcast_streams[i] = 0;
    #endif

    copyToCUDA();

    if(particle_dynamic == String("Langevin"))
        internal -> copyToCUDA(forceInternal, pos,momentum);
    else if(particle_dynamic == String("NoseHooverLangevin"))
        internal -> copyToCUDA(forceInternal, pos, momentum, random);
    else
        internal -> copyToCUDA(forceInternal, pos);

    // IMD Stuff
    void* sock = NULL;
    void* clientsock = NULL;
    int length;
    if (imd_on) 
    {
        printf("Setting up incoming socket\n");
        vmdsock_init();
        sock = vmdsock_create();
        clientsock = NULL;
        vmdsock_bind(sock, imd_port);

        printf("Waiting for IMD connection on port %d...\n", imd_port);
        vmdsock_listen(sock);
        while (!clientsock) 
        {
            if (vmdsock_selread(sock, 0) > 0) 
            {
                clientsock = vmdsock_accept(sock);
                if (imd_handshake(clientsock))
                    clientsock = NULL;
            }
        }
        sleep(1);
        if (vmdsock_selread(clientsock, 0) != 1 || imd_recv_header(clientsock, &length) != IMD_GO) 
        {
            clientsock = NULL;
        }
        imdForces = new Vector3[num*numReplicas];
        for (size_t i = 0; i < num; ++i) // clear old forces
            imdForces[i] = Vector3(0.0f);

    } // endif (imd_on)

    // Start timers
    wkf_timer_start(timer0);
    wkf_timer_start(timerS);

    //// Happens at step 1 later anyway!
    // if (fullLongRange == 0)
    // {
    //     // cudaSetDevice(0);
    //     internal->decompose();
    //     gpuErrchk(cudaDeviceSynchronize());
    //     #ifdef _OPENMP
    //     omp_set_num_threads(4);
    //     #endif
    //     #pragma omp parallel for
    //     for(int i = 0; i < numReplicas; ++i)
    //         RBC[i]->updateParticleLists( (internal->getPos_d()[0])+i*(num+conf.num_rb_attached_particles), sys_d);
    //     gpuErrchk(cudaDeviceSynchronize());
    // }

    float t; // simulation time

    int numBlocks = ((num+num_rb_attached_particles) * numReplicas) / NUM_THREADS + (((num+num_rb_attached_particles) * numReplicas) % NUM_THREADS == 0 ? 0 : 1);
    int tl = temperatureGridFile.length();

    printf("Configuration: %d particles | %d replicas\n", num, numReplicas);
    for (int i=0; i< gpuman.gpus.size(); ++i) {
	gpuman.use(i);
	gpuErrchk( cudaProfilerStart() );
    }
    gpuman.use(0);

    //float total_energy = 0.f;
    // Main loop over Brownian dynamics steps

    float kinetic_energy = 0.0f;
    float potential_energy = 0.0f;
    
    bool get_energy = false;
    bool print_energy = false;
    
    for (long int s = 1; s < steps; s++)
    {
      PUSH_NVTX("Main loop timestep",0)
        //At the very first time step, the force is computed
        if(s == 1)
        {
            // 'interparticleForce' - determines whether particles interact with each other
	    internal->clear_force();
	    internal->clear_energy();
	    const std::vector<Vector3*>& _pos = internal->getPos_d();

	    if (num_rb_attached_particles > 0) {
		#pragma omp parallel for
		for(int i = 0; i < numReplicas; ++i) {
		    RBC[i]->update_attached_particle_positions(
			internal->getPos_d()[0]+num+i*(num+num_rb_attached_particles),
			internal->getForceInternal_d()[0]+num+i*(num+num_rb_attached_particles),
			internal->getEnergy()+num+i*(num+num_rb_attached_particles),
			sys_d, num, num_rb_attached_particles, numReplicas);
		}
	    }

	    gpuman.sync();
	    if (numGroupSites > 0) updateGroupSites<<<(numGroupSites*numReplicas/32+1),32>>>(_pos[0], groupSiteData_d, num + num_rb_attached_particles, numGroupSites, numReplicas);
	    gpuman.sync();

	    #ifdef USE_NCCL
	    if (gpuman.gpus.size() > 1) {
		gpuman.nccl_broadcast(0, _pos, _pos, (num+num_rb_attached_particles+numGroupSites)*numReplicas, -1);
		gpuman.sync();
	    }
	    #endif

            #ifdef _OPENMP
            omp_set_num_threads(4);
            #endif
            #pragma omp parallel for
            for(int i = 0; i < numReplicas; ++i)
                RBC[i]->clearForceAndTorque(); //Han-Yi Chou
            
            if (interparticleForce)
            {
                if (tabulatedPotential)
                {
                    switch (fullLongRange)
                    {
                        case 0: // [ N*log(N) ] interactions, + cutoff | decomposition
                            {
                                // cudaSetDevice(0);
                                 internal -> decompose();
                                #ifdef _OPENMP
                                omp_set_num_threads(4);
                                #endif
                                #pragma omp parallel for
                                for(int i = 0; i < numReplicas; ++i)
                                    RBC[i]->updateParticleLists( (internal->getPos_d()[0])+i*(num+num_rb_attached_particles), sys_d);
                            }
                            internal -> computeTabulated(get_energy);
                            break;
                        default:
                            internal->computeTabulatedFull(get_energy);
                            break;
                    }
                }
                else
                {
                    // Not using tabulated potentials.
                    switch (fullLongRange)
                    {
                        case 0: // Use cutoff | cell decomposition.
                            if (s % decompPeriod == 0)
                            {
                                // cudaSetDevice(0);
                                internal->decompose();
                                #ifdef _OPENMP
                                omp_set_num_threads(4);
                                #endif
                                #pragma omp parallel for
                                for(int i = 0; i < numReplicas; ++i)
                                    RBC[i]->updateParticleLists( (internal->getPos_d()[0])+i*(num+num_rb_attached_particles), sys_d);
                            }
                            internal->compute(get_energy);
                            break;

                        case 1: // Do not use cutoff
                            internal->computeFull(get_energy);
                            break;

                        case 2: // Compute only softcore forces.
                            internal->computeSoftcoreFull(get_energy);
                            break;

                        case 3: // Compute only electrostatic forces.
                            internal->computeElecFull(get_energy);
                            break;
                    }
                }
            }//if inter-particle force

	    if (num_rb_attached_particles > 0) {
		if (get_energy) {
		    compute_position_dependent_force_for_rb_attached_particles
			<<< numBlocks, NUM_THREADS >>> (
			    internal -> getPos_d()[0],
			    internal -> getForceInternal_d()[0], internal -> getEnergy(),
			    internal -> getType_d(), part_d, electricField, num, num_rb_attached_particles, numReplicas, ParticleInterpolationType);
		} else {
		    compute_position_dependent_force_for_rb_attached_particles
			<<< numBlocks, NUM_THREADS >>> (
			    internal -> getPos_d()[0], internal -> getForceInternal_d()[0],
			    internal -> getType_d(), part_d, electricField, num, num_rb_attached_particles, numReplicas, ParticleInterpolationType);
		}
	    }

            #ifdef _OPENMP
            omp_set_num_threads(4);
            #endif
            #pragma omp parallel for
            for(int i = 0; i < numReplicas; ++i)
                RBC[i]->updateForces(internal->getPos_d()[0]+i*(num+num_rb_attached_particles),
				     internal->getForceInternal_d()[0]+i*(num+num_rb_attached_particles),
				     s,
				     internal->getEnergy()+i*(num+num_rb_attached_particles),
				     get_energy,
				     RigidBodyInterpolationType, sys, sys_d, num, num_rb_attached_particles);
            if(rigidbody_dynamic == String("Langevin"))
            {
                #ifdef _OPENMP
                omp_set_num_threads(4);
                #endif
                #pragma omp parallel for
                for(int i = 0; i < numReplicas; ++i)
                {
                    RBC[i]->SetRandomTorques();
                    RBC[i]->AddLangevin();
                }
            }
	    #ifdef USE_NCCL
	    if (gpuman.gpus.size() > 1) {
		const std::vector<Vector3*>& _f = internal->getForceInternal_d();
		gpuman.nccl_reduce(0, _f, _f, (num+num_rb_attached_particles+numGroupSites)*numReplicas, -1);
	    }
	    #endif

	    if (numGroupSites > 0) distributeGroupSiteForces<false><<<(numGroupSites*numReplicas/32+1),32>>>(internal->getForceInternal_d()[0], groupSiteData_d, num+num_rb_attached_particles, numGroupSites, numReplicas);

        }//if step == 1

	gpuman.sync();

	PUSH_NVTX("Integrate particles",2)
        if(particle_dynamic == String("Langevin"))
	    if (s == 1) {
		updateKernelBAOAB<true><<< numBlocks, NUM_THREADS >>>(internal->getPos_d()[0], internal->getMom_d(), internal->getForceInternal_d()[0], internal->getType_d(), part_d, kT, kTGrid_d, electricField, tl, timestep, num, num_rb_attached_particles, sys_d, randoGen_d, numReplicas, ParticleInterpolationType);
	    } else {
		updateKernelBAOAB<false><<< numBlocks, NUM_THREADS >>>(internal->getPos_d()[0], internal->getMom_d(), internal->getForceInternal_d()[0], internal->getType_d(), part_d, kT, kTGrid_d, electricField, tl, timestep, num, num_rb_attached_particles, sys_d, randoGen_d, numReplicas, ParticleInterpolationType);
	    }
        else if(particle_dynamic == String("NoseHooverLangevin"))
            //kernel for Nose-Hoover Langevin dynamic
	    if (s == 1) {
		updateKernelNoseHooverLangevin<true><<< numBlocks, NUM_THREADS >>>(internal -> getPos_d()[0], internal -> getMom_d(), internal -> getRan_d(), internal -> getForceInternal_d()[0], internal -> getType_d(), part_d, kT, kTGrid_d, electricField, tl, timestep, num, num_rb_attached_particles, sys_d, randoGen_d, numReplicas, ParticleInterpolationType);
	    } else {
		updateKernelNoseHooverLangevin<false><<< numBlocks, NUM_THREADS >>>(internal -> getPos_d()[0], internal -> getMom_d(), internal -> getRan_d(), internal -> getForceInternal_d()[0], internal -> getType_d(), part_d, kT, kTGrid_d, electricField, tl, timestep, num, num_rb_attached_particles, sys_d, randoGen_d, numReplicas, ParticleInterpolationType);
	    }
        ////For Brownian motion
        else
            updateKernel<<< numBlocks, NUM_THREADS >>>(internal -> getPos_d()[0], internal -> getForceInternal_d()[0], internal -> getType_d(),
                                                       part_d, kT, kTGrid_d, electricField, tl, timestep, num, num_rb_attached_particles, sys_d, randoGen_d, numReplicas,
                                                       internal->getEnergy(), get_energy, ParticleInterpolationType);

	POP_NVTX

	PUSH_NVTX("Integrate rigid bodies",2)
        if(rigidbody_dynamic == String("Langevin"))
        {
            #ifdef _OPENMP
            omp_set_num_threads(4);
            #endif
            #pragma omp parallel for
            for(int i = 0; i < numReplicas; ++i)
            {
                RBC[i]->integrateDLM(sys, 0);
                RBC[i]->integrateDLM(sys, 1);
            }
        }
        else
	{
            #ifdef _OPENMP
            omp_set_num_threads(4);
            #endif
            #pragma omp parallel for ordered
            for(int i = 0; i < numReplicas; ++i)
            {
                RBC[i]->integrate(sys, s);
                #pragma omp ordered
                RBC[i]->print(s);
            }
        }
	POP_NVTX

	PUSH_NVTX("Update rigid body attached particle positions",3)
	if (num_rb_attached_particles > 0) {
	    #pragma omp parallel for
	    for(int i = 0; i < numReplicas; ++i) {
		RBC[i]->update_attached_particle_positions(
		    internal->getPos_d()[0]+num+i*(num+num_rb_attached_particles),
		    internal->getForceInternal_d()[0]+num+i*(num+num_rb_attached_particles),
		    internal->getEnergy()+num+i*(num+num_rb_attached_particles),
		    sys_d, num, num_rb_attached_particles, numReplicas);
	    }
	}
	POP_NVTX

        if (s % outputPeriod == 0) {
	    PUSH_NVTX("Copy particle positions to host for output",7)
            // Copy particle positions back to CPU
	    gpuErrchk(cudaDeviceSynchronize());
            gpuErrchk(cudaMemcpy(pos, internal ->getPos_d()[0], sizeof(Vector3) * (num+num_rb_attached_particles) * numReplicas, cudaMemcpyDeviceToHost));
	    POP_NVTX
	}
        if (imd_on && clientsock && s % outputPeriod == 0)
        {
	    assert(gpuman.gpus.size()==1); // TODO: implement IMD with multiple gpus
	    gpuErrchk(cudaDeviceSynchronize());
            float* coords = new float[num*3]; // TODO: move allocation out of run loop
            int* atomIds = new int[num]; // TODO: move allocation out of run loop
            int length;

            bool paused = false;
            while (vmdsock_selread(clientsock, 0) > 0 || paused)
            {
                switch (imd_recv_header(clientsock, &length))
                {
                        case IMD_DISCONNECT:
                            printf("[IMD] Disconnecting...\n");
                            imd_disconnect(clientsock);
                            clientsock = NULL;
                            sleep(5);
                            break;
                        case IMD_KILL:
                            printf("[IMD] Killing...\n");
                            imd_disconnect(clientsock);
                            clientsock = NULL;
                            steps = s; // Stop the simulation at this step
                            sleep(5);
                            break;
                        case IMD_PAUSE:
                            paused = !paused;
                            break;
                        case IMD_GO:
                            printf("[IMD] Caught IMD_GO\n");
                            break;
                        case IMD_MDCOMM:
                            for (size_t i = 0; i < num; ++i) // clear old forces
                                imdForces[i] = Vector3(0.0f);

                            if (imd_recv_mdcomm(clientsock, length, atomIds, coords))
                            {
                                printf("[IMD] Error receiving forces\n");
                            }
                            else
                            {
                                for (size_t j = 0; j < length; ++j)
                                {
                                    int i = atomIds[j];
                                    imdForces[i] = Vector3(coords[j*3], coords[j*3+1], coords[j*3+2]) * conf.imdForceScale;
                                }
                            }
                            break;
                        default:
                            printf("[IMD] Something weird happened. Disconnecting..\n");
                            break;
                }
            }
            if (clientsock)
            {
                    // float* coords = new float[num*3]; // TODO: move allocation out of run loop
                    for (size_t i = 0; i < num; i++)
                    {
                        const Vector3& p = pos[i];
                        coords[3*i] = p.x;
                        coords[3*i+1] = p.y;
                        coords[3*i+2] = p.z;
                    }
                    imd_send_fcoords(clientsock, num, coords);
            }
                delete[] coords;
                delete[] atomIds;
        }

        #ifdef _OPENMP
        omp_set_num_threads(4);
        #endif

	PUSH_NVTX("Clear rigid body forces",2)
        #pragma omp parallel for
        for(int i = 0; i < numReplicas; ++i) 
            RBC[i]->clearForceAndTorque();
	POP_NVTX

	
	if (numGroupSites > 0) {
 	  PUSH_NVTX("Update collective coordinates",2)
	    gpuman.sync();
	    updateGroupSites<<<(numGroupSites*numReplicas/32+1),32>>>(internal->getPos_d()[0], groupSiteData_d, num + num_rb_attached_particles, numGroupSites, numReplicas);
	    gpuman.sync();
	  POP_NVTX
	}

	if (get_energy) {
	    PUSH_NVTX("Fetch particle energy data",1)
	    thrust::device_ptr<float> en_d(internal->getEnergy());
	    potential_energy = (thrust::reduce(en_d, en_d+(num+num_rb_attached_particles)*numReplicas)) / numReplicas;
	    if(particle_dynamic == String("Langevin") || particle_dynamic == String("NoseHooverLangevin"))
	    {
		gpuErrchk(cudaMemcpy(momentum, internal->getMom_d(), sizeof(Vector3)*num*numReplicas,cudaMemcpyDeviceToHost));
		kinetic_energy = KineticEnergy();
	    }
	    gpuman.sync();
	    POP_NVTX	
	}
	
	get_energy = (((s+1) % outputEnergyPeriod) == 0) || (((s+1) % outputPeriod) == 0);
	print_energy = ((s % outputEnergyPeriod) == 0);

	if (get_energy) {
	    PUSH_NVTX("Clear particle energy data",1)
	    internal->clear_energy();
	    POP_NVTX	
	}
	    
        if (imd_on && clientsock)
            internal->setForceInternalOnDevice(imdForces); // TODO ensure replicas are mutually exclusive with IMD // TODO add multigpu support with IMD
	else {
	    PUSH_NVTX("Clear particle forces",2)
            internal->clear_force();
	    #ifdef USE_NCCL
	    if (gpuman.gpus.size() > 1) {
		const std::vector<Vector3*>& _p = internal->getPos_d();
		nccl_broadcast_streams[0] = gpuman.gpus[0].get_next_stream();
		gpuman.nccl_broadcast(0, _p, _p, (num+num_rb_attached_particles+numGroupSites)*numReplicas, nccl_broadcast_streams);
	    }
	    #endif
	    POP_NVTX
    	}

        if (interparticleForce)
        {
            // 'tabulatedPotential' - determines whether interaction is described with tabulated potentials or formulas
            if (tabulatedPotential)
            {
                switch (fullLongRange)
                {
                    case 0: // [ N*log(N) ] interactions, + cutoff | decomposition
                        if (s % decompPeriod == 0)
                        {
			  PUSH_NVTX("Decompose particles",5)
                            internal -> decompose();
			  POP_NVTX
                            #ifdef _OPENMP
                            omp_set_num_threads(4);
                            #endif
			    PUSH_NVTX("Update rigid body particle lists",6)
                            #pragma omp parallel for
                            for(int i = 0; i < numReplicas; ++i)
                                RBC[i]->updateParticleLists( (internal->getPos_d()[0])+i*(num+num_rb_attached_particles), sys_d);
			    POP_NVTX
                        }
			PUSH_NVTX("Calculate particle-particle forces",7)
                        internal -> computeTabulated(get_energy);
			POP_NVTX
			#ifdef USE_NCCL
			if (gpuman.gpus.size() > 1) {
			  PUSH_NVTX("Reduce particle forces",6)
			    const std::vector<Vector3*>& _f = internal->getForceInternal_d();
			    gpuman.nccl_reduce(0, _f, _f, (num+num_rb_attached_particles)*numReplicas, -1);
			  POP_NVTX
			}
			#endif
                        break;
                    default: // [ N^2 ] interactions, no cutoff | decompositions
                        internal->computeTabulatedFull(get_energy);
                        break;
                }
            }
            else
            { 
                // Not using tabulated potentials.
                switch (fullLongRange)
                {
                        case 0: // Use cutoff | cell decomposition.
                            if (s % decompPeriod == 0)
                            {
                               internal->decompose();
                               #ifdef _OPENMP
                               omp_set_num_threads(4);
                               #endif
                               #pragma omp parallel for
                               for(int i = 0; i < numReplicas; ++i)
                                   RBC[i]->updateParticleLists( (internal->getPos_d()[0])+i*num, sys_d);
                            }
                            internal->compute(get_energy);
                            break;
                        case 1: // Do not use cutoff
                            internal->computeFull(get_energy);
                            break;

                        case 2: // Compute only softcore forces.
                            internal->computeSoftcoreFull(get_energy);
                            break;

                        case 3: // Compute only electrostatic forces.
                            internal->computeElecFull(get_energy);
                            break;
                }
            }
        }

	PUSH_NVTX("Compute RB attached particle forces",4)
	if (get_energy) {
	    compute_position_dependent_force_for_rb_attached_particles
		<<< numBlocks, NUM_THREADS >>> (
		    internal -> getPos_d()[0],
		    internal -> getForceInternal_d()[0], internal -> getEnergy(),
		    internal -> getType_d(), part_d, electricField, num, num_rb_attached_particles, numReplicas, ParticleInterpolationType);
	} else {
	    compute_position_dependent_force_for_rb_attached_particles
		<<< numBlocks, NUM_THREADS >>> (
		    internal -> getPos_d()[0], internal -> getForceInternal_d()[0],
		    internal -> getType_d(), part_d, electricField, num, num_rb_attached_particles, numReplicas, ParticleInterpolationType);
	}
	POP_NVTX


        //compute the force for rigid bodies
        #ifdef _OPENMP
        omp_set_num_threads(4);
        #endif
	PUSH_NVTX("Compute RB-RB forces forces",5)
        #pragma omp parallel for
        for(int i = 0; i < numReplicas; ++i) // TODO: Use different buffer for RB particle forces to avoid race condition
            RBC[i]->updateForces((internal->getPos_d()[0])+i*(num+num_rb_attached_particles), (internal->getForceInternal_d()[0])+i*(num+num_rb_attached_particles), s, (internal->getEnergy())+i*(num+num_rb_attached_particles), get_energy,
				 RigidBodyInterpolationType, sys, sys_d, num, num_rb_attached_particles);
	POP_NVTX

	if (numGroupSites > 0) {
	  PUSH_NVTX("Spread collective coordinate forces to constituent particles",4)
	    gpuman.sync();
	    // if ((s%100) == 0) {
	    distributeGroupSiteForces<false><<<(numGroupSites*numReplicas/32+1),32>>>(internal->getForceInternal_d()[0], groupSiteData_d, num+num_rb_attached_particles, numGroupSites, numReplicas);
	// } else {
	//     distributeGroupSiteForces<false><<<(numGroupSites*numReplicas/32+1),32>>>(internal->getForceInternal_d()[0], groupSiteData_d, num+num_rb_attached_particles, numGroupSites, numReplicas);
	// }
	    gpuman.sync();
	  POP_NVTX
	}

	PUSH_NVTX("Update particle coordinates",2)
        if(particle_dynamic == String("Langevin") || particle_dynamic == String("NoseHooverLangevin"))
            LastUpdateKernelBAOAB<<< numBlocks, NUM_THREADS >>>(internal -> getPos_d()[0], internal -> getMom_d(), internal -> getForceInternal_d()[0], 
            internal -> getType_d(), part_d, kT, kTGrid_d, electricField, tl, timestep, num, num_rb_attached_particles, sys_d, randoGen_d, numReplicas, internal->getEnergy(), get_energy,
            ParticleInterpolationType);
            //gpuErrchk(cudaDeviceSynchronize());
	POP_NVTX
  
	PUSH_NVTX("Update RB coordinates",3)
        if(rigidbody_dynamic == String("Langevin"))
        {
            #ifdef _OPENMP
            omp_set_num_threads(4);
            #endif
            #pragma omp parallel for ordered
            for(int i = 0; i < numReplicas; ++i)
            {
                RBC[i]->SetRandomTorques();
                RBC[i]->AddLangevin();
                RBC[i]->integrateDLM(sys, 2);
                #pragma omp ordered
                RBC[i]->print(s);
            }
        }
	POP_NVTX

        if (s % outputPeriod == 0)
        {
	  PUSH_NVTX("Copy and write particle and RB coordinates for output",3)
            if(particle_dynamic == String("Langevin") || particle_dynamic == String("NoseHooverLangevin"))
            {
                gpuErrchk(cudaMemcpy(momentum, internal ->  getMom_d(), sizeof(Vector3) * (num) * numReplicas, cudaMemcpyDeviceToHost));
            }
            t = s*timestep;
            // Loop over all replicas
            for (int repID = 0; repID < numReplicas; ++repID)
            {

                if (numberFluct == 1)
                    updateNameList(); // no need for it here if particles stay the same

                // Write the trajectory.
                writers[repID]->append(pos + repID*(num+num_rb_attached_particles), name, serial, t, num);

                if(particle_dynamic == String("Langevin") || particle_dynamic == String("NoseHooverLangevin"))
                {
                    momentum_writers[repID]->append(momentum + repID * (num+num_rb_attached_particles), name, serial, t, num);
                    //force_writers[repID]->append(force + repID * (num+num_rb_attached_particles), name, serial, t, num);
                }
            }
            // TODO: Currently, not compatible with replicas. Needs a fix.
            if (numberFluct == 1)
                updateReservoirs();

           remember(t);
        }
        if (print_energy)
        {
                wkf_timer_stop(timerS);
                t = s * timestep;
                // Simulation progress and statistics.
                float percent = (100.0f * s) / steps;
                float msPerStep = wkf_timer_time(timerS) * 1000.0f / outputEnergyPeriod;
                float nsPerDay = numReplicas * timestep / msPerStep * 864E5f;

                // Nice thousand separator
                setlocale(LC_NUMERIC, "");

                // Do the output
                printf("\rStep %ld [%.2f%% complete | %.3f ms/step | %.3f ns/day]",s, percent, msPerStep, nsPerDay);

                // Copy positions from GPU to CPU.
                //gpuErrchk(cudaMemcpy(pos, internal->getPos_d(), sizeof(Vector3)*num*numReplicas,cudaMemcpyDeviceToHost));
                std::fstream energy_file;
                energy_file.open( (outFilePrefixes[0]+".energy.dat").c_str(), std::fstream::out | std::fstream::app);
                if(energy_file.is_open())
                {
                    energy_file << "Kinetic Energy: " << kinetic_energy*num*0.5f*(2.388458509e-1) << " (kT) "<< std::endl;
                    energy_file << "Potential Energy: " << potential_energy << " (kcal/mol) " << std::endl;
                    energy_file.close();
                }
                else
                {
                    std::cout << "Error in opening energy files\n";
                }
                
                if(rigidbody_dynamic == String("Langevin"))
                {
                    #ifdef _OPENMP
                    omp_set_num_threads(4);
                    #endif
                    #pragma omp parallel for
                    for(int i = 0; i < numReplicas; ++i)
                        RBC[i]->KineticEnergy();
                }
                std::fstream rb_energy_file;
                rb_energy_file.open( (outFilePrefixes[0]+".rb_energy.dat").c_str(), std::fstream::out | std::fstream::app);
                if(rb_energy_file.is_open())
                {
                    float k_tol = 0.f;
                    float v_tol = 0.f;
                    float (RigidBody::*func_ptr)();
                    #ifdef _OPENMP
                    omp_set_num_threads(4);
                    #endif
                    #pragma omp parallel for private(func_ptr) reduction(+:k_tol,v_tol)
                    for(int i = 0; i < numReplicas; ++i)
                    {
                        func_ptr = &RigidBody::getKinetic;
                        k_tol += RBC[i]->getEnergy(func_ptr);
                        func_ptr = &RigidBody::getEnergy;
                        v_tol += RBC[i]->getEnergy(func_ptr);
                    }
                    rb_energy_file << "Kinetic Energy "   << k_tol/numReplicas << " (kT)" << std::endl;
                    rb_energy_file << "Potential Energy " << v_tol/numReplicas << " (kcal/mol)" << std::endl;
                    rb_energy_file.close();
                }
                else
                {
                    std::cout << "Error in opening rb energy files\n"; 
                }

                // Write restart files for each replica.
                for (int repID = 0; repID < numReplicas; ++repID)
                    writeRestart(repID);

                wkf_timer_start(timerS);
		POP_NVTX
         } // s % outputEnergyPeriod
     POP_NVTX
     } // done with all Brownian dynamics steps

     if (imd_on and clientsock)
     {
            if (vmdsock_selread(clientsock, 0) == 1)
            {
                int length;
                switch (imd_recv_header(clientsock, &length))
                {
                    case IMD_DISCONNECT:
                        printf("\n[IMD] Disconnecting...\n");
                        imd_disconnect(clientsock);
                        clientsock = NULL;
                        sleep(5);
                        break;
                    case IMD_KILL:
                        printf("\n[IMD] Killing...\n");
                        imd_disconnect(clientsock);
                        clientsock = NULL;
                        sleep(5);
                        break;
                    default:
                        printf("\n[IMD] Something weird happened. Disconnecting..\n");
                        break;
                }
            }
     }
     // Stop the main timer.
     wkf_timer_stop(timer0);

     // Compute performance data.
     const float elapsed = wkf_timer_time(timer0); // seconds
     int tot_hrs = (int) std::fmod(elapsed / 3600.0f, 60.0f);
     int tot_min = (int) std::fmod(elapsed / 60.0f, 60.0f);
     float tot_sec   = std::fmod(elapsed, 60.0f);

     printf("\nFinal Step: %d\n", (int) steps);

     printf("Total Run Time: ");
     if (tot_hrs > 0) printf("%dh%dm%.1fs\n", tot_hrs, tot_min, tot_sec);
     else if (tot_min > 0) printf("%dm%.1fs\n", tot_min, tot_sec);
     else printf("%.2fs\n", tot_sec);

} // GrandBrownTown::run()

// --------------------------------------------
// Populate lists of types and serial numbers.
//
void GrandBrownTown::populate() {
	for (int repID = 0; repID < numReplicas; repID++) {
	    const int offset = repID * (num+num_rb_attached_particles);
		int pn = 0;
		int p = 0;
		for (int i = 0; i < num; i++) {
			type[i + offset] = p;
			serial[i + offset] = currSerial++;

			if (++pn >= part[p].num) {
				p++;
				pn = 0;
			}
		}
	}
}



void GrandBrownTown::writeRestart(int repID) const 
{
    FILE* out   = fopen(restartFiles[repID].c_str(), "w");
    const int offset = repID * (num+num_rb_attached_particles);

    for (int i = 0; i < num; ++i) 
    {
        const int ind = i + offset;
        const Vector3& p = pos[ind];
	fprintf(out, "%d %.10g %.10g %.10g\n", type[ind], p.x, p.y, p.z); 
    }
    fclose(out);

    if(particle_dynamic == String("Langevin") || particle_dynamic == String("NoseHooverLangevin"))
    {
        out = fopen(restartMomentumFiles[repID].c_str(), "w");
        
        for (int i = 0; i < num; ++i) 
        {
            const int ind = i + offset;
            const Vector3& p = momentum[ind];
            fprintf(out, "%d %.10g %.10g %.10g\n", type[ind], p.x, p.y, p.z); 
        }
        fclose(out);
    }
   
}

/*the center is defined by the first pmf*/
void GrandBrownTown::initialCondCen() {
	for (int i = 0; i < num; i++)
		pos[i] = part[ type[i] ].pmf[0]->getCenter();
}


// Set random initial positions for all particles and replicas
void GrandBrownTown::initialCond() {
	for (int repID = 0; repID < numReplicas; repID++) {
	    const int offset = repID * (num+num_rb_attached_particles);
		for (int i = 0; i < num; i++) {
			pos[i + offset] = findPos(type[i + offset]);
		}
	}
}

// A couple old routines for getting particle positions.
Vector3 GrandBrownTown::findPos(int typ) {
    // TODO: sum over grids
	Vector3 r;
	const BrownianParticleType& pt = part[typ];
	do {
		const float rx = sysDim.x * randoGen->uniform(); 
		const float ry = sysDim.y * randoGen->uniform();
		const float rz = sysDim.z * randoGen->uniform();
		r = sys->wrap( Vector3(rx, ry, rz) );
	} while (pt.pmf[0]->interpolatePotential(r) > *pt.meanPmf);
	return r;
}


Vector3 GrandBrownTown::findPos(int typ, float minZ) {
	Vector3 r;
	const BrownianParticleType& pt = part[typ];
	do {
		const float rx = sysDim.x * randoGen->uniform();
		const float ry = sysDim.y * randoGen->uniform();
		const float rz = sysDim.z * randoGen->uniform();
		r = sys->wrap( Vector3(rx, ry, rz) );
	} while (pt.pmf[0]->interpolatePotential(r) > *pt.meanPmf and fabs(r.z) > minZ);
	return r;
}

//Compute the kinetic energy of particle and rigid body Han-Yi Chou
float GrandBrownTown::KineticEnergy()
{
    float *vec_red, *energy;
    float particle_energy;

    gpuErrchk(cudaMalloc((void**)&vec_red, 512*sizeof(float)));
    gpuErrchk(cudaMalloc((void**)&energy, sizeof(float)));
    gpuErrchk(cudaMemset((void*)energy,0,sizeof(float)));

    BrownParticlesKineticEnergy<64><<<dim3(512),dim3(64)>>>(internal->getMom_d(), internal -> getType_d(), part_d, vec_red, num, num_rb_attached_particles, numReplicas);
    gpuErrchk(cudaDeviceSynchronize());

    Reduction<64><<<dim3(1),dim3(64)>>>(vec_red, energy, 512);

    gpuErrchk(cudaMemcpy((void*)&particle_energy, energy, sizeof(float), cudaMemcpyDeviceToHost));

    gpuErrchk(cudaFree(vec_red));
    gpuErrchk(cudaFree(energy));

    return 2. * particle_energy / kT / num / numReplicas; //In the unit of 0.5kT
}
/*
void GrandBrownTown::RotKineticEnergy()
{
    RBC.KineticEnergy();

    return 2. * e / numReplicas / kT; //In the unit of 0.5kT
}
*/
void GrandBrownTown::InitNoseHooverBath(int N)
{
    printf("Entering Nose-Hoover Langevin\n");
    int count = 0;

    for(int i = 0; i < N; ++i)
    {
        int typ = type[i];
        double mu = part[typ].mu;
        
        double sigma = sqrt(kT / mu);

        float tmp = sigma * randoGen->gaussian();
        random[(size_t)count] = tmp;
        ++count;
    }
    printf("Done in nose-hoover bath\n");
}

void GrandBrownTown::init_cuda_group_sites()
{
    // Count the number of particles that form groups
    int num_particles = 0;
    for (auto it = conf.groupSiteData.begin(); it != conf.groupSiteData.end(); ++it) {
	num_particles += it->size();
    }

    // Create GPU-friendly data structure
    int* tmp = new int[numGroupSites+1+num_particles];
    num_particles = 0;
    int i = 0;
    for (auto it = conf.groupSiteData.begin(); it != conf.groupSiteData.end(); ++it) {
	tmp[i] = num_particles+numGroupSites+1;
	// printf("DEBUG: tmp[%d] = %d\n", i, tmp[i]);
	for (auto it2 = it->begin(); it2 != it->end(); ++it2) {
	    tmp[num_particles+numGroupSites+1] = *it2;
	    // printf("DEBUG: tmp[%d] = %d\n", num_particles+numGroupSites+1, *it2);
	    num_particles++;
	}
	i++;
    }
    assert(i == numGroupSites);
    tmp[i] = num_particles+numGroupSites+1;
    // printf("DEBUG: tmp[%d] = %d\n", i, tmp[i]);

    // printf("DEBUG: Finally:\n");
    // for (int j = 0; j < numGroupSites+1+num_particles; j++) {
    //         printf("DEBUG: tmp[%d] = %d\n", j, tmp[j]);
    // }

    // Copy data structure to GPU
    gpuErrchk(cudaMalloc((void**) &groupSiteData_d, sizeof(int)*(numGroupSites+1+num_particles)));
    gpuErrchk(cudaMemcpy(groupSiteData_d, tmp, sizeof(int)*(numGroupSites+1+num_particles), cudaMemcpyHostToDevice));
    // TODO deallocate CUDA
    delete[] tmp;

}

// -----------------------------------------------------------------------------
// Initialize file for recording ionic current
void GrandBrownTown::newCurrent(int repID) const {
    /*
	FILE* out = fopen(outCurrFiles[repID].c_str(), "w");
	fclose(out);
    */
}


// -----------------------------------------------------------------------------
// Record the ionic current flowing through the entire system
void GrandBrownTown::writeCurrent(int repID, float t) const {
    return;
    /*
	FILE* out = fopen(outCurrFiles[repID].c_str(), "a");
	fprintf(out, "%.10g %.10g %d\n", 0.5f*(t+timeLast), current(t), num);
	fclose(out);
    */
}


// -----------------------------------------------------------------------------
// Record the ionic current in a segment -segZ < z < segZ
void GrandBrownTown::writeCurrentSegment(int repID, float t, float segZ) const {
    return;
    /*
	FILE* out = fopen(outCurrFiles[repID].c_str(), "a");
	int i;
	fprintf(out, "%.10g ", 0.5f * (t + timeLast));
	for (i = -1; i < numParts; i++)
		fprintf(out, "%.10g ", currentSegment(t,segZ,i));
	fprintf(out, "%d\n", num);
	fclose(out);
    */
}


// ----------------------------------------------------
// Compute the current in nanoamperes for entire system
//
float GrandBrownTown::current(float t) const {
	float curr = 0.0f;
	float dt = timeLast - t;

	for (int i = 0; i < num; i++) {
		Vector3 d = sys->wrapDiff(pos[i]-posLast[i]);
		curr += part[type[i]].charge*d.z/(sysDim.z*dt)*1.60217733e-1f;
	}
	return curr;
}


// -----------------------------------------------------
// Compute the current in nanoamperes for a restricted segment (-segZ < z < segZ).
//
float GrandBrownTown::currentSegment(float t, float segZ, int carrier) const {
	float curr = 0.0f;
	float dt = t - timeLast;

	for (int i = 0; i < num; i++) {
		float z0 = posLast[i].z;
		float z1 = pos[i].z;

		// Ignore carriers outside the range for both times.
		if (fabs(z0) > segZ && fabs(z1) > segZ) continue;

		// Cut the pieces outside the range.
		if (z0 < -segZ) z0 = -segZ;
		if (z1 < -segZ) z1 = -segZ;
		if (z0 > segZ) z0 = segZ;
		if (z1 > segZ) z1 = segZ;

		float dz = sys->wrapDiff(z1 - z0, sysDim.z);
		if ( carrier == type[i] || carrier == -1) {
			curr += part[type[i]].charge*dz/(2.0f*segZ*dt)*1.60217733e-1f;
		}
	}
	return curr;
}


int GrandBrownTown::getReservoirCount(int partInd, int resInd) const {
	int count = 0;
	const Reservoir* res = part[partInd].reservoir;
	for (int i = 0; i < num; ++i)
		if (type[i] == partInd and res->inside(i, pos[i]))
			count++;
	return count;
}

IndexList GrandBrownTown::getReservoirList(int partInd, int resInd) const {
	IndexList ret;
	const Reservoir* res = part[partInd].reservoir;
	for (int i = 0; i < num; ++i)
		if (type[i] == partInd and res->inside(resInd, pos[i]))
			ret.add(i);
	return ret;
}

Vector3 GrandBrownTown::freePosition(Vector3 r0, Vector3 r1, float minDist) {
	const int maxTries = 1000;
	bool tooClose = true;
	Vector3 r;
	Vector3 d = r1 - r0;
	float minDist2 = minDist*minDist;

	const CellDecomposition& decomp = internal->getDecomp();
	const CellDecomposition::cell_t *cells = decomp.getCells();

	int tries = 0;
	while (tooClose) {
		r.x = r0.x + d.x*randoGen->uniform();
		r.y = r0.y + d.y*randoGen->uniform();
		r.z = r0.z + d.z*randoGen->uniform();

		tooClose = false;
		// Check to make sure we are not too near another particle.
		const CellDecomposition::cell_t cell = decomp.getCell(decomp.getCellID(r));
		for (int i = -1; i <= 1; ++i) {
			for (int j = -1; j <= 1; ++j) {
				for (int k = -1; k <= 1; ++k) {
					int nID = decomp.getNeighborID(cell, i, j, k);
					// TODO: Determine which replica to use to look for free position.
					const CellDecomposition::range_t range = decomp.getRange(nID, 0);
					for (int n = range.first; n < range.last; ++n) {
						Vector3 dj = pos[ cells[n].particle ];
						if (dj.length2() < minDist2) {
							tooClose = true;
							break;
						}
					}
				}
			}
		}

		// Don't try too many times.
		if (++tries > maxTries) {
			printf("WARNING: freePosition too many tries to find free position.\n");
			break;
		}
	}

	return r;
}

// -----------------------------------------------------------------------------
// Update the list of particle names[] for simulations with varying number of
// particles
void GrandBrownTown::updateNameList() {
	if (outputFormat == TrajectoryWriter::formatTraj) {
		char typeNum[64];
		for (int i = 0; i < num; ++i) {
			sprintf(typeNum, "%d", type[i]);
			name[i] = typeNum;
		}
	} else {
		for (int i = 0; i < num; ++i)
			name[i] = part[ type[i] ].name;
	}
}

// -----------------------------------------------------------------------------
// Save particle positions for analysis purposes.
// TODO: Fix for multiple replicas.
void GrandBrownTown::remember(float t) {
	timeLast = t;
	std::copy(pos, pos + num * numReplicas, posLast);
        if(particle_dynamic == String("Langevin") || particle_dynamic == String("NoseHooverLangevin"))
            std::copy(momentum, momentum + num * numReplicas, momLast);
}

// -----------------------------------------------------------------------------
// Delete particles listed in the list 'p'
void GrandBrownTown::deleteParticles(IndexList& p) {
	int n = 0;
	for (int i = 0; i < num; ++i) {
		pos[n] = pos[i];
		type[n] = type[i];
		serial[n] = serial[i];
		if (p.find(i) == -1) n++;
	}
	num = n;
}


// -----------------------------------------------------------------------------
// Add particles, obey numCap limit
void GrandBrownTown::addParticles(int n, int typ) {
	if (num + n > numCap) n = numCap - num;

	for (int i = num; i < num + n; i++) {
		pos[i] = findPos(typ, initialZ);
		type[i] = typ;
		serial[i] = currSerial;
		currSerial++;
	}
	num += n;
}

// -----------------------------------------------------------------------------
// Add particles randomly within a region between r0 and r1.
// TODO: Fix for CUDA.
void GrandBrownTown::addParticles(int n, int typ, Vector3 r0, Vector3 r1) {
	if (num + n > numCap) n = numCap - num;

	Vector3 d = r1 - r0;
	for (int i = num; i < num + n; ++i) {
		Vector3 r;
		r.x = r0.x + d.x * randoGen->uniform();
		r.y = r0.y + d.y * randoGen->uniform();
		r.z = r0.z + d.z * randoGen->uniform();

		pos[i] = r;
		type[i] = typ;
		serial[i] = currSerial++;
	}
	num += n;
}

// -----------------------------------------------------------------------------
// Add particles randomly within the region defined by r0 and r1. Maintain a
// minimum distance of minDist between particles.
// TODO: Fix for CUDA.
void GrandBrownTown::addParticles(int n, int typ, Vector3 r0, Vector3 r1, float minDist) {
	if (num + n > numCap) n = numCap - num;
	const int n0 = num;
	for (int i = n0; i < n0 + n; i++) {
		// Generate a position for the new particle.
		pos[i] = freePosition(r0, r1, minDist);
		type[i] = typ;
		num++;
		// Update the cell decomposition
		internal->updateNumber(num); /* RBTODO: unsure if type arg is ok */
	}
}

// -----------------------------------------------------------------------------
// Add or delete particles in the reservoirs. Reservoirs are not wrapped.
void GrandBrownTown::updateReservoirs() {
	bool numberChange = false;
	for (int p = 0; p < numParts; ++p) {
		if (part[p].reservoir == NULL) continue;

		const int n = part[p].reservoir->length();

		for (int res = 0; res < n; res++) {
			// Get the current number of particles in the reservoir.
			IndexList resPart = getReservoirList(p, res);
			int numberCurr = resPart.length();

			// Determine the new number for this particle from a Poisson distribution.
			float number0 = part[p].reservoir->getMeanNumber(res);
			int number = randoGen->poisson(number0);

			// If the number is the same nothing needs to be done.
			if (number == numberCurr) continue;

			if (number < numberCurr) {
				int dn = numberCurr - number;

				// We need to delete particles.  Choose them at random.
				IndexList delPart;
				int pick = static_cast<int>(randoGen->uniform() *numberCurr) % numberCurr;
				if (pick + dn >= numberCurr) {
					int dn0 = dn - (numberCurr - pick);
					delPart = resPart.range(pick, numberCurr-1);
					delPart.add(resPart.range(0, dn0-1));
				} else {
					delPart = resPart.range(pick, pick + dn-1);
				}

				deleteParticles(delPart);
				numberChange = true;
			} else {
				// We need to add particles.
				Vector3 r0 = part[p].reservoir->getOrigin(res);
				Vector3 r1 = part[p].reservoir->getDestination(res);
				addParticles(number - numberCurr, p, r0, r1, minimumSep);
				numberChange = true;
			}
		} // end reservoir loop
	} // end particle loop

	if (numberChange)
		internal->updateNumber(num);
}

void GrandBrownTown::copyRandToCUDA() {
	gpuErrchk(cudaMalloc((void**)&randoGen_d, sizeof(Random)));
        gpuErrchk(cudaMemcpy(&(randoGen_d->states), &(randoGen->states), sizeof(curandState_t*),cudaMemcpyHostToDevice));
}


// -----------------------------------------------------------------------------
// Allocate memory on GPU(s) and copy to device
void GrandBrownTown::copyToCUDA() {
	/* const size_t tot_num = num * numReplicas;
	gpuErrchk(cudaMalloc(&pos_d, sizeof(Vector3) * tot_num));
	gpuErrchk(cudaMemcpyAsync(pos_d, pos, sizeof(Vector3) * tot_num,
														cudaMemcpyHostToDevice));

	gpuErrchk(cudaMalloc(&forceInternal_d, sizeof(Vector3) * num * numReplicas));
	gpuErrchk(cudaMemcpyAsync(forceInternal_d, forceInternal, sizeof(Vector3) * tot_num,
														cudaMemcpyHostToDevice));*/

	// gpuErrchk(cudaDeviceSynchronize());
}

/*void GrandBrownTown::createBondList()
{
	size_t size = (numBonds / 2) * numReplicas * sizeof(int3);
	gpuErrchk( cudaMalloc( &bondList_d, size ) );
	gpuErrchk( cudaMemcpyAsync( bondList_d, bondList, size, cudaMemcpyHostToDevice) );

	for(int i = 0 ; i < (numBonds / 2) * numReplicas ; i++)
	{
		cout << "Displaying: bondList_d["<< i <<"].x = " << bondList[i].x << ".\n"
			<< "Displaying: bondList_d["<< i <<"].y = " << bondList[i].y << ".\n"
			<< "Displaying: bondList_d["<< i <<"].z = " << bondList[i].z << ".\n";

	}
}*/
