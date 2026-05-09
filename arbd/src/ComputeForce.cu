///////////////////////////////////////////////////////////////////////
// Brownian dynamics base class
// Author: Jeff Comer <jcomer2@illinois.edu>

#include "ComputeForce.h"
#include "ComputeForce.cuh"
#include "Configuration.h"
#include <cuda_profiler_api.h>
#include <fstream>
#include <iostream>

#ifndef gpuErrchk
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true) {
   if (code != cudaSuccess) {
      fprintf(stderr,"CUDA Error: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}
#endif 

#define gpuKernelCheck() {kernelCheck( __FILE__, __LINE__); }
inline void kernelCheck(const char* file, int line)
{
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        std::fprintf(stderr,"Error: %s in %s %d\n", cudaGetErrorString(err),file, line);
        assert(1==2);
    }
    //gpuErrchk(cudaDeviceSynchronize());
}

cudaEvent_t start, stop;

GPUManager ComputeForce::gpuman = GPUManager();

void runSort(int2 *d1, int *d2, float *key,
				int2 *scratch1, int  *scratch2, float *scratchKey,
				unsigned int count);

ComputeForce::ComputeForce(const Configuration& c, const int numReplicas = 1) :
    num(c.num), numParts(c.numParts), num_rb_attached_particles(c.num_rb_attached_particles),
    sys(c.sys), switchStart(c.switchStart),
    switchLen(c.switchLen), electricConst(c.coulombConst),
    cutoff2((c.switchLen + c.switchStart) * (c.switchLen + c.switchStart)),
    decomp(c.sys->getBox(), c.sys->getOrigin(), c.switchStart + c.switchLen + c.pairlistDistance, numReplicas),
    numBonds(c.numBonds), numTabBondFiles(c.numTabBondFiles),
    numExcludes(c.numExcludes), numAngles(c.numAngles),
    numTabAngleFiles(c.numTabAngleFiles), numDihedrals(c.numDihedrals),
    numTabDihedralFiles(c.numTabDihedralFiles), numRestraints(c.numRestraints),
    numBondAngles(c.numBondAngles), numProductPotentials(c.numProductPotentials),
    numGroupSites(c.numGroupSites),
    numReplicas(numReplicas) {

	// Grow vectors for per-gpu device pointers
	for (int i = 0; i < gpuman.gpus.size(); ++i) {
	    int s = gpuman.gpus.size();
	    sys_d	= std::vector<BaseGrid*>(s);
	    tablePot_addr = std::vector<TabulatedPotential**>(s);
	    tablePot_d	= std::vector<TabulatedPotential**>(s);
	    pairLists_d = std::vector<int2*>(s);
	    pairLists_tex = std::vector<cudaTextureObject_t>(s);
	    pairTabPotType_d = std::vector<int*>(s);
	    pairTabPotType_tex = std::vector<cudaTextureObject_t>(s);
	    numPairs_d = std::vector<int*>(s);
	    pos_d = std::vector<Vector3*>(s);
	    pos_tex = std::vector<cudaTextureObject_t>(s);
	    forceInternal_d = std::vector<Vector3*>(s);
	}

	// Allocate the parameter tables.
	decomp_d = NULL;

	pairlistdist2 = (sqrt(cutoff2) + c.pairlistDistance);
	pairlistdist2 *= pairlistdist2;

	int np2     = numParts*numParts;
	tableEps    = new float[np2];
	tableRad6   = new float[np2];
	tableAlpha  = new float[np2];

	const size_t tableSize = sizeof(float) * np2;
	gpuErrchk(cudaMalloc(&tableEps_d, tableSize));
	gpuErrchk(cudaMalloc(&tableRad6_d, tableSize));
	gpuErrchk(cudaMalloc(&tableAlpha_d, tableSize));
	for (std::size_t i = 0; i < gpuman.gpus.size(); ++i) {
	    gpuman.use(i);
	    gpuErrchk(cudaMalloc(&sys_d[i], sizeof(BaseGrid)));
	    gpuErrchk(cudaMemcpyAsync(sys_d[i], sys, sizeof(BaseGrid), cudaMemcpyHostToDevice));
	}
	gpuman.use(0);

	// Build the parameter tables.
	makeTables(c.part);

	gpuErrchk(cudaMemcpyAsync(tableAlpha_d, tableAlpha, tableSize, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpyAsync(tableEps_d, tableEps, tableSize, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpyAsync(tableRad6_d, tableRad6, tableSize, cudaMemcpyHostToDevice));

	// Create the potential table
	tablePot = new TabulatedPotential*[np2];
	for (int i = 0; i < np2; ++i) tablePot[i] = NULL;

	for (std::size_t i = 0; i < gpuman.gpus.size(); ++i) {
	    tablePot_addr[i] = new TabulatedPotential*[np2];
	    for (int j = 0; j < np2; ++j) tablePot_addr[i][j] = NULL;
	    gpuman.use(i);
	    gpuErrchk(cudaMalloc(&tablePot_d[i], sizeof(TabulatedPotential*) * np2));
	}
	gpuman.use(0);

	// Create the bond table
	tableBond = new TabulatedPotential*[numTabBondFiles];
	tableBond_addr = new TabulatedPotential*[numTabBondFiles];
	bondList_d = NULL;
	tableBond_d = NULL;
	for (int i = 0; i < numTabBondFiles; i++) {
		tableBond_addr[i] = NULL;
		tableBond[i] = NULL;
	}
	gpuErrchk(cudaMalloc(&tableBond_d, sizeof(TabulatedPotential*) * numTabBondFiles));

	// Create the angle table
	tableAngle = new TabulatedAnglePotential*[numTabAngleFiles];
	tableAngle_addr = new TabulatedAnglePotential*[numTabAngleFiles];
	angleList_d = NULL;
	tableAngle_d = NULL;
	for (int i = 0; i < numTabAngleFiles; i++) {
		tableAngle_addr[i] = NULL;
		tableAngle[i] = NULL;
	}
	gpuErrchk(cudaMalloc(&tableAngle_d, sizeof(TabulatedAnglePotential*) * numTabAngleFiles));

	// Create the dihedral table
	tableDihedral = new TabulatedDihedralPotential*[numTabDihedralFiles];
	tableDihedral_addr = new TabulatedDihedralPotential*[numTabDihedralFiles];
	dihedralList_d = NULL;
	tableDihedral_d = NULL;
	for (int i = 0; i < numTabDihedralFiles; i++) {
		tableDihedral_addr[i] = NULL;
		tableDihedral[i] = NULL;
	}
	gpuErrchk(cudaMalloc(&tableDihedral_d, sizeof(TabulatedDihedralPotential*) * numTabDihedralFiles));

	{	// allocate device for pairlists
		// RBTODO: select maxpairs in better way; add assertion in kernel to avoid going past this
		const int maxPairs = MAX_NLIST_PAIRS;
		for (std::size_t i = 0; i < gpuman.gpus.size(); ++i) {
		    gpuman.use(i);
		    gpuErrchk(cudaMalloc(&numPairs_d[i],       sizeof(int)));
		    gpuErrchk(cudaMalloc(&pairLists_d[i],      sizeof(int2)*maxPairs));
		    // gpuErrchk(cudaBindTexture(0, pairListsTex, pairLists_d[i], sizeof(int2)*maxPairs)); //Han-Yi
		    gpuErrchk(cudaMalloc(&pairTabPotType_d[i], sizeof(int)*maxPairs));
		}

		// create texture object
		for (std::size_t i = 0; i < gpuman.gpus.size(); ++i) {
		    gpuman.use(i);
		    cudaResourceDesc resDesc;
		    memset(&resDesc, 0, sizeof(resDesc));
		    resDesc.resType = cudaResourceTypeLinear;
		    resDesc.res.linear.devPtr = pairLists_d[i];
		    resDesc.res.linear.desc.f = cudaChannelFormatKindSigned;
		    resDesc.res.linear.desc.x = 32; // bits per channel
		    resDesc.res.linear.desc.y = 32; // bits per channel
		    resDesc.res.linear.sizeInBytes = maxPairs*sizeof(int2);

		    cudaTextureDesc texDesc;
		    memset(&texDesc, 0, sizeof(texDesc));
		    texDesc.readMode = cudaReadModeElementType;

		    // create texture object: we only have to do this once!
		    pairLists_tex[i]=0;
		    cudaCreateTextureObject(&pairLists_tex[i], &resDesc, &texDesc, NULL);
		}

		// create texture object
		for (std::size_t i = 0; i < gpuman.gpus.size(); ++i) {
		    gpuman.use(i);
		    cudaResourceDesc resDesc;
		    memset(&resDesc, 0, sizeof(resDesc));
		    resDesc.resType = cudaResourceTypeLinear;
		    resDesc.res.linear.devPtr = pairTabPotType_d[i];
		    resDesc.res.linear.desc.f = cudaChannelFormatKindSigned;
		    resDesc.res.linear.desc.x = 32; // bits per channel
		    resDesc.res.linear.sizeInBytes = maxPairs*sizeof(int);

		    cudaTextureDesc texDesc;
		    memset(&texDesc, 0, sizeof(texDesc));
		    texDesc.readMode = cudaReadModeElementType;

		    // create texture object: we only have to do this once!
		    pairTabPotType_tex[i] = 0;
		    cudaCreateTextureObject(&pairTabPotType_tex[i], &resDesc, &texDesc, NULL);

		}
		gpuman.use(0);


                //Han-Yi Chou
                int nCells = decomp.nCells.x * decomp.nCells.y * decomp.nCells.z;
                //int* nCells_dev;
		if (nCells < MAX_CELLS_FOR_CELLNEIGHBORLIST) {
		    int3 *Cells_dev;
		    size_t sz = 27*nCells*sizeof(int);
		    gpuErrchk(cudaMalloc(&CellNeighborsList, sz));
		    //gpuErrchk(cudaMalloc(&nCells_dev,sizeof(int)));
		    gpuErrchk(cudaMalloc(&Cells_dev,sizeof(int3)));
		    //gpuErrchk(cudaMemcpy(nCells_dev,&nCells,1,cudaMemcpyHostToDevice);
		    gpuErrchk(cudaMemcpy(Cells_dev,&(decomp.nCells),sizeof(int3),cudaMemcpyHostToDevice));
		    createNeighborsList<<<256,256>>>(Cells_dev,CellNeighborsList);
		    gpuErrchk(cudaFree(Cells_dev));

		    // create texture object
		    {
			cudaResourceDesc resDesc;
			memset(&resDesc, 0, sizeof(resDesc));
			resDesc.resType = cudaResourceTypeLinear;
			resDesc.res.linear.devPtr = CellNeighborsList;
			resDesc.res.linear.desc.f = cudaChannelFormatKindSigned;
			resDesc.res.linear.desc.x = 32; // bits per channel
			resDesc.res.linear.sizeInBytes = sz;

			cudaTextureDesc texDesc;
			memset(&texDesc, 0, sizeof(texDesc));
			texDesc.readMode = cudaReadModeElementType;

			// create texture object: we only have to do this once!
			neighbors_tex=0;
			cudaCreateTextureObject(&neighbors_tex, &resDesc, &texDesc, NULL);
		    }
		}
	}
	
	restraintIds_d = NULL;
	bondAngleList_d = NULL;
	product_potential_list_d = NULL;

	//Calculate the number of blocks the grid should contain
	gridSize =  (num+num_rb_attached_particles) / NUM_THREADS + 1;

	// Create and allocate the energy arrays
	gpuErrchk(cudaMalloc(&energies_d, sizeof(float) * (num+num_rb_attached_particles+numGroupSites) * numReplicas));
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
}

ComputeForce::~ComputeForce() {
	delete[] tableEps;
	delete[] tableRad6;
	delete[] tableAlpha;
	gpuErrchk(cudaFree(tableEps_d));
	gpuErrchk(cudaFree(tableAlpha_d));
	gpuErrchk(cudaFree(tableRad6_d));
	
	for (int i = 0; i < numParts; ++i) {
	    for (int j = i; j < numParts; ++j) {
		int ind = i+j*numParts;
		if (tablePot[ind] != NULL) {
		    for (std::size_t g = 0; g < gpuman.gpus.size(); ++g) {
			gpuman.use(g);
			tablePot_addr[g][ind]->free_from_cuda(tablePot_addr[g][ind]);
			tablePot_addr[g][ind] = NULL;
		    }
		    delete tablePot[ind];
		}
	    }
	}
	delete[] tablePot;
	for (auto& tpa : tablePot_addr) delete[] tpa;

	for (int j = 0; j < numTabBondFiles; ++j)
		delete tableBond[j];
	delete[] tableBond;
	delete[] tableBond_addr;
	gpuErrchk(cudaFree(tableBond_d));

	for (int j = 0; j < numTabAngleFiles; ++j)
		if (tableAngle[j] != NULL)
			delete tableAngle[j];
	delete[] tableAngle;
	delete[] tableAngle_addr;

	if(type_d != NULL)
	{
		gpuErrchk(cudaFree(tableAngle_d));

		gpuErrchk(cudaFree(energies_d));

		gpuErrchk( cudaFree(type_d) );
		if (numBonds > 0) {
			gpuErrchk( cudaFree(bonds_d) );
			gpuErrchk( cudaFree(bondMap_d) );
			gpuErrchk( cudaFree(bondList_d) );
		}
		if (numAngles > 0) {
			gpuErrchk( cudaFree(angles_d) );
			gpuErrchk( cudaFree(angleList_d) );
		}
		if (numDihedrals > 0) {
			gpuErrchk( cudaFree(dihedrals_d) );
			gpuErrchk( cudaFree(dihedralList_d) );
			gpuErrchk( cudaFree(dihedralPotList_d) );
		}
		if (numExcludes > 0) {
			gpuErrchk( cudaFree(excludes_d) );
			gpuErrchk( cudaFree(excludeMap_d) );
		}
		if (numRestraints > 0) {
			gpuErrchk( cudaFree(restraintIds_d) );
			gpuErrchk( cudaFree(restraintLocs_d) );
			gpuErrchk( cudaFree(restraintSprings_d) );
		}
	}

	for (std::size_t i = 0; i < gpuman.gpus.size(); ++i) {
	    gpuErrchk(cudaFree(forceInternal_d[i]) );
	    gpuErrchk(cudaFree(sys_d[i]));
	    gpuErrchk(cudaDestroyTextureObject(pos_tex[i]));
	    gpuErrchk(cudaFree(pos_d[i]) );
	    gpuErrchk(cudaFree(numPairs_d[i]));
	    gpuErrchk(cudaDestroyTextureObject(pairLists_tex[i]));
	    gpuErrchk(cudaFree(pairLists_d[i]));
	    gpuErrchk(cudaDestroyTextureObject(pairTabPotType_tex[i]));
	    gpuErrchk(cudaFree(pairTabPotType_d[i]));
	}
        gpuErrchk(cudaDestroyTextureObject(neighbors_tex));
        gpuErrchk(cudaFree( CellNeighborsList));

}

void ComputeForce::updateNumber(int newNum) {
	if (newNum == num or newNum < 0) return;

	// Set the new number.
	num = newNum;

	// Reallocate the neighbor list.
	//delete[] neigh;
	//neigh = new IndexList[num];
	decompose();

	printf("updateNumber() called\n");
	// Reallocate CUDA arrays

	// Recalculate the number of blocks in the grid
	gridSize = 0;
	while ((int)sqrt(NUM_THREADS) * gridSize < num+num_rb_attached_particles)
		++gridSize;

	gpuErrchk(cudaFree(energies_d));
	gpuErrchk(cudaMalloc(&energies_d, sizeof(float) * gridSize));
}

void ComputeForce::makeTables(const BrownianParticleType part[]) {
	for (int i = 0; i < numParts; ++i) {
		const BrownianParticleType& pi = part[i];
		for (int j = 0; j < numParts; ++j) {
			const BrownianParticleType& pj = part[j];
			int ind = i * numParts + j;
			tableEps[ind] = sqrtf(pi.eps * pj.eps);
			float r = pi.radius + pj.radius;
			tableRad6[ind] = r * r * r * r * r * r;
			tableAlpha[ind] = electricConst * pi.charge * pj.charge;
		}
	}
}

bool ComputeForce::addTabulatedPotential(String fileName, int type0, int type1) {
	if (type0 < 0 or type0 >= numParts) return false;
	if (type1 < 0 or type1 >= numParts) return false;

	int ind = type0 + type1 * numParts;
	int ind1 = type1 + type0 * numParts;

	// If an entry already exists for this particle type, delete it
	if (tablePot[ind] != NULL) {
	    for (std::size_t i = 0; i < gpuman.gpus.size(); ++i) {
		gpuman.use(i);
		tablePot_addr[i][ind]->free_from_cuda(tablePot_addr[i][ind]);
		tablePot_addr[i][ind] = NULL;
	    }
	    gpuman.use(0);
	    delete tablePot[ind];
	}
	// if (tablePot[ind1] != NULL) {
	//     // gpuErrchk(cudaFree(tablePot_addr[ind1]));
	// 	delete tablePot[ind1];
	// 	// tablePot[ind1] = NULL;
	// 	// tablePot_addr[ind1] = NULL;
	// }

	tablePot[ind] = tablePot[ind1] = new TabulatedPotential(fileName);
	tablePot[ind]->truncate(switchStart, sqrtf(cutoff2), 0.0f);

	for (std::size_t i = 0; i < gpuman.gpus.size(); ++i) {
	    gpuman.use(i);
	    tablePot_addr[i][ind] = tablePot_addr[i][ind1] = tablePot[ind]->copy_to_cuda();
	    gpuErrchk(cudaMemcpy(tablePot_d[i], tablePot_addr[i],
				 sizeof(TabulatedPotential*) * numParts * numParts, cudaMemcpyHostToDevice));
	}
	gpuman.use(0);
	return true;
}

bool ComputeForce::addBondPotential(String fileName, int ind, Bond bonds[], BondAngle bondAngles[])
{
    // TODO: see if tableBond_addr can be removed
    if (tableBond[ind] != NULL) {
	delete tableBond[ind];
	// gpuErrchk(cudaFree(tableBond_addr[ind])); //TODO free this a little more cleanly
    }

    tableBond[ind] = new TabulatedPotential(fileName);

	for (int i = 0; i < numBonds; ++i)
		if (bonds[i].fileName == fileName)
			bonds[i].tabFileIndex = ind;

	for (int i = 0; i < numBondAngles; i++)
	{
	    if (bondAngles[i].bondFileName == fileName)
		bondAngles[i].tabFileIndex2 = ind;
	}

	gpuErrchk(cudaMemcpyAsync(bonds_d, bonds, sizeof(Bond) * numBonds, cudaMemcpyHostToDevice));

	tableBond_addr[ind] = tableBond[ind]->copy_to_cuda();
	gpuErrchk(cudaMemcpy(tableBond_d, tableBond_addr,
			     sizeof(TabulatedPotential*) * numTabBondFiles, cudaMemcpyHostToDevice));
	return true;
}

bool ComputeForce::addAnglePotential(String fileName, int ind, Angle* angles, BondAngle* bondAngles) {
	if (tableAngle[ind] != NULL) {
		delete tableAngle[ind];
		gpuErrchk(cudaFree(tableAngle_addr[ind]));
		tableAngle[ind] = NULL;
		tableAngle_addr[ind] = NULL;
	}

	tableAngle[ind] = new TabulatedAnglePotential(fileName);
	TabulatedAnglePotential *t = new TabulatedAnglePotential(*tableAngle[ind]);

	// Copy tableAngle[ind] to the device
	float *pot;
	int size = tableAngle[ind]->size;
	gpuErrchk(cudaMalloc(&pot, sizeof(float) * size));
	gpuErrchk(cudaMemcpyAsync(pot, tableAngle[ind]->pot, sizeof(float) * size, cudaMemcpyHostToDevice));
	t->pot = pot;
	gpuErrchk(cudaMalloc(&tableAngle_addr[ind], sizeof(TabulatedAnglePotential)));
	gpuErrchk(cudaMemcpy(tableAngle_addr[ind], t, sizeof(TabulatedAnglePotential), cudaMemcpyHostToDevice));
	t->pot = NULL;
	delete t;

	gpuErrchk(cudaMemcpyAsync(tableAngle_d, tableAngle_addr,
			sizeof(TabulatedAnglePotential*) * numTabAngleFiles, cudaMemcpyHostToDevice));

	for (int i = 0; i < numAngles; i++)
		if (angles[i].fileName == fileName)
			angles[i].tabFileIndex = ind;

	for (int i = 0; i < numBondAngles; i++) {
	    if (bondAngles[i].angleFileName1 == fileName)
		bondAngles[i].tabFileIndex1 = ind;
	    if (bondAngles[i].angleFileName2 == fileName)
		bondAngles[i].tabFileIndex3 = ind;
	}
	gpuErrchk(cudaMemcpy(angles_d, angles, sizeof(Angle) * numAngles,
			cudaMemcpyHostToDevice));
	return true;
}

bool ComputeForce::addDihedralPotential(String fileName, int ind, Dihedral dihedrals[])
{
	for (int i = 0; i < numDihedrals; i++)
		if (dihedrals[i].fileName == fileName)
			dihedrals[i].tabFileIndex = ind;

	gpuErrchk(cudaMemcpyAsync(dihedrals_d, dihedrals, sizeof(Dihedral) * numDihedrals,
			cudaMemcpyHostToDevice));

	if (tableDihedral[ind] != NULL) {
		delete tableDihedral[ind];
		gpuErrchk(cudaFree(tableDihedral_addr[ind]));
		tableDihedral[ind] = NULL;
		tableDihedral_addr[ind] = NULL;
	}

	tableDihedral[ind] = new TabulatedDihedralPotential(fileName);
	TabulatedDihedralPotential t = TabulatedDihedralPotential(*tableDihedral[ind]);

	// Copy tableAngle[ind] to the device
	float *pot;
	int size = tableDihedral[ind]->size;
	gpuErrchk(cudaMalloc(&pot, sizeof(float) * size));
	gpuErrchk(cudaMemcpyAsync(pot, tableDihedral[ind]->pot,
			sizeof(float) * size, cudaMemcpyHostToDevice));
	t.pot = pot;

	gpuErrchk(cudaMalloc(&tableDihedral_addr[ind], sizeof(TabulatedDihedralPotential)));
	gpuErrchk(cudaMemcpyAsync(tableDihedral_addr[ind], &t,
			sizeof(TabulatedDihedralPotential), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(tableDihedral_d, tableDihedral_addr,
			sizeof(TabulatedDihedralPotential*) * numTabDihedralFiles, cudaMemcpyHostToDevice));
	t.pot = NULL;
	return true;
}

void ComputeForce::decompose() {
	//gpuErrchk( cudaProfilerStart() );

	// Reset the cell decomposition.
	if (decomp_d != NULL)
        {
            cudaFree(decomp_d);
            decomp_d = NULL;
	}	
	decomp.decompose_d(pos_d[0], num+num_rb_attached_particles);
	decomp_d = decomp.copyToCUDA();

	// Update pairlists using cell decomposition (not sure this is really needed or good) 
	//RBTODO updatePairlists<<< nBlocks, NUM_THREADS >>>(pos_d[0], num, numReplicas, sys_d[0], decomp_d);	

	/* size_t free, total; */
	/* { */
	/* 	cuMemGetInfo(&free,&total); */
	/* 	printf("Free memory: %zu / %zu\n", free, total); */
	/* } */
	
	// initializePairlistArrays
	int nCells = decomp.nCells.x * decomp.nCells.y * decomp.nCells.z;
	
	/* cuMemGetInfo(&free,&total); */
	/* printf("Free memory: %zu / %zu\n", free, total); */
	
	int tmp = 0;
	gpuErrchk(cudaMemcpyAsync(numPairs_d[0], &tmp,	sizeof(int), cudaMemcpyHostToDevice));
	gpuErrchk(cudaDeviceSynchronize());

#ifdef DEBUGEXCLUSIONS
	initExSum();
	gpuErrchk(cudaDeviceSynchronize()); /* RBTODO: sync needed here? */
#endif
      
      #if __CUDA_ARCH__ >= 520
      createPairlists<64,64,8><<<dim3(128,128,numReplicas),dim3(64,1,1)>>>(pos_d[0], num+num_rb_attached_particles, numReplicas, sys_d[0], decomp_d, nCells, numPairs_d[0],
                                                                             pairLists_d[0], numParts, type_d, pairTabPotType_d[0], excludes_d,
									   excludeMap_d, numExcludes, pairlistdist2, pos_tex[0], neighbors_tex);
      #else //__CUDA_ARCH__ == 300
      createPairlists<64,64,8><<<dim3(256,256,numReplicas),dim3(64,1,1)>>>(pos_d[0], num+num_rb_attached_particles, numReplicas, sys_d[0], decomp_d, nCells, numPairs_d[0],
                                                                           pairLists_d[0], numParts, type_d, pairTabPotType_d[0], excludes_d, 
                                                                           excludeMap_d, numExcludes, pairlistdist2, pos_tex[0], neighbors_tex);
      #endif
       
      gpuKernelCheck();
      gpuErrchk(cudaDeviceSynchronize()); /* RBTODO: sync needed here? */

      #ifdef USE_NCCL
      if (gpuman.gpus.size() > 1) {
	  // Currently we don't use numPairs_d[i] for i > 0... might be able to reduce data transfer with some kind nccl scatter, and in that case we'd prefer to use all numPairs_d[i]
	  gpuErrchk(cudaMemcpy(&numPairs, numPairs_d[0], sizeof(int), cudaMemcpyDeviceToHost));
	  gpuman.nccl_broadcast(0, pairTabPotType_d, pairTabPotType_d, numPairs, -1);
	  gpuman.nccl_broadcast(0, pairLists_d, pairLists_d, numPairs, -1);
      }
      gpuman.sync();
      #endif
}

IndexList ComputeForce::decompDim() const {
	IndexList ret;
	ret.add(decomp.getNx());
	ret.add(decomp.getNy());
	ret.add(decomp.getNz());
	return ret;
}

CellDecomposition ComputeForce::getDecomp() { return decomp; }

float ComputeForce::decompCutoff() { return decomp.getCutoff(); }

float ComputeForce::computeFull(bool get_energy) {
	float energy = 0.0f;
	gridSize = ((num+num_rb_attached_particles) * numReplicas) / NUM_THREADS + 1;
	dim3 numBlocks(gridSize, 1, 1);
	dim3 numThreads(NUM_THREADS, 1, 1);

	// Call the kernel to calculate forces
	computeFullKernel<<< numBlocks, numThreads >>>(forceInternal_d[0], pos_d[0], type_d, tableAlpha_d,
		tableEps_d, tableRad6_d, num+num_rb_attached_particles, numParts, sys_d[0], energies_d, gridSize,
		numReplicas, get_energy);

	// Calculate energy based on the array created by the kernel
	if (get_energy) {
		gpuErrchk(cudaDeviceSynchronize());
		thrust::device_ptr<float> en_d(energies_d);
		energy = thrust::reduce(en_d, en_d + num + num_rb_attached_particles + numGroupSites);
	}

	return energy;
}

float ComputeForce::computeSoftcoreFull(bool get_energy) {
	float energy = 0.0f;
	gridSize = ((num+num_rb_attached_particles) * numReplicas) / NUM_THREADS + 1;
	dim3 numBlocks(gridSize, 1, 1);
	dim3 numThreads(NUM_THREADS, 1, 1);

	// Call the kernel to calculate forces
	computeSoftcoreFullKernel<<<numBlocks, numThreads>>>(forceInternal_d[0], pos_d[0], type_d,
			tableEps_d, tableRad6_d, num+num_rb_attached_particles, numParts, sys_d[0], energies_d, gridSize,
			numReplicas, get_energy);

	// Calculate energy based on the array created by the kernel
	if (get_energy) {
		cudaDeviceSynchronize();
		thrust::device_ptr<float> en_d(energies_d);
		energy = thrust::reduce(en_d, en_d + num + num_rb_attached_particles);
	}

	return energy;
}

float ComputeForce::computeElecFull(bool get_energy) {
	float energy = 0.0f;

	gridSize = num/NUM_THREADS + 1;
	dim3 numBlocks(gridSize, 1, 1);
	dim3 numThreads(NUM_THREADS, 1, 1);

	// Call the kernel to calculate forces
	computeElecFullKernel<<<numBlocks, numThreads>>>(forceInternal_d[0], pos_d[0], type_d,
			tableAlpha_d, num, numParts, sys_d[0], energies_d, gridSize, numReplicas,
			get_energy);

	// Calculate energy based on the array created by the kernel
	if (get_energy) {
		gpuErrchk(cudaDeviceSynchronize());
		thrust::device_ptr<float> en_d(energies_d);
		energy = thrust::reduce(en_d, en_d + num);
	}

	return energy;
}


float ComputeForce::compute(bool get_energy) {
	float energy = 0.0f;

	gridSize = ((num+num_rb_attached_particles) * numReplicas) / NUM_THREADS + 1;
	dim3 numBlocks(gridSize, 1, 1);
	dim3 numThreads(NUM_THREADS, 1, 1);

	// Call the kernel to calculate forces
	computeKernel<<<numBlocks, numThreads>>>(forceInternal_d[0], pos_d[0], type_d,
			tableAlpha_d, tableEps_d, tableRad6_d, num+num_rb_attached_particles, numParts, sys_d[0],
			decomp_d, energies_d, switchStart, switchLen, gridSize, numReplicas,
			get_energy);

	gpuErrchk(cudaDeviceSynchronize());
	// Calculate the energy based on the array created by the kernel
	if (get_energy) {
		gpuErrchk(cudaDeviceSynchronize());
		thrust::device_ptr<float> en_d(energies_d);
		energy = thrust::reduce(en_d, en_d + num + num_rb_attached_particles + numGroupSites);
	}

	return energy;
}

//MLog: added Bond* bondList to the list of passed in variables.
/*float ComputeForce::computeTabulated(Vector3* force, Vector3* pos, int* type,
		Bond* bonds, int2* bondMap, Exclude* excludes, int2* excludeMap,
		Angle* angles, Dihedral* dihedrals, bool get_energy, Bond* bondList) {*/
float ComputeForce::computeTabulated(bool get_energy) {
	float energy = 0.0f;

	gridSize = ((num+num_rb_attached_particles) * numReplicas) / NUM_THREADS + 1;
	dim3 numBlocks(gridSize, 1, 1);
	dim3 numThreads(NUM_THREADS, 1, 1);
	
	// Call the kernel to calculate the forces
	// int nb = (decomp.nCells.x * decomp.nCells.y * decomp.nCells.z);
	// int nb = (1+(decomp.nCells.x * decomp.nCells.y * decomp.nCells.z)) * 75; /* RBTODO: number of pairLists */
	const int nb = 800;
	// printf("ComputeTabulated\n");

	// RBTODO: get_energy
	if (get_energy)
	//if (false) 
	{
		//clearEnergies<<< nb, numThreads >>>(energies_d,num);
		//gpuErrchk(cudaDeviceSynchronize());
	        cudaMemset((void*)energies_d, 0, sizeof(float)*(num+num_rb_attached_particles+numGroupSites)*numReplicas);
		computeTabulatedEnergyKernel<<< nb, numThreads >>>(forceInternal_d[0], pos_d[0], sys_d[0],
						cutoff2, numPairs_d[0], pairLists_d[0], pairTabPotType_d[0], tablePot_d[0], energies_d);
	}
	
	else
	{
	    // Copy positions from device 0 to all others

                //gpuErrchk(cudaBindTexture(0,  PosTex, pos_d[0],sizeof(Vector3)*num*numReplicas));
		//computeTabulatedKernel<<< nb, numThreads >>>(forceInternal_d[0], pos_d[0], sys_d[0],

	    int ngpu = gpuman.gpus.size();
	    if (ngpu == 1) {
		int i = 0;
		computeTabulatedKernel<64><<< dim3(2048,1,1), dim3(64,1,1), 0, gpuman.gpus[i].get_next_stream() >>>
		    (forceInternal_d[i], sys_d[i], cutoff2, numPairs_d[i], pairLists_d[i], pairTabPotType_d[i], tablePot_d[i], pairLists_tex[i], pos_tex[i], pairTabPotType_tex[i]);

	    } else {
	    for (size_t i = 0; i < ngpu; ++i) {
		gpuman.use(i);
		int start =            floor( ((float) numPairs*i    )/ngpu );
		int end   = i < ngpu-1 ? floor( ((float) numPairs*(i+1))/ngpu ) : numPairs;
		
		if (i == ngpu-1) assert(end == numPairs);
		computeTabulatedKernel<64><<< dim3(2048,1,1), dim3(64,1,1), 0, gpuman.gpus[i].get_next_stream() >>>(forceInternal_d[i], sys_d[i],
														    cutoff2, pairLists_d[i], pairTabPotType_d[i], tablePot_d[i], pairLists_tex[i], pos_tex[i], pairTabPotType_tex[i], start, end-start);
                  gpuKernelCheck();
	    }
	    gpuman.use(0);
	    }
                //gpuErrchk(cudaUnbindTexture(PosTex));
	}
	/* printPairForceCounter<<<1,32>>>(); */

	//Mlog: the commented function doesn't use bondList, uncomment for testing.
	//if(bondMap_d != NULL && tableBond_d != NULL)

	if(product_potential_list_d != NULL && product_potentials_d != NULL)
	{
	    computeProductPotentials <<<nb, numThreads, 0, gpuman.get_next_stream()>>> ( forceInternal_d[0], pos_d[0], sys_d[0], numReplicas*numProductPotentials, product_potential_particles_d, product_potentials_d, product_potential_list_d, productCount_d, energies_d, get_energy);
	}

	if(bondAngleList_d != NULL && tableBond_d != NULL && tableAngle_d != NULL)
	{
	    computeTabulatedBondAngles <<<nb, numThreads, 0, gpuman.get_next_stream()>>> ( forceInternal_d[0], pos_d[0], sys_d[0], numReplicas*numBondAngles, bondAngleList_d, tableAngle_d, tableBond_d, energies_d, get_energy);
	}

	if(bondList_d != NULL && tableBond_d != NULL)

	{
	    //computeTabulatedBonds <<<numBlocks, numThreads>>> ( force, pos, num, numParts, sys_d[0], bonds, bondMap_d, numBonds, numReplicas, energies_d, get_energy, tableBond_d);
	//computeTabulatedBonds <<<nb, numThreads>>> ( forceInternal_d[0], pos_d[0], sys_d[0], numReplicas*numBonds/2, bondList_d, tableBond_d);
	  //if(get_energy)
              //cudaMemset(bond_energy_d, 0, sizeof(float)*num);
		computeTabulatedBonds <<<nb, numThreads, 0, gpuman.get_next_stream()>>> ( forceInternal_d[0], pos_d[0], sys_d[0], numReplicas*numBonds/2, bondList_d, tableBond_d, energies_d, get_energy);
	}

	if (angleList_d != NULL && tableAngle_d != NULL)
        {
            //if(get_energy)
		//computeTabulatedAngles<<<nb, numThreads>>>(forceInternal_d[0], pos_d[0], sys_d[0], numAngles*numReplicas, angleList_d, tableAngle_d);
	    computeTabulatedAngles<<<nb, numThreads, 0, gpuman.get_next_stream()>>>(forceInternal_d[0], pos_d[0], sys_d[0], numAngles*numReplicas, angleList_d, tableAngle_d, energies_d, get_energy);
        }
	if (dihedralList_d != NULL && tableDihedral_d != NULL)
        {
            //if(get_energy)
		//computeTabulatedDihedrals<<<nb, numThreads>>>(forceInternal_d[0], pos_d[0], sys_d[0], numDihedrals*numReplicas, dihedralList_d, dihedralPotList_d, tableDihedral_d);
	    computeTabulatedDihedrals<<<nb, numThreads, 0, gpuman.get_next_stream()>>>(forceInternal_d[0], pos_d[0], sys_d[0], numDihedrals*numReplicas, 
                dihedralList_d, dihedralPotList_d, tableDihedral_d, energies_d, get_energy);
        }

	// TODO: Sum energy
	if (restraintIds_d != NULL )
	    computeHarmonicRestraints<<<1, numThreads, 0, gpuman.get_next_stream()>>>(forceInternal_d[0], pos_d[0], sys_d[0], numRestraints*numReplicas, restraintList_d, restraintLocs_d, restraintSprings_d);
	

	// Calculate the energy based on the array created by the kernel
	// TODO: return energy
	/*if (get_energy) 
        {
            float e = 0.f;
	    gpuErrchk(cudaDeviceSynchronize());
	    thrust::device_ptr<float> en_d(energies_d);
	    e = (thrust::reduce(en_d, en_d+num*numReplicas)) / numReplicas;
            std::fstream energy_file;
            energy_file.open("energy_config.txt", std::fstream::out | std::fstream::app);
            if(energy_file.is_open())
            {
                energy_file << "Configuation Energy: "  << e << " kcal/mol " << std::endl;
                energy_file.close();
            }
            else
            {
                std::cout << "Error in opening energ files\n";
            }
            energy = e;
        }*/
	return energy;
}

float ComputeForce::computeTabulatedFull(bool get_energy) {
	energy = 0.0f;

	gridSize = ((num+num_rb_attached_particles) * numReplicas) / NUM_THREADS + 1;
	dim3 numBlocks(gridSize, 1, 1);
	dim3 numThreads(NUM_THREADS, 1, 1);

	// Call the kernel to calculate forces
	computeTabulatedFullKernel<<< numBlocks, numThreads >>>(forceInternal_d[0], pos_d[0], type_d,	tablePot_d[0], tableBond_d, num+num_rb_attached_particles, numParts, sys_d[0], bonds_d, bondMap_d, numBonds, excludes_d, excludeMap_d, numExcludes, energies_d, gridSize, numReplicas, get_energy, angles_d);
	gpuErrchk(cudaDeviceSynchronize());

	computeAngles<<< numBlocks, numThreads >>>(forceInternal_d[0], pos_d[0], angles_d, tableAngle_d,
																						 numAngles, num+num_rb_attached_particles, sys_d[0], energies_d,
																						 get_energy);
	gpuErrchk(cudaDeviceSynchronize());
	computeDihedrals<<< numBlocks, numThreads >>>(forceInternal_d[0], pos_d[0], dihedrals_d,
																							  tableDihedral_d, numDihedrals,
																								num+num_rb_attached_particles, sys_d[0], energies_d,
																								get_energy);
	// Calculate the energy based on the array created by the kernel
	if (get_energy) {
		gpuErrchk(cudaDeviceSynchronize());
		thrust::device_ptr<float> en_d(energies_d);
		energy = thrust::reduce(en_d, en_d + num + num_rb_attached_particles);
	}

	return energy;
}

void ComputeForce::copyToCUDA(Vector3* forceInternal, Vector3* pos)
{
    const size_t tot_num = (num+num_rb_attached_particles+numGroupSites) * numReplicas;

	for (std::size_t i = 0; i < gpuman.gpus.size(); ++i) {
	    gpuman.use(i);
	    gpuErrchk(cudaMalloc(&pos_d[i], sizeof(Vector3) * tot_num));
	    //Han-Yi bind to the texture
	    cudaResourceDesc resDesc;
	    memset(&resDesc, 0, sizeof(resDesc));
	    resDesc.resType = cudaResourceTypeLinear;
	    resDesc.res.linear.devPtr = pos_d[i];
	    resDesc.res.linear.desc.f = cudaChannelFormatKindFloat;
	    resDesc.res.linear.desc.x = 32; // bits per channel
	    resDesc.res.linear.desc.y = 32; // bits per channel
	    resDesc.res.linear.desc.z = 32; // bits per channel
	    resDesc.res.linear.desc.w = 32; // bits per channel
	    resDesc.res.linear.sizeInBytes = tot_num*sizeof(float4);
	    
	    cudaTextureDesc texDesc;
	    memset(&texDesc, 0, sizeof(texDesc));
	    texDesc.readMode = cudaReadModeElementType;
	    
	    // create texture object: we only have to do this once!
	    pos_tex[i] = 0;
	    cudaCreateTextureObject(&pos_tex[i], &resDesc, &texDesc, NULL);
	    gpuErrchk(cudaDeviceSynchronize());
	}
	gpuman.use(0);

	gpuErrchk(cudaMemcpyAsync(pos_d[0], pos, sizeof(Vector3) * tot_num, cudaMemcpyHostToDevice));

	for (std::size_t i = 0; i < gpuman.gpus.size(); ++i) {
	    gpuman.use(i);
	    gpuErrchk(cudaMalloc(&forceInternal_d[i], sizeof(Vector3) * tot_num));
	}
	gpuman.use(0);
	gpuErrchk(cudaMemcpyAsync(forceInternal_d[0], forceInternal, sizeof(Vector3) * tot_num, cudaMemcpyHostToDevice));

	gpuErrchk(cudaDeviceSynchronize());
}
void ComputeForce::copyToCUDA(Vector3* forceInternal, Vector3* pos, Vector3* mom)
{
    const size_t tot_num = num * numReplicas;

        gpuErrchk(cudaMalloc(&mom_d, sizeof(Vector3) * tot_num));
        gpuErrchk(cudaMemcpyAsync(mom_d, mom, sizeof(Vector3) * tot_num, cudaMemcpyHostToDevice));

	copyToCUDA(forceInternal,pos);
        gpuErrchk(cudaDeviceSynchronize());
}
void ComputeForce::copyToCUDA(Vector3* forceInternal, Vector3* pos, Vector3* mom, float* random)
{
    const size_t tot_num = num * numReplicas;

        gpuErrchk(cudaMalloc(&ran_d, sizeof(float) * tot_num));
        gpuErrchk(cudaMemcpyAsync(ran_d, random, sizeof(float) * tot_num, cudaMemcpyHostToDevice));

	copyToCUDA(forceInternal, pos, mom);
        gpuErrchk(cudaDeviceSynchronize());
}

void ComputeForce::setForceInternalOnDevice(Vector3* f) {
    // const size_t tot_num = (num+numGroupSites) * numReplicas;
    assert(numGroupSites == 0); // IMD, the only feature using this function, is currently incompatible with group sites
    const size_t tot_num = num * numReplicas;
	gpuErrchk(cudaMemcpy(forceInternal_d[0], f, sizeof(Vector3) * tot_num, cudaMemcpyHostToDevice));
}

void ComputeForce::copyToCUDA(int simNum, int *type, Bond* bonds, int2* bondMap, Exclude* excludes, int2* excludeMap, Angle* angles, Dihedral* dihedrals, const Restraint* const restraints, const BondAngle* const bondAngles, const XpotMap simple_potential_map, const std::vector<SimplePotential> simple_potentials, const ProductPotentialConf* const product_potential_confs)
{
    assert(simNum == numReplicas); // Not sure why we have both of these things
    int tot_num_with_rb = (num+num_rb_attached_particles) * simNum;
    int tot_num_with_rb_group = (num+num_rb_attached_particles+numGroupSites) * simNum;
	// type_d
	gpuErrchk(cudaMalloc(&type_d, sizeof(int) * tot_num_with_rb));
	gpuErrchk(cudaMemcpyAsync(type_d, type, sizeof(int) * tot_num_with_rb, cudaMemcpyHostToDevice));
	
	if (numBonds > 0)
	{
		// bonds_d
		gpuErrchk(cudaMalloc(&bonds_d, sizeof(Bond) * numBonds));
		gpuErrchk(cudaMemcpyAsync(bonds_d, bonds, sizeof(Bond) * numBonds, cudaMemcpyHostToDevice));
		
		// bondMap_d
		gpuErrchk(cudaMalloc(&bondMap_d, sizeof(int2) * tot_num_with_rb_group));
		gpuErrchk(cudaMemcpyAsync(bondMap_d, bondMap, sizeof(int2) * tot_num_with_rb_group, cudaMemcpyHostToDevice));
	}

	if (numExcludes > 0) {
	    // printf("Copying %d exclusions to the GPU\n", numExcludes);
	    
		// excludes_d
		gpuErrchk(cudaMalloc(&excludes_d, sizeof(Exclude) * numExcludes));
		gpuErrchk(cudaMemcpyAsync(excludes_d, excludes, sizeof(Exclude) * numExcludes,
				cudaMemcpyHostToDevice));
		
		// excludeMap_d
		gpuErrchk(cudaMalloc(&excludeMap_d, sizeof(int2) * tot_num_with_rb));
		gpuErrchk(cudaMemcpyAsync(excludeMap_d, excludeMap, sizeof(int2) * tot_num_with_rb,
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
	}

	if (numRestraints > 0) {
	    int restraintIds[numRestraints];
	    Vector3 restraintLocs[numRestraints];
	    float restraintSprings[numRestraints];
	    for (int i = 0; i < numRestraints; ++i) {
		restraintIds[i]     = restraints[i].id;
		restraintLocs[i]    = restraints[i].r0;
		restraintSprings[i] = restraints[i].k;
	    }

	    gpuErrchk(cudaMalloc(&restraintIds_d, sizeof(int) * numRestraints));
	    gpuErrchk(cudaMalloc(&restraintLocs_d, sizeof(Vector3) * numRestraints));
	    gpuErrchk(cudaMalloc(&restraintSprings_d, sizeof(float) * numRestraints));
	    
	    gpuErrchk(cudaMemcpyAsync(restraintIds_d, restraintIds,
				      sizeof(int)     * numRestraints, cudaMemcpyHostToDevice));
	    gpuErrchk(cudaMemcpyAsync(restraintLocs_d, restraintLocs,
				      sizeof(Vector3) * numRestraints, cudaMemcpyHostToDevice));
	    gpuErrchk(cudaMemcpyAsync(restraintSprings_d, restraintSprings,
				      sizeof(float)   * numRestraints, cudaMemcpyHostToDevice));
	}	    

	if (numBondAngles > 0) {
		gpuErrchk(cudaMalloc(&bondAngles_d, sizeof(BondAngle) * numBondAngles));
		gpuErrchk(cudaMemcpyAsync(bondAngles_d, bondAngles, sizeof(BondAngle) * numBondAngles,
				cudaMemcpyHostToDevice));
	}

	if (simple_potentials.size() > 0) {
	    float **val = simple_potential_pots_d = new float*[simple_potentials.size()];
	    // float **tmp = new float*[simple_potentials.size()];
	    for (int i=0; i < simple_potentials.size(); ++i) {
		const SimplePotential sp = simple_potentials[i];
		gpuErrchk(cudaMalloc(&val[i], sizeof(float)*sp.size));
		gpuErrchk(cudaMemcpyAsync(val[i], sp.pot, sizeof(float)*sp.size, cudaMemcpyHostToDevice));
		// tmp[i] = sp.pot;
		// // sp.pot = val[i];
	    }

	    // size_t sz =  sizeof(SimplePotential) * simple_potentials.size();
	    // gpuErrchk(cudaMalloc(&simple_potentials_d, sz));
	    // gpuErrchk(cudaMemcpyAsync(simple_potentials_d, &simple_potentials[0], sz,
	    // 				  cudaMemcpyHostToDevice));
	    
	    // for (int i=0; i < simple_potentials.size(); ++i) { // Restore host pointers on host object
	    // 	SimplePotential &sp = simple_potentials[i];
	    // 	sp.pot = tmp[i];
	    // }
	    // // delete[] val;
	    // delete[] tmp;

	}
	
	if (numProductPotentials > 0) {
	    // Count particles
	    int n_pots = 0;
	    int n_particles = 0;
	    for (int i=0; i < numProductPotentials; ++i) {
		const ProductPotentialConf& c = product_potential_confs[i];
		n_pots += c.indices.size();
		for (int j=0; j < c.indices.size(); ++j) {
		    n_particles += c.indices[j].size();
		}
	    }
	    // printf("DEBUG: Found %d particles participating in %d potentials forming %d productPotentials\n",
	    // 	   n_particles, n_pots, numProductPotentials);

	    // Build productPotentialLists on host
	    int *particle_list = new int[n_particles*numReplicas];
	    SimplePotential *product_potentials = new SimplePotential[n_pots];
	    uint2 *product_potential_list = new uint2[numProductPotentials*numReplicas];
	    unsigned short *productCount = new unsigned short[numProductPotentials*numReplicas];

	    n_particles = 0;
	    
	    for (unsigned int r=0; r < numReplicas; ++r) {
		n_pots = 0;
		for (int i=0; i < numProductPotentials; ++i) {
		    const ProductPotentialConf& c = product_potential_confs[i];
		    product_potential_list[i+r*numProductPotentials] = make_uint2( n_pots, n_particles );

		    for (int j=0; j < c.indices.size(); ++j) {
			if (r == 0) {
			    const unsigned int sp_i = simple_potential_map.find( std::string( c.potential_names[j].val() ) )->second;
			    product_potentials[n_pots] = simple_potentials[sp_i];
			    product_potentials[n_pots].pot = simple_potential_pots_d[sp_i];
			}
			++n_pots;
			for (int k=0; k < c.indices[j].size(); ++k) {
			    particle_list[n_particles++] = c.indices[j][k]+r*num;
			}
		    }
		    productCount[i+r*numProductPotentials] = c.indices.size();
		}
	    }

	    // Copy to device
	    size_t sz = n_particles*numReplicas * sizeof(int);
	    gpuErrchk(cudaMalloc(&product_potential_particles_d, sz));
	    gpuErrchk(cudaMemcpyAsync(product_potential_particles_d, particle_list, sz,
	    				  cudaMemcpyHostToDevice));
	    sz = n_pots * sizeof(SimplePotential);
	    gpuErrchk(cudaMalloc(&product_potentials_d, sz));
	    gpuErrchk(cudaMemcpyAsync(product_potentials_d, product_potentials, sz,
	    				  cudaMemcpyHostToDevice));
	    sz = numProductPotentials*numReplicas * sizeof(uint2);
	    gpuErrchk(cudaMalloc(&product_potential_list_d, sz));
	    gpuErrchk(cudaMemcpyAsync(product_potential_list_d, product_potential_list, sz,
	    				  cudaMemcpyHostToDevice));
	    sz = numProductPotentials*numReplicas * sizeof(unsigned short);
	    gpuErrchk(cudaMalloc(&productCount_d, sz));
	    gpuErrchk(cudaMemcpyAsync(productCount_d, productCount, sz,
	    				  cudaMemcpyHostToDevice));

	    // Clean up
	    delete[] particle_list;
	    delete[] product_potentials;
	    delete[] product_potential_list;
	    delete[] productCount;
	}

	gpuErrchk(cudaDeviceSynchronize());
}

// void ComputeForce::createBondList(int3 *bondList)
// {
// 	size_t size = (numBonds / 2) * numReplicas * sizeof(int3);
// 	gpuErrchk( cudaMalloc( &bondList_d, size ) );
// 	gpuErrchk( cudaMemcpyAsync( bondList_d, bondList, size, cudaMemcpyHostToDevice) );

// 	for(int i = 0 ; i < (numBonds / 2) * numReplicas ; i++)
// 	{
// 		cout << "Displaying: bondList_d["<< i <<"].x = " << bondList[i].x << ".\n"
// 			<< "Displaying: bondList_d["<< i <<"].y = " << bondList[i].y << ".\n"
// 			<< "Displaying: bondList_d["<< i <<"].z = " << bondList[i].z << ".\n";

// 	}
// }

void ComputeForce::copyBondedListsToGPU(int3 *bondList, int4 *angleList, int4 *dihedralList, int *dihedralPotList, int4* bondAngleList, int2* restraintList) {
	size_t size;

	if (numBonds > 0) {
	size = (numBonds / 2) * numReplicas * sizeof(int3);
	gpuErrchk( cudaMalloc( &bondList_d, size ) );
	gpuErrchk( cudaMemcpyAsync( bondList_d, bondList, size, cudaMemcpyHostToDevice) );
	}
	
	if (numAngles > 0) {
    size = numAngles * numReplicas * sizeof(int4);
    gpuErrchk( cudaMalloc( &angleList_d, size ) );
    gpuErrchk( cudaMemcpyAsync( angleList_d, angleList, size, cudaMemcpyHostToDevice) );
	}
	
	if (numDihedrals > 0) {
    size = numDihedrals * numReplicas * sizeof(int4);
    gpuErrchk( cudaMalloc( &dihedralList_d, size ) );
    gpuErrchk( cudaMemcpyAsync( dihedralList_d, dihedralList, size, cudaMemcpyHostToDevice) );

    size = numDihedrals * numReplicas * sizeof(int);
    gpuErrchk( cudaMalloc( &dihedralPotList_d, size ) );
    gpuErrchk( cudaMemcpyAsync( dihedralPotList_d, dihedralPotList, size, cudaMemcpyHostToDevice) );
	}

	if (numBondAngles > 0) {
	    size = 2*numBondAngles * numReplicas * sizeof(int4);
	    gpuErrchk( cudaMalloc( &bondAngleList_d, size ) );
	    gpuErrchk( cudaMemcpyAsync( bondAngleList_d, bondAngleList, size, cudaMemcpyHostToDevice) );
	}

	if (numRestraints > 0) {
	    size = numRestraints * numReplicas * sizeof(int2);
	    gpuErrchk( cudaMalloc( &restraintList_d, size ) );
	    gpuErrchk( cudaMemcpyAsync( restraintList_d, restraintList, size, cudaMemcpyHostToDevice) );
	}

}
