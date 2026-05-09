#include "BrownianParticleType.h"

//////////////////////////////////////////
// BrownianParticleType Implementations //
//////////////////////////////////////////
void BrownianParticleType::clear() {
	if (pmf != NULL) delete [] pmf;
	if (pmf_scale != NULL) delete [] pmf_scale;
	if (diffusionGrid != NULL) delete diffusionGrid;
	if (forceXGrid != NULL) delete forceXGrid;
	if (forceYGrid != NULL) delete forceYGrid;
	if (forceZGrid != NULL) delete forceZGrid;
	if (reservoir != NULL) delete reservoir;
        if (meanPmf != NULL) delete []  meanPmf;
	pmf = NULL, diffusionGrid = NULL;
	pmf_scale = NULL;
	pmf_boundary_conditions = NULL;
	forceXGrid = NULL, forceYGrid = NULL, forceZGrid = NULL;
	reservoir = NULL, meanPmf = NULL;
}

void BrownianParticleType::copy(const BrownianParticleType& src) {
	name = src.name;
	num = src.num;
	diffusion = src.diffusion;
        mass      = src.mass;
	charge = src.charge;
	radius = src.radius;
	eps = src.eps;
        pmf = src.pmf;
        pmf_scale = src.pmf_scale;
        pmf_boundary_conditions = src.pmf_boundary_conditions;
	meanPmf = src.meanPmf;
        numPartGridFiles = src.numPartGridFiles;
        //Han-Yi Chou
        transDamping = src.transDamping;
        mu = src.mu;
        diffusionGrid = NULL;
        forceXGrid = NULL, forceYGrid = NULL, forceZGrid = NULL;
        reservoir = NULL;
        //if (src.pmf != NULL) pmf = new BaseGrid(*src.pmf);
        if (src.diffusionGrid != NULL) diffusionGrid = new BaseGrid(*src.diffusionGrid);
        if (src.forceXGrid != NULL) forceXGrid = new BaseGrid(*src.forceXGrid);
        if (src.forceYGrid != NULL) forceYGrid = new BaseGrid(*src.forceYGrid);
        if (src.forceZGrid != NULL) forceZGrid = new BaseGrid(*src.forceZGrid);
        if (src.reservoir != NULL) reservoir = new Reservoir(*src.reservoir);

        /*
	pmf = NULL, diffusionGrid = NULL;
	forceXGrid = NULL, forceYGrid = NULL, forceZGrid = NULL;
	reservoir = NULL;
	if (src.pmf != NULL) pmf = new BaseGrid(*src.pmf);
	if (src.diffusionGrid != NULL) diffusionGrid = new BaseGrid(*src.diffusionGrid);
	if (src.forceXGrid != NULL) forceXGrid = new BaseGrid(*src.forceXGrid);
	if (src.forceYGrid != NULL) forceYGrid = new BaseGrid(*src.forceYGrid);
	if (src.forceZGrid != NULL) forceZGrid = new BaseGrid(*src.forceZGrid);
	if (src.reservoir != NULL) reservoir = new Reservoir(*src.reservoir);*/
}

BrownianParticleType& BrownianParticleType::operator=(const BrownianParticleType& src) {
        if(&src != this)
        {
	    clear();
	    copy(src);
        }
	return *this;
}
/*
bool BrownianParticleType::crop(int x0, int y0, int z0,
																int x1, int y1, int z1, bool keep_origin) {
	bool success = true;
	
	// Try cropping
	BaseGrid *new_pmf(NULL), *new_diffusionGrid(NULL);
	BaseGrid *new_forceXGrid(NULL), *new_forceYGrid(NULL), *new_forceZGrid(NULL);
	if (pmf != NULL) {
		new_pmf = new BaseGrid(*pmf);
		success = new_pmf->crop(x0, y0, z0, x1, y1, z1, keep_origin);
	}
	if (success && diffusionGrid != NULL) {
		new_diffusionGrid = new BaseGrid(*diffusionGrid);
		success = new_diffusionGrid->crop(x0, y0, z0, x1, y1, z1, keep_origin);
	}
	if (success && forceXGrid != NULL) {
		new_forceXGrid = new BaseGrid(*forceXGrid);
		success = new_forceXGrid->crop(x0, y0, z0, x1, y1, z1, keep_origin);
	}
	if (success && forceYGrid != NULL) {
		new_forceYGrid = new BaseGrid(*forceYGrid);
		success = new_forceYGrid->crop(x0, y0, z0, x1, y1, z1, keep_origin); 
	}
	if (success && forceZGrid != NULL) {
		new_forceZGrid = new BaseGrid(*forceZGrid);
		success = new_forceZGrid->crop(x0, y0, z0, x1, y1, z1, keep_origin); 
	}
	
	// Save results
	if (success) {
		if (pmf != NULL) {
			delete pmf;
			pmf = new_pmf;
		}
		if (diffusionGrid != NULL) {
			delete diffusionGrid;
			diffusionGrid = new_diffusionGrid;
		}
		if (forceXGrid != NULL) {
			delete forceXGrid;
			forceXGrid = new_forceXGrid;
		}
		if (forceYGrid != NULL) {
			delete forceYGrid;
			forceYGrid = new_forceYGrid;
		}
		if (forceZGrid != NULL) {
			delete forceZGrid;
			forceZGrid = new_forceZGrid;
		}
	} else {
		if (new_pmf != NULL) delete new_pmf;
		if (new_diffusionGrid != NULL) delete new_diffusionGrid;
		if (new_forceXGrid != NULL) delete new_forceXGrid;
		if (new_forceYGrid != NULL) delete new_forceYGrid;
		if (new_forceZGrid != NULL) delete new_forceZGrid;
	}
		
	return success;
}
*/
/*
///////////////////////////////////////
// TypeDecomposition Implementations //
///////////////////////////////////////
TypeDecomposition::TypeDecomposition(const CellDecomposition &decomp,
		const BrownianParticleType *parts, size_t num_parts) :
		num_cells_(decomp.size()), num_parts_(num_parts), parts_(decomp.size()) {
	int cutoff = (int) decomp.getCutoff();
	
	for (size_t c = 0; c < num_cells_; c++) {
		parts_[c] = new BrownianParticleType[num_parts_];
		int3 pos = decomp.getCellPos(c);
	
		printf("pos[%lu] (%d, %d, %d)\n", c, pos.x, pos.y, pos.z);
		int x0 = cutoff * pos.x;
		int y0 = cutoff * pos.y;
		int z0 = cutoff * pos.z;
		for (size_t type = 0; type < num_parts; type++) {
			parts_[c][type] = BrownianParticleType(parts[type]);
			bool success = parts_[c][type].crop(x0, y0, z0,
					x0 + cutoff, y0 + cutoff, z0 + cutoff, false);
			if (!success)
				printf("WARNING: parts[%lu][%lu] was not cropped, %s %d\n",
						c, type, __FILE__, __LINE__);
		}
	}
}

TypeDecomposition::~TypeDecomposition() {
	for (size_t c = 0; c < num_cells_; c++) {
		delete[] parts_[c];
	}
}


const BrownianParticleType* TypeDecomposition::at(size_t i) const {
	if (i >= num_cells_) {
		printf("ERROR: out of bounds [%lu] %s %d\n", i, __FILE__, __LINE__);
		return NULL;
	}
	return parts_[i];
}
*/
