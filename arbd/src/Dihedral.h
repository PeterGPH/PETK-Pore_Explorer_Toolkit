// Dihedral.h
// Authors: Justin Dufresne and Terrance Howard, 2013

#ifndef DIHEDRAL_H
#define DIHEDRAL_H

#include "useful.h"
#include <cuda.h>

class Dihedral {
public:
	int ind1, ind2, ind3, ind4;
	int tabFileIndex;
		// This will be assigned after ComputeForce.cu loads the TabulatedDihedralPotential objects.
		// The tabFileIndex is used by ComputeForce to discern which TabDiPot this Dihedral object uses.
	String fileName;
	Dihedral() : ind1(-1), ind2(-1), ind3(-1), ind4(-1), tabFileIndex(-1) {}
	Dihedral(int ind1, int ind2, int ind3, int ind4, String fileName);
	Dihedral(const Dihedral& d);
	HOST DEVICE inline int getIndex(int index) const {
		if (index == ind1) return 1;
		if (index == ind2) return 2;
		if (index == ind3) return 3;
		if (index == ind4) return 4;
		return -1;
	}
	String toString();
	void print();
};

#endif
