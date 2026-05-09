// Dihedral.cu 
// Authors: Justin Dufresne and Terrance Howard, 2013

#include "Dihedral.h"

Dihedral::Dihedral(int ind1, int ind2, int ind3, int ind4, String fileName) : 
		ind1(ind1), ind2(ind2), ind3(ind3), ind4(ind4), fileName(fileName) {}

Dihedral::Dihedral(const Dihedral& d) : ind1(d.ind1), ind2(d.ind2), ind3(d.ind3), ind4(d.ind4),
		tabFileIndex(d.tabFileIndex), fileName(d.fileName) {}

String Dihedral::toString() {
	return String("DIHEDRAL ") + ind1 + " " + ind2 + " " + ind3 + " " + ind4 + " " + fileName;
}

void Dihedral::print() {
	printf("DIHEDRAL (%d %d %d %d) %s\n", ind1, ind2, ind3, ind4, fileName.val());
}
