/* Bond, James Bond.h Copyright Justin Dufresne and Terrance Howard, 2013.
 * We prefer our code shaken, not stirred.
 */

#ifndef BOND_H
#define BOND_H

#ifdef __CUDACC__
    #define HOST __host__
    #define DEVICE __device__
#else
    #define HOST
    #define DEVICE
#endif

#include "useful.h"
#include "TabulatedPotential.h"
#include <cuda.h>

const String flags[] = { "DEFAULT", "REPLACE", "ADD" };

class Bond {
public:
	enum {
		DEFAULT = 1,
		REPLACE = 1,
		ADD = 2
	};

	Bond() : flag(DEFAULT), ind1(-1), ind2(-1) { }

	Bond(int flag, int ind1, int ind2, String fileName) :
			flag(flag),
			ind1(ind1), ind2(ind2),
			fileName(fileName) { }

	Bond(String strflag, int ind1, int ind2, String fileName);

	void print();

	String toString();

public:
	int flag;
	int ind1;
	int ind2;
	int tabFileIndex;
	String fileName;
};

#endif
