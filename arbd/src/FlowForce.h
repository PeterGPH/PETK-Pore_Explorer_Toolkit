///////////////////////////////////////////////////////////////////////
// Author: Jeff Comer <jcomer2@illinois.edu>
#ifndef FLOWFORCE_H
#define FLOWFORCE_H

#include <cmath>
#include "useful.h"
// using namespace std;

class FlowForce {
public:
	FlowForce(float v);

	Vector3 force(Vector3 r, float diffusion) const;

private:
	float chanHalfLen;
	float chanHalfWidth;
	float chanVel0;
	float buffVel;
};

#endif
