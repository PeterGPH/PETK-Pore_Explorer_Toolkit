///////////////////////////////////////////////////////////////////////
// An array of positions.
// Author: Jeff Comer <jcomer2@illinois.edu>
#ifndef SCATTER_H
#define SCATTER_H

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include "useful.h"
// using namespace std;

class Scatter {
public:
	Scatter(const char* coordFile);
	Scatter(const char* coordFile, float cutTime);
	Scatter(const char* coordFile, float cutTime0, float cutTime1);

	~Scatter();

	Matrix3 topMatrix() const;
	
	Vector3 get(int i) const;
	
	int length() const;

	Vector3 minBound() const;
	Vector3 maxBound() const;

	static int countCoordinates(const char* fileName);
	static int countTrajectory(const char* fileName, float cutTime);
	static int countTrajectory(const char* fileName, float cutTime0, float cutTime1);
private:
	int n;
	Vector3* r;

	Scatter(const Scatter&){}

	// Read coordinates into a Vector array.
	void readCoordinates(const char* fileName, int num, Vector3* r);
	void readTrajectory(const char* fileName, int num, Vector3* r, float cutTime);
	void readTrajectory(const char* fileName, int num, Vector3* r, float cutTime0, float cutTime1);
};
#endif
