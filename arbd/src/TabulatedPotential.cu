///////////////////////////////////////////////////////////////////////
// Author: Jeff Comer <jcomer2@illinois.edu>

#include "useful.h"
#include "TabulatedPotential.h"
#include <cuda.h>

TabulatedPotential::TabulatedPotential() {
    n = 2;
    drInv = 1.0f;
    r0 = 0.0f;
    v0 = new float[n];
    for (int i = 0; i < n; i++) v0[i] = 0.0f;
}

TabulatedPotential::TabulatedPotential(const float* dist, const float* pot, int n0) {
	n = abs(n0);
	drInv = 1.0f/(dist[1]-dist[0]);
	r0 = dist[0];

	v0 = new float[n];
	for (int i = 0; i < n; i++) v0[i] = pot[i];
}
TabulatedPotential::TabulatedPotential(const TabulatedPotential& tab) {
    n = tab.n;
    drInv = tab.drInv;
    r0 = tab.r0;

    v0 = new float[n];
    for (int i = 0; i < n; i++) v0[i] = tab.v0[i];
}

TabulatedPotential::~TabulatedPotential() {
    if (v0 != NULL) delete [] v0;
}

void TabulatedPotential::truncate(float cutoff) {
	int home = int(floor((cutoff - r0) * drInv));
	if (home > n) return;

	float v = v0[home];
	for (int i = home; i < n; i++) v0[i] = v;
	// interpolate();
}

bool TabulatedPotential::truncate(float switchDist, float cutoff, float value) {
	int indOff = int(floor((cutoff - r0) * drInv));
	int indSwitch = int(floor((switchDist - r0) * drInv));

	if (indSwitch > n) return false;

	// Set everything after the cutoff to "value".
	for (int i = indOff; i < n; i++) v0[i] = value;
    
	// Apply a linear switch.
	float v = v0[indSwitch];
	float m = (value - v)/(indOff - indSwitch);
	for (int i = indSwitch; i < indOff; i++) v0[i] = m*(i - indSwitch) + v;

	// interpolate();
	return true;
}

// Vector3 TabulatedPotential::computeForce(Vector3 r) {
// 	float d = r.length();
// 	Vector3 rUnit = -r/d;
// 	int home = int(floor((d - r0)*drInv));

// 	if (home < 0) return Vector3(0.0f);
// 	if (home >= n) return Vector3(0.0f);
        
// 	float homeR = home*dr + r0;
// 	float w = (d - homeR)/dr;
   
// 	// Interpolate.
// 	Vector3 force = -(3.0f*v3[home]*w*w + 2.0f*v2[home]*w + v1[home])*rUnit/dr;
// 	return force;
// }
 
  
// void TabulatedPotential::interpolate() { for cubic interpolation
// 	v1 = new float[n];
// 	v2 = new float[n];
// 	v3 = new float[n];

// 	for (int i = 0; i < n; i++) {
// 		int i0 = i - 1;
// 		int i1 = i;
// 		int i2 = i + 1;
// 		int i3 = i + 2;

// 		if (i0 < 0) i0 = 0;
// 		if (i2 >= n) i2 = n-1;
// 		if (i3 >= n) i3 = n-1;

// 		v3[i] = 0.5f*(-v0[i0] + 3.0f*v0[i1] - 3.0f*v0[i2] + v0[i3]);
// 		v2[i] = 0.5f*(2.0f*v0[i0] - 5.0f*v0[i1] + 4.0f*v0[i2] - v0[i3]);
// 		v1[i] = 0.5f*(-v0[i0] + v0[i2]);
// 	}
// 	e0 = v3[n-1] + v2[n-1] + v1[n-1] + v0[n-1];
// }

// void TabulatedPotential::init(const float* dist, const float* pot, int n0) {
// 	n = abs(n0);
// 	dr = dist[1]-dist[0];
// 	r0 = dist[0];
// 	r1 = r0 + n*dr;

// 	v0 = new float[n];
// 	for (int i = 0; i < n; i++) v0[i] = pot[i];
// }


FullTabulatedPotential::FullTabulatedPotential(const char* fileName) : fileName(fileName) {
	// printf("File: %s\n", fileName);
	FILE* inp = fopen(fileName, "r");
	if (inp == NULL) {
		printf("TabulatedPotential:TabulatedPotential Could not open file '%s'\n", fileName);
		exit(-1);
	}
	
	char line[256];
	
	numLines = countValueLines(fileName);
	float* r = new float[numLines];
	float* v = new float[numLines];
	
	int count = 0;
	while (fgets(line, 256, inp) != NULL) {
		// Ignore comments.
		int len = strlen(line);
		if (line[0] == '#') continue;
		if (len < 2) continue;
		
		String s(line);
		int numTokens = s.tokenCount();
		if (numTokens != 2) {
			printf("TabulatedPotential:TabulatedPotential Invalid tabulated potential file line: %s\n", line);
			exit(-1);
		}
		
		String* tokenList = new String[numTokens];
		s.tokenize(tokenList);
		if (tokenList == NULL) {
			printf("TabulatedPotential:TabulatedPotential Invalid tabulated potential file line: %s\n", line);
			exit(-1);
		}
		r[count] = (float) strtod(tokenList[0], NULL);
		v[count] = (float) strtod(tokenList[1], NULL);
		count++;
		
		delete[] tokenList;
	}
	fclose(inp);
	pot = new TabulatedPotential(r,v,count);
	// init(r, v, count);
	// interpolate();
	delete[] r;
	delete[] v;
}

FullTabulatedPotential::FullTabulatedPotential(const FullTabulatedPotential& tab) {
    pot = new TabulatedPotential(*tab.pot);
    numLines = tab.numLines;
    fileName = String(tab.fileName);
}

FullTabulatedPotential::~FullTabulatedPotential() {
	delete pot;
}

int FullTabulatedPotential::countValueLines(const char* fileName) {
	FILE* inp = fopen(fileName, "r");
	if (inp == NULL) {
		printf("TabulatedPotential::countValueLines Could not open file '%s'\n", fileName);
		exit(-1);
	}
	char line[256];
	int count = 0;

	while (fgets(line, 256, inp) != NULL) {
		// Ignore comments.
		int len = strlen(line);
		if (line[0] == '#') continue;
		if (len < 2) continue;
      
		count++;
	}
	fclose(inp);

	return count;
}

int countValueLines(const char* fileName) {
	FILE* inp = fopen(fileName, "r");
	if (inp == NULL) {
		printf("SimplePotential::countValueLines Could not open file '%s'\n", fileName);
		exit(-1);
	}
	char line[256];
	int count = 0;

	while (fgets(line, 256, inp) != NULL) {
		// Ignore comments.
		int len = strlen(line);
		if (line[0] == '#') continue;
		if (len < 2) continue;
		count++;
	}
	fclose(inp);
	return count;
}
 
SimplePotential::SimplePotential(const char* filename, SimplePotentialType type) : type(type) {
	FILE* inp = fopen(filename, "r");
	if (inp == NULL) {
		printf("SimplePotential::SimplePotential Could not open file '%s'\n", filename);
		exit(-1);
	}
	
	char line[256];
	
	size = (unsigned int) countValueLines(filename);
	float* r = new float[size];
	pot = new float[size];
	
	int count = 0;
	while (fgets(line, 256, inp) != NULL) {
		// Ignore comments.
		int len = strlen(line);
		if (line[0] == '#') continue;
		if (len < 2) continue;
		
		String s(line);
		int numTokens = s.tokenCount();
		if (numTokens != 2) {
			printf("SimplePotential::SimplePotential Invalid tabulated potential file line: %s\n", line);
			exit(-1);
		}
		
		String* tokenList = new String[numTokens];
		s.tokenize(tokenList);
		if (tokenList == NULL) {
			printf("SimplePotential::SimplePotential Invalid tabulated potential file line: %s\n", line);
			exit(-1);
		}
		r[count] = (float) strtod(tokenList[0], NULL);
		pot[count] = (float) strtod(tokenList[1], NULL);
		count++;
		
		delete[] tokenList;
	}
	fclose(inp);

	if (type == BOND) {
	    step_inv = (size-1.0f) / (r[size-1]-r[0]);
	} else {
	    step_inv = 57.29578f * (size-1.0f) / (r[size-1]-r[0]);
	}
	delete[] r;
}
