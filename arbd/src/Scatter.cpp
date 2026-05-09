///////////////////////////////////////////////////////////////////////
// An array of positions.
// Author: Jeff Comer <jcomer2@illinois.edu>

#include "Scatter.h"

Scatter::Scatter(const char* coordFile) {
	// Count the number of points.
	n = countCoordinates(coordFile);
    
	// Load the coordinates.
	r = new Vector3[n];
	readCoordinates(coordFile, n, r);
}

Scatter::Scatter(const char* coordFile, float cutTime) {
	// Count the number of points.
	n = countTrajectory(coordFile, cutTime);
    
	// Load the coordinates.
	r = new Vector3[n];
	readTrajectory(coordFile, n, r, cutTime);
}
  
Scatter::Scatter(const char* coordFile, float cutTime0, float cutTime1) {
	// Count the number of points.
	n = countTrajectory(coordFile, cutTime0, cutTime1);
    
	// Load the coordinates.
	r = new Vector3[n];
	readTrajectory(coordFile, n, r, cutTime0, cutTime1);
}

Scatter::~Scatter() {
	delete[] r;
}

Matrix3 Scatter::topMatrix() const {
	if (n < 3) return Matrix3(1.0f);
	return Matrix3(r[0], r[1], r[2]);
}

Vector3 Scatter::get(int i) const {
#ifdef DEBUG 
	if (i < 0 || i >= n) {
		printf("Warning! Scatter::get out of bounds.\n");
		return Vector3(0.0f);
	}
#endif
	return r[i];
}
int Scatter::length() const {
	return n;
}

Vector3 Scatter::minBound() const {
	Vector3 ret = r[0];
	for (int i = 1; i < n; i++) {
		if (r[i].x < ret.x) ret.x = r[i].x;
		if (r[i].y < ret.y) ret.y = r[i].y;
		if (r[i].z < ret.z) ret.z = r[i].z;
	}
	return ret;
}

Vector3 Scatter::maxBound() const {
	Vector3 ret = r[0];
	for (int i = 1; i < n; i++) {
		if (r[i].x > ret.x) ret.x = r[i].x;
		if (r[i].y > ret.y) ret.y = r[i].y;
		if (r[i].z > ret.z) ret.z = r[i].z;
	}
	return ret;
}

int Scatter::countCoordinates(const char* fileName) {
	int nRead;
	int n = 0;
	float x, y, z;
	char line[256];

	// Open the file.
	FILE* inp = fopen(fileName,"r");
	if (inp == NULL) {
		printf("Scatter:countCoordinates Couldn't open file %s\n.",fileName);
		exit(-1);
	}

	while (fgets(line, 256, inp) != NULL) {
		// Ignore comments.
		int len = strlen(line);
		if (line[0] == '#') continue;
		if (len < 2) continue;
    
		// Read values.
		nRead = sscanf(line, "%f %f %f", &x, &y, &z);
		if (nRead >= 3) n++;
	}
    
	fclose(inp);
	return n;
}

int Scatter::countTrajectory(const char* fileName, float cutTime) {
	int nRead;
	int n = 0;
	float t, x, y, z;
	char line[256];

	// Open the file.
	FILE* inp = fopen(fileName,"r");
	if (inp == NULL) {
		printf("Scatter:countCoordinates Couldn't open file %s\n.",fileName);
		exit(-1);
	}

	while (fgets(line, 256, inp) != NULL) {
		// Ignore comments.
		int len = strlen(line);
		if (line[0] == '#') continue;
		if (len < 2) continue;
    
		// Read values.
		nRead = sscanf(line, "%f %f %f %f", &t, &x, &y, &z);
		if (nRead >= 4 && t >= cutTime) n++;
	}
    
	fclose(inp);
	return n;
}

int Scatter::countTrajectory(const char* fileName, float cutTime0, float cutTime1) {
	int nRead;
	int n = 0;
	float t, x, y, z;
	char line[256];

	// Open the file.
	FILE* inp = fopen(fileName,"r");
	if (inp == NULL) {
		printf("Scatter:countCoordinates Couldn't open file %s\n.",fileName);
		exit(-1);
	}

	while (fgets(line, 256, inp) != NULL) {
		// Ignore comments.
		int len = strlen(line);
		if (line[0] == '#') continue;
		if (len < 2) continue;
    
		// Read values.
		nRead = sscanf(line, "%f %f %f %f", &t, &x, &y, &z);
		if (nRead >= 4 && t >= cutTime0 && t< cutTime1) n++;
	}
    
	fclose(inp);
	return n;
}

void Scatter::readCoordinates(const char* fileName, int num, Vector3* r) {
	int nRead;
	int n = 0;
	float x, y, z;
	char line[256];

	// Open the file.
	FILE* inp = fopen(fileName,"r");
	if (inp == NULL) {
		printf("Scatter:countCoordinates Couldn't open file %s\n.",fileName);
		exit(-1);
	}

	while (fgets(line, 256, inp) != NULL) {
		// Ignore comments.
		int len = strlen(line);
		if (line[0] == '#') continue;
		if (len < 2) continue;
    
		// Read values.
		nRead = sscanf(line, "%f %f %f", &x, &y, &z);
		if (nRead >= 3) {
			r[n].x = x;
			r[n].y = y;
			r[n].z = z;
			n++;
			if (n >= num) break;
		}
	}
    
	fclose(inp);
}

void Scatter::readTrajectory(const char* fileName, int num, Vector3* r, float cutTime) {
	int nRead;
	int n = 0;
	float t, x, y, z;
	char line[256];

	// Open the file.
	FILE* inp = fopen(fileName,"r");
	if (inp == NULL) {
		printf("Scatter:countCoordinates Couldn't open file %s\n.",fileName);
		exit(-1);
	}

	while (fgets(line, 256, inp) != NULL) {
		// Ignore comments.
		int len = strlen(line);
		if (line[0] == '#') continue;
		if (len < 2) continue;
    
		// Read values.
		nRead = sscanf(line, "%f %f %f %f", &t, &x, &y, &z);
		if (nRead >= 4 && t >= cutTime) {
			r[n].x = x;
			r[n].y = y;
			r[n].z = z;
			n++;
			if (n >= num) break;
		}
	}
    
	fclose(inp);
}

void Scatter::readTrajectory(const char* fileName, int num, Vector3* r, float cutTime0, float cutTime1) {
	int nRead;
	int n = 0;
	float t, x, y, z;
	char line[256];

	// Open the file.
	FILE* inp = fopen(fileName,"r");
	if (inp == NULL) {
		printf("Scatter:countCoordinates Couldn't open file %s\n.",fileName);
		exit(-1);
	}

	while (fgets(line, 256, inp) != NULL) {
		// Ignore comments.
		int len = strlen(line);
		if (line[0] == '#') continue;
		if (len < 2) continue;
    
		// Read values.
		nRead = sscanf(line, "%f %f %f %f", &t, &x, &y, &z);
		if (nRead >= 4 && t >= cutTime0 && t < cutTime1) {
			r[n].x = x;
			r[n].y = y;
			r[n].z = z;
			n++;
			if (n >= num) break;
		}
	}
    
	fclose(inp);
}

