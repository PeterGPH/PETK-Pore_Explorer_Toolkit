// Configuration file reader
// Author: Jeff Comer <jcomer2@illinois.edu>

#include "Reservoir.h"
#include <cuda.h>


Reservoir::Reservoir(const char* reservoirFile) {
	reservoirs = countReservoirs(reservoirFile);
	r0 = new Vector3[reservoirs];
	r1 = new Vector3[reservoirs];
	num = new float[reservoirs];

	readReservoirs(reservoirFile);
	validateRegions();
}

Reservoir::~Reservoir() {
	delete[] r0;
	delete[] r1;
	delete[] num;
}

int Reservoir::countReservoirs(const char* reservoirFile) {
	// Open the file.
	FILE* inp = fopen(reservoirFile, "r");
	if (inp == NULL) {
		printf("Reservoir:Reservoir Couldn't open file `%s'.\n", reservoirFile);
		exit(-1);
	}

	int count = 0;
	float x0, y0, z0, x1, y1, z1;
	float n;
	char line[STRLEN];
	int nRead;
	while (fgets(line, STRLEN, inp) != NULL) {
		// Ignore comments.
		int len = strlen(line);
		if (line[0] == '#') continue;
		if (len < 2) continue;
      
		// Read definition lines.
		nRead = sscanf(line, "%f %f %f %f %f %f %f", &x0, &y0, &z0, &x1, &y1, &z1, &n);
		if (nRead < 7) {
			printf("Reservoir:Reservoir Improperly formatted line `%s'\n", line);
			fclose(inp);
			exit(-1);
		}
		count++;
	}
	return count;
}

void Reservoir::readReservoirs(const char* reservoirFile) {
	// Open the file.
	FILE* inp = fopen(reservoirFile, "r");
	if (inp == NULL) {
		printf("Reservoir:Reservoir Couldn't open file `%s'.\n", reservoirFile);
		exit(-1);
	}

	int count = 0;
	float x0, y0, z0, x1, y1, z1;
	float n;
	char line[STRLEN];
	int nRead;
	while (fgets(line, STRLEN, inp) != NULL) {
		// Ignore comments.
		int len = strlen(line);
		if (line[0] == '#') continue;
		if (len < 2) continue;
      
		// Read definition lines.
		nRead = sscanf(line, "%f %f %f %f %f %f %f", &x0, &y0, &z0, &x1, &y1, &z1, &n);
		if (nRead < 7) {
			printf("Reservoir:Reservoir Improperly formatted line `%s'\n", line);
			fclose(inp);
			exit(-1);
		}

		r0[count] = Vector3(x0, y0, z0);
		r1[count] = Vector3(x1, y1, z1);
		num[count] = n;

		count++;
	}
}

void Reservoir::validateRegions() {
	for (int i = 0; i < reservoirs; i++) {
		Vector3 a = r0[i];
		Vector3 b = r1[i];

		if (a.x > b.x) {r0[i].x = b.x; r1[i].x = a.x;}
		if (a.y > b.y) {r0[i].y = b.y; r1[i].y = a.y;}
		if (a.z > b.z) {r0[i].z = b.z; r1[i].z = a.z;}
	}
}

Vector3 Reservoir::getOrigin(int i) const {
	if (i < 0 || i >= reservoirs) return Vector3(0.0f);
	return r0[i];
}
Vector3 Reservoir::getDestination(int i) const {
	if (i < 0 || i >= reservoirs) return Vector3(0.0f);
	return r1[i];
}
Vector3 Reservoir::getDifference(int i) const {
	if (i < 0 || i >= reservoirs) return Vector3(0.0f);
	return r1[i] - r0[i];
}

//TODO: check getMeanNumber function
float Reservoir::getMeanNumber(int i) const {
	if (i < 0 || i >= reservoirs) return 0.0f;
	return num[i];
}
int Reservoir::length() const {
	return reservoirs;
}

bool Reservoir::inside(int i, Vector3 r) const {
	if (i < 0 || i >= reservoirs) return false;
	if (r.x < r0[i].x || r.x >= r1[i].x) return false;
	if (r.y < r0[i].y || r.y >= r1[i].y) return false;
	if (r.z < r0[i].z || r.z >= r1[i].z) return false;
	return true;
}

