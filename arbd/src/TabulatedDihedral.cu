// TabulatedDihedral.cu 
// Authors: Justin Dufresne and Terrance Howard, 2013

#include "TabulatedDihedral.h"
#include <cassert>

TabulatedDihedralPotential::TabulatedDihedralPotential() :
		pot(NULL), size(0), fileName("") {}

TabulatedDihedralPotential::TabulatedDihedralPotential(String fileName) : fileName(fileName), size(0) {
	FILE* inp = fopen(fileName.val(), "r");
	if (inp == NULL) {
		printf("TabulatedDihedralPotential: could not open file '%s'\n", fileName.val());
		exit(-1);
	}
	char line[256];
	int capacity = 256;
	float* angle = new float[capacity];
	pot = new float[capacity];
	while(fgets(line, 256, inp)) {
		String s(line);
		int numTokens = s.tokenCount();
		String* tokenList = new String[numTokens];
		s.tokenize(tokenList);
		
		// Legitimate TABULATED Dihedral inputs have 2 tokens
		// ANGLE | VALUE
		// Any angle input line without exactly 2 tokens should be discarded
		if (numTokens != 2) {
			printf("Invalid Dihedral input line: %s\n", line);
			continue;
		}		
		
		// Discard any empty line
		if (tokenList == NULL) {
			printf("Empty Dihedral input line: %s\n", line);
			continue;
		}
		
		if (size >= capacity) {
			float* temp = angle;
			float* temp2 = pot;
			capacity *= 2;
			angle = new float[capacity];
			pot = new float[capacity];
			for (int i = 0; i < size; i++) {
				angle[i] = temp[i];
				pot[i] = temp2[i];
			}
			delete[] temp;
			delete[] temp2;
		}	
		angle[size] = atof(tokenList[0].val());
		pot[size++] = atof(tokenList[1].val());
	}
	// units "1/deg" "1/radian"  *57.29578
	float deltaAngle = (angle[size-1]-angle[0])/(size-1); 
	assert( deltaAngle > 0 );
	assert( size*deltaAngle >= 360 );

	float tmp[size];
	for (int j = 0; j < size; ++j) {
	    // reorganize data so that the angle goes from [-Pi,Pi+delta)
	    float a = -180.0f - angle[0] + j*deltaAngle;
	    while (a < 0) a += 360.0f;
	    while (a >= size*deltaAngle) a -= 360.0f;
	    int i = round( a / deltaAngle );
	    assert(i >= 0);
	    assert(i < size);
	    tmp[j] = pot[i];
	}
	for (int i = 0; i < size; ++i) pot[i] = tmp[i];

	angle_step_inv = 57.29578f / deltaAngle;
		 
	delete[] angle;
	fclose(inp);
}

TabulatedDihedralPotential::TabulatedDihedralPotential(const TabulatedDihedralPotential &src) :
		size(src.size), fileName(src.fileName), angle_step_inv(src.angle_step_inv) {
	pot = new float[size];
	for (int i = 0; i < size; i++)
		pot[i] = src.pot[i];
}

TabulatedDihedralPotential::~TabulatedDihedralPotential() {
	delete[] pot;
}
