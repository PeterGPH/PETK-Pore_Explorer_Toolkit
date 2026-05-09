// tabulatedAngle.cu
// Authors: Justin Dufresne and Terrance Howard, 2013

#include "TabulatedAngle.h"

TabulatedAnglePotential::TabulatedAnglePotential()
{
	pot = NULL;
	size = 0;
	fileName = "";
}

TabulatedAnglePotential::TabulatedAnglePotential(String fileName) : fileName(fileName)
{
	FILE* inp = fopen(fileName.val(), "r");
	if (inp == NULL) {
		printf("TabulatedAnglePotential: could not open file '%s'\n", fileName.val());
		exit(-1);
	}
	char line[256];
	int capacity = 256;
	float* angle = new float[capacity];
	pot = new float[capacity];
	size = 0;
	while(fgets(line, 256, inp)) {
		String s(line);
		int numTokens = s.tokenCount();
		String* tokenList = new String[numTokens];
		s.tokenize(tokenList);
		
		// Legitimate TABULATED ANGLE inputs have 2 tokens
		// ANGLE | VALUE
		// Any angle input line without exactly 2 tokens should be discarded
		if (numTokens != 2) {
			printf("Invalid angle input line: %s\n", line);
			continue;
		}		
		
		// Discard any empty line
		if (tokenList == NULL) {	
			printf("Empty angle input line: %s\n", line);
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
	angle_step_inv = 57.29578f * (size-1) / (angle[size-1]-angle[0]); 
	delete[] angle;
	fclose(inp);
}

TabulatedAnglePotential::TabulatedAnglePotential(const TabulatedAnglePotential &tab)
{
	size = tab.size;
	fileName = tab.fileName;
	pot = new float[size];
	for (int i = 0; i < size; i++)
		pot[i] = tab.pot[i];
	angle_step_inv = tab.angle_step_inv;
}

TabulatedAnglePotential::~TabulatedAnglePotential() {
	delete[] pot;
}
