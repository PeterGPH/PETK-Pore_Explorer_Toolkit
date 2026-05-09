// Angle.cu
// Copyright Justin Dufresne and Terrance Howard, 2013

#include "Angle.h"

void Angle::print()
{
//	printf("about to print fileName %p\n", fileName.val());
//	fileName.print();
	printf("ANGLE (%d %d %d) %s\n", ind1, ind2, ind3, fileName.val());
}

String Angle::toString()
{
	return String("ANGLE ") + ind1 + " " + ind2 + " " + ind3 + " " + fileName;
}

// void BondAngle::print()
// {
// //	printf("about to print fileName %p\n", fileName.val());
// //	fileName.print();
//     printf("BONDANGLE (%d %d %d) %s; %s; %s\n", ind1, ind2, ind3, angleFileName.val(), bondFileName1.val(), bondFileName2.val());
// }

// String BondAngle::toString()
// {
//     return String("BONDANGLE ") + ind1 + " " + ind2 + " " + ind3 + " " + angleFileName + " " + bondFileName1 + " " + bondFileName2;
// }

