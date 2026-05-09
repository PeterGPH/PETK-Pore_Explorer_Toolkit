// tabulatedAngle.h
// Authors: Justin Dufresne and Terrance Howard, 2013
#ifndef TABULATEDANGLE_H
#define TABULATEDANGLE_H

#include "useful.h"
#include "Angle.h"
#include "TabulatedPotential.h"
#include "BaseGrid.h"

__device__ void atomicAdd( Vector3* address, Vector3 val);


class TabulatedAnglePotential
{
public:
	TabulatedAnglePotential();
	TabulatedAnglePotential(String fileName);
	TabulatedAnglePotential(const TabulatedAnglePotential &tab);
	~TabulatedAnglePotential();
	float* pot;			// actual potential values
	float angle_step_inv;	// '1/step' angle in potential file. potential file might not go 1, 2, 3,...,360, it could be in steps of .5 or something smaller 
	int size;			// The number of data points in the file
	String fileName;

	HOST DEVICE inline EnergyForce computeOLD(Angle* a, Vector3* pos, BaseGrid* sys, int index) {
		// First, we must find the actual angle we're working with. 
		// Grab the positions of each particle in the angle
		const Vector3 posa = pos[a->ind1];
		const Vector3 posb = pos[a->ind2];
		const Vector3 posc = pos[a->ind3];

		// The vectors between each pair of particles
		const Vector3 ab = sys->wrapDiff(posa - posb);
		const Vector3 bc = sys->wrapDiff(posb - posc);
		const Vector3 ac = sys->wrapDiff(posc - posa);
  
		// Find the distance between each pair of particles
		const float distab = ab.length();
		const float distbc = bc.length();
		const float distac = ac.length();
  
		// Find the cosine of the angle we want - <ABC	
		float cos = (distbc * distbc + distab * distab - distac * distac) / (2.0f * distbc * distab);
  
		// If the cosine is illegitimate, set it to 1 or -1 so that acos won't fail
		if (cos < -1.0f) cos = -1.0f;
		if (cos > 1.0f) cos = 1.0f;

		// Find the sine while we're at it.
		float sin = sqrtf(1.0f - cos*cos);

		// Now we can use the cosine to find the actual angle (in radians)		
		float angle = acos(cos);

		// tableAngle is divided into units of angle_step length
		// 'convertedAngle' is the angle, represented in these units
		float convertedAngle = angle * angle_step_inv;

		// tableAngle[0] stores the potential at angle_step
		// tableAngle[1] stores the potential at angle_step * 2, etc.
		// 'home' is the index after which 'convertedAngle' would appear if it were stored in the table	

		int home = int(floor(convertedAngle));

		// diffHome is the distance between the convertedAngle and the home index
		float diffHome = convertedAngle - home;

		// Linear interpolation for the potential
		float pot0 = pot[home];
		float delta_pot = pot[(home+1) % size] - pot0;
		float energy = (delta_pot * angle_step_inv) * diffHome + pot0;
		float diff = -delta_pot * angle_step_inv;
		diff /= sin;

		// Don't know what these are for, so I didn't bother giving them better names. 
		// Sorry, future person.
		float c1 = diff / distab;
		float c2 = diff / distbc;

		// Calculate the forces
		Vector3 force1 = c1 * (ab * (cos / distab) - bc / distbc); // force on particle 1
		Vector3 force3 = c2 * (bc * (cos / distbc) - ab / distab); // force on particle 3
		Vector3 force2 = -(force1 + force3); // the force on particle 2 (the central particle)

		EnergyForce ret;
		if (index == 1)
			ret = EnergyForce(energy, force1);
		if (index == 2)
			ret = EnergyForce(energy, force2);
		if (index == 3)
			ret = EnergyForce(energy, force3);
		return ret;
	}
};

#endif
