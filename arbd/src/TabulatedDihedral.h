// TabulatedDihedral.h
// Authors: Justin Dufresne and Terrance Howard, 2013

#ifndef TABULATEDDIHEDRAL_H
#define TABULATEDDIHEDRAL_H

#include "useful.h"

#include "Dihedral.h"
#include "TabulatedPotential.h"
#include "BaseGrid.h"

// #include <math.h>
// #define _USING_MATH_DEFINES

class TabulatedDihedralPotential {
public:
	TabulatedDihedralPotential();
	TabulatedDihedralPotential(String fileName);
	TabulatedDihedralPotential(const TabulatedDihedralPotential &src);
	~TabulatedDihedralPotential();

	float* pot;				// actual potential values
	float angle_step_inv;	// angular increments of potential file
	int size;					// number of data points in the file
	String fileName;

	// RBTODO: deprecate
	HOST DEVICE inline EnergyForce computeOLD(Dihedral* d, Vector3* pos, BaseGrid* sys, int index) { 
		const Vector3 posa = d->ind1;
		const Vector3 posb = d->ind2;
		const Vector3 posc = d->ind3;
		const Vector3 posd = d->ind4;
		
		const Vector3 ab = sys->wrapDiff(posa - posb);
		const Vector3 bc = sys->wrapDiff(posb - posc);
		const Vector3 cd = sys->wrapDiff(posc - posd);
		
		//const float distab = ab.length();
		const float distbc = bc.length();
		//const float distcd = cd.length();
	
		Vector3 crossABC = ab.cross(bc);
		Vector3 crossBCD = bc.cross(cd);
		Vector3 crossX = bc.cross(crossABC);

		const float cos_phi = crossABC.dot(crossBCD) / (crossABC.length() * crossBCD.length());
		const float sin_phi = crossX.dot(crossBCD) / (crossX.length() * crossBCD.length());
		
		const float angle = -atan2(sin_phi, cos_phi);

		float energy = 0.0f;
		float force = 0.0f;
	
		Vector3 f1, f2, f3; // forces
		f1 = -distbc * crossABC.rLength2() * crossABC;
		f3 = -distbc * crossBCD.rLength2() * crossBCD;
		f2 = -(ab.dot(bc) * bc.rLength2()) * f1 - (bc.dot(cd) * bc.rLength2()) * f3;
	
		// Shift "angle" by "PI" since    -PI < dihedral < PI
		// And our tabulated potential data: 0 < angle < 2 PI
		float t = (angle + BD_PI) * angle_step_inv;
		int home = (int) floorf(t);
		t = t - home;

		home = home % size;
		int home1 = (home + 1) % size;

		//================================================
		// Linear interpolation
		float U0 = pot[home];       // Potential
		float dU = pot[home1] - U0; // Change in potential
		
		energy = dU * t + U0;
		force = -dU * angle_step_inv;
		//================================================
		// TODO: add an option for cubic interpolation

		if (crossABC.rLength() > 1.0f || crossBCD.rLength() > 1.0f)
			// avoid singularity when one angle is straight 
			force = 0.0f;

		f1 *= force;
		f2 *= force;
		f3 *= force;

		switch (index) {
			// Return energy and forces to appropriate particles 
			case 1: return EnergyForce(energy, f1);       // a
			case 2: return EnergyForce(energy, f2 - f1);  // b
			case 3: return EnergyForce(energy, f3 - f2);  // c
			case 4: return EnergyForce(energy, -f3);      // d
			default: return EnergyForce(0.0f, Vector3(0.0f));
		}
	}
};

#endif
