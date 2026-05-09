#pragma once

// Defined elsewhere: constexpr float BD_PI = 3.1415927f;

struct AngleForce {
    __host__ __device__
    AngleForce(Vector3 f1, Vector3 f3, float e) : f1(f1), f3(f3), e(e) { }
    Vector3 f1;
    Vector3 f3;
    float e;
};

__device__ inline void computeAngle(const TabulatedAnglePotential* __restrict__ a, const BaseGrid* __restrict__ sys, Vector3* force, const Vector3* __restrict__ pos,
				const int& i, const int& j, const int& k, float* energy, bool get_energy) {
	    
	    
	// Particle's type and position
	Vector3 posa = pos[i];
	Vector3 posb = pos[j];
	Vector3 posc = pos[k];
		
	// The vectors between each pair of particles
	const Vector3 ab = sys->wrapDiff(posa - posb);
	const Vector3 bc = sys->wrapDiff(posb - posc);
	const Vector3 ac = sys->wrapDiff(posc - posa);
  
	// Find the distance between each pair of particles
	float distab = ab.length2();
	float distbc = bc.length2();
	const float distac2 = ac.length2();
  
	// Find the cosine of the angle we want - <ABC	
	float cos = (distab + distbc - distac2);

	distab = 1.0f/sqrt(distab); //TODO: test other functiosn
	distbc = 1.0f/sqrt(distbc);
	cos *= 0.5f * distbc * distab;
  
	// If the cosine is illegitimate, set it to 1 or -1 so that acos won't fail
	if (cos < -1.0f) cos = -1.0f;
	if (cos > 1.0f) cos = 1.0f;

	// Find the sine while we're at it.

	// Now we can use the cosine to find the actual angle (in radians)		
	float angle = acos(cos);

	// transform angle to units of tabulated array index
	angle *= a->angle_step_inv;

	// tableAngle[0] stores the potential at angle_step
	// tableAngle[1] stores the potential at angle_step * 2, etc.
	// 'home' is the index after which 'convertedAngle' would appear if it were stored in the table	
	int home = int(floorf(angle));
        home =  (home >= a->size) ? (a->size)-1 : home; 
	//assert(home >= 0);
	//assert(home+1 < a->size);

	// // Make angle the distance from [0,1) from the first index in the potential array index
	// angle -= home;
		
	// Linearly interpolate the potential	
	float U0 = a->pot[home];
	float dUdx = (a->pot[(((home+1)==(a->size)) ? (a->size)-1 : home+1)] - U0) * a->angle_step_inv;
        if(get_energy)
        {
	    float e = ((dUdx * (angle-home)) + U0)*0.3333333333;
            atomicAdd( &energy[i], e);
            atomicAdd( &energy[j], e);
            atomicAdd( &energy[k], e);
            
        }
	float sin = sqrtf(1.0f - cos*cos);
	dUdx /= abs(sin) > 1e-3 ? sin : 1e-3; // avoid singularity 

	// Calculate the forces
	Vector3 force1 = -(dUdx*distab) * (ab * (cos*distab) + bc * distbc); // force on particle 1
	Vector3 force3 = (dUdx*distbc) * (bc * (cos*distbc) + ab * distab); // force on particle 3

	// assert( force1.length() < 10000.0f );
	// assert( force3.length() < 10000.0f );
	
	atomicAdd( &force[i], force1 );
	atomicAdd( &force[j], -(force1 + force3) );
	atomicAdd( &force[k], force3 );
}

__device__ inline AngleForce calcAngle(const TabulatedAnglePotential* __restrict__ a, const Vector3 ab, const Vector3 bc, const Vector3 ac) {
	// // The vectors between each pair of particles
	// const Vector3 ab = sys->wrapDiff(posa - posb);
	// const Vector3 bc = sys->wrapDiff(posb - posc);
	// const Vector3 ac = sys->wrapDiff(posc - posa);
 
	// Find the distance between each pair of particles
	float distab = ab.length2();
	float distbc = bc.length2();
	const float distac2 = ac.length2();
  
	// Find the cosine of the angle we want - <ABC	
	float cos = (distab + distbc - distac2);

	distab = 1.0f/sqrt(distab); //TODO: test other functiosn
	distbc = 1.0f/sqrt(distbc);
	cos *= 0.5f * distbc * distab;
  
	// If the cosine is illegitimate, set it to 1 or -1 so that acos won't fail
	if (cos < -1.0f) cos = -1.0f;
	if (cos > 1.0f) cos = 1.0f;

	// Find the sine while we're at it.

	// Now we can use the cosine to find the actual angle (in radians)		
	float angle = acos(cos);

	// transform angle to units of tabulated array index
	angle *= a->angle_step_inv;

	// tableAngle[0] stores the potential at angle_step
	// tableAngle[1] stores the potential at angle_step * 2, etc.
	// 'home' is the index after which 'convertedAngle' would appear if it were stored in the table	
	int home = int(floorf(angle));
        home =  (home >= a->size) ? (a->size)-1 : home; 
	//assert(home >= 0);
	//assert(home+1 < a->size);

	// // Make angle the distance from [0,1) from the first index in the potential array index
	// angle -= home;
		
	// Linearly interpolate the potential	
	float U0 = a->pot[home];
	float dUdx = (a->pot[(((home+1)==(a->size)) ? (a->size)-1 : home+1)] - U0) * a->angle_step_inv;
	float e = ((dUdx * (angle-home)) + U0);

	float sin = sqrtf(1.0f - cos*cos);
	dUdx /= abs(sin) > 1e-3 ? sin : 1e-3; // avoid singularity 

	// Calculate the forces
	Vector3 force1 = -(dUdx*distab) * (ab * (cos*distab) + bc * distbc); // force on particle 1
	Vector3 force3 = (dUdx*distbc) * (bc * (cos*distbc) + ab * distab); // force on particle 3

	return AngleForce(force1,force3,e);
}

__device__ inline void computeBondAngle(const TabulatedAnglePotential* __restrict__ a1,
					const TabulatedPotential* __restrict__ b, const TabulatedAnglePotential* __restrict__ a2,
					const BaseGrid* __restrict__ sys, Vector3* force, const Vector3* __restrict__ pos,
					const int& i, const int& j, const int& k, const int& l, float* energy, bool get_energy) {

	// Particle's type and position
	Vector3 posa = pos[i];
	Vector3 posb = pos[j];
	Vector3 posc = pos[k];
	Vector3 posd = pos[l];

	// The vectors between each pair of particles
	const Vector3 ab = sys->wrapDiff(posb - posa);
	const Vector3 bc = sys->wrapDiff(posc - posb);
	const Vector3 ca = sys->wrapDiff(posc - posa);
	AngleForce fe_a1 = calcAngle(a1, -ab,-bc,ca);

	float distbc = bc.length2();
	EnergyForce fe_b = b->compute(bc,distbc);

	const Vector3 cd = sys->wrapDiff(posd - posc);
	const Vector3 db = sys->wrapDiff(posd - posb);
	AngleForce fe_a2 = calcAngle(a2, -bc,-cd,db);

        if(get_energy)
        {
	    float e =  fe_a1.e * fe_b.e * fe_a2.e * 0.25f;
            atomicAdd( &energy[i], e);
            atomicAdd( &energy[j], e);
            atomicAdd( &energy[k], e);
            atomicAdd( &energy[l], e);
        }
	atomicAdd( &force[i], fe_a1.f1 * fe_b.e * fe_a2.e );
	atomicAdd( &force[j], 
		   -(fe_a1.f1 + fe_a1.f3) * fe_b.e * fe_a2.e
		   + fe_b.f * fe_a1.e * fe_a2.e
		   + fe_a2.f1 * fe_b.e * fe_a1.e 
	    );
	atomicAdd( &force[k], 
		   fe_a1.f3 * fe_b.e * fe_a2.e
		   - fe_b.f * fe_a1.e * fe_a2.e 
		   - (fe_a2.f1 + fe_a2.f3) * fe_b.e * fe_a1.e
	    );
	atomicAdd( &force[l], fe_a2.f3 * fe_b.e * fe_a1.e );
}


__device__ inline void computeDihedral(const TabulatedDihedralPotential* __restrict__ d,
				const BaseGrid* __restrict__ sys, Vector3* forces, const Vector3* __restrict__ pos,
				const int& i, const int& j, const int& k, const int& l, float* energy, bool get_energy) {
	const Vector3 posa = pos[i];
	const Vector3 posb = pos[j];
	const Vector3 posc = pos[k];
	const Vector3 posd = pos[l];
	// if (i >= 8738)  printf("Dihedral posa: (%f,%f,%f)\nDihedral posb: (%f,%f,%f)\nDihedral posc: (%f,%f,%f)\nDihedral posd: (%f,%f,%f)\n",
	// 		       posa.x,posa.y,posa.z,
	// 		       posb.x,posb.y,posb.z,
	// 		       posc.x,posc.y,posc.z,
	// 		       posd.x,posd.y,posd.z);

		
	const Vector3 ab = sys->wrapDiff(posa - posb);
	const Vector3 bc = sys->wrapDiff(posb - posc);
	const Vector3 cd = sys->wrapDiff(posc - posd);
		
	//const float distab = ab.length();
	const float distbc = bc.length();
	//const float distcd = cd.length();
	
	Vector3 crossABC = ab.cross(bc);
	Vector3 crossBCD = bc.cross(cd);
	Vector3 crossX = bc.cross(crossABC);
	// assert( crossABC.rLength2() <= 1.0f );
	// assert( crossBCD.rLength2() <= 1.0f );

	
	const float cos_phi = crossABC.dot(crossBCD) / (crossABC.length() * crossBCD.length());
	const float sin_phi = crossX.dot(crossBCD) / (crossX.length() * crossBCD.length());
		
	const float angle = -atan2(sin_phi, cos_phi);

	// float energy = 0.0f;
	float force = 0.0f;
	
	Vector3 f1, f2, f3; // forces
	f1 = -distbc * crossABC.rLength2() * crossABC;
	f3 = -distbc * crossBCD.rLength2() * crossBCD;
	f2 = -(ab.dot(bc) * bc.rLength2()) * f1 - (bc.dot(cd) * bc.rLength2()) * f3;
	
	// Shift "angle" by "PI" since    -PI < dihedral < PI
	// And our tabulated potential data: 0 < angle < 2 PI
	float t = (angle + BD_PI) * d->angle_step_inv;
	int home = (int) floorf(t);
	t = t - home;
        //home = home % (d->size);
        home = (home < d->size) ? home : d->size-1;
        int home1 = (home + 1) >= d->size ? (home+1-d->size) : home+1;

	//assert(home >= 0);
	//assert(home < d->size);
	// home = home % size;
	//int home1 = (home + 1) >= d->size ? (home+1-d->size) : home+1;

	//assert(home1 >= 0);
	//assert(home1 < d->size);

	//================================================
	// Linear interpolation
	float U0 = d->pot[home];       // Potential
	float dU = d->pot[home1] - U0; // Change in potential
	if(get_energy)
        {	
	    float e_local = (dU * t + U0)*0.25f;
            atomicAdd( &energy[i], e_local );
            atomicAdd( &energy[j], e_local );
            atomicAdd( &energy[k], e_local );
            atomicAdd( &energy[l], e_local );
        }
	force = -dU * d->angle_step_inv;
	// if (i >= 8738)  printf("Dihedral (angle,U0,dUdT): (%f,%f,%f)\n", angle*180.0f/BD_PI, U0, force);

	// avoid singularity when one angle is straight 
	// force = (distbc*distbc*crossABC.rLength2() > 1000.0f || distbc*distbc*crossBCD.rLength2() > 1000.0f) ? 0.0f : force;
	force = (ab.length2()*bc.length2()*crossABC.rLength2() > 100.0f || bc.length2()*cd.length2()*crossBCD.rLength2() > 100.0f) ? 0.0f : force;

	// if ( force > 1000.0f )
	//     printf("%f %d %d (%.4f %.4f) %.2f %f\n",force,home,home1, d->pot[home], d->pot[home1], dU, d->angle_step_inv);	    
	//assert( force < 10000.0f )

	//if( force != force ) 
            //force = 0.f;
        assert(force == force);
	if ( force > 1000.0f ) 
	    force = 1000.0f;
	if ( force < -1000.0f ) 
	    force = -1000.0f;
	//assert( force < 10000.0f );

	f1 *= force;
	f2 *= force;
	f3 *= force;

	// assert( f1.length() < 10000.0f );
	// assert( f2.length() < 10000.0f );
	// assert( f3.length() < 10000.0f );

	atomicAdd( &forces[i], f1 );
	atomicAdd( &forces[j], f2-f1 );
	atomicAdd( &forces[k], f3-f2 );
	atomicAdd( &forces[l], -f3 );
}
