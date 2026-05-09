///////////////////////////////////////////////////////////////////////
// Author: Jeff Comer <jcomer2@illinois.edu>
#ifndef TABULATEDPOTENTIAL_H
#define TABULATEDPOTENTIAL_H

#ifdef __CUDACC__
    #define HOST __host__
    #define DEVICE __device__
#else
    #define HOST
    #define DEVICE
#endif

#include "useful.h"
#include "BaseGrid.h"

#ifdef __CUDA_ARCH__
#include "CudaUtil.cuh"
#endif

#include <cuda.h>

#ifndef gpuErrchk
#define delgpuErrchk
#define gpuErrchk(code) { if ((code) != cudaSuccess) {			                            \
	    fprintf(stderr,"CUDA Error: %s %s %d\n", cudaGetErrorString(code), __FILE__, __LINE__); \
	}}
#endif
constexpr float BD_PI = 3.1415927f;

class EnergyForce {
public:
  HOST DEVICE
	inline EnergyForce(float energy=0.0f, const Vector3& force=Vector3(0.0f)) :
			e(energy), f(force) { }
  HOST DEVICE
	inline EnergyForce operator+=(const EnergyForce& ef) {
		e += ef.e;
		f += ef.f;
		return *this;
	}
  float e;
  Vector3 f;
};

class TabulatedPotential;

class FullTabulatedPotential {
public:
  FullTabulatedPotential();
  FullTabulatedPotential(const char* fileName);
  FullTabulatedPotential(const FullTabulatedPotential& tab);
  ~FullTabulatedPotential();

  static int countValueLines(const char* fileName);

  /* HOST DEVICE inline EnergyForce computeOLD(Vector3 r) { */
  /* 		float d = r.length(); */
  /* 		Vector3 rUnit = -r/d; */
  /* 		int home = int(floorf((d - r0)/dr)); */
  /* 		if (home < 0) return EnergyForce(v0[0], Vector3(0.0f)); */
  /* 		if (home >= n) return EnergyForce(e0, Vector3(0.0f)); */
  /* 		float homeR = home*dr + r0; */
  /* 		float w = (d - homeR)/dr; */
		
  /* 		// Interpolate. */
  /* 		float energy = v3[home]*w*w*w + v2[home]*w*w + v1[home]*w + v0[home]; */
  /* 		Vector3 force = -(3.0f*v3[home] * w * w */
  /* 										+ 2.0f*v2[home] * w */
  /* 										+ v1[home]) * rUnit/dr; */
  /* 		return EnergyForce(energy,force); */
  /* 	} */

  TabulatedPotential* pot;

private:
  int numLines;
  String fileName;
};

class TabulatedPotential {
public:
  TabulatedPotential();
  TabulatedPotential(const TabulatedPotential& tab);
  TabulatedPotential(const float* dist, const float* pot, int n0);
    TabulatedPotential(const FullTabulatedPotential& tab) : TabulatedPotential(*tab.pot) {}
    TabulatedPotential(const char* filename) : TabulatedPotential(FullTabulatedPotential(filename)) {}
  ~TabulatedPotential();

  void truncate(float cutoff);
  bool truncate(float switchDist, float cutoff, float value);

  Vector3 computeForce(Vector3 r);

  int size() const { return n; }

    TabulatedPotential* copy_to_cuda() const {
	// Allocate data for array 
	TabulatedPotential* dev_ptr;
	TabulatedPotential tmp(*this); // TODO consider avoiding allocating v0

	float *v;
	{
	    size_t sz = sizeof(float) * n;
	    gpuErrchk(cudaMalloc(&v, sz));
	    gpuErrchk(cudaMemcpy(v, v0, sz, cudaMemcpyHostToDevice));
	}
	delete [] tmp.v0;
	tmp.v0 = v;

	size_t sz = sizeof(TabulatedPotential);
	gpuErrchk(cudaMalloc(&dev_ptr, sz));
	gpuErrchk(cudaMemcpy(dev_ptr, &tmp, sz, cudaMemcpyHostToDevice));
	tmp.v0 = NULL;
	return dev_ptr;
    }
    void free_from_cuda(TabulatedPotential* dev_ptr) const {
	TabulatedPotential tmp = TabulatedPotential();
	delete [] tmp.v0;
	gpuErrchk(cudaMemcpy(&tmp, dev_ptr, sizeof(TabulatedPotential), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaFree(dev_ptr));
	gpuErrchk(cudaFree(tmp.v0));
	tmp.v0 = NULL;
    }


  HOST DEVICE inline EnergyForce compute(Vector3 r) {
		float d = r.length();
		float w = (d - r0) * drInv;
		int home = int( floorf(w) );
		w = w - home;
		// if (home < 0) return EnergyForce(v0[0], Vector3(0.0f));
		home = home < 0 ? 0 : home;
		if (home >= n) return EnergyForce(v0[n-1], Vector3(0.0f));
		
		float u0 = v0[home];
		float du = home+1 < n ? v0[home+1]-u0 : 0;
				
		// Interpolate.
		float energy = du*w+u0;
		Vector3 force = (du*drInv/d)*r;
		return EnergyForce(energy,force);
	}

  HOST DEVICE inline EnergyForce compute(const Vector3 r, float d) const {
		d = sqrt(d);
		// float d = r.length();
		float w = (d - r0) * drInv;
		int home = int( floorf(w) );
		w = w - home;
		// if (home < 0) return EnergyForce(v0[0], Vector3(0.0f));
		home = home < 0 ? 0 : home;
		if (home >= n) return EnergyForce(v0[n-1], Vector3(0.0f));
		
		float u0 = v0[home];
		float du = home+1 < n ? v0[home+1]-u0 : 0;
				
		// Interpolate.
		float energy = du*w+u0;
		Vector3 force = (du*drInv/d)*r;
		return EnergyForce(energy,force);
	}
  HOST DEVICE inline Vector3 computef(const Vector3 r, float d) const {
		d = sqrt(d);
		// float d = r.length();
		// RBTODO: precompute so that initial blocks are zero; reduce computation here
		float w = (d - r0)*drInv;
		int home = int( floorf(w) );
		w = w - home;
		// if (home < 0) return EnergyForce(v0[0], Vector3(0.0f));
		home = home < 0 ? 0 : home;
		if (home >= n) return Vector3(0.0f);
		
		if (home+1 < n) 
		    return ((v0[home+1]-v0[home])*drInv/d)*r;
		else
		    return Vector3(0.0f);
	}

// private:
private:
  float* v0;
  int n;
  float drInv; //TODO replace with drInv
  float r0;
};

/* // New unified/simplified classes for working with potentials */
/* <template int num_indices, int max_integer> */
/* class BitMaskInts { */
/*     BitMaskInts(); */

/* private: */
/*     static_assert( ceil(log2(max_integer)) <= CHAR_BIT ); */

/*     char data[ ceil(num_indices * ceil(log2(max_integer)) / CHAR_BIT) ]; */

/*     HOST DEVICE inline unsigned short int get_int(i) const { */
/* 	unsigned int first_bit = i * ceil(log2(max_integer)); */
/* 	unsigned int last_bit  = (i+1) * ceil(log2(max_integer))-1; */
/* 	char c0 = data[floor(first_bit/CHAR_BIT)]; */
/* 	char c1 = data[floor(last_bit/CHAR_BIT)]; */

/* 	unsigned short int ret = c0 << (first_bit % CHAR_BIT) /\* shift left *\/ */
/*     };     */

/* } */


enum SimplePotentialType { UNSET, BOND, ANGLE, DIHEDRAL, VECANGLE };
// enum PotentialTypeAtoms { bond=2, angle=3, dihedral=4 };


class SimplePotential {
public:
    SimplePotential() { }
    SimplePotential(const char* filename, SimplePotentialType type);
    SimplePotential(float* pot, float step_inv, unsigned int size, SimplePotentialType type) :
	pot(pot), step_inv(step_inv), size(size), type(type) { }
    

    float* pot;	     // actual potential values
    float  step_inv; // angular increments of potential file
    unsigned int size;     // number of data points in the file

    SimplePotentialType type;

    /* float start = 0;  */
    /* bool is_periodic = false; */

    /* HOST void copy_to_device(SimplePotential* device_addr_p, unsigned int offset=0) { */
    /* 	/\* Assumes device_addr_p is already allocated, allocates space for pot *\/ */
    /* 	float* val, tmp; */
    /* 	gpuErrchk(cudaMalloc(&val, sizeof(float)*size)); // TODO equivalent cudaFree */
    /* 	gpuErrchk(cudaMemcpyAsync(val, pot, sizeof(float)*size, cudaMemcpyHostToDevice)); */
    /* 	tmp = pot; */
    /* 	pot = val; */
    /* 	gpuErrchk(cudaMemcpyAsync(device_addr_p+offset, this, sizeof(SimplePotential), cudaMemcpyHostToDevice)); */
    /* 	pot = tmp; */
    /* 	// return val; */
    /* } */

    HOST DEVICE inline float compute_value(const Vector3* __restrict__ pos,
					   const BaseGrid* __restrict__ sys,
					   const int* __restrict__ particles) const {
	float val;
	if (type == BOND)
	    val = compute_bond(pos, sys, particles[0], particles[1]);
	else if (type == ANGLE)
	    val = compute_angle(pos, sys, particles[0], particles[1], particles[2]);
	else if (type == DIHEDRAL)
	    val = compute_dihedral(pos, sys, particles[0], particles[1], particles[2], particles[3]);
	else if (type == VECANGLE)
	    val = compute_vecangle(pos, sys, particles[0], particles[1], particles[2], particles[3]);
	return val;
    }

    HOST DEVICE inline float2 compute_energy_and_deriv(float value) {
	float2 ret;
	if (type == DIHEDRAL) {
	    ret = linearly_interpolate<true>(value, -BD_PI);
	} else {
	    ret = linearly_interpolate<false>(value);
	}
	return ret;
    }

    HOST DEVICE inline float compute_bond(const Vector3* __restrict__ pos,
					      const BaseGrid* __restrict__ sys,
					      int i, int j) const {
	return sys->wrapDiff( pos[j] - pos[i] ).length();
    }

    HOST DEVICE inline float compute_angle(const Vector3* __restrict__ pos,
					   const BaseGrid* __restrict__ sys,
					   int i, int j, int k) const {
	const Vector3 ab = sys->wrapDiff(pos[j] - pos[i]);
	const Vector3 bc = sys->wrapDiff(pos[k] - pos[j]);
	const Vector3 ac = sys->wrapDiff(pos[k] - pos[i]);
	return compute_angle( ab.length2(), bc.length2(), ac.length2() );
    }

    HOST DEVICE inline float compute_vecangle(const Vector3* __restrict__ pos,
					      const BaseGrid* __restrict__ sys,
					      int i, int j, int k, int l) const {
	const Vector3 ab = sys->wrapDiff(pos[j] - pos[i]);
	const Vector3 bc = sys->wrapDiff(pos[l] - pos[k]);
	const Vector3 ac = bc+ab;
	return compute_angle( ab.length2(), bc.length2(), ac.length2() );
    }

    HOST DEVICE inline float compute_angle(float distab2, float distbc2, float distac2) const {
	// Find the cosine of the angle we want - <ABC
	float cos = (distab2 + distbc2 - distac2);

	distab2 = 1.0f/sqrt(distab2); //TODO: test other functions
	distbc2 = 1.0f/sqrt(distbc2);
	cos *= 0.5f * distbc2 * distab2;

	// If the cosine is illegitimate, set it to 1 or -1 so that acos won't fail
	if (cos < -1.0f) cos = -1.0f;
	if (cos > 1.0f) cos = 1.0f;

	return acos(cos);
    }

    HOST DEVICE inline float compute_dihedral(const Vector3* __restrict__ pos,
					      const BaseGrid* __restrict__ sys,
					      int i, int j, int k, int l) const {
	const Vector3 ab = -sys->wrapDiff(pos[j] - pos[i]);
	const Vector3 bc = -sys->wrapDiff(pos[k] - pos[j]);
	const Vector3 cd = -sys->wrapDiff(pos[l] - pos[k]);

	const Vector3 crossABC = ab.cross(bc);
	const Vector3 crossBCD = bc.cross(cd);
	const Vector3 crossX = bc.cross(crossABC);

	const float cos_phi = crossABC.dot(crossBCD) / (crossABC.length() * crossBCD.length());
	const float sin_phi = crossX.dot(crossBCD) / (crossX.length() * crossBCD.length());

	return -atan2(sin_phi,cos_phi);
    }

    template <bool is_periodic>
	HOST DEVICE inline float2 linearly_interpolate(float x, float start=0.0f) const {
	float w = (x - start) * step_inv;
	int home = int( floorf(w) );
	w = w - home;
	// if (home < 0) return EnergyForce(v0[0], Vector3(0.0f));
	if (home < 0) {
	    if (is_periodic) home += size;
	    else return make_float2(pot[0],0.0f);
	}
	else if (home >= size) {
	    if (is_periodic) home -= size;
	    else return make_float2(pot[size-1],0.0f);
	}

	float u0 = pot[home];
	float du = home+1 < size ? pot[home+1]-u0 : is_periodic ? pot[0]-u0 : 0;

	return make_float2(du*w+u0, du*step_inv);
    }

    DEVICE inline void apply_force(const Vector3* __restrict__ pos,
				   const BaseGrid* __restrict__ sys,
				   Vector3* __restrict__ forces,
				   int* particles, float energy_deriv) const {
	if (type == BOND)
	    apply_bond_force(pos, sys, forces, particles[0], particles[1], energy_deriv);
	else if (type == ANGLE)
	    apply_angle_force(pos, sys, forces, particles[0], particles[1],
			     particles[2], energy_deriv);
	else if (type == DIHEDRAL)
	    apply_dihedral_force(pos, sys, forces, particles[0], particles[1],
				 particles[2], particles[3], energy_deriv);
	else if (type == VECANGLE)
	    apply_vecangle_force(pos, sys, forces, particles[0], particles[1],
				 particles[2], particles[3], energy_deriv);
    }

    __device__ inline void apply_bond_force(const Vector3* __restrict__ pos,
					const BaseGrid* __restrict__ sys,
					Vector3* __restrict__ force,
					int i, int j, float energy_deriv) const {
#ifdef __CUDA_ARCH__
	Vector3 f = sys->wrapDiff( pos[j] - pos[i] );
	f = f * energy_deriv / f.length();
	atomicAdd(&force[i], f);
	atomicAdd(&force[j], -f);
#endif
    }

    struct TwoVector3 {
	Vector3 v1;
	Vector3 v2;
    };

    DEVICE inline TwoVector3 get_angle_force(const Vector3& ab,
					     const Vector3& bc,
					     float energy_deriv) const {
	// Find the distance between each pair of particles
	float distab = ab.length2();
	float distbc = bc.length2();
	const float distac2 = (ab+bc).length2();

	// Find the cosine of the angle we want - <ABC
	float cos = (distab + distbc - distac2);

	distab = 1.0f/sqrt(distab); //TODO: test other functions
	distbc = 1.0f/sqrt(distbc);
	cos *= 0.5f * distbc * distab;

	// If the cosine is illegitimate, set it to 1 or -1 so that acos won't fail
	if (cos < -1.0f) cos = -1.0f;
	if (cos > 1.0f) cos = 1.0f;

	float sin = sqrtf(1.0f - cos*cos);
	energy_deriv /= abs(sin) > 1e-3 ? sin : 1e-3; // avoid singularity
	if (abs(sin) < 1e-3) {
	    printf("BAD ANGLE: sin, cos, energy_deriv, distab, distbc, distac2: (%f %f %f %f %f)\n",
		   sin,cos,energy_deriv,distab,distbc);	}

	// Calculate the forces
	TwoVector3 force;
	force.v1 = (energy_deriv*distab) * (ab * (cos*distab) + bc * distbc); // force on 1st particle
	force.v2 = -(energy_deriv*distbc) * (bc * (cos*distbc) + ab * distab); // force on last particle
	return force;
    }

    // DEVICE inline TwoVector3 get_angle_force(const Vector3& ab,
    // 					     const Vector3& bc,
    // 					     float energy_deriv) const {
    // 	// Find the distance between each pair of particles
    // 	float distab = ab.length2();
    // 	float distbc = bc.length2();

    // 	float pre = distab*distbc - pow(ab.dot(bc),2);
    // 	// if (pre < 1e-6) {
    // 	//     pre = 1e-3;
    // 	//     printf("BAD ANGLE: pre, energy_deriv, distab, distbc, ab.dot(bc): (%f %f %f %f %f)\n",
    // 	// 	   pre,energy_deriv,distab,distbc,ab.dot(bc));	
    // 	// } else pre = sqrt(pre);
    // 	// if (distab == distbc) {
    // 	//     printf("GOOD ANGLE: pre, energy_deriv, distab, distbc, ab.dot(bc): (%f %f %f %f %f)\n",
    // 	// 	   pre,energy_deriv,distab,distbc,ab.dot(bc));	
    // 	// }
	    
    // 	pre = pre > 1e-6 ? sqrt(pre) : 1e-3;
    // 	energy_deriv /= pre;

    // 	TwoVector3 force;
    // 	//force.v1 = energy_deriv * Vector3::element_mult( 1-Vector3::element_mult(ab,ab)/distab, bc);
    // 	//force.v2 = energy_deriv * Vector3::element_mult( 1-Vector3::element_mult(bc,bc)/distbc, ab);

    // 	Vector3 abbc = Vector3::element_mult(-ab,bc);
    // 	force.v1 = -energy_deriv * (bc-Vector3::element_mult(abbc, -ab/distab));
    // 	force.v2 = -energy_deriv * (-ab-Vector3::element_mult(abbc, bc/distbc));
    // 	return force;
    // }

    DEVICE inline void apply_angle_force(const Vector3* __restrict__ pos,
					 const BaseGrid* __restrict__ sys,
					 Vector3* __restrict__ force,
					 int i, int j, int k, float energy_deriv) const {

#ifdef __CUDA_ARCH__
	const Vector3 ab = sys->wrapDiff(pos[j] - pos[i]);
	const Vector3 bc = sys->wrapDiff(pos[k] - pos[j]);
	// const Vector3 ac = sys->wrapDiff(pos[k] - pos[i]);

	TwoVector3 f = get_angle_force(ab,bc, energy_deriv);

	atomicAdd( &force[i], f.v1 );
	atomicAdd( &force[j], -(f.v1 + f.v2) );
	atomicAdd( &force[k], f.v2 );
#endif
    }

    DEVICE inline void apply_dihedral_force(const Vector3* __restrict__ pos,
					    const BaseGrid* __restrict__ sys,
					    Vector3* __restrict__ force,
					    int i, int j, int k, int l, float energy_deriv) const {
#ifdef __CUDA_ARCH__
	const Vector3 ab = -sys->wrapDiff(pos[j] - pos[i]);
	const Vector3 bc = -sys->wrapDiff(pos[k] - pos[j]);
	const Vector3 cd = -sys->wrapDiff(pos[l] - pos[k]);

	const Vector3 crossABC = ab.cross(bc);
	const Vector3 crossBCD = bc.cross(cd);
	const Vector3 crossX = bc.cross(crossABC);

	const float cos_phi = crossABC.dot(crossBCD) / (crossABC.length() * crossBCD.length());
	const float sin_phi = crossX.dot(crossBCD) / (crossX.length() * crossBCD.length());

	// return -atan2(sin_phi,cos_phi);
	Vector3 f1, f2, f3; // forces
	float distbc = bc.length();

	f1 = -distbc * crossABC.rLength2() * crossABC;
	f3 = -distbc * crossBCD.rLength2() * crossBCD;
	f2 = -(ab.dot(bc) * bc.rLength2()) * f1 - (bc.dot(cd) * bc.rLength2()) * f3;

	// energy_deriv = (ab.length2()*bc.length2()*crossABC.rLength2() > 100.0f || bc.length2()*cd.length2()*crossBCD.rLength2() > 100.0f) ? 0.0f : energy_deriv;
	/* if ( energy_deriv > 1000.0f ) */
	/*     energy_deriv = 1000.0f; */
	/* if ( energy_deriv < -1000.0f ) */
	/*     energy_deriv = -1000.0f; */

	f1 *= energy_deriv;
	f2 *= energy_deriv;
	f3 *= energy_deriv;

	atomicAdd( &force[i], f1 );
	atomicAdd( &force[j], f2-f1 );
	atomicAdd( &force[k], f3-f2 );
	atomicAdd( &force[l], -f3 );
#endif
    }
    DEVICE inline void apply_vecangle_force(const Vector3* __restrict__ pos,
					    const BaseGrid* __restrict__ sys,
					    Vector3* __restrict__ force,
					    int i, int j, int k, int l, float energy_deriv) const {

#ifdef __CUDA_ARCH__

	const Vector3 ab = -sys->wrapDiff(pos[j] - pos[i]);
	const Vector3 bc = -sys->wrapDiff(pos[k] - pos[j]);
	// const Vector3 ac = sys->wrapDiff(pos[k] - pos[i]);

	TwoVector3 f = get_angle_force(ab,bc, energy_deriv);

	atomicAdd( &force[i], f.v1 );
	atomicAdd( &force[j], -f.v1 );
	atomicAdd( &force[k], -f.v2 );
	atomicAdd( &force[l], f.v2 );
#endif
    }

};

#ifndef delgpuErrchk
#undef  delgpuErrchk
#undef  gpuErrchk(code)
#endif

#endif
