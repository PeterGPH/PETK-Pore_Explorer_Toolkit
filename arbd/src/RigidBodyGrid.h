//////////////////////////////////////////////////////////////////////
// Copy of BaseGrid with some modificaitons
// 
#ifndef RBBASEGRID_H
#define RBBASEGRID_H
// #pragma once

#ifdef __CUDACC__
    #define HOST __host__
    #define DEVICE __device__
#else
    #define HOST 
    #define DEVICE 
#endif

#include "BaseGrid.h"
#include "useful.h"
#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cuda.h>

// using namespace std;

#define STRLEN 512

/* class ForceEnergy { */
/* public: */
/* 	DEVICE ForceEnergy(Vector3 &f, float &e) : */
/* 		f(f), e(e) {}; */
/* 	Vector3 f; */
/* 	float e; */
/* }; */

class RigidBodyGrid { 
	friend class SparseGrid;
	
public:
	/*                               \
	| CONSTRUCTORS, DESTRUCTORS, I/O |
	\===============================*/
	
	// RBTODO Fix?
	RigidBodyGrid(); // cmaffeo2 (2015) moved this out of protected, cause I wanted RigidBodyGrid in a struct
  // The most obvious of constructors.
	RigidBodyGrid(int nx0, int ny0, int nz0);

  // Make a copy of a BaseGrid grid.
  RigidBodyGrid(const BaseGrid& g);

  // Make an exact copy of a grid.
  RigidBodyGrid(const RigidBodyGrid& g);

  RigidBodyGrid mult(const RigidBodyGrid& g);

  RigidBodyGrid& operator=(const RigidBodyGrid& g);
  
	virtual ~RigidBodyGrid();

	/*             \
	| DATA METHODS |
	\=============*/
		
	void zero();
  
  bool setValue(int j, float v);

  bool setValue(int ix, int iy, int iz, float v);

  virtual float getValue(int j) const;

  HOST DEVICE float getValue(int ix, int iy, int iz) const;

  HOST DEVICE Vector3 getPosition(int j) const;
	HOST DEVICE Vector3 getPosition(int j, Matrix3 basis, Vector3 origin) const;
		
  IndexList index(int j) const;
  int indexX(int j) const;
  int indexY(int j) const;
  int indexZ(int j) const;
  int index(int ix, int iy, int iz) const;
  
  /* int index(Vector3 r) const; */
  /* int nearestIndex(Vector3 r) const; */

  HOST DEVICE inline int length() const { return nx*ny*nz; }

  HOST DEVICE inline int getNx() const {return nx;}
  HOST DEVICE inline int getNy() const {return ny;}
  HOST DEVICE inline int getNz() const {return nz;}
  HOST DEVICE inline int getSize() const {return nx*ny*nz;}

  HOST DEVICE inline int getRadius(Matrix3 basis) const {
	  // return radius of smallest sphere circumscribing grid
	  float radius = basis.transform(Vector3(nx,ny,nz)).length2();

	  float tmp = basis.transform(Vector3(-nx,ny,nz)).length2();
	  radius = tmp > radius ? tmp : radius;

	  tmp = basis.transform(Vector3(nx,-ny,nz)).length2();
	  radius = tmp > radius ? tmp : radius;

	  tmp = basis.transform(Vector3(nx,ny,-nz)).length2();
	  radius = tmp > radius ? tmp : radius;

	  return 0.5 * sqrt(radius);
  }
  DEVICE ForceEnergy interpolateForceDnamd(const Vector3& l) const;
  DEVICE float compute_V(float *a, float *x, float *y, float *z) const;
  DEVICE Vector3 compute_dV(float *a, float *x, float *y, float *z) const;
  DEVICE void compute_a(float *a, float *b) const;
  DEVICE void compute_b(float * __restrict__ b, int * __restrict__ inds) const; 
  // Add a fixed value to the grid.
  void shift(float s);

  // Multiply the grid by a fixed value.
  void scale(float s);
	
	DEVICE ForceEnergy interpolateForceDLinearly(const Vector3& l) const;
	DEVICE ForceEnergy interpolateForceD(Vector3 l) const;
  
public:
  int nx, ny, nz;
  int size;
  float* val;
};

#endif
