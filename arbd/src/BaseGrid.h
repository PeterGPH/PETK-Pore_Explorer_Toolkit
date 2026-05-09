//////////////////////////////////////////////////////////////////////
// Grid base class that does just the basics.
// Author: Jeff Comer <jcomer2@illinois.edu>
#ifndef BASEGRID_H
#define BASEGRID_H


#ifdef __CUDACC__
    #define HOST __host__
    #define DEVICE __device__
#else
    #define HOST 
    #define DEVICE 
#endif

#include "useful.h"
#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <ctime>
// #include <cuda.h>

#ifndef gpuErrchk
#define delgpuErrchk
#define gpuErrchk(code) { if ((code) != cudaSuccess) {			                            \
	    fprintf(stderr,"CUDA Error: %s %s %d\n", cudaGetErrorString(code), __FILE__, __LINE__); \
	}}
#endif


enum BoundaryCondition { dirichlet, neumann, periodic };
enum InterpolationOrder { linear = 1, cubic = 3 };

#define INTERPOLATE_FORCE(result, function, boundary_condition, args) \
switch (boundary_condition) { \
case dirichlet: \
    result = function<dirichlet>(args); break; \
case neumann: \
    result = function<neumann>(args); break; \
case periodic: \
    result = function<periodic>(args); break; \
}

// using namespace std;

#define STRLEN 512

DEVICE float __fdividef(float x, float y);

class NeighborList {
public:
  float v[3][3][3];
};

class BaseGrid {
  friend class SparseGrid;
 
private:
  // Initialize the variables that get used a lot.
  // Also, allocate the main value array.
  void init();

public:

	/*                               \
	| CONSTRUCTORS, DESTRUCTORS, I/O |
	\===============================*/
	
	// RBTODO Fix?
	BaseGrid(); // cmaffeo2 (2015) moved this out of protected, cause I wanted BaseGrid in a struct
  // The most obvious of constructors.
		BaseGrid(Matrix3 basis0, Vector3 origin0, int nx0, int ny0, int nz0);

  // Make an orthogonal grid given the box dimensions and resolution.
  BaseGrid(Vector3 box, float dx);

  // The box gives the system geometry.
  // The grid point numbers define the resolution.
  BaseGrid(Matrix3 box, int nx0, int ny0, int nz0);

  // The box gives the system geometry.
  // dx is the approx. resolution.
  // The grid spacing is always a bit larger than dx.
  BaseGrid(Matrix3 box, Vector3 origin0, float dx);

  // The box gives the system geometry.
  // dx is the approx. resolution.
  // The grid spacing is always a bit smaller than dx.
  BaseGrid(Matrix3 box, float dx);

  // Make an exact copy of a grid.
  BaseGrid(const BaseGrid& g);

  BaseGrid mult(const BaseGrid& g);

  BaseGrid& operator=(const BaseGrid& g);

  // Make a copy of a grid, but at a different resolution.
  BaseGrid(const BaseGrid& g, int nx0, int ny0, int nz0);

  // Read a grid from a file.
  BaseGrid(const char* fileName);
  
  // Write without comments.
  virtual void write(const char* fileName) const;

  // Writes the grid as a file in the dx format.
  virtual void write(const char* fileName, const char* comments) const;

  // Writes the grid data as a single column in the order:
  // nx ny nz ox oy oz dxx dyx dzx dxy dyy dzy dxz dyz dzz val0 val1 val2 ...
  virtual void writeData(const char* fileName);
 
  // Write the valies in a single column.
  virtual void writePotential(const char* fileName) const;
  
	virtual ~BaseGrid();

	/*             \
	| DATA METHODS |
	\=============*/
		
	void zero();
  
  bool setValue(int j, float v);

  bool setValue(int ix, int iy, int iz, float v);

  virtual float getValue(int j) const;

  //virtual float getValue(int ix, int iy, int iz) const;

  Vector3 getPosition(int ix, int iy, int iz) const;

  Vector3 getPosition(int j) const;

  // Does the point r fall in the grid?
  // Obviously this is without periodic boundary conditions.
  bool inGrid(Vector3 r) const;

  bool inGridInterp(Vector3 r) const;

  Vector3 transformTo(Vector3 r) const;

  Vector3 transformFrom(Vector3 l) const;

  IndexList index(int j) const;
  int indexX(int j) const;
  int indexY(int j) const;
  int indexZ(int j) const;
  int index(int ix, int iy, int iz) const;
  
  int index(Vector3 r) const;

  int nearestIndex(Vector3 r) const;

  HOST DEVICE inline int length() const {
		return size;
	}
  void setBasis(const Matrix3& b);
  void setOrigin(const Vector3& o);

  HOST DEVICE inline Vector3 getOrigin() const {return origin;}
  HOST DEVICE inline Matrix3 getBasis() const {return basis;}
  HOST DEVICE inline Matrix3 getInverseBasis() const {return basisInv;}
  HOST DEVICE inline int getNx() const {return nx;}
  HOST DEVICE inline int getNy() const {return ny;}
  HOST DEVICE inline int getNz() const {return nz;}
  HOST DEVICE inline int getSize() const {return size;}

  
  // A matrix defining the basis for the entire system.
  Matrix3 getBox() const;
  // The diagonal (nx,ny,nz) of the system.
  Vector3 getExtent() const;
  // The length of diagonal (nx,ny,nz) of the system.
  float getDiagonal() const;

  HOST DEVICE inline int getRadius() const {
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

  // The position farthest from the origin.
  Vector3 getDestination() const;
  // The center of the grid.
  Vector3 getCenter() const;
  // The volume of a single cell.
  float getCellVolume() const;
  // The volume of the entire system.
  float getVolume() const;
  Vector3 getCellDiagonal() const;

  // Add a fixed value to the grid.
  void shift(float s);

  // Multiply the grid by a fixed value.
  void scale(float s);

	/*         \
	| COMPUTED |
	\=========*/
	
  // Get the mean of the entire grid.
  float mean() const;
	
  // Compute the average profile along an axis.
  // Assumes that the grid axis with index "axis" is aligned with the world axis of index "axis".
  void averageProfile(const char* fileName, int axis);

  // Get the potential at the closest node.
  virtual float getPotential(Vector3 pos) const;

	// crop
	// Cuts the grid down
	// @param		boundries to crop to (x0, y0, z0) -> (x1, y1, z1);
	//					whether to change the origin
	// @return	success of the function
	bool crop(int x0, int y0, int z0, int x1, int y1, int z1, bool keep_origin);

  // Added by Rogan for times when simpler calculations are required.
  // virtual float interpolatePotentialLinearly(Vector3 pos) const;

  HOST DEVICE inline float interpolateDiffX(Vector3 pos, float w[3], float g1[4][4][4]) const {
    float a0, a1, a2, a3;

		// RBTODO parallelize loops?
		
    // Mix along x, taking the derivative.
    float g2[4][4];
    for (int iy = 0; iy < 4; iy++) {
      for (int iz = 0; iz < 4; iz++) {
				a3 = 0.5f*(-g1[0][iy][iz] + 3.0f*g1[1][iy][iz] - 3.0f*g1[2][iy][iz] + g1[3][iy][iz]);
				a2 = 0.5f*(2.0f*g1[0][iy][iz] - 5.0f*g1[1][iy][iz] + 4.0f*g1[2][iy][iz] - g1[3][iy][iz]);
				a1 = 0.5f*(-g1[0][iy][iz] + g1[2][iy][iz]);
				a0 = g1[1][iy][iz];

				//g2[iy][iz] = a3*w[0]*w[0]*w[0] + a2*w[0]*w[0] + a1*w[0] + a0;
				g2[iy][iz] = 3.0f*a3*w[0]*w[0] + 2.0f*a2*w[0] + a1;
      }
    }


    // Mix along y.
    float g3[4];
    for (int iz = 0; iz < 4; iz++) {
      a3 = 0.5f*(-g2[0][iz] + 3.0f*g2[1][iz] - 3.0f*g2[2][iz] + g2[3][iz]);
      a2 = 0.5f*(2.0f*g2[0][iz] - 5.0f*g2[1][iz] + 4.0f*g2[2][iz] - g2[3][iz]);
      a1 = 0.5f*(-g2[0][iz] + g2[2][iz]);
      a0 = g2[1][iz];
   
      g3[iz] = a3*w[1]*w[1]*w[1] + a2*w[1]*w[1] + a1*w[1] + a0;
    }

    // Mix along z.
    a3 = 0.5f*(-g3[0] + 3.0f*g3[1] - 3.0f*g3[2] + g3[3]);
    a2 = 0.5f*(2.0f*g3[0] - 5.0f*g3[1] + 4.0f*g3[2] - g3[3]);
    a1 = 0.5f*(-g3[0] + g3[2]);
    a0 = g3[1];
 
    float retval = -(a3*w[2]*w[2]*w[2] + a2*w[2]*w[2] + a1*w[2] + a0);
    return retval;
  }

  HOST DEVICE inline float interpolateDiffY(Vector3 pos, float w[3], float g1[4][4][4]) const {
    float a0, a1, a2, a3;
  
    // Mix along x, taking the derivative.
    float g2[4][4];
    for (int iy = 0; iy < 4; iy++) {
      for (int iz = 0; iz < 4; iz++) {
				a3 = 0.5f*(-g1[0][iy][iz] + 3.0f*g1[1][iy][iz] - 3.0f*g1[2][iy][iz] + g1[3][iy][iz]);
				a2 = 0.5f*(2.0f*g1[0][iy][iz] - 5.0f*g1[1][iy][iz] + 4.0f*g1[2][iy][iz] - g1[3][iy][iz]);
				a1 = 0.5f*(-g1[0][iy][iz] + g1[2][iy][iz]);
				a0 = g1[1][iy][iz];

				g2[iy][iz] = a3*w[0]*w[0]*w[0] + a2*w[0]*w[0] + a1*w[0] + a0;
      }
    }

    // Mix along y.
    float g3[4];
    for (int iz = 0; iz < 4; iz++) {
      a3 = 0.5f*(-g2[0][iz] + 3.0f*g2[1][iz] - 3.0f*g2[2][iz] + g2[3][iz]);
      a2 = 0.5f*(2.0f*g2[0][iz] - 5.0f*g2[1][iz] + 4.0f*g2[2][iz] - g2[3][iz]);
      a1 = 0.5f*(-g2[0][iz] + g2[2][iz]);
      a0 = g2[1][iz];
   
      //g3[iz] = a3*w[1]*w[1]*w[1] + a2*w[1]*w[1] + a1*w[1] + a0;
      g3[iz] = 3.0f*a3*w[1]*w[1] + 2.0f*a2*w[1] + a1;
    }

    // Mix along z.
    a3 = 0.5f*(-g3[0] + 3.0f*g3[1] - 3.0f*g3[2] + g3[3]);
    a2 = 0.5f*(2.0f*g3[0] - 5.0f*g3[1] + 4.0f*g3[2] - g3[3]);
    a1 = 0.5f*(-g3[0] + g3[2]);
    a0 = g3[1];

    return -(a3*w[2]*w[2]*w[2] + a2*w[2]*w[2] + a1*w[2] + a0);
  }

  HOST DEVICE inline float interpolateDiffZ(Vector3 pos, float w[3], float g1[4][4][4]) const {
    float a0, a1, a2, a3;
  
    // Mix along x, taking the derivative.
    float g2[4][4];
    for (int iy = 0; iy < 4; iy++) {
      for (int iz = 0; iz < 4; iz++) {
				a3 = 0.5f*(-g1[0][iy][iz] + 3.0f*g1[1][iy][iz] - 3.0f*g1[2][iy][iz] + g1[3][iy][iz]);
				a2 = 0.5f*(2.0f*g1[0][iy][iz] - 5.0f*g1[1][iy][iz] + 4.0f*g1[2][iy][iz] - g1[3][iy][iz]);
				a1 = 0.5f*(-g1[0][iy][iz] + g1[2][iy][iz]);
				a0 = g1[1][iy][iz];

				g2[iy][iz] = a3*w[0]*w[0]*w[0] + a2*w[0]*w[0] + a1*w[0] + a0;
      }
    }

    // Mix along y.
    float g3[4];
    for (int iz = 0; iz < 4; iz++) {
      a3 = 0.5f*(-g2[0][iz] + 3.0f*g2[1][iz] - 3.0f*g2[2][iz] + g2[3][iz]);
      a2 = 0.5f*(2.0f*g2[0][iz] - 5.0f*g2[1][iz] + 4.0f*g2[2][iz] - g2[3][iz]);
      a1 = 0.5f*(-g2[0][iz] + g2[2][iz]);
      a0 = g2[1][iz];
   
      g3[iz] = a3*w[1]*w[1]*w[1] + a2*w[1]*w[1] + a1*w[1] + a0;
    }

    // Mix along z.
    a3 = 0.5f*(-g3[0] + 3.0f*g3[1] - 3.0f*g3[2] + g3[3]);
    a2 = 0.5f*(2.0f*g3[0] - 5.0f*g3[1] + 4.0f*g3[2] - g3[3]);
    a1 = 0.5f*(-g3[0] + g3[2]);
    a0 = g3[1];

    return -(3.0f*a3*w[2]*w[2] + 2.0f*a2*w[2] + a1);
  }

  HOST DEVICE inline float interpolatePotential(const Vector3& pos) const {
    // Find the home node.
    Vector3 l = basisInv.transform(pos - origin);

		const int homeX = int(floor(l.x));
		const int homeY = int(floor(l.y));
		const int homeZ = int(floor(l.z));
		const float wx = l.x - homeX;
		const float wy = l.y - homeY;
		const float wz = l.z - homeZ;
		const float wx2 = wx*wx;
		const float wy2 = wy*wy;
		const float wz2 = wz*wz;

		float g3[4];
		for (int iz = 0; iz < 4; iz++) {
			float g2[4];
			const int jz = (iz + homeZ - 1);
			for (int iy = 0; iy < 4; iy++) {
				float v[4];
				const int jy = (iy + homeY - 1);
				for (int ix = 0; ix < 4; ix++) {
					const int jx = (ix + homeX - 1);
					const int ind = jz + jy*nz + jx*nz*ny;
					v[ix] = jz < 0 || jz >= nz || jy < 0 || jy >= ny || jx < 0 || jx >= nx ?
						0 : val[ind];
				}
				g2[iy] = 0.5f*(-v[0] + 3.0f*v[1] - 3.0f*v[2] + v[3])*wx2*wx +
					0.5f*(2.0f*v[0] - 5.0f*v[1] + 4.0f*v[2] - v[3])   *wx2  +
					0.5f*(-v[0] + v[2])                               *wx +
					v[1];
			}

			// Mix along y.
			g3[iz] = 0.5f*(-g2[0] + 3.0f*g2[1] - 3.0f*g2[2] + g2[3])*wy2*wy +
				0.5f*(2.0f*g2[0] - 5.0f*g2[1] + 4.0f*g2[2] - g2[3])   *wy2  +
				0.5f*(-g2[0] + g2[2])                                 *wy +
				g2[1];
		}
		// Mix along z.
		const float e = 0.5f*(-g3[0] + 3.0f*g3[1] - 3.0f*g3[2] + g3[3])*wz2*wz +
			0.5f*(2.0f*g3[0] - 5.0f*g3[1] + 4.0f*g3[2] - g3[3])          *wz2  +
			0.5f*(-g3[0] + g3[2])                                        *wz +
			g3[1];
    return e;
  }

  HOST DEVICE inline static int wrap(int i, int n) {
		while (i < 0) {
			//i %= n;
			i += n;
		}
		// The portion above allows i == n, so no else keyword.
		if (i >= n) {
			i %= n;
		} 
		return i;
	}
	HOST DEVICE float interpolatePotentialLinearly(const Vector3& pos) const {
		Vector3 f;
 		const Vector3 l = basisInv.transform(pos - origin);

		// Find the home node.
		const int homeX = int(floor(l.x));
		const int homeY = int(floor(l.y));
		const int homeZ = int(floor(l.z));

		const float wx = l.x - homeX;
		const float wy = l.y - homeY;	
		const float wz = l.z - homeZ;

		float v[2][2][2];
		for (int iz = 0; iz < 2; iz++) {
			int jz = (iz + homeZ);
			for (int iy = 0; iy < 2; iy++) {
				int jy = (iy + homeY);
				for (int ix = 0; ix < 2; ix++) {
					int jx = (ix + homeX);
					int ind = jz + jy*nz + jx*nz*ny;
					v[ix][iy][iz] = jz < 0 || jz >= nz || jy < 0 || jy >= ny || jx < 0 || jx >= nx ?
						0 : val[ind];
				}
			}
		}

		float g3[2];
		for (int iz = 0; iz < 2; iz++) {
			float g2[2];
			for (int iy = 0; iy < 2; iy++) {
				g2[iy] = wx * (v[1][iy][iz] - v[0][iy][iz]) + v[0][iy][iz];
			}
			// Mix along y.
			g3[iz] = wy * (g2[1] - g2[0]) + g2[0];
		}
		// Mix along z.
		float e = wz * (g3[1] - g3[0]) + g3[0];
		return e;
	}


	/** interpolateForce() to be used on CUDA Device **/
	DEVICE inline ForceEnergy interpolateForceD(const Vector3& pos) const {
		Vector3 f;
 		const Vector3 l = basisInv.transform(pos - origin);

		const int homeX = int(floor(l.x));
		const int homeY = int(floor(l.y));
		const int homeZ = int(floor(l.z));
		const float wx = l.x - homeX;
		const float wy = l.y - homeY;
		const float wz = l.z - homeZ;
		const float wx2 = wx*wx;

	/* f.x */
	float g3[3][4];
	for (int iz = 0; iz < 4; iz++) {
		float g2[2][4];
		const int jz = (iz + homeZ - 1);
		for (int iy = 0; iy < 4; iy++) {
			float v[4];
			const int jy = (iy + homeY - 1);
			for (int ix = 0; ix < 4; ix++) {
				const int jx = (ix + homeX - 1);
				const int ind = jz + jy*nz + jx*nz*ny;
				v[ix] = jz < 0 || jz >= nz || jy < 0 || jy >= ny || jx < 0 || jx >= nx ?
					0 : val[ind];
			}
			const float a3 = 0.5f*(-v[0] + 3.0f*v[1] - 3.0f*v[2] + v[3])*wx2;
			const float a2 = 0.5f*(2.0f*v[0] - 5.0f*v[1] + 4.0f*v[2] - v[3])*wx;
			const float a1 = 0.5f*(-v[0] + v[2]);
			g2[0][iy] = 3.0f*a3 + 2.0f*a2 + a1;				/* f.x (derivative) */
			g2[1][iy] = a3*wx + a2*wx + a1*wx + v[1]; /* f.y & f.z */
		}

		// Mix along y.
		{
			g3[0][iz] = 0.5f*(-g2[0][0] + 3.0f*g2[0][1] - 3.0f*g2[0][2] + g2[0][3])*wy*wy*wy +
				0.5f*(2.0f*g2[0][0] - 5.0f*g2[0][1] + 4.0f*g2[0][2] - g2[0][3])      *wy*wy +
				0.5f*(-g2[0][0] + g2[0][2])                                          *wy +
				g2[0][1];
		}

		{
			const float a3 = 0.5f*(-g2[1][0] + 3.0f*g2[1][1] - 3.0f*g2[1][2] + g2[1][3])*wy*wy;
			const float a2 = 0.5f*(2.0f*g2[1][0] - 5.0f*g2[1][1] + 4.0f*g2[1][2] - g2[1][3])*wy;
			const float a1 = 0.5f*(-g2[1][0] + g2[1][2]);
			g3[1][iz] = 3.0f*a3 + 2.0f*a2 + a1;						/* f.y */
			g3[2][iz] = a3*wy + a2*wy + a1*wy + g2[1][1]; /* f.z */
		}
	}

	// Mix along z.
	f.x = -0.5f*(-g3[0][0] + 3.0f*g3[0][1] - 3.0f*g3[0][2] + g3[0][3])*wz*wz*wz +
		-0.5f*(2.0f*g3[0][0] - 5.0f*g3[0][1] + 4.0f*g3[0][2] - g3[0][3])*wz*wz +
		-0.5f*(-g3[0][0] + g3[0][2])                                    *wz -
		g3[0][1];
	f.y = -0.5f*(-g3[1][0] + 3.0f*g3[1][1] - 3.0f*g3[1][2] + g3[1][3])*wz*wz*wz +
		-0.5f*(2.0f*g3[1][0] - 5.0f*g3[1][1] + 4.0f*g3[1][2] - g3[1][3])*wz*wz +
		-0.5f*(-g3[1][0] + g3[1][2])                                    *wz -
		g3[1][1];
	f.z = -1.5f*(-g3[2][0] + 3.0f*g3[2][1] - 3.0f*g3[2][2] + g3[2][3])*wz*wz -
		(2.0f*g3[2][0] - 5.0f*g3[2][1] + 4.0f*g3[2][2] - g3[2][3])      *wz -
		0.5f*(-g3[2][0] + g3[2][2]);
	float e = 0.5f*(-g3[2][0] + 3.0f*g3[2][1] - 3.0f*g3[2][2] + g3[2][3])*wz*wz*wz +
		0.5f*(2.0f*g3[2][0] - 5.0f*g3[2][1] + 4.0f*g3[2][2] - g3[2][3])    *wz*wz +
		0.5f*(-g3[2][0] + g3[2][2])                                        *wz +
		g3[2][1];

	f = basisInv.transpose().transform(f);
	return ForceEnergy(f,e);
	
	}

	template <BoundaryCondition bc>
	DEVICE inline ForceEnergy interpolateForceDLinearly(const Vector3& pos) const {
 		const Vector3 l = basisInv.transform(pos - origin);

		// Find the home node.
		const int homeX = int(floor(l.x));
		const int homeY = int(floor(l.y));
		const int homeZ = int(floor(l.z));

		const float wx = l.x - homeX;
		const float wy = l.y - homeY;	
		const float wz = l.z - homeZ;

		if (bc == neumann) {
		    if (homeX < -1 || homeX >= nx+1 ||
		    	homeY < -1 || homeY >= ny+1 ||
		    	homeZ < -1 || homeZ >= nz+1)
		    	return ForceEnergy();
		}

		float v[2][2][2];
		for (int iz = 0; iz < 2; iz++) {
			int jz = (iz + homeZ);
			for (int iy = 0; iy < 2; iy++) {
				int jy = (iy + homeY);
				for (int ix = 0; ix < 2; ix++) {
					int jx = (ix + homeX);
					int ind;
					switch (bc) {
					case dirichlet:
					    ind = jz + jy*nz + jx*nz*ny;
					    v[ix][iy][iz] = 
						jz < 0 || jz >= nz ||
						jy < 0 || jy >= ny ||
						jx < 0 || jx >= nx ?
						0 : val[ind];
					    break;
					case neumann:
					    ind =
						(jz < 0 ? 0 : jz >= nz ? nz-1 : jz) +
						(jy < 0 ? 0 : jy >= ny ? ny-1 : jy)*nz + 
						(jx < 0 ? 0 : jx >= nx ? nx-1 : jx)*nz*ny;
					    v[ix][iy][iz] = val[ind];
					    break;
					case periodic:
					    ind =
						(jz < 0 ? nz-1 : jz >= nz ? 0 : jz) +
						(jy < 0 ? ny-1 : jy >= ny ? 0 : jy)*nz + 
						(jx < 0 ? nx-1 : jx >= nx ? 0 : jx)*nz*ny;
					    v[ix][iy][iz] = val[ind];
					    break;
					}
				}
			}
		}

		float g3[3][2];
		for (int iz = 0; iz < 2; iz++) {
			float g2[2][2];
			for (int iy = 0; iy < 2; iy++) {
				g2[0][iy] = (v[1][iy][iz] - v[0][iy][iz]); /* f.x */
				g2[1][iy] = wx * (v[1][iy][iz] - v[0][iy][iz]) + v[0][iy][iz]; /* f.y & f.z */
			}
			// Mix along y.
			g3[0][iz] = wy * (g2[0][1] - g2[0][0]) + g2[0][0];
			g3[1][iz] = (g2[1][1] - g2[1][0]);
			g3[2][iz] = wy * (g2[1][1] - g2[1][0]) + g2[1][0];
		}
		// Mix along z.
		Vector3 f;
		f.x = -(wz * (g3[0][1] - g3[0][0]) + g3[0][0]);
		f.y = -(wz * (g3[1][1] - g3[1][0]) + g3[1][0]);
		f.z = -      (g3[2][1] - g3[2][0]);

		f = basisInv.transpose().transform(f);
		float e = wz * (g3[2][1] - g3[2][0]) + g3[2][0];
		return ForceEnergy(f,e);
	}
        DEVICE inline ForceEnergy interpolateForceDnamd(const Vector3& pos) const {
                Vector3 f;
                const Vector3 l = basisInv.transform(pos - origin);

                const int homeX = int(floor(l.x));
                const int homeY = int(floor(l.y));
                const int homeZ = int(floor(l.z));
                const float wx = l.x - homeX;
                const float wy = l.y - homeY;
                const float wz = l.z - homeZ;

                Vector3 dg = Vector3(wx,wy,wz);

                int inds[3];
                inds[0] = homeX;
                inds[1] = homeY;
                inds[2] = homeZ;

                // TODO: handle edges

                // Compute b
                                   float b[64];    // Matrix of values at 8 box corners
                compute_b(b, inds);

                // Compute a
                                   float a[64];
                compute_a(a, b);

                // Calculate powers of x, y, z for later use
                                   // e.g. x[2] = x^2
                                                      float x[4], y[4], z[4];
                x[0] = 1; y[0] = 1; z[0] = 1;
                for (int j = 1; j < 4; j++) {
                    x[j] = x[j-1] * dg.x;
                    y[j] = y[j-1] * dg.y;
                    z[j] = z[j-1] * dg.z;
                }
                float e = compute_V(a, x, y, z);
                f = compute_dV(a, x, y, z);

                f = basisInv.transpose().transform(f);
                return ForceEnergy(f,e);
        }
        DEVICE inline float compute_V(float *a, float *x, float *y, float *z) const
        {
            float V = 0.0;
            long int ind = 0;
            for (int l = 0; l < 4; l++) {
                for (int k = 0; k < 4; k++) {
                    for (int j = 0; j < 4; j++) {
                        V += a[ind] * x[j] * y[k] * z[l];
                        ind++;
                    }
                }
            }
            return V;
        }
        DEVICE inline Vector3 compute_dV(float *a, float *x, float *y, float *z) const
        {
            Vector3 dV = Vector3(0.0f);
            long int ind = 0;
            for (int l = 0; l < 4; l++) {
                for (int k = 0; k < 4; k++) {
                    for (int j = 0; j < 4; j++) {
                        if (j > 0) dV.x += a[ind] * j * x[j-1] * y[k]   * z[l];         // dV/dx
                        if (k > 0) dV.y += a[ind] * k * x[j]   * y[k-1] * z[l];         // dV/dy
                        if (l > 0) dV.z += a[ind] * l * x[j]   * y[k]   * z[l-1];       // dV/dz
                        ind++;
                    }
                }
            }
            return dV*(-1.f);
        }
        DEVICE inline void compute_a(float *a, float *b) const
        {
            // Static sparse 64x64 matrix times vector ... nicer looking way than this?
                           a[0] = b[0];
            a[1] = b[8];
            a[2] = -3*b[0] + 3*b[1] - 2*b[8] - b[9];
            a[3] = 2*b[0] - 2*b[1] + b[8] + b[9];
            a[4] = b[16];
            a[5] = b[32];
            a[6] = -3*b[16] + 3*b[17] - 2*b[32] - b[33];
            a[7] = 2*b[16] - 2*b[17] + b[32] + b[33];
            a[8] = -3*b[0] + 3*b[2] - 2*b[16] - b[18];
            a[9] = -3*b[8] + 3*b[10] - 2*b[32] - b[34];
            a[10] = 9*b[0] - 9*b[1] - 9*b[2] + 9*b[3] + 6*b[8] + 3*b[9] - 6*b[10] - 3*b[11]
                + 6*b[16] - 6*b[17] + 3*b[18] - 3*b[19] + 4*b[32] + 2*b[33] + 2*b[34] + b[35];
            a[11] = -6*b[0] + 6*b[1] + 6*b[2] - 6*b[3] - 3*b[8] - 3*b[9] + 3*b[10] + 3*b[11]
                - 4*b[16] + 4*b[17] - 2*b[18] + 2*b[19] - 2*b[32] - 2*b[33] - b[34] - b[35];
            a[12] = 2*b[0] - 2*b[2] + b[16] + b[18];
            a[13] = 2*b[8] - 2*b[10] + b[32] + b[34];
            a[14] = -6*b[0] + 6*b[1] + 6*b[2] - 6*b[3] - 4*b[8] - 2*b[9] + 4*b[10] + 2*b[11]
                - 3*b[16] + 3*b[17] - 3*b[18] + 3*b[19] - 2*b[32] - b[33] - 2*b[34] - b[35];
            a[15] = 4*b[0] - 4*b[1] - 4*b[2] + 4*b[3] + 2*b[8] + 2*b[9] - 2*b[10] - 2*b[11]
                + 2*b[16] - 2*b[17] + 2*b[18] - 2*b[19] + b[32] + b[33] + b[34] + b[35];
            a[16] = b[24];
            a[17] = b[40];
            a[18] = -3*b[24] + 3*b[25] - 2*b[40] - b[41];
            a[19] = 2*b[24] - 2*b[25] + b[40] + b[41];
            a[20] = b[48];
            a[21] = b[56];
            a[22] = -3*b[48] + 3*b[49] - 2*b[56] - b[57];
            a[23] = 2*b[48] - 2*b[49] + b[56] + b[57];
            a[24] = -3*b[24] + 3*b[26] - 2*b[48] - b[50];
            a[25] = -3*b[40] + 3*b[42] - 2*b[56] - b[58];
            a[26] = 9*b[24] - 9*b[25] - 9*b[26] + 9*b[27] + 6*b[40] + 3*b[41] - 6*b[42] - 3*b[43]
                + 6*b[48] - 6*b[49] + 3*b[50] - 3*b[51] + 4*b[56] + 2*b[57] + 2*b[58] + b[59];
            a[27] = -6*b[24] + 6*b[25] + 6*b[26] - 6*b[27] - 3*b[40] - 3*b[41] + 3*b[42] + 3*b[43]
                - 4*b[48] + 4*b[49] - 2*b[50] + 2*b[51] - 2*b[56] - 2*b[57] - b[58] - b[59];
            a[28] = 2*b[24] - 2*b[26] + b[48] + b[50];
            a[29] = 2*b[40] - 2*b[42] + b[56] + b[58];
            a[30] = -6*b[24] + 6*b[25] + 6*b[26] - 6*b[27] - 4*b[40] - 2*b[41] + 4*b[42] + 2*b[43]
                - 3*b[48] + 3*b[49] - 3*b[50] + 3*b[51] - 2*b[56] - b[57] - 2*b[58] - b[59];
            a[31] = 4*b[24] - 4*b[25] - 4*b[26] + 4*b[27] + 2*b[40] + 2*b[41] - 2*b[42] - 2*b[43]
                + 2*b[48] - 2*b[49] + 2*b[50] - 2*b[51] + b[56] + b[57] + b[58] + b[59];
            a[32] = -3*b[0] + 3*b[4] - 2*b[24] - b[28];
            a[33] = -3*b[8] + 3*b[12] - 2*b[40] - b[44];
            a[34] = 9*b[0] - 9*b[1] - 9*b[4] + 9*b[5] + 6*b[8] + 3*b[9] - 6*b[12] - 3*b[13]
                + 6*b[24] - 6*b[25] + 3*b[28] - 3*b[29] + 4*b[40] + 2*b[41] + 2*b[44] + b[45];
            a[35] = -6*b[0] + 6*b[1] + 6*b[4] - 6*b[5] - 3*b[8] - 3*b[9] + 3*b[12] + 3*b[13]
                - 4*b[24] + 4*b[25] - 2*b[28] + 2*b[29] - 2*b[40] - 2*b[41] - b[44] - b[45];
            a[36] = -3*b[16] + 3*b[20] - 2*b[48] - b[52];
            a[37] = -3*b[32] + 3*b[36] - 2*b[56] - b[60];
            a[38] = 9*b[16] - 9*b[17] - 9*b[20] + 9*b[21] + 6*b[32] + 3*b[33] - 6*b[36] - 3*b[37]
                + 6*b[48] - 6*b[49] + 3*b[52] - 3*b[53] + 4*b[56] + 2*b[57] + 2*b[60] + b[61];
            a[39] = -6*b[16] + 6*b[17] + 6*b[20] - 6*b[21] - 3*b[32] - 3*b[33] + 3*b[36] + 3*b[37]
                - 4*b[48] + 4*b[49] - 2*b[52] + 2*b[53] - 2*b[56] - 2*b[57] - b[60] - b[61];
            a[40] = 9*b[0] - 9*b[2] - 9*b[4] + 9*b[6] + 6*b[16] + 3*b[18] - 6*b[20] - 3*b[22]
                + 6*b[24] - 6*b[26] + 3*b[28] - 3*b[30] + 4*b[48] + 2*b[50] + 2*b[52] + b[54];
            a[41] = 9*b[8] - 9*b[10] - 9*b[12] + 9*b[14] + 6*b[32] + 3*b[34] - 6*b[36] - 3*b[38]
                + 6*b[40] - 6*b[42] + 3*b[44] - 3*b[46] + 4*b[56] + 2*b[58] + 2*b[60] + b[62];
            a[42] = -27*b[0] + 27*b[1] + 27*b[2] - 27*b[3] + 27*b[4] - 27*b[5] - 27*b[6] + 27*b[7]
                - 18*b[8] - 9*b[9] + 18*b[10] + 9*b[11] + 18*b[12] + 9*b[13] - 18*b[14] - 9*b[15]
                - 18*b[16] + 18*b[17] - 9*b[18] + 9*b[19] + 18*b[20] - 18*b[21] + 9*b[22] - 9*b[23]
                - 18*b[24] + 18*b[25] + 18*b[26] - 18*b[27] - 9*b[28] + 9*b[29] + 9*b[30] - 9*b[31]
                - 12*b[32] - 6*b[33] - 6*b[34] - 3*b[35] + 12*b[36] + 6*b[37] + 6*b[38] + 3*b[39]
                - 12*b[40] - 6*b[41] + 12*b[42] + 6*b[43] - 6*b[44] - 3*b[45] + 6*b[46] + 3*b[47]
                - 12*b[48] + 12*b[49] - 6*b[50] + 6*b[51] - 6*b[52] + 6*b[53] - 3*b[54] + 3*b[55]
                - 8*b[56] - 4*b[57] - 4*b[58] - 2*b[59] - 4*b[60] - 2*b[61] - 2*b[62] - b[63];
            a[43] = 18*b[0] - 18*b[1] - 18*b[2] + 18*b[3] - 18*b[4] + 18*b[5] + 18*b[6] - 18*b[7]
                + 9*b[8] + 9*b[9] - 9*b[10] - 9*b[11] - 9*b[12] - 9*b[13] + 9*b[14] + 9*b[15]
                + 12*b[16] - 12*b[17] + 6*b[18] - 6*b[19] - 12*b[20] + 12*b[21] - 6*b[22] + 6*b[23]
                + 12*b[24] - 12*b[25] - 12*b[26] + 12*b[27] + 6*b[28] - 6*b[29] - 6*b[30] + 6*b[31]
                + 6*b[32] + 6*b[33] + 3*b[34] + 3*b[35] - 6*b[36] - 6*b[37] - 3*b[38] - 3*b[39]
                + 6*b[40] + 6*b[41] - 6*b[42] - 6*b[43] + 3*b[44] + 3*b[45] - 3*b[46] - 3*b[47]
                + 8*b[48] - 8*b[49] + 4*b[50] - 4*b[51] + 4*b[52] - 4*b[53] + 2*b[54] - 2*b[55]
                + 4*b[56] + 4*b[57] + 2*b[58] + 2*b[59] + 2*b[60] + 2*b[61] + b[62] + b[63];
            a[44] = -6*b[0] + 6*b[2] + 6*b[4] - 6*b[6] - 3*b[16] - 3*b[18] + 3*b[20] + 3*b[22]
                - 4*b[24] + 4*b[26] - 2*b[28] + 2*b[30] - 2*b[48] - 2*b[50] - b[52] - b[54];
            a[45] = -6*b[8] + 6*b[10] + 6*b[12] - 6*b[14] - 3*b[32] - 3*b[34] + 3*b[36] + 3*b[38]
                - 4*b[40] + 4*b[42] - 2*b[44] + 2*b[46] - 2*b[56] - 2*b[58] - b[60] - b[62];
            a[46] = 18*b[0] - 18*b[1] - 18*b[2] + 18*b[3] - 18*b[4] + 18*b[5] + 18*b[6] - 18*b[7]
                + 12*b[8] + 6*b[9] - 12*b[10] - 6*b[11] - 12*b[12] - 6*b[13] + 12*b[14] + 6*b[15]
                + 9*b[16] - 9*b[17] + 9*b[18] - 9*b[19] - 9*b[20] + 9*b[21] - 9*b[22] + 9*b[23]
                + 12*b[24] - 12*b[25] - 12*b[26] + 12*b[27] + 6*b[28] - 6*b[29] - 6*b[30] + 6*b[31]
                + 6*b[32] + 3*b[33] + 6*b[34] + 3*b[35] - 6*b[36] - 3*b[37] - 6*b[38] - 3*b[39]
                + 8*b[40] + 4*b[41] - 8*b[42] - 4*b[43] + 4*b[44] + 2*b[45] - 4*b[46] - 2*b[47]
                + 6*b[48] - 6*b[49] + 6*b[50] - 6*b[51] + 3*b[52] - 3*b[53] + 3*b[54] - 3*b[55]
                + 4*b[56] + 2*b[57] + 4*b[58] + 2*b[59] + 2*b[60] + b[61] + 2*b[62] + b[63];
            a[47] = -12*b[0] + 12*b[1] + 12*b[2] - 12*b[3] + 12*b[4] - 12*b[5] - 12*b[6] + 12*b[7]
                - 6*b[8] - 6*b[9] + 6*b[10] + 6*b[11] + 6*b[12] + 6*b[13] - 6*b[14] - 6*b[15]
                - 6*b[16] + 6*b[17] - 6*b[18] + 6*b[19] + 6*b[20] - 6*b[21] + 6*b[22] - 6*b[23]
                - 8*b[24] + 8*b[25] + 8*b[26] - 8*b[27] - 4*b[28] + 4*b[29] + 4*b[30] - 4*b[31]
                - 3*b[32] - 3*b[33] - 3*b[34] - 3*b[35] + 3*b[36] + 3*b[37] + 3*b[38] + 3*b[39]
                - 4*b[40] - 4*b[41] + 4*b[42] + 4*b[43] - 2*b[44] - 2*b[45] + 2*b[46] + 2*b[47]
                - 4*b[48] + 4*b[49] - 4*b[50] + 4*b[51] - 2*b[52] + 2*b[53] - 2*b[54] + 2*b[55]
                - 2*b[56] - 2*b[57] - 2*b[58] - 2*b[59] - b[60] - b[61] - b[62] - b[63];
            a[48] = 2*b[0] - 2*b[4] + b[24] + b[28];
            a[49] = 2*b[8] - 2*b[12] + b[40] + b[44];
            a[50] = -6*b[0] + 6*b[1] + 6*b[4] - 6*b[5] - 4*b[8] - 2*b[9] + 4*b[12] + 2*b[13]
                - 3*b[24] + 3*b[25] - 3*b[28] + 3*b[29] - 2*b[40] - b[41] - 2*b[44] - b[45];
            a[51] = 4*b[0] - 4*b[1] - 4*b[4] + 4*b[5] + 2*b[8] + 2*b[9] - 2*b[12] - 2*b[13]
                + 2*b[24] - 2*b[25] + 2*b[28] - 2*b[29] + b[40] + b[41] + b[44] + b[45];
            a[52] = 2*b[16] - 2*b[20] + b[48] + b[52];
            a[53] = 2*b[32] - 2*b[36] + b[56] + b[60];
            a[54] = -6*b[16] + 6*b[17] + 6*b[20] - 6*b[21] - 4*b[32] - 2*b[33] + 4*b[36] + 2*b[37]
                - 3*b[48] + 3*b[49] - 3*b[52] + 3*b[53] - 2*b[56] - b[57] - 2*b[60] - b[61];
            a[55] = 4*b[16] - 4*b[17] - 4*b[20] + 4*b[21] + 2*b[32] + 2*b[33] - 2*b[36] - 2*b[37]
                + 2*b[48] - 2*b[49] + 2*b[52] - 2*b[53] + b[56] + b[57] + b[60] + b[61];
            a[56] = -6*b[0] + 6*b[2] + 6*b[4] - 6*b[6] - 4*b[16] - 2*b[18] + 4*b[20] + 2*b[22]
                - 3*b[24] + 3*b[26] - 3*b[28] + 3*b[30] - 2*b[48] - b[50] - 2*b[52] - b[54];
            a[57] = -6*b[8] + 6*b[10] + 6*b[12] - 6*b[14] - 4*b[32] - 2*b[34] + 4*b[36] + 2*b[38]
                - 3*b[40] + 3*b[42] - 3*b[44] + 3*b[46] - 2*b[56] - b[58] - 2*b[60] - b[62];
            a[58] = 18*b[0] - 18*b[1] - 18*b[2] + 18*b[3] - 18*b[4] + 18*b[5] + 18*b[6] - 18*b[7]
                + 12*b[8] + 6*b[9] - 12*b[10] - 6*b[11] - 12*b[12] - 6*b[13] + 12*b[14] + 6*b[15]
                + 12*b[16] - 12*b[17] + 6*b[18] - 6*b[19] - 12*b[20] + 12*b[21] - 6*b[22] + 6*b[23]
                + 9*b[24] - 9*b[25] - 9*b[26] + 9*b[27] + 9*b[28] - 9*b[29] - 9*b[30] + 9*b[31]
                + 8*b[32] + 4*b[33] + 4*b[34] + 2*b[35] - 8*b[36] - 4*b[37] - 4*b[38] - 2*b[39]
                + 6*b[40] + 3*b[41] - 6*b[42] - 3*b[43] + 6*b[44] + 3*b[45] - 6*b[46] - 3*b[47]
                + 6*b[48] - 6*b[49] + 3*b[50] - 3*b[51] + 6*b[52] - 6*b[53] + 3*b[54] - 3*b[55]
                + 4*b[56] + 2*b[57] + 2*b[58] + b[59] + 4*b[60] + 2*b[61] + 2*b[62] + b[63];
            a[59] = -12*b[0] + 12*b[1] + 12*b[2] - 12*b[3] + 12*b[4] - 12*b[5] - 12*b[6] + 12*b[7]
                - 6*b[8] - 6*b[9] + 6*b[10] + 6*b[11] + 6*b[12] + 6*b[13] - 6*b[14] - 6*b[15]
                - 8*b[16] + 8*b[17] - 4*b[18] + 4*b[19] + 8*b[20] - 8*b[21] + 4*b[22] - 4*b[23]
                - 6*b[24] + 6*b[25] + 6*b[26] - 6*b[27] - 6*b[28] + 6*b[29] + 6*b[30] - 6*b[31]
                - 4*b[32] - 4*b[33] - 2*b[34] - 2*b[35] + 4*b[36] + 4*b[37] + 2*b[38] + 2*b[39]
                - 3*b[40] - 3*b[41] + 3*b[42] + 3*b[43] - 3*b[44] - 3*b[45] + 3*b[46] + 3*b[47]
                - 4*b[48] + 4*b[49] - 2*b[50] + 2*b[51] - 4*b[52] + 4*b[53] - 2*b[54] + 2*b[55]
                - 2*b[56] - 2*b[57] - b[58] - b[59] - 2*b[60] - 2*b[61] - b[62] - b[63];
            a[60] = 4*b[0] - 4*b[2] - 4*b[4] + 4*b[6] + 2*b[16] + 2*b[18] - 2*b[20] - 2*b[22]
                + 2*b[24] - 2*b[26] + 2*b[28] - 2*b[30] + b[48] + b[50] + b[52] + b[54];
            a[61] = 4*b[8] - 4*b[10] - 4*b[12] + 4*b[14] + 2*b[32] + 2*b[34] - 2*b[36] - 2*b[38]
                + 2*b[40] - 2*b[42] + 2*b[44] - 2*b[46] + b[56] + b[58] + b[60] + b[62];
            a[62] = -12*b[0] + 12*b[1] + 12*b[2] - 12*b[3] + 12*b[4] - 12*b[5] - 12*b[6] + 12*b[7]
                - 8*b[8] - 4*b[9] + 8*b[10] + 4*b[11] + 8*b[12] + 4*b[13] - 8*b[14] - 4*b[15]
                - 6*b[16] + 6*b[17] - 6*b[18] + 6*b[19] + 6*b[20] - 6*b[21] + 6*b[22] - 6*b[23]
                - 6*b[24] + 6*b[25] + 6*b[26] - 6*b[27] - 6*b[28] + 6*b[29] + 6*b[30] - 6*b[31]
                - 4*b[32] - 2*b[33] - 4*b[34] - 2*b[35] + 4*b[36] + 2*b[37] + 4*b[38] + 2*b[39]
                - 4*b[40] - 2*b[41] + 4*b[42] + 2*b[43] - 4*b[44] - 2*b[45] + 4*b[46] + 2*b[47]
                - 3*b[48] + 3*b[49] - 3*b[50] + 3*b[51] - 3*b[52] + 3*b[53] - 3*b[54] + 3*b[55]
                - 2*b[56] - b[57] - 2*b[58] - b[59] - 2*b[60] - b[61] - 2*b[62] - b[63];
            a[63] = 8*b[0] - 8*b[1] - 8*b[2] + 8*b[3] - 8*b[4] + 8*b[5] + 8*b[6] - 8*b[7]
                + 4*b[8] + 4*b[9] - 4*b[10] - 4*b[11] - 4*b[12] - 4*b[13] + 4*b[14] + 4*b[15]
                + 4*b[16] - 4*b[17] + 4*b[18] - 4*b[19] - 4*b[20] + 4*b[21] - 4*b[22] + 4*b[23]
                + 4*b[24] - 4*b[25] - 4*b[26] + 4*b[27] + 4*b[28] - 4*b[29] - 4*b[30] + 4*b[31]
                + 2*b[32] + 2*b[33] + 2*b[34] + 2*b[35] - 2*b[36] - 2*b[37] - 2*b[38] - 2*b[39]
                + 2*b[40] + 2*b[41] - 2*b[42] - 2*b[43] + 2*b[44] + 2*b[45] - 2*b[46] - 2*b[47]
                + 2*b[48] - 2*b[49] + 2*b[50] - 2*b[51] + 2*b[52] - 2*b[53] + 2*b[54] - 2*b[55]
                + b[56] + b[57] + b[58] + b[59] + b[60] + b[61] + b[62] + b[63];
        }
        DEVICE void compute_b(float * __restrict__ b, int * __restrict__ inds) const
        {
            int k[3];
            k[0] = nx;
            k[1] = ny;
            k[2] = nz;

            int inds2[3] = {0,0,0};

            for (int i0 = 0; i0 < 8; i0++) {
                inds2[0] = 0;
                inds2[1] = 0;
                inds2[2] = 0;

                /* printf("%d\n", inds2[0]); */
                /* printf("%d\n", inds2[1]); */
                /* printf("%d\n", inds2[2]); */

                bool zero_derivs = false;

                int bit = 1;    // bit = 2^i1 in the below loop
                for (int i1 = 0; i1 < 3; i1++) {
                    inds2[i1] = (inds[i1] + ((i0 & bit) ? 1 : 0)) % k[i1];
                    bit <<= 1;  // i.e. multiply by 2
                }
                int d_lo[3] = {1, 1, 1};
                float voffs[3] = {0.0f, 0.0f, 0.0f};
                float dscales[3] = {0.5, 0.5, 0.5};

                for (int i1 = 0; i1 < 3; i1++) {
                    if (inds2[i1] == 0) {
                        zero_derivs = true;
                    }
                    else if (inds2[i1] == k[i1]-1) {
                        zero_derivs = true;
                    }
                    else {
                        // printf("%d\n",i1);
                        voffs[i1] = 0.0;
                    }
                }

                // V
                b[i0] = getValue(inds2[0],inds2[1],inds2[2]);

                if (zero_derivs) {
                    b[8+i0] = 0.0;
                    b[16+i0] = 0.0;
                    b[24+i0] = 0.0;
                    b[32+i0] = 0.0;
                    b[40+i0] = 0.0;
                    b[48+i0] = 0.0;
                    b[56+i0] = 0.0;
                } else {
                    b[8+i0]  = dscales[0] * (getValue(inds2[0]+1,inds2[1],inds2[2]) - getValue(inds2[0]-d_lo[0],inds2[1],inds2[2]) + voffs[0]);
                    b[16+i0] = dscales[1] * (getValue(inds2[0],inds2[1]+1,inds2[2]) - getValue(inds2[0],inds2[1]-d_lo[1],inds2[2]) + voffs[1]);
                    b[24+i0] = dscales[2] * (getValue(inds2[0],inds2[1],inds2[2]+1) - getValue(inds2[0],inds2[1],inds2[2]-d_lo[2]) + voffs[2]);
                    b[32+i0] = dscales[0] * dscales[1]
                        * (getValue(inds2[0]+1,inds2[1]+1,inds2[2]) - getValue(inds2[0]-d_lo[0],inds2[1]+1,inds2[2]) -
                           getValue(inds2[0]+1,inds2[1]-d_lo[1],inds2[2]) + getValue(inds2[0]-d_lo[0],inds2[1]-d_lo[1],inds2[2]));
                    b[40+i0] = dscales[0] * dscales[2]
                        * (getValue(inds2[0]+1,inds2[1],inds2[2]+1) - getValue(inds2[0]-d_lo[0],inds2[1],inds2[2]+1) -
                           getValue(inds2[0]+1,inds2[1],inds2[2]-d_lo[2]) + getValue(inds2[0]-d_lo[0],inds2[1],inds2[2]-d_lo[2]));
                    b[48+i0] = dscales[1] * dscales[2]
                        * (getValue(inds2[0],inds2[1]+1,inds2[2]+1) - getValue(inds2[0],inds2[1]-d_lo[1],inds2[2]+1) -
                           getValue(inds2[0],inds2[1]+1,inds2[2]-d_lo[2]) + getValue(inds2[0],inds2[1]-d_lo[1],inds2[2]-d_lo[2]));
                    b[56+i0] = dscales[0] * dscales[1] * dscales[2]
                        * (getValue(inds2[0]+1,inds2[1]+1,inds2[2]+1) - getValue(inds2[0]+1,inds2[1]+1,inds2[2]-d_lo[2]) -
                           getValue(inds2[0]+1,inds2[1]-d_lo[1],inds2[2]+1) - getValue(inds2[0]-d_lo[0],inds2[1]+1,inds2[2]+1) +
                           getValue(inds2[0]+1,inds2[1]-d_lo[1],inds2[2]-d_lo[2]) + getValue(inds2[0]-d_lo[0],inds2[1]+1,inds2[2]-d_lo[2]) +
                           getValue(inds2[0]-d_lo[0],inds2[1]-d_lo[1],inds2[2]+1) - getValue(inds2[0]-d_lo[0],inds2[1]-d_lo[1],inds2[2]-d_lo[2]));
                }
            }
        }
        HOST DEVICE inline float getValue(int ix, int iy, int iz) const {
        
            if (ix < 0 || ix >= nx) return 0.0f;
            if (iy < 0 || iy >= ny) return 0.0f;
            if (iz < 0 || iz >= nz) return 0.0f;

            int j = iz + iy*nz + ix*ny*nz;
            return val[j];
        /*
           if(ix < 0) ix = 0;
           else if(ix >= nx) ix = nx -1;

           if(iy < 0) iy = 0;
           else if(iy >= ny) iy = ny-1;

           if(iz < 0) iz = 0;
           else if(iz >= nz) iz = nz-1;

           int j = iz + nz * (iy + ny * ix);
           return val[j];*/
        }

  inline virtual Vector3 interpolateForce(Vector3 pos) const {
		Vector3 f;
 		Vector3 l = basisInv.transform(pos - origin);
		int homeX = int(floor(l.x));
		int homeY = int(floor(l.y));
		int homeZ = int(floor(l.z));
		// Get the array jumps with shifted indices.
		int jump[3];
		jump[0] = nz*ny;
		jump[1] = nz;
		jump[2] = 1;
		// Shift the indices in the home array.
		int home[3];
		home[0] = homeX;
		home[1] = homeY;
		home[2] = homeZ;

		// Shift the indices in the grid dimensions.
		int g[3];
		g[0] = nx;
		g[1] = ny;
		g[2] = nz;

		// Get the interpolation coordinates.
		float w[3];
		w[0] = l.x - homeX;
		w[1] = l.y - homeY;
		w[2] = l.z - homeZ;
		// Find the values at the neighbors.
		float g1[4][4][4];
		//RBTODO parallelize?
		for (int ix = 0; ix < 4; ix++) {
			for (int iy = 0; iy < 4; iy++) {
				for (int iz = 0; iz < 4; iz++) {
	  			// Wrap around the periodic boundaries. 
					int jx = ix-1 + home[0];
					jx = wrap(jx, g[0]);
					int jy = iy-1 + home[1];
					jy = wrap(jy, g[1]);
					int jz = iz-1 + home[2];
					jz = wrap(jz, g[2]);
					int ind = jz*jump[2] + jy*jump[1] + jx*jump[0];
					g1[ix][iy][iz] = val[ind];
				}
			}
		}  
		f.x = interpolateDiffX(pos, w, g1);
		f.y = interpolateDiffY(pos, w, g1);
		f.z = interpolateDiffZ(pos, w, g1);
		Vector3 f1 = basisInv.transpose().transform(f);
		return f1;
	}

  // Wrap coordinate: 0 <= x < l
  HOST DEVICE   inline int quotient(float x, float l) const {
#if __CUDA_ARCH__ > 0
	  return int(floorf( __fdividef(x,l) ));
#else
	  return int(floor(x/l));
#endif
  }

  HOST DEVICE inline float wrapFloat(float x, float l) const {
		int image = int(floor(x/l));
		x -= image*l;
		return x;
  }
  
  // Wrap distance: -0.5*l <= x < 0.5*l
  HOST DEVICE static inline float wrapDiff(float x, float l) {
		int image = int(floor(x/l));
		x -= image*l;
		if (x >= 0.5f * l)
			x -= l;
		return x;
  }

  // TODO: implement simpler approach when basis is diagonal
  // TODO: implement device version using __fdividef and floorf
  // TODO: make BaseGrid an abstract class that diagGrid and nonDiagGrid inherit from 
  // TODO: test wrap and wrapDiff for non-diagonal basis
  // Wrap vector, 0 <= x < lx  &&  0 <= y < ly  &&  0 <= z < lz
  HOST DEVICE inline Vector3 wrap(Vector3 r) const {
	Vector3 l = basisInv.transform(r - origin);
	if ( basis.isDiagonal() ) {
		r = r - Vector3(quotient(l.x,nx) * nx*basis.exx,
						quotient(l.y,ny) * ny*basis.eyy,
						quotient(l.z,nz) * nz*basis.ezz);
	} else {
		r = r - quotient(l.x,nx) * nx*basis.ex();
		r = r - quotient(l.y,ny) * ny*basis.ey();
		r = r - quotient(l.z,nz) * nz*basis.ez();
	}
	return r;
  }

  HOST DEVICE inline Vector3 wrapDiff(Vector3 r) const {
	Vector3 l = basisInv.transform(r);
	if ( basis.isDiagonal() ) {
		r = r - Vector3(quotient(l.x+0.5f*nx,nx) * nx*basis.exx,
						quotient(l.y+0.5f*ny,ny) * ny*basis.eyy,
						quotient(l.z+0.5f*nz,nz) * nz*basis.ezz);
	} else {
		r = r - quotient(l.x+0.5f*nx,nx) * nx*basis.ex();
		r = r - quotient(l.y+0.5f*ny,ny) * ny*basis.ey();
		r = r - quotient(l.z+0.5f*nz,nz) * nz*basis.ez();
	}
	return r;
  }
  
  // Wrap vector distance, -0.5*l <= x < 0.5*l  && ...
  /* HOST DEVICE inline Vector3 wrapDiff(Vector3 r) const { */
  /*   Vector3 l = basisInv.transform(r); */
  /*   l.x = wrapDiff(l.x, nx); */
  /*   l.y = wrapDiff(l.y, ny); */
  /*   l.z = wrapDiff(l.z, nz); */
  /*   return basis.transform(l); */
  /* } */
  HOST DEVICE inline Vector3 wrapDiffOrig(Vector3 r) const {
    Vector3 l = basisInv.transform(r);
    l.x = wrapDiff(l.x, nx);
    l.y = wrapDiff(l.y, ny);
    l.z = wrapDiff(l.z, nz);
    return basis.transform(l);
  }
  Vector3 wrapDiffNearest(Vector3 r) const;

  // Includes the home node.
  // indexBuffer must have a size of at least 27.
  void getNeighbors(int j, int* indexBuffer) const;
  
  // Get the values at the neighbors of a node.
  // Note that homeX, homeY, and homeZ do not need to be wrapped,
  // since we do it here.
  void getNeighborValues(NeighborList* neigh, int homeX, int homeY, int homeZ) const;
  inline void setVal(float* v) { val = v; }

    BaseGrid* copy_to_cuda() const {
	BaseGrid* g_d = NULL;
	BaseGrid g_tmp;
	float* val_d = NULL;
	size_t sz = sizeof(float) * size;
	gpuErrchk(cudaMalloc(&g_d, sizeof(BaseGrid)));
	gpuErrchk(cudaMalloc(&val_d, sz));
	gpuErrchk(cudaMemcpy(val_d, val, sz, cudaMemcpyHostToDevice));
	g_tmp.origin = origin;
	g_tmp.basis = basis;
	g_tmp.nx = nx;
	g_tmp.ny = ny;
	g_tmp.nz = nz;
	g_tmp.size= size;
	g_tmp.basisInv = basisInv;
	g_tmp.val = val_d;
	gpuErrchk(cudaMemcpy(g_d, &g_tmp, sizeof(BaseGrid), cudaMemcpyHostToDevice));
	g_tmp.val = NULL;
	return g_d;
    }

    static void remove_from_cuda(BaseGrid* g_d) {
	BaseGrid g_tmp;
	gpuErrchk(cudaMemcpy(&g_tmp, g_d, sizeof(BaseGrid), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaFree(&(g_tmp.val)));
	g_tmp.val = NULL;
	gpuErrchk(cudaMemcpy(g_d, &g_tmp, sizeof(BaseGrid), cudaMemcpyHostToDevice)); // copy NULL back to device
	gpuErrchk(cudaFree(&g_d));
    }

public:
  Vector3 origin;
  Matrix3 basis;
  int nx, ny, nz;
  int size;
  Matrix3 basisInv;
public:
  float* val;
};

#ifndef delgpuErrchk
#undef  delgpuErrchk
#undef  gpuErrchk
#endif

#endif
