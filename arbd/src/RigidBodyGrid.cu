//////////////////////////////////////////////////////////////////////
// Grid base class that does just the basics.
// Author: Jeff Comer <jcomer2@illinois.edu>

#include "RigidBodyGrid.h"
#include <cuda.h>

#define STRLEN 512

	/*                               \
	| CONSTRUCTORS, DESTRUCTORS, I/O |
	\===============================*/

RigidBodyGrid::RigidBodyGrid() {
	RigidBodyGrid tmp(1,1,1);
	val = new float[1];
	*this = tmp;									// TODO: verify that this is OK
}

// The most obvious of constructors.
RigidBodyGrid::RigidBodyGrid(int nx0, int ny0, int nz0) {
	nx = abs(nx0);
	ny = abs(ny0);
	nz = abs(nz0);
	
	val = new float[nx*ny*nz];
	zero();
}

RigidBodyGrid::RigidBodyGrid(const BaseGrid& g) {
	nx = g.nx;
	ny = g.ny;
	nz = g.nz;
	
	val = new float[nx*ny*nz];
	for (int i = 0; i < nx*ny*nz; i++) val[i] = g.val[i];
}

// Make an exact copy of a grid.
RigidBodyGrid::RigidBodyGrid(const RigidBodyGrid& g) {
	nx = g.nx;
	ny = g.ny;
	nz = g.nz;
	
	val = new float[nx*ny*nz];
	for (int i = 0; i < nx*ny*nz; i++) val[i] = g.val[i];
}

RigidBodyGrid RigidBodyGrid::mult(const RigidBodyGrid& g) {
	for (int i = 0; i < nx*ny*nz; i++) val[i] *= g.val[i];
	return *this;
}

RigidBodyGrid& RigidBodyGrid::operator=(const RigidBodyGrid& g) {
	if(val!=NULL) 
            delete[] val;
	val = NULL;
	nx = g.nx;
	ny = g.ny;
	nz = g.nz;
	
	val = new float[nx*ny*nz];
	for (int i = 0; i < nx*ny*nz; i++) val[i] = g.val[i];

	return *this;
}

RigidBodyGrid::~RigidBodyGrid() {
	if (val != NULL)
        {
		delete[] val;
                val = NULL;
        }
}

void RigidBodyGrid::zero() {
	for (int i = 0; i < nx*ny*nz; i++) val[i] = 0.0f;
}

bool RigidBodyGrid::setValue(int j, float v) {
	if (j < 0 || j >= nx*ny*nz) return false;
	val[j] = v;
	return true;
}

bool RigidBodyGrid::setValue(int ix, int iy, int iz, float v) {
	if (ix < 0 || ix >= nx) return false;
	if (iy < 0 || iy >= ny) return false;
	if (iz < 0 || iz >= nz) return false;
	int j = iz + iy*nz + ix*ny*nz;

	val[j] = v;
	return true;
}

float RigidBodyGrid::getValue(int j) const {

	if (j < 0 || j >= nx*ny*nz) return 0.0f;
	return val[j];
/*
    Vector3 idx = getPosition(j)
    return getValue(idx.x,idx.y,idx.z);
*/
}

HOST DEVICE float RigidBodyGrid::getValue(int ix, int iy, int iz) const {
/*
           if(ix < 0) ix = 0;
           else if(ix >= nx) ix = nx -1;

           if(iy < 0) iy = 0;
           else if(iy >= ny) iy = ny-1;

           if(iz < 0) iz = 0;
           else if(iz >= nz) iz = nz-1;

           int j = iz + nz * (iy + ny * ix);
           return val[j];
*/

	if (ix < 0 || ix >= nx) return 0.0f;
	if (iy < 0 || iy >= ny) return 0.0f;
	if (iz < 0 || iz >= nz) return 0.0f;
	
	int j = iz + iy*nz + ix*ny*nz;
	return val[j];

}

Vector3 RigidBodyGrid::getPosition(const int j) const {
	/* const int iz = j%nz; */
	/* const int iy = (j/nz)%ny; */
	/* const int ix = j/(nz*ny); */
	const int jy = j/nz;
	const int jx = jy/ny;

	const int iz = j - jy*nz;
	const int iy = jy - jx*ny;
	// const int ix = jx;

	return Vector3(jx,iy,iz);
}

Vector3 RigidBodyGrid::getPosition(int j, Matrix3 basis, Vector3 origin) const {
	int iz = j%nz;
	int iy = (j/nz)%ny;
	int ix = j/(nz*ny);

	return basis.transform(Vector3(ix, iy, iz)) + origin;
}

IndexList RigidBodyGrid::index(int j) const {
	int iz = j%nz;
	int iy = (j/nz)%ny;
	int ix = j/(nz*ny);
	IndexList ret;
	ret.add(ix);
	ret.add(iy);
	ret.add(iz);
	return ret;
}
int RigidBodyGrid::indexX(int j) const { return j/(nz*ny); }
int RigidBodyGrid::indexY(int j) const { return (j/nz)%ny; }
int RigidBodyGrid::indexZ(int j) const { return j%nz; }
int RigidBodyGrid::index(int ix, int iy, int iz) const { return iz + iy*nz + ix*ny*nz; }

// Add a fixed value to the grid.
void RigidBodyGrid::shift(float s) {
	for (int i = 0; i < nx*ny*nz; i++) val[i] += s;
}

// Multiply the grid by a fixed value.
void RigidBodyGrid::scale(float s) {
	for (int i = 0; i < nx*ny*nz; i++) val[i] *= s;
}

/** interpolateForce() to be used on CUDA Device **/
DEVICE ForceEnergy RigidBodyGrid::interpolateForceD(const Vector3 l) const {
	Vector3 f;
	// Vector3 l = basisInv.transform(pos - origin);
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
	
	return ForceEnergy(f,e);
}
//#define cubic
DEVICE ForceEnergy RigidBodyGrid::interpolateForceDLinearly(const Vector3& l) const {
//#ifdef cubic
//return interpolateForceD(l);
//#elif defined(cubic_namd)
//return interpolateForceDnamd(l);
//#else
	// Find the home node.
	const int homeX = int(floor(l.x));
	const int homeY = int(floor(l.y));
	const int homeZ = int(floor(l.z));

	Vector3 f;

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
	f.x = -(wz * (g3[0][1] - g3[0][0]) + g3[0][0]);
	f.y = -(wz * (g3[1][1] - g3[1][0]) + g3[1][0]);
	f.z = -      (g3[2][1] - g3[2][0]);
	float e = wz * (g3[2][1] - g3[2][0]) + g3[2][0];
	return ForceEnergy(f,e);
//#endif
}
DEVICE ForceEnergy RigidBodyGrid::interpolateForceDnamd(const Vector3& l) const
{
                Vector3 f;
                //const Vector3 l = basisInv.transform(pos - origin);

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

                //f = basisInv.transpose().transform(f);
                return ForceEnergy(f,e);
        }

DEVICE float RigidBodyGrid::compute_V(float *a, float *x, float *y, float *z) const
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
DEVICE Vector3 RigidBodyGrid::compute_dV(float *a, float *x, float *y, float *z) const
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
DEVICE void RigidBodyGrid::compute_a(float *a, float *b) const
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
DEVICE void RigidBodyGrid::compute_b(float * __restrict__ b, int * __restrict__ inds) const
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
                //int d_hi[3] = {1, 1, 1};
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
                    b[8+i0]  = dscales[0] * (getValue(inds2[0]+1,inds2[1],inds2[2]) - getValue(inds2[0]-d_lo[0],inds2[1],inds2[2]) + voffs[0]); //  dV/dx
                    b[16+i0] = dscales[1] * (getValue(inds2[0],inds2[1]+1,inds2[2]) - getValue(inds2[0],inds2[1]-d_lo[1],inds2[2]) + voffs[1]); //  dV/dy
                    b[24+i0] = dscales[2] * (getValue(inds2[0],inds2[1],inds2[2]+1) - getValue(inds2[0],inds2[1],inds2[2]-d_lo[2]) + voffs[2]); //  dV/dz
                    b[32+i0] = dscales[0] * dscales[1] *
                        (getValue(inds2[0]+1,inds2[1]+1,inds2[2]) - getValue(inds2[0]-d_lo[0],inds2[1]+1,inds2[2])
                       - getValue(inds2[0]+1,inds2[1]-d_lo[1],inds2[2]) + getValue(inds2[0]-d_lo[0],inds2[1]-d_lo[1],inds2[2]));      //  d2V/dxdy

                    b[40+i0] = dscales[0] * dscales[2] *
                              (getValue(inds2[0]+1,inds2[1],inds2[2]+1) - getValue(inds2[0]-d_lo[0],inds2[1],inds2[2]+1)
                             - getValue(inds2[0]+1,inds2[1],inds2[2]-d_lo[2]) + getValue(inds2[0]-d_lo[0],inds2[1],inds2[2]-d_lo[2]));      //  d2V/dxdz

                    b[48+i0] = dscales[1] * dscales[2] *
                               (getValue(inds2[0],inds2[1]+1,inds2[2]+1) - getValue(inds2[0],inds2[1]-d_lo[1],inds2[2]+1)
                              - getValue(inds2[0],inds2[1]+1,inds2[2]-d_lo[2]) + getValue(inds2[0],inds2[1]-d_lo[1],inds2[2]-d_lo[2]));      //  d2V/dydz

                    b[56+i0] = dscales[0] * dscales[1] * dscales[2] *                                    // d3V/dxdydz
                       (getValue(inds2[0]+1,inds2[1]+1,inds2[2]+1) - getValue(inds2[0]+1,inds2[1]+1,inds2[2]-d_lo[2]) -
                        getValue(inds2[0]+1,inds2[1]-d_lo[1],inds2[2]+1) - getValue(inds2[0]-d_lo[0],inds2[1]+1,inds2[2]+1) +
                        getValue(inds2[0]+1,inds2[1]-d_lo[1],inds2[2]-d_lo[2]) + getValue(inds2[0]-d_lo[0],inds2[1]+1,inds2[2]-d_lo[2]) +
                        getValue(inds2[0]-d_lo[0],inds2[1]-d_lo[1],inds2[2]+1) - getValue(inds2[0]-d_lo[0],inds2[1]-d_lo[1],inds2[2]-d_lo[2]));

                        }
                    }
                }

