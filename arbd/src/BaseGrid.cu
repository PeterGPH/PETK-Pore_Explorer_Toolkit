
//////////////////////////////////////////////////////////////////////
// Grid base class that does just the basics.
// Author: Jeff Comer <jcomer2@illinois.edu>

#include "BaseGrid.h"
#include <cuda.h>


#define STRLEN 512

// Initialize the variables that get used a lot.
// Also, allocate the main value array.
void BaseGrid::init() {
	basisInv = basis.inverse();
	size = nx*ny*nz;
	val = new float[size];
}
BaseGrid::BaseGrid() {
	BaseGrid tmp(Matrix3(),Vector3(),1,1,1);
	val = new float[1];
	*this = tmp;									// TODO: verify that this is OK
	
	// basis = Matrix3();
	// origin = Vector3();
	// nx = 1;
	// ny = 1;
	// nz = 1;
	
	// init();
	// zero();
}

// The most obvious of constructors.
BaseGrid::BaseGrid(Matrix3 basis0, Vector3 origin0, int nx0, int ny0, int nz0) {
	basis = basis0;
	origin = origin0;
	nx = abs(nx0);
	ny = abs(ny0);
	nz = abs(nz0);
	
	init();
	zero();
}

// Make an orthogonal grid given the box dimensions and resolution.
BaseGrid::BaseGrid(Vector3 box, float dx) {
	dx = fabsf(dx);
	box.x = fabsf(box.x);
	box.y = fabsf(box.y);
	box.z = fabsf(box.z);

	// Tile the grid into the system box.
	// The grid spacing is always a bit smaller than dx.
	nx = int(ceilf(box.x/dx));
	ny = int(ceilf(box.y/dx));
	nz = int(ceilf(box.z/dx));
	if (nx <= 0) nx = 1;
	if (ny <= 0) ny = 1;
	if (nz <= 0) nz = 1;
	basis = Matrix3(box.x/nx, box.y/ny, box.z/nz);
	origin = -0.5f*box;

	init();
	zero();
}

// The box gives the system geometry.
// The grid point numbers define the resolution.
BaseGrid::BaseGrid(Matrix3 box, int nx0, int ny0, int nz0) {
	nx = nx0;
	ny = ny0;
	nz = nz0;

	// Tile the grid into the system box.
	if (nx <= 0) nx = 1;
	if (ny <= 0) ny = 1;
	if (nz <= 0) nz = 1;
	basis = Matrix3(box.ex()/nx, box.ey()/ny, box.ez()/nz);
	origin = -0.5f*(box.ex() + box.ey() + box.ez());

	init();
	zero();
}

// The box gives the system geometry.
// dx is the approx. resolution.
// The grid spacing is always a bit larger than dx.
BaseGrid::BaseGrid(Matrix3 box, Vector3 origin0, float dx) {
	dx = fabs(dx);
	
	// Tile the grid into the system box.
	// The grid spacing is always a bit larger than dx.
	nx = int(floor(box.ex().length()/dx))-1;
	ny = int(floor(box.ey().length()/dx))-1;
	nz = int(floor(box.ez().length()/dx))-1;
	if (nx <= 0) nx = 1;
	if (ny <= 0) ny = 1;
	if (nz <= 0) nz = 1;

	basis = Matrix3(box.ex()/nx, box.ey()/ny, box.ez()/nz);
	origin = origin0;

	init();
	zero();
}

// The box gives the system geometry.
// dx is the approx. resolution.
// The grid spacing is always a bit smaller than dx.
BaseGrid::BaseGrid(Matrix3 box, float dx) {
	dx = fabs(dx);
	
	// Tile the grid into the system box.
	// The grid spacing is always a bit smaller than dx.
	nx = int(ceilf(box.ex().length()/dx));
	ny = int(ceilf(box.ey().length()/dx));
	nz = int(ceilf(box.ez().length()/dx));
	if (nx <= 0) nx = 1;
	if (ny <= 0) ny = 1;
	if (nz <= 0) nz = 1;

	basis = Matrix3(box.ex()/nx, box.ey()/ny, box.ez()/nz);
	origin = -0.5f*(box.ex() + box.ey() + box.ez());

	init();
	zero();
}

// Make an exact copy of a grid.
BaseGrid::BaseGrid(const BaseGrid& g) {
	nx = g.nx;
	ny = g.ny;
	nz = g.nz;
	basis = g.basis;
	origin = g.origin;
	
	init();
	for (int i = 0; i < size; i++) val[i] = g.val[i];
}

BaseGrid BaseGrid::mult(const BaseGrid& g) {
	for (int i = 0; i < size; i++) val[i] *= g.val[i];
	return *this;
}

BaseGrid& BaseGrid::operator=(const BaseGrid& g) {
	if(val != NULL) delete[] val;
	val = NULL;
	nx = g.nx;
	ny = g.ny;
	nz = g.nz;
	basis = g.basis;
	origin = g.origin;
	
	init();
	for (int i = 0; i < size; i++) val[i] = g.val[i];

	return *this;
}


// Make a copy of a grid, but at a different resolution.
BaseGrid::BaseGrid(const BaseGrid& g, int nx0, int ny0, int nz0) : nx(nx0),  ny(ny0), nz(nz0) {
	if (nx <= 0) nx = 1;
	if (ny <= 0) ny = 1;
	if (nz <= 0) nz = 1;

	// Tile the grid into the box of the template grid.
	Matrix3 box = g.getBox();
	basis = Matrix3(box.ex()/nx, box.ey()/ny, box.ez()/nz);

	origin = g.origin;
	init();

	// Do an interpolation to obtain the values.
	for (int i = 0; i < size; i++) {
		Vector3 r = getPosition(i);
		val[i] = g.interpolatePotential(r);
	}
}

// Read a grid from a file.
BaseGrid::BaseGrid(const char* fileName) {
		 // Open the file.
	FILE* inp = fopen(fileName,"r");
	if (inp == NULL) {
		printf("ERROR BaseGrid::BaseGrid Couldn't open file %s.\n",fileName);
		exit(-1);
	}
	//printf("Reading dx file %s...\n", fileName);
	
	size = 0;
	nx = 0;
	ny = 0;
	nz = 0;
	basis = Matrix3(1.0f);
	origin = Vector3(0.0f);    

	int n = 0;
	float x, y, z;
	char line[STRLEN];
	int p, nRead;
	int deltaCount = 0;
	Vector3 base[3];
	while (fgets(line, STRLEN, inp) != NULL) {
		// Ignore comments.
		int len = strlen(line);
		if (line[0] == '#') continue;
		if (len < 2) continue;
	
		if (isInt(line[0]) && n < size) {
// Read grid values.
nRead = sscanf(line, "%f %f %f", &x, &y, &z);
if (size > 0) {
	switch(nRead) {
	case 1:
		val[n] = x;
		n++;
		if (n != size) {
			printf("ERROR BaseGrid::BaseGrid Improperly formatted dx file %s.\n", fileName);
			printf("line `%s'\n", line);
		}
		break;
	case 2:
		val[n] = x;
		val[n+1] = y;
		n += 2;
		if (n != size) {
			printf("ERROR BaseGrid::BaseGrid Improperly formatted dx file %s.\n", fileName);
			printf("line `%s'\n", line);
		}
		break;
	case 3:
		val[n] = x;
		val[n+1] = y;
		val[n+2] = z;
		n += 3;
		break;
	}
}
		} else if (len > 5) {
// Read the grid parameters.
char start[6];
for (int i = 0; i < 5; i++) start[i] = line[i];
start[5] = '\0';

if(strcmp("origi", start) == 0) {
	// Get an origin line.
	p = firstSpace(line, STRLEN);
	sscanf(&(line[p+1]), "%f %f %f", &x, &y, &z);
	origin = Vector3(x, y, z);
	//printf("Origin: %.12g %.12g %.12g\n", x, y, z);
} else if(strcmp("delta", start) == 0) {
	// Get a delta matrix line.
	p = firstSpace(line, STRLEN);
	sscanf(&(line[p+1]), "%f %f %f", &x, &y, &z);
	base[deltaCount] = Vector3(x, y, z);
	//printf("Delta %d: %.12g %.12g %.12g\n", deltaCount, x, y, z);
	if (deltaCount < 2) deltaCount = deltaCount + 1;
} else if(strcmp("objec", start) == 0) {
	//printf("%s", line);
	// Get the system dimensions.
	if (line[7] != '1') continue;
	int read = sscanf(line, "object 1 class gridpositions counts %d %d %d\n", &nx, &ny, &nz);
	//printf("Size: %d %d %d\n", nx, ny, nz);
	if (read == 3) {
		size = nx*ny*nz;
		val = new float[size];
		zero();
	}
}
		}
	}
	fclose(inp);

	basis = Matrix3(base[0], base[1], base[2]);
	basisInv = basis.inverse();
	if (size == 0 || n != size) {
		printf("ERROR BaseGrid::BaseGrid Improperly formatted dx file %s.\n",fileName);
		printf("declared size: %d, items: %d\n", size, n);
		printf("first value: %10g, final value: %.10g\n", val[0], val[n-1]);
		exit(-1);
	}
}  

// Write without comments.
void BaseGrid::write(const char* fileName) const {
	write(fileName, "");
}

// Writes the grid as a file in the dx format.
void BaseGrid::write(const char* fileName, const char* comments) const {
	// Open the file.
	FILE* out = fopen(fileName,"w");
	if (out == NULL) {
		printf("ERROR BaseGrid::write Couldn't open file %s.\n",fileName);
		exit(-1);
	}

	// Write the header.
	fprintf(out, "# %s\n", comments);
	fprintf(out, "object 1 class gridpositions counts %d %d %d\n", nx, ny, nz);
	fprintf(out, "origin %.12g %.12g %.12g\n", origin.x, origin.y, origin.z);
	fprintf(out, "delta %.12g %.12g %.12g\n", basis.exx, basis.eyx, basis.ezx);
	fprintf(out, "delta %.12g %.12g %.12g\n", basis.exy, basis.eyy, basis.ezy);
	fprintf(out, "delta %.12g %.12g %.12g\n", basis.exz, basis.eyz, basis.ezz);
	fprintf(out, "object 2 class gridconnections counts %d %d %d\n", nx, ny, nz);
	fprintf(out, "object 3 class array type float rank 0 items %d data follows\n", size);
	
	// Write the data.
	int penultima = 3*(size/3);
	int mod = size - penultima;

	int i;
	for (i = 0; i < penultima; i+=3) {
		fprintf(out, "%.12g %.12g %.12g\n", val[i], val[i+1], val[i+2]);
	}
	if (mod == 1) {
		fprintf(out, "%.12g\n", val[size-1]);
	} else if (mod == 2) {
		fprintf(out, "%.12g %.12g\n", val[size-2], val[size-1]);
	}
	fclose(out);
}

// Writes the grid data as a single column in the order:
// nx ny nz ox oy oz dxx dyx dzx dxy dyy dzy dxz dyz dzz val0 val1 val2 ...
void BaseGrid::writeData(const char* fileName) {
	// Open the file.
	FILE* out = fopen(fileName,"w");
	if (out == NULL) {
		printf("Couldn't open file %s.\n",fileName);
		exit(-1);
	}

	fprintf(out, "%d\n%d\n%d\n", nx, ny, nz);
	fprintf(out, "%.12g\n%.12g\n%.12g\n", origin.x, origin.y, origin.z);
	fprintf(out, "%.12g\n%.12g\n%.12g\n", basis.exx, basis.eyx, basis.ezx);
	fprintf(out, "%.12g\n%.12g\n%.12g\n", basis.exx, basis.eyx, basis.ezx);
	fprintf(out, "%.12g\n%.12g\n%.12g\n", basis.exx, basis.eyx, basis.ezx);

	for (int i = 0; i < size; i++) fprintf(out, "%.12g\n", val[i]);
	fclose(out);
}

// Write the valies in a single column.
void BaseGrid::writePotential(const char* fileName) const {
	FILE* out = fopen(fileName, "w");
	for (int i = 0; i < size; i++) fprintf(out, "%.12g\n", val[i]);
	fclose(out);
}

BaseGrid::~BaseGrid() {
	if (val != NULL)
		delete[] val;
}

void BaseGrid::zero() {
	for (int i = 0; i < size; i++) val[i] = 0.0f;
}

bool BaseGrid::setValue(int j, float v) {
	if (j < 0 || j >= size) return false;
	val[j] = v;
	return true;
}

bool BaseGrid::setValue(int ix, int iy, int iz, float v) {
	if (ix < 0 || ix >= nx) return false;
	if (iy < 0 || iy >= ny) return false;
	if (iz < 0 || iz >= nz) return false;
	int j = iz + iy*nz + ix*ny*nz;

	val[j] = v;
	return true;
}

float BaseGrid::getValue(int j) const {
	if (j < 0 || j >= size) return 0.0f;
	return val[j];
}
/*
float BaseGrid::getValue(int ix, int iy, int iz) const {
	if (ix < 0 || ix >= nx) return 0.0f;
	if (iy < 0 || iy >= ny) return 0.0f;
	if (iz < 0 || iz >= nz) return 0.0f;
	
	int j = iz + iy*nz + ix*ny*nz;
	return val[j];
}
*/
Vector3 BaseGrid::getPosition(int ix, int iy, int iz) const {
	return basis.transform(Vector3(ix, iy, iz)) + origin;
}

Vector3 BaseGrid::getPosition(int j) const {
	int iz = j%nz;
	int iy = (j/nz)%ny;
	int ix = j/(nz*ny);

	return basis.transform(Vector3(ix, iy, iz)) + origin;
}

// Does the point r fall in the grid?
// Obviously this is without periodic boundary conditions.
bool BaseGrid::inGrid(Vector3 r) const {
	Vector3 l = basisInv.transform(r-origin);

	if (l.x < 0.0f || l.x >= nx) return false;
	if (l.y < 0.0f || l.y >= ny) return false;
	if (l.z < 0.0f || l.z >= nz) return false;
	return true;
}

bool BaseGrid::inGridInterp(Vector3 r) const {
	Vector3 l = basisInv.transform(r-origin);

	if (l.x < 2.0f || l.x >= nx-3.0f) return false;
	if (l.y < 2.0f || l.y >= ny-3.0f) return false;
	if (l.z < 2.0f || l.z >= nz-3.0f) return false;
	return true;
}

Vector3 BaseGrid::transformTo(Vector3 r) const {
	return basisInv.transform(r-origin);
}
Vector3 BaseGrid::transformFrom(Vector3 l) const {
	return basis.transform(l) + origin;
}

IndexList BaseGrid::index(int j) const {
	int iz = j%nz;
	int iy = (j/nz)%ny;
	int ix = j/(nz*ny);
	IndexList ret;
	ret.add(ix);
	ret.add(iy);
	ret.add(iz);
	return ret;
}
int BaseGrid::indexX(int j) const { return j/(nz*ny); }
int BaseGrid::indexY(int j) const { return (j/nz)%ny; }
int BaseGrid::indexZ(int j) const { return j%nz; }
int BaseGrid::index(int ix, int iy, int iz) const { return iz + iy*nz + ix*ny*nz; }

int BaseGrid::index(Vector3 r) const {
	Vector3 l = basisInv.transform(r-origin);
	
	int ix = int(floor(l.x));
	int iy = int(floor(l.y));
	int iz = int(floor(l.z));

	ix = wrap(ix, nx);
	iy = wrap(iy, ny);
	iz = wrap(iz, nz);
	
	return iz + iy*nz + ix*ny*nz;
}

int BaseGrid::nearestIndex(Vector3 r) const {
	Vector3 l = basisInv.transform(r-origin);
	
	int ix = int(floorf(l.x + 0.5f));
	int iy = int(floorf(l.y + 0.5f));
	int iz = int(floorf(l.z + 0.5f));

	ix = wrap(ix, nx);
	iy = wrap(iy, ny);
	iz = wrap(iz, nz);
	
	return iz + iy*nz + ix*ny*nz;
}

// A matrix defining the basis for the entire system.
Matrix3 BaseGrid::getBox() const {
	return Matrix3(nx*basis.ex(), ny*basis.ey(), nz*basis.ez());
} 
// The longest diagonal of the system.
Vector3 BaseGrid::getExtent() const {
	return basis.transform(Vector3(nx,ny,nz));
}
// The longest diagonal of the system.
float BaseGrid::getDiagonal() const {
	return getExtent().length();
}
// The position farthest from the origin.
Vector3 BaseGrid::getDestination() const {
	return basis.transform(Vector3(nx,ny,nz)) + origin;
}
// The center of the grid.
Vector3 BaseGrid::getCenter() const {
	return basis.transform(Vector3(0.5f*nx,0.5f*ny,0.5f*nz)) + origin;
}
// The volume of a single cell.
float BaseGrid::getCellVolume() const {
	return fabs(basis.det());
}
// The volume of the entire system.
float BaseGrid::getVolume() const {
	return getCellVolume()*size;
}
Vector3 BaseGrid::getCellDiagonal() const {
	return basis.ex() + basis.ey() + basis.ez();
}

// Add a fixed value to the grid.
void BaseGrid::shift(float s) {
	for (int i = 0; i < size; i++) val[i] += s;
}

// Multiply the grid by a fixed value.
void BaseGrid::scale(float s) {
	for (int i = 0; i < size; i++) val[i] *= s;
}

// Get the mean of the entire grid.
float BaseGrid::mean() const {
	float sum = 0.0f;
	for (int i = 0; i < size; i++) sum += val[i];
	return sum/size;
}

// Compute the average profile along an axis.
// Assumes that the grid axis with index "axis" is aligned with the world axis of index "axis".
void BaseGrid::averageProfile(const char* fileName, int axis) {
	FILE* out = fopen(fileName,"w");
	if (out == NULL) {
		printf("Couldn't open file %s.\n",fileName);
		exit(-1);
	}

	int dir0 = wrap(axis, 3);
	int dir1 = (axis+1)%3;
	int dir2 = (axis+2)%3;

	int jump[3];
	jump[0] = ny*nz;
	jump[1] = nz;
	jump[2] = 1;

	int n[3];
	n[0] = nx;
	n[1] = ny;
	n[2] = nz;
 
	for (int i0 = 0; i0 < n[dir0]; i0++) {
		float sum = 0;

		for (int i1 = 0; i1 < n[dir1]; i1++) {
			for (int i2 = 0; i2 < n[dir2]; i2++) {
				int j = i0*jump[dir0] + i1*jump[dir1] + i2*jump[dir2];
				sum += val[j];
			}
		}
		
		float v = sum/(n[dir1]*n[dir2]);
		float x = 0.0f;
		switch (dir0) {
		case 0:
			x = origin.x + i0*basis.exx;
			break;
		case 1:
			x = origin.y + i0*basis.eyy;
			break;
		case 2:
			x = origin.z + i0*basis.ezz;
			break;
		}
		fprintf(out, "%0.10g %0.10g\n", x, v);
	}

	fclose(out);
}

// Get the potential at the closest node.
float BaseGrid::getPotential(Vector3 pos) const {
	// Find the nearest node.
	int j = nearestIndex(pos);

	return val[j];
}

bool BaseGrid::crop(int x0, int y0, int z0, int x1, int y1, int z1, bool keep_origin) {
	if (x0 < 0 || x0 >= 2 * nx) x0 = 0;
	if (y0 < 0 || y0 >= 2 * ny) y0 = 0;
	if (z0 < 0 || z0 >= 2 * nz) z0 = 0;
	if (x1 < 0 || x1 >= 2 * nx) x1 = 2 * nx - 1;
	if (y1 < 0 || y1 >= 2 * ny) y1 = 2 * ny - 1;
	if (z1 < 0 || z1 >= 2 * nz) z1 = 2 * nz - 1;
	printf("Cropping to (%d, %d, %d) -> (%d, %d, %d)\n", x0, y0, z0, x1, y1, z1);

	if (x0 >= x1 || y0 >= y1 || z0 >= z1)
		return false;

	int new_nx = x1 - x0 + 1;
	int new_ny = y1 - y0 + 1;
	int new_nz = z1 - z0 + 1;
	int new_size = new_nx * new_ny * new_nz;
	float *new_val = new float[new_size];

	int ind = 0;
	for (int i = x0; i < x1; i++)
		for (int j = y0; j < y1; j++)
			for (int k = z0; k < z1; k++) {
				int ind1 = k + j * nz + i * ny*nz;
				new_val[ind++] = val[ind1];
			}

	if (!keep_origin)
		origin += basis.transform(Vector3(x0, y0, z0));
	nx = new_nx;
	ny = new_ny;
	nz = new_nz;
	size = new_size;
	delete[] val;
	val = new_val;

	return true;
}

Vector3 BaseGrid::wrapDiffNearest(Vector3 r) const {
	Vector3 l = basisInv.transform(r);
	l.x = wrapDiff(l.x, nx);
	l.y = wrapDiff(l.y, ny);
	l.z = wrapDiff(l.z, nz);

	float length2 = basis.transform(l).length2();

	for (int dx = -1; dx <= 1; dx++) {
		for (int dy = -1; dy <= 1; dy++) {
			for (int dz = -1; dz <= 1; dz++) {
				//if (dx == 0 && dy == 0 && dz == 0) continue;
				Vector3 tmp = Vector3(l.x+dx*nx, l.y+dy*ny, l.z+dz*nz);
				if (basis.transform(tmp).length2() < length2) {
					l = tmp;
					length2 = basis.transform(l).length2();
				}
			}
		}
	}

	return basis.transform(l);
}


// Includes the home node.
// indexBuffer must have a size of at least 27.
void BaseGrid::getNeighbors(int j, int* indexBuffer) const {
	int jx = indexX(j);
	int jy = indexY(j);
	int jz = indexZ(j);

	int k = 0;
	for (int ix = -1; ix <= 1; ix++) {
		for (int iy = -1; iy <= 1; iy++) {
			for (int iz = -1; iz <= 1; iz++) {
				int ind = wrap(jz+iz,nz) + nz*wrap(jy+iy,ny) + ny*nz*wrap(jx+ix,nx);
				indexBuffer[k] = ind;
				k++;
			}
		}
	}
}

// Get the values at the neighbors of a node.
// Note that homeX, homeY, and homeZ do not need to be wrapped,
// since we do it here.
void BaseGrid::getNeighborValues(NeighborList* neigh, int homeX, int homeY, int homeZ) const {
	for (int ix = -1; ix <= 1; ix++) {
		for (int iy = -1; iy <= 1; iy++) {
			for (int iz = -1; iz <= 1; iz++) {
				int ind = wrap(homeZ+iz,nz) + nz*wrap(homeY+iy,ny) + ny*nz*wrap(homeX+ix,nx);
				neigh->v[ix+1][iy+1][iz+1] = val[ind];
			}
		}
	}
}  
