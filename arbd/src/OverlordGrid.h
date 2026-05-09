
///////////////////////////////////////////////////////////////////////
// Cell decomposition of points.
// Author: Jeff Comer <jcomer2@illinois.edu>
#ifndef OVERLORDGRID_H
#define OVERLORDGRID_H

#include "BaseGrid.h"
#include "useful.h"

#ifdef __CUDACC__
    #define HOST __host__
    #define DEVICE __device__
#else
    #define HOST 
    #define DEVICE 
#endif

class OverlordGrid : public BaseGrid {
public:
  OverlordGrid(const BaseGrid& grid) : BaseGrid(grid) {
    initSubgrids();
    initUniqueGrids();
  }
  OverlordGrid() {};  
  /*OverlordGrid(const char* systemDefFile) : BaseGrid(readDefFirst(systemDefFile)) {
    printf("size: %d\n", size);
    
    // Initialize stuff.
    initSubgrids();
    initUniqueGrids();

    // Load the rest of the system definition file.
    readDef(systemDefFile);
    }*/
  

  // Read a grid from a file.
  OverlordGrid(const char* rootGrid) : BaseGrid(rootGrid) {
    // Initialize stuff.
    initSubgrids();
    initUniqueGrids();
  }

private:
  void initSubgrids() {
    subgrid = new const BaseGrid*[size];
    subtrans = new Matrix3[size];
    for (int i = 0; i < size; i++) {
      subtrans[i] = Matrix3(1.0f);
      subgrid[i] = NULL;
    }
  }

  void initUniqueGrids() {
    uniqueGridNum = 0;
    uniqueGrid = new BaseGrid*[size];
    uniqueGridName = new String[size];
    for (int i = 0; i < size; i++) uniqueGrid[i] = NULL;
  }

public:
  int readDef(const char* systemDefFile) {
    // Open the file.
    FILE* inp = fopen(systemDefFile, "r");
    if (inp == NULL) {
      printf("OverlordGrid:readDef Couldn't open file `%s'.\n", systemDefFile);
      exit(-1);
    }
    
    int ind;
    char gridFile[STRLEN];
    char transform[STRLEN];
    char line[STRLEN];
    int nRead;
    int count = 0;
    while (fgets(line, STRLEN, inp) != NULL) {
      // Ignore comments.
      int len = strlen(line);
      if (line[0] == '#') continue;
      if (len < 2) continue;
      
      // Read definition lines.
      nRead = sscanf(line, "%d %s %s", &ind, gridFile, transform);
      if (nRead < 3) {
	printf("OverlordGrid:readDef Improperly formatted line `%s'\n", line);
	fclose(inp);
	exit(-1);
      }

      // Skip the root grid.
      if (ind < 0) continue;
      
      // Die for an improper index.
      if (ind >= size) {
	printf("OverlordGrid:readDef Index %d does not exist for %d nodes.\n", ind, size);
	fclose(inp);
	exit(-1);
      }

      // Find the grid to link to.
      String gridName(gridFile);
      int gridInd = -1;
      for (int i = 0; i < uniqueGridNum; i++) {
	if (gridName == uniqueGridName[i]) {
	  gridInd = i;
	  break;
	}
      }

      // This is new grid.
      // Load it.
      if (gridInd < 0) {
	if (uniqueGridNum >= size) {
	  printf("OverlordGrid:readDef Too many unique grids.\n");
	  fclose(inp);
	  exit(-1);
	}

	uniqueGrid[uniqueGridNum] = new BaseGrid(gridFile);
	uniqueGridName[uniqueGridNum] = gridFile;
	gridInd = uniqueGridNum;
	printf("New grid: %s\n", gridFile);
	uniqueGridNum++;
      }

      // Link the subgrid.
      link(ind, uniqueGrid[gridInd], parseTransform(transform));
      count++;
    }

    return count;
  }

  static String readDefFirst(const char* systemDefFile) {
    // Open the file.
    FILE* inp = fopen(systemDefFile, "r");
    if (inp == NULL) {
      printf("OverlordGrid:readDefFirst Couldn't open file `%s'.\n", systemDefFile);
      exit(-1);
    }
    
    int ind;
    char gridFile[STRLEN];
    char transform[STRLEN];
    char line[STRLEN];
    int nRead;
    while (fgets(line, STRLEN, inp) != NULL) {
      // Ignore comments.
      int len = strlen(line);
      if (line[0] == '#') continue;
      if (len < 2) continue;
      
      // Read definition lines.
      nRead = sscanf(line, "%d %s %s", &ind, gridFile, transform);
      if (nRead < 3 || ind != -1) {
	printf("OverlordGrid:readDefFirst Improperly formatted line `%s'\n", line);
	fclose(inp);
	exit(-1);
      }
      
      // Just get the root grid and return.
      return String(gridFile);
    }
    
    return String();
  }

  virtual ~OverlordGrid() {
    delete[] subgrid;

    for (int i = 0; i < uniqueGridNum; i++) delete uniqueGrid[i];
    delete[] uniqueGrid;
    delete[] uniqueGridName;
  }

  static Matrix3 parseTransform(const char* trans) {    
    if (strlen(trans) < 2) return Matrix3(1.0f);

    char sgn = trans[0];
    char axis = trans[1];

    Matrix3 ret(1.0f);
    switch(axis) {
    case 'x':
    case 'X':
      if (sgn == '-') ret = Matrix3(0.0f, 0.0f, -1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f);
      else ret = Matrix3(0.0f, 0.0f, 1.0f, 0.0f, -1.0f, 0.0f, 1.0f, 0.0f, 0.0f);
      break;

    case 'y':
    case 'Y':
      if (sgn == '-') ret = Matrix3(0.0f, -1.0f, 0.0f, 0.0f, 0.0f, -1.0f, 1.0f, 0.0f, 0.0f);
      else ret = Matrix3(0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f);
      break;

    case 'z':
    case 'Z':
      if (sgn == '-') ret = Matrix3(-1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, -1.0f);
      else ret = Matrix3(1.0f);
      break;

    }

    return ret.transpose();
  }

  // Link a grid node to a subgrid.
  bool link(int j, const BaseGrid* p, Matrix3 trans) {
    if (j < 0 || j >= size) return false;
    subgrid[j] = p;
    subtrans[j] = trans;

    return true;
  }
  bool link(int j, const BaseGrid& g, Matrix3 trans) {
    return link(j, &g, trans);
  }

  virtual float getPotential(Vector3 pos) const {
    // Find the nearest node.
    int j = nearestIndex(pos);
    
    // Return nothing for a null subgrid.
    if (subgrid[j] == NULL) return 0.0f;

    // Shift the point to get into node j's space.
    Vector3 r = subtrans[j].transform(pos - getPosition(j));
    // Do a getPotential on the subgrid.
    return subgrid[j]->getPotential(r);
  }

  DEVICE virtual float interpolatePotential(Vector3 pos) const {
    // Find the nearest node.
    int j = nearestIndex(pos);
    
    // Return nothing for a null subgrid.
    if (subgrid[j] == NULL) return 0.0f;

    // Shift the point to get into node j's space.
    Vector3 r = subtrans[j].transform(pos - getPosition(j));
    // Do a getPotential on the subgrid.
    return subgrid[j]->interpolatePotential(r);
  }

  DEVICE virtual float interpolatePotentialLinearly(Vector3 pos) const {
    // Find the nearest node.
    int j = nearestIndex(pos);
    
    // Return nothing for a null subgrid.
    if (subgrid[j] == NULL) return 0.0f;

    // Shift the point to get into node j's space.
    Vector3 r = subtrans[j].transform(pos - getPosition(j));
    // Do a getPotential on the subgrid.
    return subgrid[j]->interpolatePotentialLinearly(r);
  }
  
  Vector3 interpolateForce(Vector3 pos) const {
    // Find the nearest node.
    int j = nearestIndex(pos);
    
    // Return nothing for a null subgrid.
    if (subgrid[j] == NULL) return Vector3(0.0f);
    // Shift the point to get into node j's space.
    Vector3 r = subtrans[j].transform(pos - getPosition(j));

    Vector3 f;
 	Vector3 l = subgrid[j]->getInverseBasis().transform(r - subgrid[j]->getOrigin());
    	int homeX = int(floor(l.x));
    	int homeY = int(floor(l.y));
    	int homeZ = int(floor(l.z));
       	 // Get the array jumps with shifted indices.
   	 int jump[3];
    	jump[0] = subgrid[j]->getNz()*subgrid[j]->getNy();
    	jump[1] = subgrid[j]->getNz();
    	jump[2] = 1;
   	// Shift the indices in the home array.
   	int home[3];
    	home[0] = homeX;
   	home[1] = homeY;
    	home[2] = homeZ;

    	// Shift the indices in the grid dimensions.
    	int g[3];
	g[0] = subgrid[j]->getNx();
	g[1] = subgrid[j]->getNy();
	g[2] = subgrid[j]->getNz();

	// Get the interpolation coordinates.
	   float w[3];
	w[0] = l.x - homeX;
	w[1] = l.y - homeY;
	w[2] = l.z - homeZ;
	// Find the values at the neighbors.
	float g1[4][4][4];
	for (int ix = 0; ix < 4; ix++) {
	      	for (int iy = 0; iy < 4; iy++) {
			for (int iz = 0; iz < 4; iz++) {
	  		// Wrap around the periodic boundaries. 
				int jx = ix-1 + home[0];
		 		 jx = subgrid[j]->wrap(jx, g[0]);
		  		int jy = iy-1 + home[1];
		 		 jy = subgrid[j]->wrap(jy, g[1]);
		 		 int jz = iz-1 + home[2];
		  		jz = subgrid[j]->wrap(jz, g[2]);
		  
				 int ind = jz*jump[2] + jy*jump[1] + jx*jump[0];
				  g1[ix][iy][iz] = subgrid[j]->val[ind];
			}
	      	}
	}  

    f.x = subgrid[j]->interpolateDiffX(r, w, g1);
    f.y = subgrid[j]->interpolateDiffY(r, w, g1);
    f.z = subgrid[j]->interpolateDiffZ(r, w, g1);
    Matrix3 m = subgrid[j]->getInverseBasis();
    Vector3 f1 = m.transpose().transform(f);
    Vector3 f2 = subtrans[j].transpose().transform(f1);
    return f2;
  }
 
  int getUniqueGridNum() const { return uniqueGridNum; }

  bool writeSubgrids(const char* fileName) const {
    FILE* out = fopen(fileName, "w");
    if (out == NULL) return false;

    for (int i = 0; i < size; i++) {
      if (subgrid[i] != NULL)
	fprintf(out, "%d %g %s\n", i, subgrid[i]->mean(), subtrans[i].toString1().val());
    }
    fclose(out);

    return true;
  }

private:  
  const BaseGrid** subgrid;
  Matrix3* subtrans;
  BaseGrid** uniqueGrid;
  String* uniqueGridName;
  int uniqueGridNum;
};
#endif
