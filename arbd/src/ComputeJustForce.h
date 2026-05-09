///////////////////////////////////////////////////////////////////////
// Brownian dynamics base class
// Author: Jeff Comer <jcomer2@illinois.edu>
#ifndef COMPUTEFORCE_H
#define COMPUTEFORCE_H

#include "BaseGrid.h"
#include "BrownianParticleType.h"
#include "CellDecomposition.h"
#include "TabulatedPotential.h"

class ComputeForce {
public:
  ComputeForce(int num0, const BrownianParticleType* part, int numParts0, const BaseGrid* g, float switchStart0, float switchLen0, float electricConst0) : 
    num(num0),
    numParts(numParts0),
    sys(g), switchStart(switchStart0),
    switchLen(switchLen0), electricConst(electricConst0),
    cutoff2((switchLen0+switchStart0)*(switchLen0+switchStart0)), 
    decomp(g->getBox(), g->getOrigin(), switchStart0+switchLen0) {
    
    // Allocate the parameter tables.
    tableEps = new float*[numParts];
    tableRad6 = new float*[numParts];
    tableAlpha = new float*[numParts];
    for (int i = 0; i < numParts; i++) {
      tableEps[i] = new float[numParts];
      tableRad6[i] = new float[numParts];
      tableAlpha[i] = new float[numParts];
    }    

    // Form the parameter tables.
    makeTables(part);
    tablePot = new TabulatedPotential*[numParts*numParts];
    for (int i = 0; i < numParts*numParts; i++) tablePot[i] = NULL;
    
    // Make the cell decomposition.
    neigh = new IndexList[num];
  }

  ~ComputeForce() {
    for (int i = 0; i < numParts; i++) {
      delete[] tableEps[i];
      delete[] tableRad6[i];
      delete[] tableAlpha[i];
    }
    delete[] tableEps;
    delete[] tableRad6;
    delete[] tableAlpha;

    for (int j = 0; j < numParts*numParts; j++) {
      if (tablePot[j] != NULL) {
	delete tablePot[j];
	tablePot[j] = NULL;
      }
    }
    delete[] tablePot;

    delete[] neigh;
  }

  void updateNumber(const Vector3* pos, int newNum) {
    if (newNum == num || newNum < 0) return;

    // Set the new number.
    num = newNum;

    // Reallocate the neighbor list.
    delete[] neigh;
    neigh = new IndexList[num];
    decompose(pos);
  }

  void makeTables(const BrownianParticleType* part) {
    for (int i = 0; i < numParts; i++) {
      for (int j = 0; j < numParts; j++) {
	tableEps[i][j] = sqrtf(part[i].eps*part[j].eps);
	float r = part[i].radius + part[j].radius;
	tableRad6[i][j] = r*r*r*r*r*r;
	tableAlpha[i][j] = electricConst*part[i].charge*part[j].charge;
      }
    }
  }

  bool addTabulatedPotential(String fileName, int type0, int type1) {
    if (type0 < 0 || type0 >= numParts) return false;
    if (type1 < 0 || type1 >= numParts) return false;

    int ind = type0 + type1*numParts;
    int ind1 = type1 + type0*numParts;

    if (tablePot[ind] != NULL) {
      delete tablePot[ind];
      tablePot[ind] = NULL;
    }
    if (tablePot[ind1] != NULL) delete tablePot[ind1];
    
    tablePot[ind] = new TabulatedPotential(fileName);
    tablePot[ind1] = new TabulatedPotential(*tablePot[ind]);

    return true;
  }

  void decompose(const Vector3* pos) {
    // Reset the cell decomposition.
    decomp.clearCells();
    decomp.decompose(pos, num);

    // Regenerate the neighbor lists.
    for (int i = 0; i < num; i++) neigh[i] = decomp.neighborhood(pos[i]);
  }

  IndexList decompDim() const {
    IndexList ret;
    ret.add(decomp.getNx());
    ret.add(decomp.getNy());
    ret.add(decomp.getNz());
    return ret;
  }

  float decompCutoff() const {
    return decomp.getCutoff();
  }

  IndexList neighborhood(Vector3 r) const {
    return decomp.neighborhood(r);
  }

  void computeFull(Vector3* force, const Vector3* pos, const int* type) const {
    // Zero the force.
    for (int i = 0; i < num; i++) force[i] = Vector3(0.0);
    
    // Compute the force for all pairs.
    for (int i = 0; i < num-1; i++) {
      for (int j = i + 1; j < num; j++) {
	float alpha = tableAlpha[type[i]][type[j]];
	float eps = tableEps[type[i]][type[j]];
	float rad6 = tableRad6[type[i]][type[j]];
	Vector3 dr = sys->wrapDiff(pos[j] - pos[i]);

	Vector3 fc = coulombForceFull(dr, alpha);
	Vector3 fh = softcoreForce(dr, eps, rad6);

	force[i] += fc + fh;
	force[j] -= fc + fh;
      }
    }
  }

  void compute(Vector3* force, const Vector3* pos, const int* type) const {
    for (int i = 0; i < num; i++) {
      // Zero the force.
      force[i] = Vector3(0.0);

      // Loop through the neighbors.
      for (int n = 0; n < neigh[i].length(); n++) {
	int j = neigh[i].get(n);
	if (j == i) continue;

	float alpha = tableAlpha[type[i]][type[j]];
	float eps = tableEps[type[i]][type[j]];
	float rad6 = tableRad6[type[i]][type[j]];
	Vector3 dr = sys->wrapDiff(pos[j] - pos[i]);

	Vector3 fc = coulombForceFull(dr, alpha);
	Vector3 fh = softcoreForce(dr, eps, rad6);

	force[i] += fc + fh;
      }
    }
  }

  void computeTabulated(Vector3* force, const Vector3* pos, const int* type) {
    for (int i = 0; i < num; i++) {
      // Zero the force.
      force[i] = Vector3(0.0);
      
      // Loop through the neighbors.
      for (int n = 0; n < neigh[i].length(); n++) {
	int j = neigh[i].get(n);
	if (j == i) continue;
	int ind = type[i] + type[j]*numParts;
	if (tablePot[ind] == NULL) continue;
	Vector3 dr = sys->wrapDiff(pos[j] - pos[i]);
	
	if (dr.length2() > cutoff2) continue;
	Vector3 ft = tablePot[ind]->computeForce(dr);

	force[i] += ft;
      }
    }
  }

  void computeTabulatedFull(Vector3* force, const Vector3* pos, const int* type) {
    // Zero the force.
    for (int i = 0; i < num; i++) force[i] = Vector3(0.0);
    
    // Compute the force for all pairs.
    for (int i = 0; i < num-1; i++) {
      for (int j = i + 1; j < num; j++) {
	int ind = type[i] + type[j]*numParts;
	if (tablePot[ind] == NULL) continue;
	Vector3 dr = sys->wrapDiff(pos[j] - pos[i]);

	Vector3 ft = tablePot[ind]->computeForce(dr);

	force[i] += ft;
	force[j] -= ft;
      }
    }
  }

  static Vector3 coulombForce(Vector3 r, float alpha, float start, float len) {
    float d = r.length();

    if (d >= start + len) return Vector3(0.0);
    if (d <= start) {
      Vector3 force = -alpha/(d*d*d)*r;
      return force;
    }

    // Switching.
    float c = alpha/(start*start);
    Vector3 force = -c*(1.0 - (d - start)/len)/d*r;
    return force;
  }

  static Vector3 coulombForceFull(Vector3 r, float alpha) {
    float d = r.length();
    
    return -alpha/(d*d*d)*r;
  }

  static Vector3 softcoreForce(Vector3 r, float eps, float rad6) {
    const float d2 = r.length2();
    const float d6 = d2*d2*d2;
  
    if (d6 < rad6) return (-12*eps*(rad6*rad6/(d6*d6*d2) - rad6/(d6*d2)))*r;
    return Vector3(0.0);
  }

private:
  int num;
  int numParts;
  float** tableEps;
  float** tableRad6;
  float** tableAlpha;
  const BaseGrid* sys;
  IndexList* neigh;
  float switchStart, switchLen, electricConst, cutoff2;
  CellDecomposition decomp;
  TabulatedPotential** tablePot;
  int numTablePots;
  float energy;
};
#endif
