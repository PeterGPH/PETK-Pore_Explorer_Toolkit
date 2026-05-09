// Configuration file reader
// Author: Jeff Comer <jcomer2@illinois.edu>
#ifndef RESERVOIR_H
#define RESERVOIR_H

#define STRLEN 512 

#include "useful.h"

class Reservoir {
public:
  Reservoir(const char* reservoirFile);
  ~Reservoir();
	
  static int countReservoirs(const char* reservoirFile);
	
  Vector3 getOrigin(int i) const;
  Vector3 getDestination(int i) const;
  Vector3 getDifference(int i) const;

  float getMeanNumber(int i) const;
  int length() const;

  bool inside(int i, Vector3 r) const;

private:
  int reservoirs;
  Vector3* r0;
  Vector3* r1;
  float* num;
	
  void readReservoirs(const char* reservoirFile);
  void validateRegions();
};
#endif
