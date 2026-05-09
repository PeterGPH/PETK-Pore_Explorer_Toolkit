///////////////////////////////////////////////////////////////////////  
// Author: Jeff Comer <jcomer2@illinois.edu>

#include <cstdio>
#include <cstdlib>

#include "useful.h"
#include "BaseGrid.h"
#include "OverlordGrid.h"

//using namespace std;

int main(int argc, char* argv[]) {
  if ( argc != 4 ) {
    printf("Usage: %s systemDefinitionFile spacing outputFile\n", argv[0]);
    printf("You entered %i arguments.\n", argc-1);
    return 0;
  }
  const char* systemDefFile = argv[1];
  const float spacing = strtod(argv[2], NULL);

  String rootGrid = OverlordGrid::readDefFirst(systemDefFile);
  OverlordGrid* over = new OverlordGrid(rootGrid.val());
  int count = over->readDef(systemDefFile);
  printf("Found %d unique grids.\n", over->getUniqueGridNum());
  printf("Linked %d subgrids.\n", count);

  over->writeSubgrids("subgrid.txt");

  Matrix3 box(over->getBox());
  Vector3 org(over->getOrigin());
  BaseGrid sample(box, org, spacing);
  const int n = sample.length();

  printf("Sampling...\n");
  for (int i = 0; i < n; i++) {
    Vector3 r = sample.getPosition(i);
    float v = over->interpolatePotential(r);
    //Vector3 f = over->interpolateForce(r);
    sample.setValue(i, v);
  }
  sample.write(argv[argc-1]);

  delete over;
  return 0;
}
