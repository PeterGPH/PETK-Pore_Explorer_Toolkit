///////////////////////////////////////////////////////////////////////  
// Author: Jeff Comer <jcomer2@illinois.edu>

#include <cstdio>
#include <cstdlib>

#include "useful.h"
#include "BaseGrid.h"
#include "FlowForce.h"

// using namespace std;

int main(int argc, char* argv[]) {
  if ( argc != 4 ) {
    printf("Usage: %s inGridFile diffusion outGridFile\n", argv[0]);
    printf("You entered %i arguments.\n", argc-1);
    return 0;
  }

  const char* inGrid = argv[1];
  const float diffusion = strtod(argv[2], NULL);
  const char* outGrid = argv[argc-1];

  BaseGrid sample(inGrid);
  FlowForce flow;
  const int n = sample.length();

  printf("Sampling...\n");
 for (int i = 0; i < n; i++) {
    Vector3 r = sample.getPosition(i);
    Vector3 f = flow.force(r,diffusion);
    sample.setValue(i, f.x);
  }
  sample.write(outGrid);

  return 0;
}
