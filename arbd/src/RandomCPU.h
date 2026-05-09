/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
 * Copyright (c) 1993 Martin Birgmeier
 * All rights reserved.
 *
 * You may redistribute unmodified or modified versions of this source
 * code provided that the above copyright notice and this and the
 * following conditions are retained.
 *
 * This software is provided ``as is'', and comes with no warranties
 * of any kind. I shall in no event be liable for anything that happens
 * to anyone/anything when using this software.
 */

/*––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––.
| Redundant with Random.h but with RandomCPU class since Random class is also |
| declared in RandomCUDA.h, and I want to use both for rigidBody code         |
`––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––*/
// RBTODO: make Random.h include RandomCPU.h and "typedef" this class

#ifndef RANDOM_CPU_H
#define RANDOM_CPU_H

#ifdef __CUDACC__
#define HOST __host__
#define DEVICE __device__
#else
#define HOST 
#define DEVICE 
#endif

#include "namd_common.h"
#include "useful.h"

#ifdef _MSC_VER
#define INT64_LITERAL(X) X ## i64
#else
#define INT64_LITERAL(X) X ## LL
#endif

#define	RAND48_SEED   INT64_LITERAL(0x00001234abcd330e)
#define	RAND48_MULT   INT64_LITERAL(0x00000005deece66d)
#define	RAND48_ADD    INT64_LITERAL(0x000000000000000b)
#define RAND48_MASK   INT64_LITERAL(0x0000ffffffffffff)

class RandomCPU {

private:

	float second_gaussian;
  bool second_gaussian_waiting;
  int64 rand48_seed;
  int64 rand48_mult;
  int64 rand48_add;

public:

  // default constructor
  RandomCPU(void) {
    init(0);
    rand48_seed = RAND48_SEED;
  }

  // constructor with seed
  RandomCPU(unsigned long seed) {
    init(seed);
  }

  // reinitialize with seed
  HOST DEVICE inline void init(unsigned long seed) {
    second_gaussian = 0;
    second_gaussian_waiting = false;
    rand48_seed = seed & INT64_LITERAL(0x00000000ffffffff);
    rand48_seed = rand48_seed << 16;
    rand48_seed |= RAND48_SEED & INT64_LITERAL(0x0000ffff);
    rand48_mult = RAND48_MULT;
    rand48_add = RAND48_ADD;
  }

  HOST DEVICE inline void print(const char* string)
		{
			printf(string);
			printf("RAND48_SEED = %lld, RAND48_MULT = %lld, RAND48_ADD = %lld, RAND48_MASK = %lld\n", RAND48_SEED, RAND48_MULT, RAND48_ADD, RAND48_MASK);
			printf("rand48_seed = %ld, rand48_mult = %ld, rand48_add = %ld\n", rand48_seed, rand48_mult, rand48_add);
		}

  // advance generator by one (seed = seed * mult + add, to 48 bits)
  HOST DEVICE inline void skip(void) {
    rand48_seed = ( rand48_seed * rand48_mult + rand48_add ) & RAND48_MASK;
  }

  // split into numStreams different steams and take stream iStream
  void split(int iStream, int numStreams) {
    int i;

    // make sure that numStreams is odd to ensure maximum period
    numStreams |= 1;

    // iterate to get to the correct stream
    for ( i = 0; i < iStream; ++i ) skip();

    // save seed and add so we can use skip() for our calculations
    int64 save_seed = rand48_seed;

    // calculate c *= ( 1 + a + ... + a^(numStreams-1) )
    rand48_seed = rand48_add;
    for ( i = 1; i < numStreams; ++i ) skip();
    int64 new_add = rand48_seed;

    // calculate a = a^numStreams
    rand48_seed = rand48_mult;
    rand48_add  = 0;
    for ( i = 1; i < numStreams; ++i ) skip();
    rand48_mult = rand48_seed;

    rand48_add  = new_add;
    rand48_seed = save_seed;

    second_gaussian = 0;
    second_gaussian_waiting = false;
    print("END SPLIT\n");
  }

  // return a number uniformly distributed between 0 and 1
  HOST DEVICE inline BigReal uniform(void) {
    skip();
    const float exp48 = ( 1.0 / (float)(INT64_LITERAL(1) << 48) );
    return ( (float) rand48_seed * exp48 );
  }

  long poisson(BigReal lambda) {
    const BigReal l = exp(-lambda);
    long k = 0;
    BigReal p = uniform();
    
    while (p >= l) {
      p *= uniform();
      k = k + 1;
    }

    return k;
  }

  // return a number from a standard gaussian distribution
  HOST DEVICE inline BigReal gaussian(void) {
    BigReal fac, r, v1, v2;

    if (second_gaussian_waiting) {
      second_gaussian_waiting = false;
      return second_gaussian;
    } else {
      r = 2.;                 // r >= 1.523e-8 ensures abs result < 6
      while (r >=1. || r < 1.523e-8) { // make sure we are within unit circle
        v1 = 2.0 * uniform() - 1.0;
        v2 = 2.0 * uniform() - 1.0;
        r = v1*v1 + v2*v2;
				// printf("r %f\n", r);
      }
      fac = sqrt(-2.0 * log(r)/r);
      // now make the Box-Muller transformation to get two normally
      // distributed random numbers. Save one and return the other.
      second_gaussian_waiting = true;
      second_gaussian = v1 * fac;
      return v2 * fac;
    }

  }

  // return a vector of gaussian random numbers
  HOST DEVICE inline Vector3 gaussian_vector(void) {
    return Vector3( gaussian(), gaussian(), gaussian() );
  }

  // return a random long
  HOST DEVICE inline long integer(void) {
    skip();
    return ( ( rand48_seed >> 17 ) & INT64_LITERAL(0x000000007fffffff) );
  }

  // randomly order an array of whatever
  template <class Elem> void reorder(Elem *a, int n) {
    for ( int i = 0; i < (n-1); ++i ) {
      int ie = i + ( integer() % (n-i) );
      if ( ie == i ) continue;
      const Elem e = a[ie];
      a[ie] = a[i];
      a[i] = e;
    }
  }

};

#endif  // RANDOM_CPU_H

