#pragma once

#ifdef __CUDACC__
    #define HOST __host__
    #define DEVICE __device__
#else
    #define HOST 
    #define DEVICE 
#endif

#include "useful.h"
#include <cuda.h>


class ProductPotentialConf {
public:
    ProductPotentialConf() {}
    ProductPotentialConf( std::vector< std::vector<int> > indices, std::vector<String> potential_names ) :
	indices(indices), potential_names(potential_names) { }

    std::vector< std::vector<int> > indices; /* indices of particles */
    std::vector<String> potential_names;

    inline ProductPotentialConf(const ProductPotentialConf& a) : indices(a.indices), potential_names(a.potential_names) { }

    
	/* String toString(); */
	/* void print(); */
};

