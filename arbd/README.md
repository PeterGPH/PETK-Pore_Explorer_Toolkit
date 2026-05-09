# Atomic Resolution Brownian Dynamics (ARBD) - May 24

Brownian dynamics (BD) simulation is method for studying biomolecules,
ions, and nanomaterials that balances detail with computational
efficiency.

ARBD supports tabulated non-bonded and bonded interactions between BD
particles that can also be influenced by grid-specified
potentials. Uniquely, ARBD also allows grid-specified densities and
potentials to be associated with rigid body particles that rotate and
translate to represent larger molecules. Most importantly, the code is
designed to run quickly on modern NVIDIA GPUs.

ARBD is a rewrite of the BrownianMover code, moving almost all
computations to the GPU and enabling grid-specified particle
models. Please be aware that ARBD is being actively developed and is
offered without warranty.


## Building

### Dependencies

Linux workstation with CUDA-compatible GPU (minimum 3.5 compute capability)
  - CMake >= 3.9
  - gcc >= 4.9
  - cuda >= 9.0  (> 11.5 recommended)

### Build process

From the root arbd directory (where this README is found), run:
```
## Determine the compute capability of your CUDA-enabled graphics card
export CMAKE_CUDA_ARCHITECTURES="35;50;75;80"   ;# especially important for CMake < 3.24.0
## export CUDA_INCLUDE_DIRS="$CUDA_HOME/include" ;# optionally be explicit about cuda include paths; usually not needed
cmake -S . -B build &&
(
  cd build
  make -j
)
```

If your CUDA toolkit is installed in a nonstandard location that CMake
is unable to find, you may provide use the environement variable
`CMAKE_CUDA_COMPILER` to specify the path to nvcc. You may also find
it neccesary to set the environment variable `CUDA_INCLUDE_DIRS` if
compilation fails due to the compiler being unable to find <cuda.h>.

Note that ARBD has been developed using CUDA-9.0 and targets NVIDIA
GPUs featuring 6.0 compute capability. The code should work with
devices with compute capability >=2.0, but there are no guarantees.

## Usage

Please explore the examples in the 'tests' directory.

For example, try the following commands:

cd tests/argon-small
mkdir output
../../src/arbd BrownDyn.bd output/BrownDyn > output/BrownDyn.log

You may use the '-g n' option to specify the n-th GPU on your machine,
counting from 0.

## Citing

If you publish results obtained using ARBD, please cite the
following manuscripts:

"Predicting the DNA sequence dependence of nanopore ion current using atomic-resolution Brownian dynamics"
Jeffrey Comer and Aleksei Aksimentiev.
J Phys Chem C Nanomater Interfaces 116:3376-3393 (2012).

## Authors

ARBD is developed by the Aksimentiev group
(http://bionano.physics.illinois.edu).

Please direct questions or problems to Chris.

- Christopher Maffeo <cmaffeo2@illinois.edu>
- Han-yi Chao
- Jeffrey Comer
- Max Belkin
- Emmanual Guzman
- Justin Dufresne
- Terrance Howard

## Outstanding issues

* There are no checks to ensure that pairlists are recalculated before
  particles further than the pairlist distance move to within the
  cutoff

* A large amount of GPU memory for pairlists is allocated statically,
  which may cause out-of-memory crashes in older hardware
