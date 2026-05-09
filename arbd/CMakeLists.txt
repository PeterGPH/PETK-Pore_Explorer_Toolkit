## Specify the project
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)

# option(USE_CUDA "Use CUDA" ON)
set(USE_CUDA ON)

if(${USE_CUDA} AND NOT CMAKE_CUDA_ARCHITECTURES)
  if(${CMAKE_VERSION} VERSION_LESS "3.18.0")
    message("CMake support for CUDA prior to CMake 3.18 may not work well; consider switching to a more recent version of CMake if you encounter errors at compile- or run-time.")
  endif()
  if(${CMAKE_VERSION} VERSION_LESS "3.24.0") 
    message("CMAKE_CUDA_ARCHITECTURES not specified; checking GPUs on your system for CUDA architecture...")
    execute_process(COMMAND nvidia-smi --query-gpu=compute_cap --format=csv,noheader
    COMMAND sort
    COMMAND uniq
    COMMAND tr "\\n" " "
    COMMAND sed "s/ \\(.\\)/,\\1/g"
    TIMEOUT 1
    RESULT_VARIABLE GPU_QUERY_RESULT
    OUTPUT_VARIABLE GPU_QUERY_OUTPUT
    )
    if (GPU_QUERY_RESULT EQUAL 0)
      set(CMAKE_CUDA_ARCHITECTURES ${GPU_QUERY_OUTPUT})
    else()
      message("CUDA architecture could not be found! Guessing SM 7.5")  
      set(CMAKE_CUDA_ARCHITECTURES 7.5)
    endif()
  else()
    set(CMAKE_CUDA_ARCHITECTURES native)
  endif()
endif()

## specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set a default build type if none was specified
set(default_build_type "Release")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
    STRING "Choose the type of build." FORCE)
  ## Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "RelWithDebInfo" "MinSizeRel")
endif()


if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  ## Some CMake versions use -O2 for release; replace with -O3
  if(CMAKE_CXX_FLAGS_RELEASE)
    string(REPLACE "-O2" "-O3" CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
  endif()
  if(CMAKE_CUDA_FLAGS_RELEASE)
    string(REPLACE "-O2" "-O3" CMAKE_CUDA_FLAGS_RELEASE ${CMAKE_CUDA_FLAGS_RELEASE})
  endif()
endif()


## OPTIONS
option(DEBUG "Build with debug flags" False)
option(USE_NCCL "Use NCCL for single node GPU peer communication" False)
option(USE_NVTX "Build with NVTX profiling ranges" False)

# (not yet optional) message(STATUS "USE_CUDA: ${USE_CUDA}")
message(STATUS "DEBUG: ${DEBUG}")
if(${USE_CUDA})
  message(STATUS "CMAKE_CUDA_ARCHITECTURES: ${CMAKE_CUDA_ARCHITECTURES}")
  message(STATUS "USE_NCCL: ${USE_NCCL}")
endif()
message(STATUS "USE_NVTX: ${USE_NVTX}")


## Set flags before adding executable 
# configure_file(TutorialConfig.h.in TutorialConfig.h)
if(USE_CUDA)
  add_definitions(-DUSE_CUDA)
  ## CUDA_INCLUDE_DIRS wasn't getting set on my system with cmake 3.14.1, so check if in env
  if(DEFINED ENV{CUDA_INCLUDE_DIRS})
    set(CUDA_INCLUDE_DIRS $ENV{CUDA_INCLUDE_DIRS})
  endif()
  # message(STATUS "CUDA_INC: ${CUDA_INCLUDE_DIRS}")
  include_directories(${CUDA_INCLUDE_DIRS})
endif()

if(DEBUG)
  set(CMAKE_BUILD_TYPE Debug)
  set(PROJECT_NAME arbd_dbg)
else(DEBUG)
  set(PROJECT_NAME arbd)
endif()
if(USE_NVTX)
  add_definitions(-DUSE_NVTX)
endif()
if(USE_NCCL)
  add_definitions(-DUSE_NCCL)
  target_link_libraries("${PROJECT_NAME}" PRIVATE nccl)
endif()

## Two lines below needed?
set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(USE_CUDA)
  project("${PROJECT_NAME}" VERSION 1.2 LANGUAGES CXX CUDA)
else(USE_CUDA)
  project("${PROJECT_NAME}" VERSION 1.2 LANGUAGES CXX)
endif()


## Print all variables by uncommenting block below 
# get_cmake_property(_variableNames VARIABLES)
# list (SORT _variableNames)
# foreach (_variableName ${_variableNames})
#     message(STATUS "${_variableName}=${${_variableName}}")
# endforeach()

## Set rpath
set(CMAKE_MACOSX_RPATH 1)	# Unsure if this works for CMAKE_BUIlD_RPATH, or just CMAKE_INSTALL_RPATH
set(CMAKE_BUILD_RPATH "${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES}")

# set(CMAKE_VERBOSE_MAKEFILE True)
add_executable("${PROJECT_NAME}" src/arbd.cpp
  src/Configuration.cpp
  src/FlowForce.cpp
  src/GPUManager.cpp
  src/Scatter.cpp
  src/SignalManager.cpp
  src/WKFUtils.cpp
  src/Angle.cu
  src/BaseGrid.cu
  src/BrownianParticle.cu
  src/BrownianParticleType.cpp
  src/RigidBodyController.cu
  src/RigidBody.cu
  src/RigidBodyGrid.cu
  src/RigidBodyType.cu
  src/CellDecomposition.cu
  src/ComputeForce.cuh
  src/ComputeForce.cu
  src/ComputeGridGrid.cuh
  src/ComputeGridGrid.cu
  src/CudaUtil.cu
  src/CudaUtil.cuh
  src/Dihedral.cu
  src/Exclude.cu
  src/GrandBrownTown.cu
  src/GrandBrownTown.cuh
  src/imd.cpp
  src/vmdsock.cpp
  src/JamesBond.cu
  src/RandomCUDA.cu
  src/Reservoir.cu
  src/TabulatedAngle.cu
  src/TabulatedDihedral.cu
  src/TabulatedMethods.cuh
  src/TabulatedPotential.cu
  src/useful.cu
  )

## Add optional libraries
if(USE_CUDA)
  target_link_libraries("${PROJECT_NAME}" PRIVATE curand)
endif()
if(USE_NCCL)
  target_link_libraries("${PROJECT_NAME}" PRIVATE nccl)
endif()
if(USE_NVTX)
  target_link_libraries("${PROJECT_NAME}" PRIVATE nvToolsExt)
endif()

install(TARGETS "${PROJECT_NAME}")
