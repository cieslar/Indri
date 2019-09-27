#top place to be

project( PulsarsPopulation )

#SET(BUILD_SHARED_LIBRARIES FALSE)

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake")
#set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake.modules")

#Zacznijmy od skompilowania jakkolwiek.
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Release Debug"
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

#set(BUILD_SHARED_LIBS ON)

#add_definitions(-DGLM_COMPILER=GLM_COMPILER_GCC44)
#add_definitions(-DGLM_COMPILER=GLM_COMPILER_CUDA30)
#add_definitions(-D__CUDACC__)
#add_definitions(-DGLM_FORCE_CUDA )
#add_definitions(-DGLM_MESSAGES)

#find_package(CUDA REQUIRED)

#if(WIN32 AND NOT UNIX)
#set(ADDITIONAL_FIND_PATHS "d:/libs/fftw/" )
#to powinno isc przez zmienne:
#CMAKE_INCLUDE_PATH 
#CMAKE_LIBRARY_PATH 
#endif(WIN32 AND NOT UNIX)


FIND_PACKAGE(HDF5 1.8.6 REQUIRED COMPONENTS CXX)
#find_package(hdf5 NAMES h5c++ COMPONENTS CXX static)

#-static-libgcc -static -s



#set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS} -O2 -ansi --pedantic -std=c++0x -pg  -static -mtune=amdfam10")
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O2 -ansi --pedantic -std=c++0x -pg -static -mtune=amdfam10")
#set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS} -O2 -ansi --pedantic -std=c11 -mtune=native -march=native ")
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O2 -ansi --pedantic -std=c++11 -mtune=native -march=native ")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS} -O2 -ansi --pedantic -std=c11 -mtune=native -march=native ")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O2 -ansi --pedantic -std=c++14 -mtune=native -march=native ")
#set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS} -O1 -ansi --pedantic -std=c++0x")
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O1 -ansi --pedantic -std=c++0x")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS} -O0 -ansi --pedantic -std=c11 -g -ggdb")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -O0 -ansi --pedantic -std=c++11 -g -ggdb")
#if(UNIX AND NOT WIN32)
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
#endif(UNIX AND NOT WIN32)
#set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -arch sm_20 -ccbin g++" )
#nie sprawdzilem jeszcze ktore dziala
#set(CUDA_NVCC_FLAGS "-arch=sm_20" CACHE STRING "nvcc flags" FORCE)

#co to robi?
#mark_as_advanced(CLEAR CUDA_64_BIT_DEVICE_CODE)



#-Wno-write-strings


#Fortran stuff:
enable_language (Fortran)
# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  # set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -std=legacy -ffixed-form -ffixed-line-length-none -pg")
  set (CMAKE_Fortran_FLAGS_RELEASE "-fno-f2c -O2 -std=legacy -ffixed-form -ffixed-line-length-none -pg")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0  -std=legacy -ffixed-form -ffixed-line-length-none  -g -ggdb")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try
-O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")
















set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")


include_directories( ${HDF5_INCLUDE_DIRS} )

include_directories( ${CMAKE_SOURCE_DIR}/DataStructures )
include_directories( ${CMAKE_SOURCE_DIR}/Shared )
include_directories( ${CMAKE_SOURCE_DIR}/Tools )
include_directories( ${CMAKE_SOURCE_DIR}/Randomizer )
include_directories( ${CMAKE_SOURCE_DIR}/ConfigContainer )
include_directories( ${CMAKE_SOURCE_DIR}/Galaxy )
include_directories( ${CMAKE_SOURCE_DIR}/Projections )
include_directories( ${CMAKE_SOURCE_DIR}/Executors )
include_directories( ${CMAKE_SOURCE_DIR}/NeutronStar )
include_directories( ${CMAKE_SOURCE_DIR}/NE2001 )
include_directories( ${CMAKE_SOURCE_DIR}/RadioDetector )
include_directories( ${CMAKE_SOURCE_DIR}/Population )
include_directories( ${CMAKE_SOURCE_DIR}/Catalogue )
include_directories( ${CMAKE_SOURCE_DIR}/MCMC )
include_directories( ${CMAKE_SOURCE_DIR}/Pugixml )
include_directories( ${CMAKE_SOURCE_DIR}/HEALPix )


link_directories( ${CMAKE_SOURCE_DIR}/DataStructures )
link_directories( ${CMAKE_SOURCE_DIR}/Shared )
link_directories( ${CMAKE_SOURCE_DIR}/Randomizer )
link_directories( ${CMAKE_SOURCE_DIR}/ConfigContainer )
link_directories( ${CMAKE_SOURCE_DIR}/Galaxy )
link_directories( ${CMAKE_SOURCE_DIR}/Projections )
link_directories( ${CMAKE_SOURCE_DIR}/Executors )
link_directories( ${CMAKE_SOURCE_DIR}/NeutronStar )
link_directories( ${CMAKE_SOURCE_DIR}/NE2001 )
link_directories( ${CMAKE_SOURCE_DIR}/RadioDetector )
link_directories( ${CMAKE_SOURCE_DIR}/Population )
link_directories( ${CMAKE_SOURCE_DIR}/Tools )
link_directories( ${CMAKE_SOURCE_DIR}/Catalogue )
link_directories( ${CMAKE_SOURCE_DIR}/MCMC )
link_directories( ${CMAKE_SOURCE_DIR}/Pugixml )
link_directories( ${CMAKE_SOURCE_DIR}/HEALPix )


add_subdirectory( ${CMAKE_SOURCE_DIR}/DataStructures )
add_subdirectory( ${CMAKE_SOURCE_DIR}/Shared )
add_subdirectory( ${CMAKE_SOURCE_DIR}/Randomizer )
add_subdirectory( ${CMAKE_SOURCE_DIR}/ConfigContainer )
add_subdirectory( ${CMAKE_SOURCE_DIR}/Galaxy )
add_subdirectory( ${CMAKE_SOURCE_DIR}/Projections )
add_subdirectory( ${CMAKE_SOURCE_DIR}/NeutronStar )
add_subdirectory( ${CMAKE_SOURCE_DIR}/NE2001 )
add_subdirectory( ${CMAKE_SOURCE_DIR}/RadioDetector )
add_subdirectory( ${CMAKE_SOURCE_DIR}/Population )
add_subdirectory( ${CMAKE_SOURCE_DIR}/Tools )
add_subdirectory( ${CMAKE_SOURCE_DIR}/Catalogue )
add_subdirectory( ${CMAKE_SOURCE_DIR}/MCMC )
add_subdirectory( ${CMAKE_SOURCE_DIR}/Pugixml )
add_subdirectory( ${CMAKE_SOURCE_DIR}/HEALPix )

set( OtherLibs  ${HDF5_CXX_LIBRARIES})

add_library(PulsarPop STATIC
	$<TARGET_OBJECTS:DataStructures>
	$<TARGET_OBJECTS:Shared>
	$<TARGET_OBJECTS:Randomizer>
	$<TARGET_OBJECTS:ConfigContainer>
	$<TARGET_OBJECTS:Galaxy>
	$<TARGET_OBJECTS:PulsarDynamics>
	$<TARGET_OBJECTS:Projections>
	$<TARGET_OBJECTS:PulsarPhysics>
	$<TARGET_OBJECTS:NE2001>
	$<TARGET_OBJECTS:NE2001_FORT>
	$<TARGET_OBJECTS:RadioDetector>
	$<TARGET_OBJECTS:Pulsar>
	$<TARGET_OBJECTS:NSNS>
	$<TARGET_OBJECTS:Population>
	$<TARGET_OBJECTS:Tools>
	$<TARGET_OBJECTS:Catalogue>
	$<TARGET_OBJECTS:MCMC>
	$<TARGET_OBJECTS:Pugixml>
	$<TARGET_OBJECTS:HEALPix>
	$<TARGET_OBJECTS:HEALPixLib>
)

#add_subdirectory( ${CMAKE_SOURCE_DIR}Tests )
add_subdirectory( ${CMAKE_SOURCE_DIR}/Executors )

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/bin/input )
#file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/bin/input/NE2001 )

file(COPY "${CMAKE_SOURCE_DIR}/NE2001/input.NE2001" DESTINATION "${CMAKE_BINARY_DIR}/bin/input")
file(COPY "${CMAKE_SOURCE_DIR}/Scripts" DESTINATION "${CMAKE_BINARY_DIR}/bin/")
file(COPY "${CMAKE_SOURCE_DIR}/RadioDetector/tsky1.ascii" DESTINATION "${CMAKE_BINARY_DIR}/bin/input")
file(COPY "${CMAKE_SOURCE_DIR}/Input/ATNF_MOD1.txt" DESTINATION "${CMAKE_BINARY_DIR}/bin/input")
file(RENAME "${CMAKE_BINARY_DIR}/bin/input/input.NE2001" "${CMAKE_BINARY_DIR}/bin/input/NE2001")

#FC     = gfortran
#FC = g77
#CXX = g++ 

#FFLAGS = -O3 -std=legacy -ffixed-form -ffixed-line-length-none -pg -g -ggdb
#FFLAGS = -O3 -std=legacy -ffixed-form -ffixed-line-length-none -fPIC -pg
#FFLAGS = -fPIC -pg

#INCL = -I../ConfigContainer/ -I../DataStructures/ -I../Galaxy/ -I../NE2001/ -I../NeutronStar/ -I../Population/ -I../Projections/ -I../Randomizer/ -I../Shared/ -I../Tools/ -I../RadioDetector/
#CXXFLAGS = -W -Wall -ansi -pedantic-errors -std=c++0x -fPIC
#CXXFLAGS = -O3 -fPIC -g -ggdb -ansi -pedantic-errors -std=c++0x -pg 
#LIBS =  -lgfortran -L../lib -lPulsar
#LIBS =  -lgfortran 
#LFLAGS = -pg -g -ggdb
#LIBS =  -L../lib -lPulsar  -lm -lc -lg2c
