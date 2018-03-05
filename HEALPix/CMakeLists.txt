#PulsarPop HEALPix

set ( HEALPIX_LIB_C
chealpix.c
)

set ( HEALPIX_LIB_H
chealpix.h
)

set ( HEALPIX_CPP
DMContainer.cpp
)

set ( HEALPIX_H
DMContainer.hpp
)

add_library( HEALPixLib OBJECT ${HEALPIX_LIB_C} ${HEALPIX_LIB_H} )
add_library( HEALPix OBJECT ${HEALPIX_CPP} ${HEALPIX_H} )
