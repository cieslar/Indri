#PulsarsPop NeutronStar

set ( PULSARPHYSICS_CPP
PulsarPhysics.cpp
)

set ( PULSARPHYSICS_H
PulsarPhysics.hpp
)

add_library( PulsarPhysics OBJECT ${PULSARPHYSICS_CPP} ${PULSARPHYSICS_H} )





