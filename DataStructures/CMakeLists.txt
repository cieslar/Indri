#PulsarPop DataStructures

set ( DATASTRUCTURES_CPP
Coordinates.cpp
ModelParameters.cpp
BasicStar.cpp
)

set ( DATASTRUCTURES_H
Coordinates.hpp
ModelParameters.hpp
BasicStar.hpp
)

add_library( DataStructures OBJECT ${DATASTRUCTURES_CPP} ${DATASTRUCTURES_H} )





