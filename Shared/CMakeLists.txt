#PulsarsPopulation Shared

set ( SHARED_CPP
Statistics.hpp
MyMath.hpp
strlcpy.hpp
)

set ( SHARED_H
Statistics.cpp
MyMath.cpp
strlcpy.cpp
)

add_library( Shared OBJECT ${SHARED_CPP} ${SHARED_H} )





