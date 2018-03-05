#CuPIC Randomizer

set ( RANDOMIZER_CPP
Randomizer.cpp
)

set ( RANDOMIZER_H
MersenneTwister.h
Randomizer.hpp
)

add_library( Randomizer OBJECT ${RANDOMIZER_CPP} ${RANDOMIZER_H} )





