#PulsarsPopulation

set ( RADIOOBSERVER_CPP
RadioDetector.cpp
Survey.cpp
)

set ( RADIOOBSERVER_H
RadioDetector.hpp
Survey.hpp
)

add_library( RadioDetector OBJECT ${RADIOOBSERVER_CPP} ${RADIOOBSERVER_H} )






