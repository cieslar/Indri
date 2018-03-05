#PulsarsPopulation

set ( POPULATION_CPP
Population.cpp
)

set ( POPULATION_H
Population.hpp
)

set ( PULSAR_CPP
Pulsar.cpp
)

set ( PULSAR_H
Pulsar.hpp
)

set ( NSNS_CPP
NSNS.cpp
)

set ( NSNS_HPP
NSNS.hpp
)

add_library( NSNS OBJECT ${NSNS_CPP} ${NSNS_H} )
#target_link_libraries( Pulsar PulsarPhysics )


add_library( Pulsar OBJECT ${PULSAR_CPP} ${PULSAR_H} )
#target_link_libraries( Pulsar PulsarPhysics )

add_library( Population OBJECT ${POPULATION_CPP} ${POPULATION_H} )
#target_link_libraries( Population Projections PulsarDynamics DataStructures
#Pulsar Galaxy Shared)



