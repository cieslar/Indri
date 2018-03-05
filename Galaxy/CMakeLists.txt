#PulsarsPop Galaxy

set ( GALAXY_CPP
Galaxy.cpp
)

set ( GALAXY_H
Galaxy.hpp
)


#add_library( Galaxy OBJECT ${GALAXY_CPP} ${GALAXY_H} )


set ( PULSARDYNAMICS_CPP
PulsarDynamics.cpp
)
set ( PULSARDYNAMICS_H
PulsarDynamics.hpp
)

add_library( Galaxy OBJECT ${GALAXY_CPP} ${GALAXY_H} )
add_library( PulsarDynamics OBJECT ${PULSARDYNAMICS_CPP} ${PULSARDYNAMICS_H})

#add_library( Galaxy OBJECT ${GALAXY_CPP} ${GALAXY_H} ${PULSARDYNAMICS_CPP} ${PULSARDYNAMICS_H})



