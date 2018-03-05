#PulsarsPop Pugixml

set ( PUGIXML_CPP
pugixml.cpp
)

set ( PUGIXML_H
pugiconfig.hpp  pugixml.hpp
)

add_library(Pugixml OBJECT ${PUGIXML_CPP} ${PUGIXML_H} )




