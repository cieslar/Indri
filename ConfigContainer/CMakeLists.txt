#PulsarsPop ConfigContainer




set ( CONFIG_CPP
ConfigContainer.cpp
)

set ( CONFIG_H
ConfigContainer.hpp
)

add_library(ConfigContainer OBJECT ${CONFIG_CPP} ${CONFIG_H} )




