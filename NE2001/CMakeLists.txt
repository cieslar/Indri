#PulsarPop NE2001

set ( NE2001_CPP
NE2001.cpp
)

set ( NE2001_H
NE2001.hpp
)

set ( NE2001_F
density.NE2001.f
dmdsm.NE2001.f
NE2001.f
neclumpN.f
neclumpN.NE2001.f
neLISM.NE2001.f
nevoidN.f
nevoidN.NE2001.f
scattering98.f
)


add_library( NE2001_FORT OBJECT ${NE2001_F} )
add_library( NE2001 OBJECT ${NE2001_CPP} ${NE2001_H} ) 





