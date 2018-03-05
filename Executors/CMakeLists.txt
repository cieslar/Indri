#PulsarsPopulation

add_executable( Model Model.cpp )
target_link_libraries( Model PulsarPop ${OtherLibs} )


add_executable( NSNSFlier NSNSFlier.cpp )
target_link_libraries(NSNSFlier PulsarPop ${OtherLibs} )


add_executable( Pop1Per100  Pop1Per100.cpp )
target_link_libraries(Pop1Per100 PulsarPop ${OtherLibs} )

add_executable( GalacticTest GalacticTests.cpp )
target_link_libraries( GalacticTest PulsarPop ${OtherLibs} )

add_executable( MCMCMaker MCMCMaker.cpp )
target_link_libraries( MCMCMaker PulsarPop ${OtherLibs} )

add_executable( MCMCTester MCMCTester.cpp )
target_link_libraries( MCMCTester PulsarPop ${OtherLibs} )

add_executable( MCMCReadWrite MCMCReadWrite.cpp )
target_link_libraries( MCMCReadWrite PulsarPop ${OtherLibs} )


add_executable( MCMCChainCombiner MCMCChainCombiner.cpp )
target_link_libraries( MCMCChainCombiner PulsarPop ${OtherLibs} )

add_executable( ModelTester ModelTester.cpp )
target_link_libraries( ModelTester PulsarPop ${OtherLibs} )


add_executable( PopulationMaker PopulationMaker.cpp )
target_link_libraries( PopulationMaker PulsarPop ${OtherLibs} )

add_executable( DMMaker DMMaker.cpp )
target_link_libraries( DMMaker PulsarPop ${OtherLibs} )

add_executable( PopMeter PopMeter.cpp )
target_link_libraries( PopMeter PulsarPop ${OtherLibs} )

add_executable( DMChanger4YMW16 DMChanger4YMW16.cpp )
target_link_libraries( DMChanger4YMW16  PulsarPop ${OtherLibs} )

add_executable( XmlReadWrite XmlReadWrite.cpp )
target_link_libraries( XmlReadWrite PulsarPop ${OtherLibs} )
