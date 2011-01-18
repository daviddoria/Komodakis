cmake_minimum_required(VERSION 2.6)

PROJECT(VTKExamples)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

ADD_EXECUTABLE(ExampleEnergy ExampleEnergy.cpp mrf.cpp ICM.cpp GCoptimization.cpp graph.cpp maxflow.cpp
       MaxProdBP.cpp LinkedBlockList.cpp regions-maxprod.cpp
       TRW-S.cpp BP-S.cpp)
TARGET_LINK_LIBRARIES(ExampleEnergy vtkHybrid)

# ADD_EXECUTABLE(example example.cpp mrf.cpp ICM.cpp GCoptimization.cpp graph.cpp maxflow.cpp
#        MaxProdBP.cpp LinkedBlockList.cpp regions-maxprod.cpp
#        TRW-S.cpp BP-S.cpp)
# TARGET_LINK_LIBRARIES(example vtkHybrid)
#
# ADD_EXECUTABLE(ExampleICM ExampleICM.cpp mrf.cpp ICM.cpp GCoptimization.cpp graph.cpp maxflow.cpp
#        LinkedBlockList.cpp regions-maxprod.cpp
#        ExampleEnergyFunctions.cpp)
# TARGET_LINK_LIBRARIES(ExampleICM vtkHybrid)
#
# ADD_EXECUTABLE(ExampleSwap ExampleSwap.cpp mrf.cpp ICM.cpp GCoptimization.cpp graph.cpp maxflow.cpp
#        MaxProdBP.cpp LinkedBlockList.cpp regions-maxprod.cpp
#        TRW-S.cpp BP-S.cpp ExampleEnergyFunctions.cpp)
# TARGET_LINK_LIBRARIES(ExampleSwap vtkHybrid)
#
# ADD_EXECUTABLE(ExampleBP ExampleBP.cpp mrf.cpp ICM.cpp GCoptimization.cpp graph.cpp maxflow.cpp
#        MaxProdBP.cpp LinkedBlockList.cpp regions-maxprod.cpp
#        TRW-S.cpp BP-S.cpp ExampleEnergyFunctions.cpp)
# TARGET_LINK_LIBRARIES(ExampleBP vtkHybrid)
#
# ADD_EXECUTABLE(ExampleTRW ExampleTRW.cpp mrf.cpp ICM.cpp GCoptimization.cpp graph.cpp maxflow.cpp
#        MaxProdBP.cpp LinkedBlockList.cpp regions-maxprod.cpp
#        TRW-S.cpp BP-S.cpp ExampleEnergyFunctions.cpp)
# TARGET_LINK_LIBRARIES(ExampleTRW vtkHybrid)
#
# ADD_EXECUTABLE(ExampleBPS ExampleBPS.cpp mrf.cpp ICM.cpp GCoptimization.cpp graph.cpp maxflow.cpp
#        MaxProdBP.cpp LinkedBlockList.cpp regions-maxprod.cpp
#        TRW-S.cpp BP-S.cpp ExampleEnergyFunctions.cpp)
# TARGET_LINK_LIBRARIES(ExampleBPS vtkHybrid)
