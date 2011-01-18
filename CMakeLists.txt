cmake_minimum_required(VERSION 2.6)

PROJECT(Komodakis)

include_directories(MRF)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

ADD_EXECUTABLE(Example Example.cxx Komodakis.cxx
MRF/mrf.cpp
#ICM.cpp
#GCoptimization.cpp graph.cpp maxflow.cpp
MRF/MaxProdBP.cpp MRF/LinkedBlockList.cpp MRF/regions-maxprod.cpp
#TRW-S.cpp BP-S.cpp
)
TARGET_LINK_LIBRARIES(Example
vtkHybrid ITKIO)

########################

ADD_EXECUTABLE(Tests Tests.cxx Komodakis.cxx
MRF/mrf.cpp
MRF/MaxProdBP.cpp MRF/LinkedBlockList.cpp MRF/regions-maxprod.cpp
)
TARGET_LINK_LIBRARIES(Tests vtkHybrid ITKIO)

