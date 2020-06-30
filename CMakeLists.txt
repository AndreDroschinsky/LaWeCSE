cmake_minimum_required(VERSION 3.1)

PROJECT(LaWeCSE CXX)

#
ADD_EXECUTABLE(LaWeCSE src/main.cpp src/auxiliary.cpp src/bbpmcsi.cpp src/graphVisualization.cpp src/hungarian.cpp src/MCIS2NK.cpp src/mwm.cpp src/TreeGen.cpp)

set(CMAKE_CXX_STANDARD 11)

set(OGDF_Path "/path_to_OGDF_release_build/OGDF")

include_directories(${CMAKE_SOURCE_DIR}/include/)
include_directories(${OGDF_Path}/include/)
include_directories(${OGDF_Path}/Release/include/)

target_link_libraries(LaWeCSE pthread)
target_link_libraries(LaWeCSE ${OGDF_Path}/Release/libOGDF.a)
target_link_libraries(LaWeCSE ${OGDF_Path}/Release/libCOIN.a)

