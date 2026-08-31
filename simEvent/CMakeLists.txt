cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(simEvent)
# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

set( ${PROJECT_NAME}_SRCS
        G4d2oGeom.cc
        ${PROJECT_NAME}.cc
)

set(${PROJECT_NAME}_HEADERS
        G4d2oGeom.h
        ${PROJECT_NAME}.h
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/include)
ROOT_GENERATE_DICTIONARY(${PROJECT_NAME}Dict
  ${${PROJECT_NAME}_HEADERS}
  MODULE ${PROJECT_NAME}
  LINKDEF ${PROJECT_NAME}_LinkDef.h
)

set(PCMFILE ${CMAKE_CURRENT_BINARY_DIR}/lib${PROJECT_NAME}_rdict.pcm)
set(RMAPFILE ${CMAKE_CURRENT_BINARY_DIR}/lib${PROJECT_NAME}.rootmap)

#---Create a shared library with geneated dictionary
add_library(${PROJECT_NAME} SHARED
  ${${PROJECT_NAME}_SRCS} ${${PROJECT_NAME}_HEADERS} ${PROJECT_NAME}Dict.cxx)

target_link_libraries(${PROJECT_NAME}
  PUBLIC ${ROOT_LIBRARIES})

############### install the library ###################
install(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
install(FILES ${${PROJECT_NAME}_HEADERS} DESTINATION include)
install(FILES ${PCMFILE} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
install(FILES ${RMAPFILE} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
