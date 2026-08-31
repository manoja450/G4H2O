cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(G4d2o)

# Find Geant4 and ROOT
find_package(Geant4 REQUIRED ui_all vis_all)
find_package(ROOT REQUIRED)

# Include Geant4 and ROOT directories
include(${Geant4_USE_FILE})
include(${ROOT_USE_FILE})

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CONFIGURATION_TYPES Release;Debug;MinSizeRel;RelWithDebInfo)
MESSAGE( STATUS "CMAKE_BUILD_TYPE:      " ${CMAKE_BUILD_TYPE} )
MESSAGE( STATUS "CMAKE_INSTALL_PREFIX:  " ${CMAKE_INSTALL_PREFIX} )

SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_MACOSX_RPATH 1)

message( STATUS "ROOT_USE_FILE:         " ${ROOT_USE_FILE})
message( STATUS "PROJECT_SOURCE_DIR:    " ${PROJECT_SOURCE_DIR})
message( STATUS "ROOT_INCLUDE_DIR:      " ${ROOT_INCLUDE_DIR})

# IMPORTANT: Use Geant4_INCLUDE_DIRS (plural) not Geant4_INCLUDE_DIR
# This is the correct variable name after find_package(Geant4)
set(INCLUDE_DIRECTORIES
        ${ROOT_INCLUDE_DIR}
        ${Geant4_INCLUDE_DIRS}
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/simEvent
        ${PROJECT_SOURCE_DIR}/include
        $ENV{CRYHOME}/src
)

include_directories( ${INCLUDE_DIRECTORIES})

add_subdirectory(simEvent)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh ${PROJECT_SOURCE_DIR}/include/*.h)

if(DEFINED ENV{CRYHOME})
  MESSAGE( STATUS "CRYHOME:  " $ENV{CRYHOME} )
  add_definitions(-DUSE_CRY)
  file(GLOB crylibs $ENV{CRYHOME}/lib/libCRY.a)
  add_executable(G4d2o G4d2o.cc ${sources} ${headers})
  target_link_libraries(G4d2o ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} simEvent ${crylibs})
else(DEFINED ENV{CRYHOME})
  message( STATUS "CRYHOME environment variable not defined")
  add_executable(G4d2o G4d2o.cc ${sources} ${headers})
  target_link_libraries(G4d2o ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} simEvent)
endif(DEFINED ENV{CRYHOME})

############# install ####################

install(TARGETS G4d2o DESTINATION ..)
