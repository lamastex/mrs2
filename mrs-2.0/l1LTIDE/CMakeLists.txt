# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( DensityTree )

set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "limited configs" FORCE) # Get around issue with CGAL not allowing config types other than build | release

set(CMAKE_VERBOSE_MAKEFILE 1)
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra -Wpedantic -g")  # enable C++11 standard

cmake_minimum_required(VERSION 2.8.11)

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()

endif()

# include helper file
include( ${CGAL_USE_FILE} )

# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()

endif()

# include for local directory
include_directories( include )
include_directories( ${MY_SOURCE_DIR}/src )
file (GLOB HEADERS include/*.h)

# include for local package


# Creating entries for target: BSP_Tree
# ############################

add_executable( BSP_Tree  hyperplane.cpp ${HEADERS} )

# Random number generation testing

# Define all the binary space files
add_library( bsp.o STATIC src/bsp.cpp include/DensityTree/bsp.hpp)
target_link_libraries( BSP_Tree bsp.o )

# General utilities for dealing with cgal points
add_library( PointUtils.o STATIC include/DensityTree/PointUtils.hpp src/PointUtils.cpp)
target_link_libraries( BSP_Tree PointUtils.o)

add_to_cached_list( CGAL_EXECUTABLE_TARGETS BSP_Tree )

# Link the executable to CGAL and third-party libraries
target_link_libraries(BSP_Tree ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )


set( CMAKE_EXPORT_COMPILE_COMMANDS 1 )

