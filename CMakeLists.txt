cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(MC CXX)

# Set module path
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)

# Make sure module directory is included when building
include_directories(${CMAKE_SOURCE_DIR}/include)

#==============================================================================
# Command line options
#==============================================================================

option(profile  "Compile with profiling flags"                   OFF)
option(debug    "Compile with debug flags"                       OFF)
option(optimize "Turn on all compiler optimization flags"        OFF)
option(test     "Include all tests"                              OFF)

#==============================================================================
# Set compile flags
#==============================================================================

list(APPEND cxxflags -std=c++11 -O2)

if(profile)
  list(APPEND cxxflags -pg)
endif()
if(optimize)
  list(REMOVE_ITEM cxxflags -O2)
  list(APPEND cxxflags -O3)
endif()
if(debug)
  list(REMOVE_ITEM cxxflags -O2)
  list(REMOVE_ITEM cxxflags -O3)
  list(APPEND cxxflags -g -O0)
endif()

# Show flags being used
message(STATUS "C++ flags: ${cxxflags}")

# Convert to string
string(REPLACE ";" " " cxxflags "${cxxflags}")

add_definitions(${cxxflags})

#==============================================================================
# pugixml library
#==============================================================================

add_library(pugixml src/pugixml/pugixml.cpp)

#==============================================================================
# HDF5 for binary output
#==============================================================================

find_package(HDF5 COMPONENTS CXX)
if(NOT HDF5_FOUND)
    message(FATAL_ERROR "Could not find HDF5")
endif()

#==============================================================================
# Define the MC program executable
#==============================================================================

file(GLOB libMC_SRC "src/*.cpp")
set(libMC_SRC
    src/Algorithm.cpp
    src/Distribution.cpp
    src/Entropy.cpp
    src/Estimator.cpp
    src/Geometry.cpp
    src/Material.cpp
    src/Nuclide.cpp
    src/Particle.cpp
    src/Random.cpp
    src/Reaction.cpp
    src/Source.cpp
    src/XSec.cpp
    src/simulator/fixed_source.cpp
    src/simulator/general.cpp
    src/simulator/handler.cpp
    src/simulator/ksearch.cpp
    src/simulator/report.cpp
    src/simulator/report.cpp
    src/simulator/setup.cpp
    src/simulator/time_dependent.cpp
)
add_library(libMC ${libMC_SRC})

set(program "MC.exe")
add_executable(${program} Main.cpp)

target_include_directories(libMC PUBLIC ${HDF5_INCLUDE_DIRS})
# Add HDF5 library directories to link line with -L
foreach(LIBDIR ${HDF5_LIBRARY_DIRS})
    list(APPEND ldflags "-L${LIBDIR}")
endforeach()
target_link_libraries(libMC ${ldflags} ${HDF5_LIBRARIES} pugixml)
target_link_libraries(${program} ${ldflags} libMC)

#==============================================================================
# Define the TRMM program executable
#==============================================================================

set(program "TRMM.exe")
add_executable(${program} TRMM.cpp)

target_include_directories(${program} PUBLIC ${HDF5_INCLUDE_DIRS})
# Add HDF5 library directories to link line with -L
target_link_libraries(${program} ${ldflags} ${HDF5_LIBRARIES})

#==============================================================================
# Copy important files
#==============================================================================

# Copy Cross-Section library (text files) and input examples
file(COPY ${CMAKE_SOURCE_DIR}/xs_library DESTINATION ${CMAKE_BINARY_DIR})
file(COPY ${CMAKE_SOURCE_DIR}/examples DESTINATION ${CMAKE_BINARY_DIR})

#==============================================================================
# Gtest
#==============================================================================

if(test)
    message(STATUS "Test: ON")
    enable_testing()
    find_package(GTest REQUIRED)
    add_executable(unit_test.exe ${CMAKE_SOURCE_DIR}/test/unit_test.cpp)
    target_link_libraries(unit_test.exe GTest::GTest libMC)
    add_test(NAME UnitTest COMMAND unit_test.exe)
else()
    message(STATUS "Test: OFF")
endif()
