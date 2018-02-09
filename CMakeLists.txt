cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(MC CXX)

# Set module path
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)

# Make sure module directory is included when building
include_directories(${CMAKE_SOURCE_DIR}/include)

#===============================================================================
# Command line options
#===============================================================================

option(profile  "Compile with profiling flags"                   OFF)
option(debug    "Compile with debug flags"                       OFF)
option(optimize "Turn on all compiler optimization flags"        OFF)

#===============================================================================
# Set compile flags
#===============================================================================

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

#===============================================================================
# pugixml library
#===============================================================================

add_library(pugixml src/pugixml/pugixml.cpp)

#===============================================================================
# HDF5 for binary output
#===============================================================================
find_package(HDF5 COMPONENTS CXX)
if(NOT HDF5_FOUND)
    message(FATAL_ERROR "Could not find HDF5")
endif()

#===============================================================================
# Define the MC program executable
#===============================================================================

file(GLOB libMC_SRC "src/*.cpp")
add_library(libMC ${libMC_SRC})

set(program "MC.exe")
add_executable(${program} Main.cpp)

#target_include_directories(libMC PUBLIC ${HDF5_INCLUDE_DIRS})
# Add HDF5 library directories to link line with -L
#foreach(LIBDIR ${HDF5_LIBRARY_DIRS})
#    list(APPEND ldflags "-L${LIBDIR}")
#endforeach()
#target_link_libraries(libMC ${ldflags} ${HDF5_LIBRARIES} pugixml)

# I don't know, I just need these
target_include_directories(libMC PUBLIC /usr/local/HDF_Group/HDF5/1.10.1/include)
set(ldflags "-I/usr/local/HDF_Group/HDF5/1.10.1/include -DBOOST_MATH_DISABLE_FLOAT128 -m64 -L/usr/local/HDF_Group/HDF5/1.10.1/lib /usr/local/HDF_Group/HDF5/1.10.1/lib/libhdf5_hl_cpp-static.a /usr/local/HDF_Group/HDF5/1.10.1/lib/libhdf5_cpp-static.a /usr/local/HDF_Group/HDF5/1.10.1/lib/libhdf5_hl-static.a /usr/local/HDF_Group/HDF5/1.10.1/lib/libhdf5-static.a -L/usr/local/HDF_Group/HDF5/1.10.1/lib -Wl,-rpath,/usr/local/HDF_Group/HDF5/1.10.1/lib -L/usr/local/HDF_Group/HDF5/1.10.1/lib -lrt -lpthread -lz -ldl -lm -Wl,-rpath -Wl,/usr/local/HDF_Group/HDF5/1.10.1/lib")
target_link_libraries(libMC ${ldflags} pugixml)

target_link_libraries(${program} ${ldflags} libMC)

#===============================================================================
# Copy important files
#===============================================================================

# Copy Cross-Section library (text files) and input examples
file(COPY ${CMAKE_SOURCE_DIR}/xs_library DESTINATION ${CMAKE_BINARY_DIR})
file(COPY ${CMAKE_SOURCE_DIR}/examples DESTINATION ${CMAKE_BINARY_DIR})

#===============================================================================
# Catch Testing
#===============================================================================

# Set up directory
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/unit_test)    
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/unit_test)

# Iterate over tests
file(GLOB TESTS "${CMAKE_SOURCE_DIR}/unit_test/test_*.cpp")
foreach(TEST ${TESTS})
    # Built tests
    get_filename_component(TEST_NAME ${TEST} NAME_WE)
    add_executable(${TEST_NAME}.tst unit_test/${TEST_NAME}.cpp)
    target_link_libraries(${TEST_NAME}.tst libMC)
    add_custom_command(
        TARGET ${TEST_NAME}.tst
        POST_BUILD
        COMMAND ${TEST_NAME}.tst
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Running unit test: ${TEST_NAME}" VERBATIM)
endforeach()
