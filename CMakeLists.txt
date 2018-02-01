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
    if(debug)
      list(REMOVE_ITEM cxxflags -O2)
      list(APPEND cxxflags -g -O0)
    endif()
    if(profile)
      list(APPEND cxxflags -pg)
    endif()
    if(optimize)
      list(REMOVE_ITEM cxxflags -O2)
      list(APPEND cxxflags -O3)
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
# Define the MC program executable
#===============================================================================

file(GLOB libMC_SRC "src/*.cpp")
add_library(libMC ${libMC_SRC})
target_link_libraries(libMC pugixml)

set(program "MC.exe")
add_executable(${program} Main.cpp)
target_link_libraries(${program} pugixml libMC)

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
