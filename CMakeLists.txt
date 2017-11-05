cmake_minimum_required (VERSION 2.6 FATAL_ERROR)

project (MC CXX)
add_definitions("-std=c++11 -O3")

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules)

enable_testing()

include_directories(Utilities)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(Utilities)

file(GLOB MC_SRC "*.cpp")
add_executable(MC.exe ${MC_SRC})
target_link_libraries(MC.exe utils)

# Copy Cross-Section library (text files) and input examples
file(COPY ${CMAKE_SOURCE_DIR}/xs_library DESTINATION ${CMAKE_BINARY_DIR})
file(COPY ${CMAKE_SOURCE_DIR}/input_example DESTINATION ${CMAKE_BINARY_DIR})
