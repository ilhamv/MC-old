file(GLOB distribution_SRC "*.cpp")
add_library(distribution ${distribution_SRC})

target_link_libraries(distribution utils)

add_subdirectory(unit_test)
