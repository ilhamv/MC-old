file(GLOB simulator_SRC "*.cpp")
add_library(simulator ${simulator_SRC})

target_link_libraries(simulator parser)

add_subdirectory(unit_test)
