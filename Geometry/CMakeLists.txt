file(GLOB geometry_SRC "*.cpp")
add_library(geometry ${geometry_SRC})

target_link_libraries(geometry utils particle estimator material)

add_subdirectory(unit_test)
