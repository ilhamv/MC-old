file(GLOB vreduction_SRC "*.cpp")
add_library(vreduction ${vreduction_SRC})

target_link_libraries(vreduction utils geometry particle)

add_subdirectory(unit_test)
