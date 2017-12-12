file(GLOB parser_SRC "*.cpp")
add_library(parser ${parser_SRC})

add_definitions("-Lyaml-cpp -Iyaml-cpp-include -lyaml-cpp")

target_link_libraries(parser material geometry particle vreduction estimator particle)

add_subdirectory(unit_test)
