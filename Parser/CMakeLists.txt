file(GLOB parser_SRC "*.cpp")
add_library(parser ${parser_SRC})

target_link_libraries(parser material geometry particle vreduction estimator particle)

add_subdirectory(unit_test)
