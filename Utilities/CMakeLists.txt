file(GLOB utils_SRC "*.cpp")
add_library(utils ${utils_SRC})

add_subdirectory(unit_test)
