file(GLOB estimator_SRC "*.cpp")
add_library(estimator ${estimator_SRC})

target_link_libraries(estimator geometry particle distribution)

add_subdirectory(unit_test)
