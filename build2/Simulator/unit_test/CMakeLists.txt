include(Catch)

file(GLOB TEST_INPUTS "${CMAKE_CURRENT_SOURCE_DIR}/*")

file(COPY ${TEST_INPUTS}                 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_SOURCE_DIR}/xs_library DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_catch_test(integral_simulation_test simulator)
