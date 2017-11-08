file(GLOB particle_SRC "*.cpp")
add_library(particle ${particle_SRC})

target_link_libraries(particle utils geometry distribution)

add_subdirectory(unit_test)
