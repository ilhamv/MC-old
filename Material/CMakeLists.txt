file(GLOB material_SRC "*.cpp")
add_library(material ${material_SRC})

target_link_libraries(material particle utils distribution)

add_subdirectory(unit_test)
