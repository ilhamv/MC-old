# add_catch_test(<target> <libs>...)
# target's source need to have name of <target>.cpp

function(add_catch_test target)
    add_executable(${target} ${target}.cpp)
    target_link_libraries(${target} ${ARGN})
    add_test(${target} ${target})

    add_custom_command(
        TARGET ${target}
        POST_BUILD
        COMMAND ${target}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Running ${target}" VERBATIM)
endfunction()
