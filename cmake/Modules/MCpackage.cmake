# Add linked libraries

function(add_MC_package dir lib)
    include_directories(${dir})
    list(APPEND MC_PCKGS_DIR ${dir})
    set(MC_PCKGS_DIR ${MC_PCKGS_DIR} PARENT_SCOPE)
    list(APPEND MC_PCKGS ${lib})
    set(MC_PCKGS ${MC_PCKGS} PARENT_SCOPE)
endfunction()

function(build_MC_packages)
    foreach(dir ${MC_PCKGS_DIR})
        add_subdirectory(${dir})
    endforeach()
endfunction()

