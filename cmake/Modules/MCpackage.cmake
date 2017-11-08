# Add linked libraries

function(add_MC_package dir lib)
    include_directories(${dir})
    add_subdirectory(${dir})
    list(APPEND MC_PCKGS ${lib})
    set(MC_PCKGS ${MC_PCKGS} PARENT_SCOPE)
endfunction()
