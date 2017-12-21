include(MCpackage)

# add_MC_package(<Folder name> <lib name>)
add_MC_package(Utilities utils)
add_MC_package(Distribution distribution)
add_MC_package(Particle particle)
add_MC_package(Geometry geometry)
add_MC_package(Estimator estimator)
add_MC_package(Parser parser)
add_MC_package(Material material)
add_MC_package(VReduction vreduction)
add_MC_package(Simulator simulator)

build_MC_packages()
