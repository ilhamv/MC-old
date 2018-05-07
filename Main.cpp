#include <iostream>
#include <cstring> 

#include "H5Cpp.h"

#include "simulator.h"

int main( int argc, char* argv[] )
{
    // I/O Directory
    if ( argc == 1 ){ 
        std::cout<< "[ERROR] Please provide input.xml directory...\n";
        std::exit(EXIT_FAILURE);
    }
    const std::string io_dir = std::string(argv[1]) + "/"; 
    
    // HDF5 output
    H5std_string output_name( io_dir + "output.h5" );
    H5::H5File output(output_name, H5F_ACC_TRUNC);
    
    //=========================================================================
    // Monte Carlo Simulation
    //=========================================================================

    Simulator MC_Simulator( io_dir );
    std::cout<<"\nSimulation setup done,\nNow running the simulation...\n\n";
    
    MC_Simulator.start();
    std::cout<<"Simulation done!\n\nReporting simulation output...\n";
    
    MC_Simulator.report( output );
    std::cout<<"Simulation output done!\n";

    return 0;
}
