#include <iostream>
#include <cstring> 

#include "Simulator.h"


int main( int argc, char* argv[] )
{
    // Input
    if ( argc == 1 ){ 
        std::cout<< "[INPUT ERROR] Please provide input.xml directory...\n";
        std::exit(EXIT_FAILURE);
    }
    std::string input_file = argv[1]; 

    // Set up simulation
    Simulator MC_Simulator( input_file );
    std::cout<<"\nSimulation setup done,\nNow running the simulation...\n\n";

    // Start simulation
    MC_Simulator.start();
    std::cout<<"Simulation done!\n\nCreating output.h5...\n";

    // Report simulation results
    MC_Simulator.report();
    std::cout<<"output.h5 Done!\n";

    return 0;
}
