#include <iostream>// cout
#include <cstring> // string, strcmp

#include "Simulator.h"


int main( int argc, char* argv[] )
{
    // Check # of arguments passed
    if ( argc == 1 ) { std::cout<< "Please provide input file\n"; return 1; }

    // Variable declarations
    std::string input_file = argv[1]; // Input file name 

    // Set up simulation
    Simulator_t MC_Simulator( input_file );
    std::cout<<"\nSimulation setup done,\nNow running the simulation...\n\n";
    std::cout.flush();

    // Start simulation
    MC_Simulator.start();
    std::cout<<"Simulation done!\n\nCreating output.h5...\n";

    // Report simulation results
    MC_Simulator.report();
    std::cout<<"output.h5 created.\n";

    return 0;
}
