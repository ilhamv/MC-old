#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file

#include <algorithm> // abs

#include "Catch.h"
#include "Simulator.h"

TEST_CASE( "Integral Simulation Tests", "" ) 
{
    SECTION ( " analytic 1D slab test " ) 
    {
        Simulator_t Sim("../examples/slab_analytic");
        Sim.start();
        Sim.report();

        double mean  = Sim.Estimator[0]->total_tally[0].mean;
        double uncer = Sim.Estimator[0]->total_tally[0].uncer;
        
        REQUIRE( std::abs( mean - 0.0149956 ) <= uncer  );
    }
    
    SECTION ( " MCNP6: detecting a sphere " ) 
    {
        Simulator_t Sim("../examples/sphere_detection");
        Sim.start();
        Sim.report();

        double mean  = Sim.Estimator[0]->total_tally[1].mean;
        double uncer = Sim.Estimator[0]->total_tally[1].uncer;
        
        REQUIRE( std::abs( mean - 6.9276e-5 ) <= uncer+1.15e-6  );
    }
}