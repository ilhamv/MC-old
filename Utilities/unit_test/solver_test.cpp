#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file

#include "Catch.h"
#include "Solver.h"

TEST_CASE( "Solvers", "Test solvers" )
{
    // Test linear interpolation
    const double x1 = 1.0;
    const double x2 = 2.0;
    const double x  = 1.5;
    const double y1 = 2.0;
    const double y2 = 4.0;
    SECTION ( "linear interpolation" )
    { 
        REQUIRE( Linterpolate(x, x1, x2, y1, y2) == 3.0 ); 
    }

    // Test quadratic solver

    // Test binary search

    // Test direction scattering
}
