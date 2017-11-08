#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file

#include "Catch.h"
#include "Distribution.h"
#include "Point.h"

TEST_CASE( "Delta distribution", "" )
{
    SECTION( "delta - double" )
    { 
        Delta_Distribution<double> Delta_double(-5.0);
        REQUIRE( Delta_double.sample() == -5.0 );
    }
    
    SECTION( "delta - point" )
    { 
        Point_t p_return(3.0, 2.0, 1.0);
        Delta_Distribution<Point_t> Delta_point( p_return );
        Point_t p_val = Delta_point.sample();
        REQUIRE( p_return == p_val);
    }
}
