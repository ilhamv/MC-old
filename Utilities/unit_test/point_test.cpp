#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file

#include "Catch.h"
#include "Solver.h"
#include "Point.h"

TEST_CASE( "Point", "Test 3-elements point object" )
{
    // Test creation (default)
    Point_t p;
    SECTION( "point default creation" )
    { 
        bool flag = true ;
        if ( p.x != 0.0 ) { flag = false; }
        if ( p.y != 0.0 ) { flag = false; }
        if ( p.z != 0.0 ) { flag = false; }
        REQUIRE( flag );
    }

    // Test creation (default)
    Point_t p2(1.0,-2.0,3.0);
    SECTION( "point creation" )
    { 
        bool flag = true ;
        if ( p2.x != 1.0 ) { flag = false; }
        if ( p2.y != -2.0 ) { flag = false; }
        if ( p2.z != 3.0 ) { flag = false; }
        REQUIRE( flag );
    }

    // Test point normalization
    p2.normalize();
    SECTION( "point normalized" )
    {
        REQUIRE( p2.x*p2.x + p2.y*p2.y + p2.z*p2.z == Approx(1.0) );
    }
}
