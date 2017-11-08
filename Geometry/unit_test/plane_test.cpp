#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file

#include <limits>

#include "Catch.h"
#include "Geometry.h"
#include "Random.h"


TEST_CASE( "Generic Plane Surface", "" ) {

    const double a =  1.0;
    const double b =  2.0;
    const double c = -3.0;
    const double d =  4.0;

    std::string plane_name = "myPlane";

    Plane_Surface planeT( plane_name, "transmission", a, b, c, d );

    // test returns appropriate name
    SECTION ( " return surface name " ) 
    {
        REQUIRE( planeT.name() == plane_name );
    }

    // test evaluation for a single point outside plane
    SECTION ( " evaluation positive side of plane " ) 
    {
        Point_t p( 3.0, 2.0, -1.0 );
        double eval_result = +6.0;
        REQUIRE( planeT.eval(p) == Approx( eval_result ) );
    }

    // test evaluation for a single point outside plane
    SECTION ( " evaluation negative side of plane " ) 
    {
        Point_t p( -3.0, -2.0, +1.0 );
        double eval_result = -14.0;
        REQUIRE( planeT.eval(p) == Approx( eval_result ) );
    }

    // test evaluation for a single point on plane
    SECTION ( " evaluation on plane " ) 
    {
        Point_t p( 3.0, 2.0, 1.0 );
        double eval_result = +0.0;
        REQUIRE( planeT.eval(p) == Approx( eval_result ) );
    }

    // test 1000 random points in box from -10 to 10
    SECTION ( " 1000 random points near plane " ) 
    {
        bool flag = true;
        for ( unsigned i = 0 ; i < 1000 ; i++ ) 
        {
            Point_t  p( 20.0 * Urand() - 10.0, 20.0 * Urand() - 10.0, 20.0 * Urand() - 10.0 );
            double eval_result = a*p.x + b*p.y + c*p.z - d;
            if ( planeT.eval(p) != eval_result ) 
            { 
                flag = false;
                break;
            };
        }
        REQUIRE( flag );
    }

/* TBD: Need to define Particle library first
    // test smallest positive distance w/ intersection
    SECTION ( " distance with intersection " ) {
      ray r( point( -3.0, -2.0, +1.0 ), point( 1.0, 0.0, 0.0 ) );
      double eval_result = 14.0;
      REQUIRE( thePlane.distance(r) == Approx( eval_result ) );
    } 

    // test smallest positive distance w/ no intersection
    SECTION ( " distance with no intersection " ) {
      ray r( point( -3.0, -2.0, +1.0 ), point( -1.0, 0.0, 0.0 ) );
      double eval_result = std::numeric_limits<double>::max();
      REQUIRE( thePlane.distance(r) == Approx( eval_result ) );
    } 
*/
}
