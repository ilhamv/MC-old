// #define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file

#include <limits>

#include "Catch.h"
#include "Geometry.h"
#include "Random.h"


TEST_CASE( "Generic Sphere ", "" ) {

    const double x0 =  1.0;
    const double y0 =  2.0;
    const double z0 = -3.0;
    const double rad_sq =  4.0;

    std::string sphere_name = "mySphere";

    Sphere_Surface sphere( sphere_name, "transmission", x0, y0, z0, rad_sq );

    // test returns appropriate name
    SECTION ( " return surface name " )
    {
        REQUIRE( sphere.name() == sphere_name );
    }

    // test evaluation for a single point outside sphere
    SECTION ( " evaluation positive side of sphere " )
    {
        Point_t p( 2.0, 0.0, -1.0 );
        double eval_result = +5.0;
        REQUIRE( sphere.eval(p) == Approx( eval_result ) );
    }

    // test evaluation for a single point outside sphere
    SECTION ( " evaluation negative side of sphere " )
    {
        Point_t p( 0.0, 3.0, -2.0 );
        double eval_result = -1.0;
        REQUIRE( sphere.eval(p) == Approx( eval_result ) );
    }

    // test evaluation for a single point on sphere
    SECTION ( " evaluation on sphere " )
    {
        Point_t p( 1.0, 2.0, -1.0 );
        double eval_result = +0.0;
        REQUIRE( sphere.eval(p) == Approx( eval_result ) );
    }

    // test 1000 random points in box from -10 to 10
    SECTION ( " 1000 random points near sphere " )
    {
        bool flag = true;
        for ( unsigned i = 0 ; i < 1000 ; i++ )
        {
            Point_t  p( 20.0 * Urand() - 10.0, 20.0 * Urand() - 10.0, 20.0 * Urand() - 10.0 );
	    double x_t = p.x - x0;
  	    double y_t = p.y - y0;
  	    double z_t = p.z - z0;
	    double eval_result = x_t*x_t + y_t*y_t + z_t*z_t - rad_sq;
            if ( sphere.eval(p) != eval_result )
            {
                flag = false;
                break;
            };
        }
        REQUIRE( flag );
    }
}
