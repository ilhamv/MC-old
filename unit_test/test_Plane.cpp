#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file

#include <limits>

#include "Catch.h"
#include "Geometry.h"
#include "Random.h"
#include "Const.h"


TEST_CASE( "Generic Plane Surface", "" ) 
{
    const double a =  1.0;
    const double b =  2.0;
    const double c = -3.0;
    const double d =  4.0;


    Plane_Surface planeT( "my_Plane", 1, "transmission", a, b, c, d );
    Plane_Surface planeR( "my_RPlane", 1, "reflective", a, b, c, d );

    // test returns appropriate name
    SECTION ( " return surface name and bc " ) 
    {
        REQUIRE( planeT.name() == "my_Plane" );
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

    // test smallest positive distance w/ intersection
    SECTION ( " distance with intersection " ) 
    {
        Particle_t P( Point_t( -3.0, -2.0, +1.0 ), Point_t( 1.0, 0.0, 0.0 ), 1.0, 1.0, 1.0 );
        double eval_result = 14.0;
        REQUIRE( planeT.distance(P) == Approx( eval_result ) );
    } 

    // test smallest positive distance w/ no intersection
    SECTION ( " distance with no intersection " ) 
    {
        Particle_t P( Point_t( -3.0, -2.0, +1.0 ), Point_t( -1.0, 0.0, 0.0 ), 1.0, 1.0, 1.0 );
        double eval_result = MAX;
        REQUIRE( planeT.distance(P) == Approx( eval_result ) );
    } 
    
    // test particle moving parallel to surface
    SECTION ( " parallel to surface " ) 
    {
        Point_t dir( 1.0, 1.0, 1.0 ); dir.normalize();
        Particle_t P( Point_t( -3.0, -2.0, +1.0 ), dir, 1.0, 1.0, 1.0 );
        double eval_result = MAX;
        REQUIRE( planeT.distance(P) == Approx( eval_result ) );
    } 
}

/*
// PlaneX test suite
TEST_CASE( "X Plane Surface", "" ) 
// PlaneY test suite
TEST_CASE( "Y Plane Surface", "" ) 
// PlaneZ test suite
TEST_CASE( "Z Plane Surface", "" ) 
*/
