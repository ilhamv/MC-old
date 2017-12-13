#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file

#include <vector>

#include "Catch.h"
#include "Solver.h"
#include "Point.h"

TEST_CASE( "Linear Interpolation", "" )
{
    const double x1 = 1.0;
    const double x2 = 2.0;
    const double x  = 1.5;
    const double y1 = 2.0;
    const double y2 = 4.0;

    SECTION( "positive gradient line" )
    { 
        REQUIRE( Linterpolate(x, x1, x2, y1, y2) == 3.0 ); 
    }

    SECTION( "negative gradient line" )
    { 
        REQUIRE( Linterpolate(x, x1, x2, y2, y1) == 3.0 ); 
    }
    
    SECTION( "zero gradient line" )
    { 
        REQUIRE( Linterpolate(x, x1, x2, 3.0, 3.0) == 3.0 ); 
    }
}

TEST_CASE( "Geometry Quadratic Solver", "" )
{
    double a,b,c;

    SECTION( "no intersection" )
    {
        b = 1;
        a = 10;
        c = 10;
        REQUIRE( solve_quad(a,b,c) == Approx(MAX) ); 
    }
    
    SECTION( "tangent to surface" )
    {
        b = 2;
        a = 1;
        c = 1;
        REQUIRE( solve_quad(a,b,c) == Approx(MAX) ); 
    }
    
    SECTION( "moving away" )
    {
        a = 1;
        b = 3;
        c = 2;
        REQUIRE( solve_quad(a,b,c) == Approx(MAX) ); 
    }
    
    SECTION( "moving in" )
    {
        a = 1;
        b = -3;
        c = 2;
        REQUIRE( solve_quad(a,b,c) == Approx(1.0) ); 
    }

    SECTION( "inside the surface" )
    {
        a = 1;
        b = -1;
        c = -2;
        REQUIRE( solve_quad(a,b,c) == Approx(2.0) );
    }
}

TEST_CASE( "Binary Search", "" )
{
    std::vector<double> vec = { 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };

    SECTION( "between points" )
    {
        REQUIRE( Binary_Search( 4.3, vec ) == 2 );
        REQUIRE( Binary_Search( 7.3, vec ) == 5 );
    }
    
    SECTION( "at point" )
    {
        REQUIRE( Binary_Search( 4.0, vec ) == 1 );
        REQUIRE( Binary_Search( 7.0, vec ) == 4 );
        REQUIRE( Binary_Search( 2.0, vec ) == -1 );
        REQUIRE( Binary_Search( 10.0, vec ) == 7 );
    }

    SECTION( "outside" )
    {
        REQUIRE( Binary_Search( 1.0, vec ) == -1 );
        REQUIRE( Binary_Search( 12.0, vec ) == 8 );
    }
}

TEST_CASE( "Direction Scattering","" )
{
    /* 
       In this test, the first four random number is considered in evaluating the expected results:
       (0.380815, 0.592166, 0.159181, 0.559543)
    */

    SECTION( "Forward scattering" )
    {
        const double mu = 0.4;
        Point_t dir( 2.0, 2.0, 1.0 );
        dir.normalize();

        dir = scatter_direction(dir,mu);

        REQUIRE( dir.x == Approx(-0.332775973037166) );
        REQUIRE( dir.y == Approx(0.549648481770809) );
        REQUIRE( dir.z == Approx(0.766254982532715) );
    }
    
    SECTION( "Backward scattering" )
    {
        const double mu = -0.4;
        Point_t dir( 2.0, 2.0, 1.0 );
        dir.normalize();

        dir = scatter_direction(dir,mu);
        REQUIRE( dir.x == Approx(-0.092800252430585) );
        REQUIRE( dir.y == Approx(-0.802140132078611) );
        REQUIRE( dir.z == Approx(0.589880769018391) );
    }
    
    SECTION( "Forward scattering, with initial particle direction (0,0,1)" )
    {
        const double mu = 0.4;
        Point_t dir( 0.0, 0.0, 1.0 );

        dir = scatter_direction(dir,mu);

        REQUIRE( dir.x == Approx(-0.771301959816257) );
        REQUIRE( dir.y == Approx(-0.495068971743938) );
        REQUIRE( dir.z == Approx(0.4) );
    }
    
    SECTION( "Backward scattering, with initial particle direction (0,0,1)" )
    {
        const double mu = -0.4;
        Point_t dir( 0.0, 0.0, 1.0 );

        dir = scatter_direction(dir,mu);

        REQUIRE( dir.x == Approx(0.334943449022552) );
        REQUIRE( dir.y == Approx(0.853119502740898) );
        REQUIRE( dir.z == Approx(-0.4) );
    }
}
