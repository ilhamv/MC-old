#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file

#include <memory>

#include "Catch.h"
#include "Particle.h"
#include "Geometry.h"
#include "Point.h"

TEST_CASE( "Particle", "" )
{
    Particle_t P(Point_t(1.0,2.0,3.0),Point_t(1.0,0.0,0.0),1E6,0.0,0.5);
    SECTION( "creation and getters" )
    {
        REQUIRE( P.pos() == Point_t(1.0,2.0,3.0) );
        REQUIRE( P.dir() == Point_t(1.0,0.0,0.0) );
        REQUIRE( P.alive() == true );
        REQUIRE( P.weight() == 0.5 );
        REQUIRE( P.time() == 0.0 );
        REQUIRE( P.energy() == 1E6 );
        REQUIRE( P.speed() == Approx(1.383160092041162E09) );
    }

    SECTION( "setters" )
    {
        P.setDirection( Point_t(0.0,1.0,0.0) );
        REQUIRE( P.dir() == Point_t(0.0,1.0,0.0) );
        
        P.setWeight( 0.7 );
        REQUIRE( P.weight() == 0.7 );
        
        P.setTime( 10.0 );
        REQUIRE( P.time() == 10.0 );
        
        P.setEnergy( 1.0 );
        REQUIRE( P.energy() == 1.0 );
        REQUIRE( P.speed() == Approx(1.383160092041162E06) );
        
        P.kill();
        REQUIRE( P.alive() == false );
        
        P.setSpeed( 1E7 );
        REQUIRE( P.energy() == Approx(52.270312948607518) );
        REQUIRE( P.speed() == 1E7 );
    }
}
