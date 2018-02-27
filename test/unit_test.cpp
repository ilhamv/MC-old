#include "gtest/gtest.h"
#include "Point.h"
#include "Particle.h"
#include "Geometry.h"


//=============================================================================
// Point
//=============================================================================

TEST( Point, Creation_Default)
{
    Point p;
    bool flag = true ;
    if ( p.x != 0.0 ) { flag = false; }
    if ( p.y != 0.0 ) { flag = false; }
    if ( p.z != 0.0 ) { flag = false; }
    ASSERT_TRUE( flag );
}

TEST( Point, Creation_Init)
{ 
    Point p(1.0,-2.0,3.0);
    bool flag = true ;
    if ( p.x != 1.0 ) { flag = false; }
    if ( p.y != -2.0 ) { flag = false; }
    if ( p.z != 3.0 ) { flag = false; }
    ASSERT_TRUE( flag );
}

//=============================================================================
// Point
//=============================================================================

TEST( Particle, Creation_Getters )
{
    Particle P( Point(1.0,2.0,3.0), Point(1.0,0.0,0.0), 1E6, 0.0, 1.0, 0, 
                NULL );
    ASSERT_TRUE( P.pos() == Point(1.0,2.0,3.0) );
    ASSERT_TRUE( P.dir() == Point(1.0,0.0,0.0) );
    ASSERT_EQ( P.energy(), 1E6 );
    ASSERT_EQ( P.time(), 0.0 );
    ASSERT_EQ( P.alive(), true );
    ASSERT_EQ( P.weight(), 1.0 );    
    ASSERT_EQ( P.speed(), std::sqrt( 1E6 * 191312955.067 ) * 100.0 );
    ASSERT_EQ( P.tdmc(), 0 );
}

//=============================================================================
// Main
//=============================================================================

int main( int argc, char** argv )
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
