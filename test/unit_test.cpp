#include "gtest/gtest.h"
#include "Point.h"


//=============================================================================
// Point
//=============================================================================

TEST( Point, Creation_Default)
{
    Point_t p;
    bool flag = true ;
    if ( p.x != 0.0 ) { flag = false; }
    if ( p.y != 0.0 ) { flag = false; }
    if ( p.z != 0.0 ) { flag = false; }
    ASSERT_TRUE( flag );
}

TEST( Point, Creation_Init)
{ 
    Point_t p(1.0,-2.0,3.0);
    bool flag = true ;
    if ( p.x != 1.0 ) { flag = false; }
    if ( p.y != -2.0 ) { flag = false; }
    if ( p.z != 3.0 ) { flag = false; }
    ASSERT_TRUE( flag );
}

//=============================================================================
// Point
//=============================================================================

//TEST( Particle

//=============================================================================
// Main
//=============================================================================

int main( int argc, char** argv )
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
