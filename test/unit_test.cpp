#include <memory>

#include "gtest/gtest.h"
#include "Point.h"
#include "Particle.h"
#include "Geometry.h"
#include "Constants.h"
#include "Algorithm.h"


//=============================================================================
// Point
//=============================================================================

TEST( Point, Point)
{
    Point p;
    bool flag = true ;
    if ( p.x != 0.0 ) { flag = false; }
    if ( p.y != 0.0 ) { flag = false; }
    if ( p.z != 0.0 ) { flag = false; }
    ASSERT_TRUE( flag );
    
    Point p2(1.0,-2.0,3.0);
    flag = true ;
    if ( p2.x != 1.0 ) { flag = false; }
    if ( p2.y != -2.0 ) { flag = false; }
    if ( p2.z != 3.0 ) { flag = false; }
    ASSERT_TRUE( flag );
}

//=============================================================================
// Particle
//=============================================================================

TEST( Particle, Particle )
{
    Particle P( Point(1.0,2.0,3.0), Point(1.0,0.0,0.0), 1E6, 0.0, 1.0, 0, 
                NULL );
    ASSERT_TRUE( point_equal( P.pos(), Point(1.0,2.0,3.0) ) );
    ASSERT_TRUE( point_equal( P.dir(), Point(1.0,0.0,0.0) ) );
    ASSERT_EQ( P.energy(), 1E6 );
    ASSERT_EQ( P.time(), 0.0 );
    ASSERT_EQ( P.alive(), true );
    ASSERT_EQ( P.weight(), 1.0 );    
    ASSERT_EQ( P.speed(), std::sqrt( 1E6 * 191312955.067 ) * 100.0 );
    ASSERT_EQ( P.tdmc(), 0 );

    Point p(0.0,0.0,1.0);
    P.set_direction(p);
    P.set_weight(0.5);
    P.set_energy(100);
    ASSERT_TRUE( point_equal( P.dir(), Point(0.0,0.0,1.0) ) );
    ASSERT_EQ( P.energy(), 1E2 );
    ASSERT_EQ( P.energy_old(), 1E6 );
    ASSERT_EQ( P.speed(), std::sqrt( 1E2 * 191312955.067 ) * 100.0 );
    ASSERT_EQ( P.weight(), 0.5 );

    P.set_speed(2200.0);
    ASSERT_EQ( P.energy(), 2200.0*2200.0 * 5.2270376e-13 );
    ASSERT_EQ( P.energy_old(), 1E2 );
    ASSERT_EQ( P.speed(), 2200.0 );

    P.move( 2.0 );
    P.increment_tdmc();P.increment_tdmc();
    ASSERT_TRUE( point_equal( P.pos(), Point(1.0,2.0,5.0) ) );
    ASSERT_EQ( P.time(), 0.0 + 2.0 / 2200.0 );
    ASSERT_EQ( P.time_old(), 0.0 );
    ASSERT_EQ( P.tdmc(), 2 );

    P.kill();
    ASSERT_EQ( P.alive(), false );
    ASSERT_EQ( P.weight(), 0.0 ); 

    std::shared_ptr<Cell> C  = std::make_shared<Cell>("cell", 1, 1.0);
    std::shared_ptr<Cell> C2 = std::make_shared<Cell>("cell2", 1, 1.0);
    Particle P2( Point(1.0,2.0,3.0), Point(1.0,0.0,0.0), 1E6, 0.0, 1.0, 0, 
                 C );
    P2.set_cell(C2);
    ASSERT_EQ( P2.cell()->name(), "cell2" );
    ASSERT_EQ( P2.cell_old()->name(), "cell" );

    std::shared_ptr<Surface> S = 
        std::make_shared<SurfacePlaneX>("surface",1,0,0.0);
    P2.set_surface_old(S);
    ASSERT_EQ( P2.surface_old()->name(), "surface" );
}

//=============================================================================
// Surface
//=============================================================================

TEST( Surface, General )
{
    std::shared_ptr<Surface> S = 
        std::make_shared<SurfacePlaneX>("surface1", 0, 0, 0.0);
    ASSERT_EQ( S->name(), "surface1" );
    ASSERT_EQ( S->ID(), 0 );
    // Estimator
}

TEST( Surface, PlaneX )
{
    std::shared_ptr<Surface> S = 
        std::make_shared<SurfacePlaneX>("surface1", 0, 0, 0.0);
    ASSERT_TRUE( S->eval( Point(1.0,0.0,0.0) ) > 0.0 );
    ASSERT_TRUE( S->eval( Point(-1.0,0.0,0.0) ) < 0.0 );
    Particle P( Point(1.0,2.0,3.0), Point(-1.0,0.0,0.0), 1E6, 0.0, 1.0, 0, 
                NULL );
    ASSERT_EQ( S->distance(P), 1.0 );
    P.set_direction( Point(1.0,0.0,0.0) );
    ASSERT_EQ( S->distance(P), MAX_float );
    P.move(-2.0);
    ASSERT_EQ( S->distance(P), 1.0 );
    P.set_direction( Point(-1.0,0.0,0.0) );
    ASSERT_EQ( S->distance(P), MAX_float );
    P.set_direction( Point(0.0,1.0,0.0) );
    ASSERT_EQ( S->distance(P), MAX_float );
    
    P.set_direction( Point(0.3,0.1,0.4) );
    S->reflect(P);
    ASSERT_TRUE( point_equal( P.dir(), Point(-0.3,0.1,0.4) ) );
}


//=============================================================================
// Cell
//=============================================================================



//=============================================================================
// Main
//=============================================================================

int main( int argc, char** argv )
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
