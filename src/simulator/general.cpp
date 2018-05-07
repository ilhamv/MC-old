#include <iostream>
#include <memory>  

#include "simulator.h"
#include "Algorithm.h"
#include "Random.h"


//=============================================================================
// Test Point
//=============================================================================

bool Simulator::test_point( const Point& p, const std::shared_ptr<Cell>& C )
{
    // Loop over surfaces of cell, if particle not on correct side return false
    for ( const auto& S : C->surfaces() ) {
    	if ( S.first->eval( p ) * S.second < 0 ) { return false; }  
    }
    return true;
}

//=============================================================================
// Search Cell
//=============================================================================

std::shared_ptr<Cell> Simulator::search_cell( const Point& p )
{
    for( const auto& C : Cells ){
        if ( test_point( p, C ) ){ return C; }
    }
    std::cout<< "[WARNING] A particle is lost:\n( x, y, z )  (" << p.x << ", " 
             << p.y << ", " << p.z << " )\n";
    std::exit(EXIT_FAILURE);
}

//=============================================================================
// Collision Distance
//=============================================================================

double Simulator::collision_distance( const Particle& P )
{
    if ( P.cell()->material() ){ 
        return exponential_sample( P.cell()->material()->SigmaT( P.energy() ) );
    }
    // Vacuum --> return sligthly less than very large number
    //            to ensure collision (kill) if no surface intersection
    else { return MAX_float_less; }
}

//=============================================================================
// Surface Intersect
//=============================================================================

std::pair<std::shared_ptr<Surface>, double> 
Simulator::surface_intersect( const Particle& P ) 
{
    double dist = MAX_float;
    std::shared_ptr< Surface > S = nullptr;
    for ( const auto& s : P.cell()->surfaces() ) {
    	double d = s.first->distance( P );
    	if ( d < dist ){ 
	    dist = d;
	    S    = s.first;
	}
    }
    return std::make_pair( S, dist ); 
}

//=============================================================================
// Move particle
//=============================================================================

void Simulator::move_particle( Particle& P, const double l )
{
    P.move( l );
    Ntrack++;

    // Tallies
    if(ksearch) { k_estimator->estimate_TL(P,l); }
    if(tally){ 
        for( auto& e : P.cell()->estimators_TL ) { e->score( P, l ); }
    }
}

//=============================================================================
// Surface Hit
//=============================================================================
        
void Simulator::surface_hit( Particle& P, const std::shared_ptr<Surface>& S )
{
    P.set_surface_old(S);

    // Transmission
    if ( S->bc() == 0 ){
	P.move( EPSILON_float );
	P.set_cell( search_cell( P.pos() ) );
    // Vacuum
    } else if ( S->bc() == -1 ){
	P.kill();
        P.set_cell( P.cell() );
    // Reflective
    } else{
	S->reflect( P );
	P.move( EPSILON_float );
        P.set_cell( P.cell() );
    }

    // Tallies
    if (tally){
	for ( auto& e : S->estimators ){ e->score( P, 0.0 ); }
    }
}

//=============================================================================
// Collision
//=============================================================================

void Simulator::collision( Particle& P )
{
    // Vacuum --> Kill particle at collision
    if( !P.cell()->material() ){ return P.kill(); }
    
    // Collision tallies
    if(tally){ 
        for( auto& e : P.cell()->estimators_C ) { e->score( P, 0 ); }
    }
   
    // Get material
    std::shared_ptr<Material> M = P.cell()->material();
    
    // Implicit Fission
    const double bank_nu = std::floor( P.weight() / k * M->nuSigmaF(P.energy())
                                       / M->SigmaT(P.energy()) + Urand() );
    
    std::shared_ptr<Nuclide> N_fission = M->nuclide_nufission( P.energy() );

    // Collision implementaiton
    if(ksearch){ collision_ksearch( P, bank_nu, N_fission ); }
    else{ collision_fixed_source( P, bank_nu, N_fission ); }

    // Implicit Absorption
    const double implicit = M->SigmaC(P.energy()) + M->SigmaF(P.energy());
    P.set_weight( P.weight() * ( M->SigmaT(P.energy()) - implicit ) 
                  / M->SigmaT(P.energy()) );

    std::shared_ptr<Nuclide> N_scatter = M->nuclide_scatter( P.energy() );
    if(!N_scatter){ return; }
    
    // The only analog reaction
    N_scatter->scatter()->sample( P );
}

//=============================================================================
// Weight Roulette
//=============================================================================

void Simulator::weight_roulette( Particle& P )
{
    if( P.weight() < wr ){
        if( Urand() < P.weight() / ws ) { P.set_weight(ws); }
        else { P.kill(); }
    }
}

//=============================================================================
// Push Particle Bank
//=============================================================================
void Simulator::push_particle_bank( const Particle& P )
{
    Pbank.push(P);
}
