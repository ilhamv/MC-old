#include "simulator.h"
#include "Random.h"


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
// Cell Importance
//=============================================================================

void Simulator::cell_importance( Particle& P, std::vector<Particle>& Pbank )
{
    // Importances
    const double Iold = P.cell_old()->importance();
    const double Inew = P.cell()->importance();

    // Same importance, do nothing
    if ( Inew == Iold ) { return; }
	
    const double rat = Inew / Iold; // Ratio of importances
	
    // Less important, Russian Roulette
    if ( rat < 1.0 ){
	if ( Urand() < rat ) { P.set_weight( P.weight() / rat ); }
	else                 { P.kill(); }
    }

    // More important, Splitting
    else{
	// Sample # of splitting
	const int n = std::floor( rat + Urand() );
		
	// Update working particle weight
	P.set_weight( P.weight() / double(n) );

	// Push n-1 identical particles into particle bank
	for ( int i = 0 ; i < n-1 ; i++ ){ Pbank.push_back( P ); }
    }
}

//=============================================================================
// Particle Comb
//=============================================================================

void Simulator::particle_comb( std::vector<Particle>& Pbank )
{
    if( Pbank.size() < bank_max ){ return; }

    // Total weight
    double W = 0.0;
    for( auto P : Pbank ) { W += P.weight(); }

    // Average weight
    const double w_avg = W / comb_teeth;

    // First comb tooth
    double tooth = Urand() * w_avg;

    // Temporary bank: combed particles
    std::vector<Particle> bank_tmp( comb_teeth, Pbank[0] );
    double sum = 0; int j = 0;
    for( int i = 0; i < Pbank.size(); i++ ){
        sum += Pbank[i].weight();
        if( sum > tooth ){
            bank_tmp[j] = Pbank[i]; 
            bank_tmp[j].set_weight( w_avg );

            tooth += w_avg; j++;
        }
    }

    // Copy bank
    Pbank = bank_tmp;
}
