#include <cmath>
#include <iostream>
#include <sstream>  // ostringstream

#include "Estimator.h"
#include "Geometry.h"
#include "Solver.h"  // Binary_Search, Linterpolate


//=============================================================================
// Scoring Kernel
//   - Supports: Neutron, Track Length, Collision
//=============================================================================

// Neutron
double ScoreKernelNeutron::score( const Particle_t& P, const double l )
{ 
    return P.weight(); 
}

// Track Length
double ScoreKernelTrackLength::score( const Particle_t& P, const double l )
{ 
    return P.weight() * l; 
}

// Collision
double ScoreKernelCollision::score( const Particle_t& P, const double l )
{ 
    return P.weight() / P.cell()->SigmaT( P.energy() ); 
}


//=============================================================================
// Score
//   - Specified Scoring Kernel
//   - Supports: Flux, Absorption, Scatter, Capture, Fission,
//               NuFission, Total
//=============================================================================

// Flux
double ScoreFlux::score( const Particle_t& P, const double l )
{ 
    return s_kernel->score(P,l); 
}

// Absorption
double ScoreAbsorption::score( const Particle_t& P, const double l )
{ 
    return P.cell()->SigmaA( P.energy() ) * s_kernel->score(P,l); 
}

// Scatter
double ScoreScatter::score( const Particle_t& P, const double l )
{ 
    return P.cell()->SigmaS( P.energy() ) * s_kernel->score(P,l); 
}

// Capture
double ScoreCapture::score( const Particle_t& P, const double l )
{ 
    return P.cell()->SigmaC( P.energy() ) * s_kernel->score(P,l); 
}

// Fission
double ScoreFission::score( const Particle_t& P, const double l )
{ 
    return P.cell()->SigmaF( P.energy() ) * s_kernel->score(P,l); 
}

// NuFission
double ScoreNuFission::score( const Particle_t& P, const double l )
{ 
    return P.cell()->nuSigmaF( P.energy() ) * s_kernel->score(P,l); 
}

// Total
double ScoreTotal::score( const Particle_t& P, const double l )
{ 
    return P.cell()->SigmaT( P.energy() ) * s_kernel->score(P,l); 
}


//==============================================================================
// Filter
//   - Supports: recently crossed Surface, Cell, Energy, Time
//==============================================================================

// Surface
std::vector<std::pair<int,double>> FilterSurface::idx_l( Particle_t& P,
                                                         const double l )
{
    std::vector<std::pair<int,double>> v_i_l;
    std::pair<int,double> i_l;

    // Index location
    i_l.first  = Binary_Search( P.surface_old()->ID(), f_grid ) + 1;
    // Track length to score
    i_l.second = l;

    v_i_l.push_back(i_l);
    return v_i_l;
}
// Cell
std::vector<std::pair<int,double>> FilterCell::idx_l( Particle_t& P,
                                                      const double l )
{
    std::vector<std::pair<int,double>> v_i_l;
    std::pair<int,double> i_l;

    // Index location
    i_l.first  = Binary_Search( P.cell()->ID(), f_grid ) + 1;
    // Track length to score
    i_l.second = l;

    v_i_l.push_back(i_l);
    return v_i_l;
}
// Energy
std::vector<std::pair<int,double>> FilterEnergy::idx_l( Particle_t& P,
                                                        const double l )
{
    std::vector<std::pair<int,double>> v_i_l;
    std::pair<int,double> i_l;
    
    // Index location
    i_l.first  = Binary_Search( P.energy(), f_grid );
    // Check if inside the grids
    if ( i_l.first < 0 && i_l.first >= f_Nbin ) { return v_i_l; }
    // Track length to score
    i_l.second = l;

    v_i_l.push_back(i_l);
    return v_i_l;
}
// Time bin
std::vector<std::pair<int,double>> FilterTime::idx_l( Particle_t& P,
                                                      const double l )
{
    std::vector<std::pair<int,double>> v_i_l;
    std::pair<int,double> i_l;
    
    // Edge bin locations
    int loc1 = Binary_Search( P.time_old(), f_grid );     // before track generation
    int loc2 = Binary_Search( P.time(), f_grid ); // after
    
    // Distribute score into spanned bins
    if ( loc1 == loc2 ) // 1 bin spanned
    {
        // 1 bin spanned
	if ( loc1 >= 0 && loc1 < f_Nbin ){ 
            i_l.first  = loc1;
            i_l.second = l;
            v_i_l.push_back(i_l);
        }else{ 
            // Outside range?
            return v_i_l; 
        }
    }else{
        // >1 bins spanned
        int    num_bin = loc2 - loc1 - 1; // # of full bins spanned
	double new_track;                 // to hold bin track
	// First partial bin
	if ( loc1 >= 0 ){
	    new_track = ( f_grid[loc1+1] - P.time_old() ) * P.speed();
	    v_i_l.push_back( std::make_pair( loc1, new_track ) );
	}
	// Intermediate full bin
	for ( int i = 1 ; i <= num_bin ; i++ ){
	    new_track = ( f_grid[loc1+i+1] - f_grid[loc1+i] ) * P.speed();
	    v_i_l.push_back( std::make_pair( loc1+i, new_track ) );
	}	
	// Last partial bin
	if ( loc2 < f_Nbin ){
	    new_track = ( P.time() - f_grid[loc2] ) * P.speed();
	    v_i_l.push_back( std::make_pair( loc2, new_track ) );
	}
	// Note that this algorithm covers the following "extreme" cases
	//   loc1 : <lowest_grid  or at grid point
	//   loc2 : >highest_grid or at grid point
    }

    return v_i_l;
}
// Time TDMC
std::vector<std::pair<int,double>> FilterTDMC::idx_l( Particle_t& P,
                                                      const double l )
{
    std::vector<std::pair<int,double>> v_i_l;
    std::pair<int,double> i_l;
    
    if( P.time() != f_grid[P.tdmc()] ) { return v_i_l; }
    // Index location
    i_l.first  = P.tdmc(); 
    P.set_tdmc(i_l.first+1);
    // Velocity to score
    i_l.second = P.speed();

    v_i_l.push_back(i_l);
    return v_i_l;
}


//==============================================================================
/// Estimator
//==============================================================================

// Initialization
void Estimator::add_score( const std::shared_ptr<Score>& S )
{ 
    e_scores.push_back( S ); 
}
void Estimator::add_filter( const std::shared_ptr<Filter>& F )
{
    e_filters.push_back( F );
    e_idx_l.push_back(std::vector<std::pair<int,double>>());
}
void Estimator::initialize_tallies()
{
    int Ntally = e_scores.size();
    for( auto& f : e_filters ){
        Ntally *= f->size();
    }
    e_tally.resize(Ntally,Tally());

    // Tally index factor
    double factor = 1.0;
    idx_factor.insert( idx_factor.begin(), factor );
    for( int i = e_filters.size()-1; i >= 0; i-- ){
        factor *= e_filters[i]->size();
        idx_factor.insert( idx_factor.begin(), factor );
    }
}

// Score
void Estimator::score( Particle_t& P, const double l )
{
    // Individual filter indexes and l
    for( int i = 0; i < e_filters.size(); i++ ){
        e_idx_l[i] = e_filters[i]->idx_l(P,l);
        if( e_idx_l[i].size() == 0 ) { return; }
    }

    // Iterate over all possible filter grid combination
    std::vector<int> idx( e_filters.size(), 0 ); // Current combination
    while( true ){
        // Find minimum length
        double l = MAX;
        for( int i = 0; i < e_filters.size(); i++ ){
            l = std::min( l, e_idx_l[i][idx[i]].second );
        }
        // Score the minimum length to the current combination
        //   (Need to convert multi-D index to 1-D)
        int idx_1D = 0;
        for( int i = 0; i < e_filters.size(); i++ ){
            idx_1D += e_idx_l[i][idx[i]].first * idx_factor[i+1];
        }
        for( int s = 0; s < e_scores.size(); s++ ){
            e_tally[idx_1D].hist += e_scores[s]->score(P,l);
            idx_1D += idx_factor[0];
        }

        // Go to next combination
        for( int i = 0; i < e_filters.size(); i++ ){
            e_idx_l[i][idx[i]].second -= l;
            if ( e_idx_l[i][idx[i]].second < EPSILON ){ 
                if ( idx[i] == e_idx_l[i].size() - 1 ) { return; }
                idx[i]++; 
            }
        }

        if( e_filters.size() == 0 ) { return; }
    }
}

// Loop closeouts
void Estimator::end_history()
{
    for( auto& tally : e_tally ){
        tally.sum     += tally.hist;
        tally.squared += tally.hist * tally.hist;
        tally.hist     = 0.0;
    }
}
void Estimator::end_cycle( const int N, const double tracks )
{
    for( auto& tally : e_tally ){
        const double mean          = tally.sum / N;
        const double uncer_squared = ( tally.squared / N - mean*mean ) 
                                     / ( N - 1.0 );
        const double uncer_rel     = std::sqrt(uncer_squared) / mean;
        const double FOM           = 1.0 / ( uncer_rel*uncer_rel * tracks );

        tally.mean  += mean;
        tally.uncer += uncer_squared;
        tally.FOM   += FOM;

        tally.sum     = 0.0;
        tally.squared = 0.0;
    }
}
void Estimator::end_simulation( const int N )
{
    for( auto& tally : e_tally ){
        tally.mean      = tally.mean / N;
        tally.uncer     = std::sqrt( tally.uncer ) / N;
        tally.uncer_rel = tally.uncer / tally.mean;
        tally.FOM       = tally.FOM / N;
    }
}

void Estimator::report( std::ostringstream& output )
{
    output << "\n\n";
    output << "Estimator report: " + e_name + "\n";
    for ( int i = 0 ; i < e_name.length()+18 ; i++ ) { output << "="; }
    output << "\n";
	
    for ( auto& tally : e_tally )
    {
        output<<tally.mean<<"  +-  "<<tally.uncer<< " " 
              << "(" << tally.uncer_rel*100<<"%)\n";
    }
}

Tally Estimator::tally( const int i )
{
    return e_tally[i];
}
