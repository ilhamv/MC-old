#include <cmath>
#include <iostream>
#include <sstream>  // ostringstream

#include "Estimator.h"
#include "Geometry.h"
#include "Solver.h"  // Binary_Search, Linterpolate
#include "H5Cpp.h"


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
// Velocity
double ScoreKernelVelocity::score( const Particle_t& P, const double l )
{ 
    return P.weight() * P.speed(); 
}
// Track Length Velocity
double ScoreKernelTrackLengthVelocity::score( const Particle_t& P, 
                                              const double l )
{ 
    return P.weight() * l * P.speed(); 
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


//=============================================================================
// Filter
//   - Supports: recently crossed Surface, Cell, Energy, Time
//=============================================================================

// Surface
std::vector<std::pair<int,double>> FilterSurface::idx_l( const Particle_t& P,
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
std::vector<std::pair<int,double>> FilterCell::idx_l( const Particle_t& P,
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
std::vector<std::pair<int,double>> FilterEnergy::idx_l( const Particle_t& P,
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
// Energy - old
std::vector<std::pair<int,double>> FilterEnergyOld::idx_l( const Particle_t& P,
                                                           const double l )
{
    std::vector<std::pair<int,double>> v_i_l;
    std::pair<int,double> i_l;
    
    // Index location
    i_l.first  = Binary_Search( P.energy_old(), f_grid );
    // Check if inside the grids
    if ( i_l.first < 0 && i_l.first >= f_Nbin ) { return v_i_l; }
    // Track length to score
    i_l.second = l;

    v_i_l.push_back(i_l);
    return v_i_l;
}
// Time bin
std::vector<std::pair<int,double>> FilterTime::idx_l( const Particle_t& P,
                                                      const double l )
{
    std::vector<std::pair<int,double>> v_i_l;
    std::pair<int,double> i_l;
    
    // Edge bin locations
    int loc1 = Binary_Search( P.time_old(), f_grid ); // before track
    int loc2 = Binary_Search( P.time(), f_grid );     // after
    
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
std::vector<std::pair<int,double>> FilterTDMC::idx_l( const Particle_t& P,
                                                      const double l )
{
    std::vector<std::pair<int,double>> v_i_l;
    std::pair<int,double> i_l;
    
    if( P.time() != f_grid[P.tdmc()] ) { return v_i_l; }
    // Index location
    i_l.first  = P.tdmc(); 
    // Velocity to score
    i_l.second = l;

    v_i_l.push_back(i_l);
    return v_i_l;
}


//=============================================================================
// Basic Estimator
//=============================================================================

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
void Estimator::score( const Particle_t& P, const double l )
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
void Estimator::end_cycle( const double tracks )
{
    for( auto& tally : e_tally ){
        const double mean          = tally.sum / e_Nsample;
        const double uncer_squared = ( tally.squared / e_Nsample - mean*mean ) 
                                     / ( e_Nsample - 1.0 );

        tally.mean  += mean;
        tally.uncer += uncer_squared;

        tally.sum     = 0.0;
        tally.squared = 0.0;
    }
}
void Estimator::end_simulation()
{
    for( auto& tally : e_tally ){
        tally.mean      = tally.mean / e_Nactive;
        tally.uncer     = std::sqrt( tally.uncer ) / e_Nactive;
    }
}
void Estimator::report( std::ostringstream& output, H5::H5File& output_H5  )
{
    output << "\n\n";
    output << "Estimator report: " + e_name + "\n";
    for ( int i = 0 ; i < e_name.length()+18 ; i++ ) { output << "="; }
    output << "\n";
	
    for ( auto& tally : e_tally )
    {
        output<<tally.mean<<"  +-  "<<tally.uncer<< " " 
              << "(" << tally.uncer/tally.mean*100<<"%)\n";
    }
    
    // Populate H5 output file
    H5::StrType type_str(0, H5T_VARIABLE);
    H5::DataSpace att_space(H5S_SCALAR);

    H5::Group e_group = output_H5.createGroup("/"+e_name);
    int idx = 0;
    std::vector<double> mean(idx_factor[0]);
    std::vector<double> uncer(idx_factor[0]);
    for( auto& score : e_scores ){
        output_H5.createGroup("/"+e_name+"/"+score->name());
        for( int i = 0; i < idx_factor[0]; i++ ){
            mean[i] = e_tally[i+idx].mean;
            uncer[i] = e_tally[i+idx].uncer;
        }
        hsize_t dims[e_filters.size()];
        for( int i = 0; i < e_filters.size(); i++ ){
            dims[i] = e_filters[i]->grid().size();
        }
        H5::DataSpace data_space_mean(e_filters.size(),dims);
        H5::DataSpace data_space_uncer(e_filters.size(),dims);
        H5::DataSet data_mean = 
            output_H5.createDataSet("/"+e_name+"/"+score->name()+"/"+"mean",
                                    H5::PredType::NATIVE_DOUBLE, 
                                    data_space_mean );
        H5::DataSet data_uncer = 
            output_H5.createDataSet("/"+e_name+"/"+score->name()+"/"
                                    +"uncertainty", 
                                    H5::PredType::NATIVE_DOUBLE, 
                                    data_space_uncer );
        data_mean.write(mean.data(), H5::PredType::NATIVE_DOUBLE);
        data_uncer.write(uncer.data(), H5::PredType::NATIVE_DOUBLE);
        idx += idx_factor[0];
    }

    for( auto& filter : e_filters ){
        hsize_t dims[1]; dims[0] = filter->grid().size();
        H5::DataSpace data_space(1,dims);
        H5::DataSet data_set = 
            output_H5.createDataSet("/"+e_name+"/"+filter->name(), 
                                    H5::PredType::NATIVE_DOUBLE, data_space );

        data_set.write(filter->grid().data(), H5::PredType::NATIVE_DOUBLE);
        H5::Attribute att = data_set.createAttribute( "unit", type_str,
                                                      att_space );
        att.write( type_str, filter->unit() );
    }
    
    H5::Attribute att = e_group.createAttribute( "indexing", type_str, 
                                                 att_space );
    std::string indexing;
    for( auto& filter : e_filters ){
        indexing += "[" + filter->name() + "]";
    }
    att.write( type_str, indexing );
}

Tally Estimator::tally( const int i )
{
    return e_tally[i];
}

// Scattering simulation estimator
//   It simulates scattering event before scoring
void EstimatorScatter::score( const Particle_t& P, const double l )
{
    Particle_t P_simulated = P;
    P.cell()->simulate_scatter(P_simulated);
    Estimator::score( P_simulated, l );
}
// Fissio simulation estimator
//   It simulates scattering event before scoring
void EstimatorFission::score( const Particle_t& P, const double l )
{
    Particle_t P_simulated = P;
    const double energy_final = P.cell()->material()
                                ->nuclide_nufission( P.energy() )
                                ->Chi( P.energy() );
    P_simulated.setEnergy(energy_final);
    Estimator::score( P_simulated, l );
}


//=============================================================================
// k-eigenvalue Estimator
//=============================================================================

EstimatorK::EstimatorK( const unsigned long long Ncycle, 
                        const unsigned long long Na,
                        const unsigned long long Ns )
{
    Nactive = Na;
    Nsample = Ns;
    k_cycle.resize(Ncycle,0.0);
    k_avg.resize(Nactive,0.0);
    k_uncer.resize(Nactive,0.0);
}

void EstimatorK::estimate_C( const Particle_t& P )
{
    k_C += P.cell()->nuSigmaF(P.energy()) * P.weight() 
           / P.cell()->SigmaT(P.energy());
}

void EstimatorK::estimate_TL( const Particle_t& P, const double l )
{
    k_TL += P.cell()->nuSigmaF(P.energy()) * P.weight() * l;
}

void EstimatorK::end_history()
{
    k_sum_C  += k_C;
    k_sum_TL += k_TL;
    k_sq_C   += k_C*k_C;
    k_sq_TL  += k_TL*k_TL;
    k_C  = 0;
    k_TL = 0;
}
void EstimatorK::report_cycle( const bool tally )
{
    const double mean_C  = k_sum_C / Nsample;
    const double mean_TL = k_sum_TL / Nsample;
    const double mean    = (mean_C + mean_TL)/2;

    k_cycle[icycle] = mean;
    k = k_cycle[icycle];
    std::cout<<icycle+1<<"   "<<k_cycle[icycle];

    if(tally){
        const double uncer_squared_C  = ( k_sq_C / Nsample - mean_C*mean_C ) 
                                        / ( Nsample - 1.0 );
        const double uncer_squared_TL = ( k_sq_TL / Nsample - mean_TL*mean_TL )
                                        / ( Nsample - 1.0 );
        const double uncer_squared    = uncer_squared_C + uncer_squared_TL;

        Navg++;
        mean_accumulator     += mean;
        uncer_sq_accumulator += uncer_squared;
        k_avg[Navg-1]   = mean_accumulator/Navg;
        k_uncer[Navg-1] = std::sqrt(uncer_sq_accumulator)/Navg/2;
        std::cout<<"   "<<k_avg[Navg-1]<<"   +/-   "<<k_uncer[Navg-1];
    }
    
    std::cout<<"\n";
    k_sum_C  = 0.0;
    k_sum_TL = 0.0;
    k_sq_C   = 0.0;
    k_sq_TL  = 0.0;
    icycle++;
}
