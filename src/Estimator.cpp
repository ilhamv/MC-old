#include <cmath>
#include <iostream>
#include <sstream>  // ostringstream

#include "Estimator.h"
#include "Geometry.h"
#include "Algorithm.h"  // binary_search, interpolate
#include "Material.h"
#include "H5Cpp.h"


//=============================================================================
// Scoring Kernel
//   - Supports: Neutron, Track Length, Collision
//=============================================================================

// Neutron
double ScoreKernelNeutron::score( const Particle& P, const double l )
{ 
    return P.weight(); 
}
// Track Length
double ScoreKernelTrackLength::score( const Particle& P, const double l )
{ 
    return P.weight() * l; 
}
// Collision
double ScoreKernelCollision::score( const Particle& P, const double l )
{ 
    return P.weight() / P.cell()->material()->SigmaT( P.energy() ); 
}
// Velocity
double ScoreKernelVelocity::score( const Particle& P, const double l )
{ 
    return P.weight() * P.speed(); 
}
// Track Length Velocity
double ScoreKernelTrackLengthVelocity::score( const Particle& P, 
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
double ScoreFlux::score( const Particle& P, const double l )
{ 
    return s_kernel->score(P,l); 
}

// Absorption
double ScoreAbsorption::score( const Particle& P, const double l )
{ 
    return P.cell()->material()->SigmaA( P.energy() ) * s_kernel->score(P,l); 
}

// Scatter
double ScoreScatter::score( const Particle& P, const double l )
{ 
    return P.cell()->material()->SigmaS( P.energy() ) * s_kernel->score(P,l); 
}

// Scatter old
double ScoreScatterOld::score( const Particle& P, const double l )
{ 
    return P.cell()->material()->SigmaS( P.energy_old() ) * s_kernel->score(P,l);
}

// Capture
double ScoreCapture::score( const Particle& P, const double l )
{ 
    return P.cell()->material()->SigmaC( P.energy() ) * s_kernel->score(P,l); 
}

// Fission
double ScoreFission::score( const Particle& P, const double l )
{ 
    return P.cell()->material()->SigmaF( P.energy() ) * s_kernel->score(P,l); 
}

// NuFission
double ScoreNuFission::score( const Particle& P, const double l )
{ 
    return P.cell()->material()->nuSigmaF( P.energy() ) * s_kernel->score(P,l); 
}

// NuFission Old
double ScoreNuFissionOld::score( const Particle& P, const double l )
{ 
    return P.cell()->material()->nuSigmaF( P.energy_old() ) * s_kernel->score(P,l); 
}

// NuFission Prompt Old
double ScoreNuFissionPromptOld::score( const Particle& P, const double l )
{ 
    return P.cell()->material()->nuSigmaF_prompt( P.energy_old() ) * s_kernel->score(P,l); 
}

// NuFission Delayed Old
double ScoreNuFissionDelayedOld::score( const Particle& P, const double l )
{ 
    return P.cell()->material()->nuSigmaF_delayed( P.energy_old(), cg ) * s_kernel->score(P,l); 
}

// Total
double ScoreTotal::score( const Particle& P, const double l )
{ 
    return P.cell()->material()->SigmaT( P.energy() ) * s_kernel->score(P,l); 
}


//=============================================================================
// Filter
//   - Supports: recently crossed Surface, Cell, Energy, Time
//=============================================================================

// Surface
std::vector<std::pair<int,double>> FilterSurface::idx_l( const Particle& P,
                                                         const double l )
{
    std::vector<std::pair<int,double>> v_i_l;
    std::pair<int,double> i_l;

    // Index location
    i_l.first  = binary_search( P.surface_old()->ID(), f_grid ) + 1;
    // Track length to score
    i_l.second = l;

    v_i_l.push_back(i_l);
    return v_i_l;
}
int FilterSurface::size() { return f_grid.size(); }
// Cell
std::vector<std::pair<int,double>> FilterCell::idx_l( const Particle& P,
                                                      const double l )
{
    std::vector<std::pair<int,double>> v_i_l;
    std::pair<int,double> i_l;

    // Index location
    i_l.first  = binary_search( P.cell()->ID(), f_grid ) + 1;
    // Track length to score
    i_l.second = l;

    v_i_l.push_back(i_l);
    return v_i_l;
}
int FilterCell::size() { return f_grid.size(); }
// Energy
std::vector<std::pair<int,double>> FilterEnergy::idx_l( const Particle& P,
                                                        const double l )
{
    std::vector<std::pair<int,double>> v_i_l;
    std::pair<int,double> i_l;
    
    // Index location
    i_l.first  = binary_search( P.energy(), f_grid );
    // Check if inside the grids
    if ( i_l.first < 0 || i_l.first >= f_Nbin ) { return v_i_l; }
    // Track length to score
    i_l.second = l;

    v_i_l.push_back(i_l);
    return v_i_l;
}
// Energy - old
std::vector<std::pair<int,double>> FilterEnergyOld::idx_l( const Particle& P,
                                                           const double l )
{
    std::vector<std::pair<int,double>> v_i_l;
    std::pair<int,double> i_l;
    
    // Index location
    i_l.first  = binary_search( P.energy_old(), f_grid );
    // Check if inside the grids
    if ( i_l.first < 0 || i_l.first >= f_Nbin ) { return v_i_l; }
    // Track length to score
    i_l.second = l;

    v_i_l.push_back(i_l);
    return v_i_l;
}
// Time bin
std::vector<std::pair<int,double>> FilterTime::idx_l( const Particle& P,
                                                      const double l )
{
    std::vector<std::pair<int,double>> v_i_l;
    std::pair<int,double> i_l;
    
    // Edge bin locations
    int loc1 = binary_search( P.time_old(), f_grid ); // before track
    int loc2 = binary_search( P.time(), f_grid );     // after
    
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
std::vector<std::pair<int,double>> FilterTDMC::idx_l( const Particle& P,
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
int FilterTDMC::size() { return f_grid.size(); }


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
void Estimator::score( const Particle& P, const double l )
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
        double l = MAX_float;
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
            if ( e_idx_l[i][idx[i]].second < EPSILON_float ){ 
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
void Estimator::end_cycle()
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
void Estimator::report( H5::H5File& output )
{
    H5::StrType type_str(0, H5T_VARIABLE);
    H5::DataSpace space_scalar(H5S_SCALAR);
    H5::DataType type_double = H5::PredType::NATIVE_DOUBLE;
    H5::Group group;
    H5::Attribute attribute;

    H5::Group group_root = output.createGroup("/"+e_name);

    // Attribute
    attribute = group_root.createAttribute( "indexing", type_str,space_scalar);
    std::string indexing;
    for( auto& filter : e_filters ){
        indexing += "[" + filter->name() + "]";
    }
    attribute.write( type_str, indexing );
    
    // Filter
    for( auto& filter : e_filters ){
        hsize_t dims[1]; dims[0] = filter->grid().size();
        H5::DataSpace space_filter(1,dims);
        H5::DataSet data_filter = group_root.createDataSet(filter->name(), 
                                                type_double, space_filter );

        data_filter.write(filter->grid().data(), type_double);
        H5::Attribute att = data_filter.createAttribute( "unit", type_str,
                                                         space_scalar );
        att.write( type_str, filter->unit() );
    }

    // Score
    hsize_t dims[e_filters.size()];
    for( int i = 0; i < e_filters.size(); i++ ){
        dims[i] = e_filters[i]->size();
    }
    H5::DataSpace space_score(e_filters.size(),dims);
    int idx = 0;
    std::vector<double> mean(idx_factor[0]);
    std::vector<double> uncer(idx_factor[0]);
    for( auto& score : e_scores ){
        H5::Group group = group_root.createGroup(score->name());
        for( int i = 0; i < idx_factor[0]; i++ ){
            mean[i] = e_tally[i+idx].mean;
            uncer[i] = e_tally[i+idx].uncer;
        }
        H5::DataSet data_mean  = group.createDataSet("mean", type_double,
                                                    space_score );
        H5::DataSet data_uncer = group.createDataSet("uncertainty",type_double,
                                                    space_score );
        data_mean.write(mean.data(), type_double);
        data_uncer.write(uncer.data(), type_double);
        idx += idx_factor[0];
    }
}

Tally Estimator::tally( const int i )
{
    return e_tally[i];
}

unsigned long long Estimator::tally_size()
{
    return e_tally.size();
}

// Scattering simulation estimator
//   It simulates scattering event before scoring
void EstimatorScatter::score( const Particle& P, const double l )
{
    Particle P_simulated = P;
    // Simulate scattering
    P.cell()->material()->nuclide_scatter(P.energy())
            ->scatter()->sample(P_simulated);
    Estimator::score( P_simulated, l );
}
// Fission simulation estimator
//   It simulates fission event before scoring
void EstimatorFission::score( const Particle& P, const double l )
{
    Particle P_simulated = P;
    const double energy_final = P.cell()->material()
                                ->nuclide_nufission( P.energy() )
                                ->fission()->Chi( P.energy() );
    P_simulated.set_energy(energy_final);
    Estimator::score( P_simulated, l );
}
// Prompt Fission simulation estimator
//   It simulates fission event before scoring
void EstimatorFissionPrompt::score( const Particle& P, const double l )
{
    Particle P_simulated = P;
    const double energy_final = P.cell()->material()
                                ->nuclide_nufission_prompt( P.energy() )
                                ->fission()->Chi( P.energy() );
    P_simulated.set_energy(energy_final);
    Estimator::score( P_simulated, l );
}
// Delayed Fission simulation estimator
//   It simulates fission event before scoring
void EstimatorFissionDelayed::score( const Particle& P, const double l )
{
    Particle P_simulated = P;
    const double energy_final = P.cell()->material()
                                ->nuclide_nufission_delayed( P.energy(), cg )
                                ->fission()->ChiD( cg, P.energy() );
    P_simulated.set_energy(energy_final);
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

void EstimatorK::estimate_C( const Particle& P )
{
    k_C += P.cell()->material()->nuSigmaF(P.energy()) * P.weight() 
           / P.cell()->material()->SigmaT(P.energy());
}

void EstimatorK::estimate_TL( const Particle& P, const double l )
{
    k_TL += P.cell()->material()->nuSigmaF(P.energy()) * P.weight() * l;
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
void EstimatorK::report( H5::H5File& output )
{
    H5::DataSet dataset;
    H5::Group group;
    H5::DataSpace space_scalar(H5S_SCALAR);
    H5::DataType type_ull    = H5::PredType::NATIVE_ULLONG;
    H5::DataType type_double = H5::PredType::NATIVE_DOUBLE;
    H5::StrType type_str(0, H5T_VARIABLE);

    group = output.createGroup("/ksearch");
    hsize_t dims[1];
    H5::DataSpace space_vector;
    
    dims[0] = k_cycle.size();
    space_vector = H5::DataSpace(1,dims);
    dataset = group.createDataSet( "k_cycle", type_double, space_vector);
    dataset.write(k_cycle.data(), type_double);
    
    dataset = group.createDataSet( "mean", type_double, space_scalar);
    dataset.write(&k_avg.back(), type_double);
    dataset = group.createDataSet( "uncertainty", type_double, space_scalar);
    dataset.write(&k_uncer.back(), type_double);

    group = output.createGroup("/ksearch/k_active");
    dims[0] = k_avg.size();
    space_vector = H5::DataSpace(1,dims);
    dataset = group.createDataSet( "mean", type_double, space_vector);
    dataset.write(k_avg.data(), type_double);
    dataset = group.createDataSet( "uncertainty", type_double, space_vector);
    dataset.write(k_uncer.data(), type_double);
}
