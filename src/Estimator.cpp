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
//   - Supports: ID(Surface, Cell), Energy, Time
//==============================================================================

// ID
std::vector<std::pair<int,double>> FilterID::idx_tl( const Particle_t& P,
                                                     const double l,
                                                     const double id )
{
    std::vector<std::pair<int,double>> v_i_l;
    std::pair<int,double> i_l;

    // Index location
    i_l.first  = Binary_Search( id, f_grid ) + 1;
    // Track length to score
    i_l.second = l;

    v_i_l.push_back(i_l);
    return v_i_l;
}

// Energy
std::vector<std::pair<int,double>> FilterEnergy::idx_tl( const Particle_t& P,
                                                     const double l,
                                                     const double id )
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

// Time
std::vector<std::pair<int,double>> FilterTime::idx_tl( const Particle_t& P,
                                                     const double l,
                                                     const double told )
{
    std::vector<std::pair<int,double>> v_i_l;
    std::pair<int,double> i_l;
    
    // Edge bin locations
    int loc1 = Binary_Search( told, f_grid );     // before track generation
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
	    new_track = ( f_grid[loc1+1] - told ) * P.speed();
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

///////////
/// Bin ///
///////////

// Energy Bin (multi-score)
void Energy_Bin::score( const Particle_t& P, const std::vector<double>& grid, const double told, const double track /*= 0.0*/ )
{
	// Search bin location	
	const int loc = Binary_Search( P.energy(), grid );

	// Check if inside the grids
	if ( loc >= 0 && loc < Nbin )
	{
		// Iterate over scores
		for ( int i = 0 ; i < Nscore ; i++ )
		{
			tally[loc][i].hist += scores[i]->score( P, track );
		}
		// Note: In case of using cross section table, 
		// might want to pass an index pointing to the XSec table location
		// to avoid XS search for each reaction score!
	}
	// Note: value exactly at the grid point constributes to bin whose upper bound is the value
}

// Time Bin (multi-score)
void Time_Bin::score( const Particle_t& P, const std::vector<double>& grid, const double told, const double track /*= 0.0*/ )
{
	// Since score migth be distributed across bins...
	
	// Pair of bin location and corresponding track to be scored
	std::vector< std::pair<int,double> > loc_track;

	// Search bin location
	int loc1    = Binary_Search( told, grid ); // before track generation
	int loc2    = Binary_Search( P.time()   , grid ); // after
	
	// 1 bin spanned [need to consider if it is outside the time grid]
	if ( loc1 == loc2 ) 
	{
		if ( loc1 >= 0 && loc1 < Nbin )
		{ loc_track.push_back( std::make_pair( loc1, track ) ); }
	}
	
	// >1 bin spanned
	else
	{				
		int    num_bin = loc2 - loc1 - 1; // # of full bins spanned
		double new_track;                 // to hold bin track
		
		// First partial bin [need to consider if it's outside]
		if ( loc1 >= 0 )
		{
			new_track = ( grid[loc1+1] - told ) * P.speed();
			loc_track.push_back( std::make_pair( loc1, new_track ) );
		}
		
		// Intermediate full bin
		for ( int i = 1 ; i <= num_bin ; i++ )
		{
			new_track = ( grid[loc1+i+1] - grid[loc1+i] ) * P.speed();
			loc_track.push_back( std::make_pair( loc1+i, new_track ) );
		}
		
		// Last partial bin [consider if it's outside]
		if ( loc2 < Nbin )
		{
			new_track = ( P.time() - grid[loc2] ) * P.speed();
			loc_track.push_back( std::make_pair( loc2, new_track ) );
		}
		// Note: this algorithm covers the following "extreme" cases
		// 	loc1 : <lowest_grid  or at grid point
		// 	loc2 : >highest_grid or at grid point
	}

	// Iterate over scores
	for ( int i = 0 ; i < Nscore ; i++ )
	{
		// Iterate over scored bins
		for ( auto& LnT : loc_track )
		{
			tally[LnT.first][i].hist += scores[i]->score( P, LnT.second );
		}	
		// Note: In case of using cross section table, 
		// might want to pass an index pointing to the location in XSec table
		// to avoid XS search for each reaction score!
	}
}




/////////////////
/// Estimator ///
/////////////////

// Generic estimator
////////////////////
		
// Add thing to be scored and push new total tally
void Generic_Estimator::addScore( const std::shared_ptr<Score>& S )
{ 
	// Push new score
	scores.push_back( S );
	Nscore++;
		
	// Push new total tally
	Tally T;
	total_tally.push_back( T );
	// Note: Index in total_tally corresponds to its score
}

// Set bin
void Generic_Estimator::setBin( const std::string type, const std::vector<double> gr )
{
	grid = gr;
	Nbin = grid.size() - 1;
	if      ( type == "energy" ) { bin = std::make_shared<Energy_Bin> ( grid, total_tally, scores ); }
	else if ( type == "time" )   { bin = std::make_shared<Time_Bin>   ( grid, total_tally, scores ); }
}

// Score at events
void Generic_Estimator::score( const Particle_t& P, const double told, const double track /*= 0.0*/ )
{
        // Total tallies
        for ( int i = 0 ; i < Nscore ; i++ )
        {
                total_tally[i].hist += scores[i]->score( P, track );
        }

	// Bin tallies, if any
	if ( Nbin != 0 )
	{ bin->score( P, grid, told, track ); }
}

// Tally operations
void Generic_Estimator::tally_endHistory( Tally& T )
{
    T.sum     += T.hist;
    T.squared += T.hist * T.hist;
    T.hist     = 0.0; 
}	
void Generic_Estimator::tally_endCycle( Tally& T, const double trackTime )
{
    const double mean          = T.sum / nhist;
    const double uncer_squared = ( T.squared / nhist - mean*mean ) / ( nhist - 1.0 );
    const double uncer_rel     = std::sqrt(uncer_squared) / mean;
   
    T.mean      += mean;
    T.uncer     += uncer_squared;
    T.FOM       = 1.0 / ( uncer_rel*uncer_rel * trackTime );

    T.sum     = 0.0;
    T.squared = 0.0;
    
};
void Generic_Estimator::tally_average( Tally& T )
{
    T.mean      = T.mean / ncycle;
    T.uncer     = sqrt( T.uncer ) / ncycle;
    T.uncer_rel = T.uncer / T.mean;
    T.FOM       = T.FOM / ncycle;
};

// Closeout
void Generic_Estimator::endHistory()
{
    nhist++;
    for ( int i = 0 ; i < Nscore ; i++ ){
	tally_endHistory( total_tally[i] );

	for ( int j = 0 ; j < Nbin ; j++ ){
	    tally_endHistory( bin->tally[j][i] );
	}
    }
}
void Generic_Estimator::endCycle( const double tracks )
{
    ncycle++;
    for ( int i = 0 ; i < Nscore ; i++ ){
	tally_endCycle( total_tally[i], tracks );
		
	for ( int j = 0 ; j < Nbin ; j++ ){
	    tally_endCycle( bin->tally[j][i], tracks );
	}
    }
    nhist = 0.0;
}

// Report results
void Generic_Estimator::report( std::ostringstream& output )
{
	// Compute mean, variance, and uncertainty of all tallies
	for ( int i = 0 ; i < Nscore ; i++ )
	{
		// Total tally
		tally_average( total_tally[i] );
		
		// Bin tally
		for ( int j = 0 ; j < Nbin ; j++ )
		{
			tally_average( bin->tally[j][i] );
		}	
	}
	
	output << "\n\n";
	output << "Estimator report: " + e_name + "\n";
	for ( int i = 0 ; i < e_name.length()+18 ; i++ ) { output << "="; }
	output << "\n";
	
	// Total tallies
	output << "Total tallies,\n";
	for ( int i = 0 ; i < Nscore ; i++ )
	{
		output << "  " + scores[i]->name() + ":\n";
		output << "  -> Mean     = " << total_tally[i].mean;
	       	output << "  +/-  " << total_tally[i].uncer;
		output << "  (" << total_tally[i].uncer_rel * 100.0;
		output << "%)\n";
		output << "  [F.O.M.: " << total_tally[i].FOM;
	        output << "]\n\n\n";
	}

	// Bin tallies
	if ( Nbin != 0 )
	{
		output << "Bin tallies," << std::endl;

		output << "bin#\t";
		output << "lower(" + bin->unit + ")\tupper(" + bin->unit + ")\t";
		for ( int i = 0 ; i < Nscore ; i++ ) 
		{ 
			output << std::setw(12) << std::left << scores[i]->name() << "\t"; 
			output << std::setw(12) << std::left << "uncertainty\t"; 
		}
		output << "\n";
	
		output << "----\t" << "------------\t" << "------------\t";
		for ( int i = 0 ; i < Nscore ; i++ ) 
		{ 
			output << "------------\t" << "------------\t"; 
		}
		output << "\n";
	
		for ( int j = 0 ; j < Nbin ; j++ )
		{
			output << j+1;
			output << "\t" << std::setw(12) << std::left << grid[j];
		        output << "\t" << std::setw(12) << std::left << grid[j+1]; 
			output << "\t";

			for ( int i = 0 ; i < Nscore ; i++ )
			{
				output << std::setw(12) << bin->tally[j][i].mean << "\t";
			        output << std::setw(12) << bin->tally[j][i].uncer << "\t";
			}
			output << "\n";
		}
	}
}



// Homogenized MG Constant Generator
////////////////////////////////////
/*
// Set bin (or group structure), and calculate Chi group constants
void MGXS_Estimator::setBin( const std::string type, const std::vector<double> gr )
{
	// Set energy grid points
	grid = gr; 

	///////////////////////////////////////////////////////
	// Set bin with multi score for basic MGXS constants //

	// Set scores = Flux, Capture, Fission, nuFission, Total, Scatter
	scores.push_back( std::make_shared<Flux_Score>()      );
	scores.push_back( std::make_shared<Capture_Score>()   );
	scores.push_back( std::make_shared<Fission_Score>()   );
	scores.push_back( std::make_shared<nuFission_Score>() );
	scores.push_back( std::make_shared<Total_Score>()     );
	scores.push_back( std::make_shared<Scatter_Score>()   );
	Nscore = scores.size();

	// Set a vector of tallies corresponding to each score
	std::vector<Tally> Tvec; 
	Tally T;
	Tvec.resize( Nscore, T );

	// Set energy bins containing copies of the vector of tallies
	bin = std::make_shared<Energy_Bin> ( grid, Tvec, scores );
	Nbin = grid.size() - 1.0;

	// Basic MGXS constants bin done//
	//////////////////////////////////


	///////////////////////////////////////////////
	// Set Legendre scattering components tensor //

	// Set scores = Scatter	
	std::vector<std::shared_ptr<Score>> temp_scores;
	temp_scores.push_back( std::make_shared<Scatter_Score>() );
	
	// Set a temporary vector of tallies corresponding to each score
	std::vector<Tally> temp_Tvec; 
	temp_Tvec.push_back( T ); // Previous single tally T is employed

	// Construct the legendre scattering components tensor and group Chi
	for ( int i = 0 ; i < Nbin ; i++ )
	{
		std::shared_ptr<Bin_t> temp_bin = std::make_shared<Energy_Bin> ( grid, temp_Tvec, temp_scores ); // Bin of initial energy with only one tally for scattering score
	}
	
	for ( int j = 0 ; j <= N ; j++ )
	{
		std::vector<std::shared_ptr<Bin_t>> temp_matrix; // Temporary matrix ( final energy bin x initial energy bin )
		                                                 // to be pushed to the legendre scattering components tensor
		for ( int i = 0 ; i < Nbin ; i++ )
		{
			// Temporary bin of initial energy with only one tally for scattering score
			std::shared_ptr<Bin_t> t_bin = std::make_shared<Energy_Bin> ( grid, temp_Tvec, temp_scores ); 
			// Push temporary bin into temporary matrix
			temp_matrix.push_back( t_bin ); 
		}
		tensor_bin.push_back( temp_matrix ); // Push temmporary matrix into the tensor
	}
	
	// Set legendre polynomials storage
	Pl.push_back( 1.0 ); // P0
	Pl.push_back( 1.0 ); // P1 (always constructed by default, even if N=0)
	for ( int i = 2 ; i <= N ; i++ ) { Pl.push_back( 1.0 ); }

	// Legendre Scattering components tensor done//
	///////////////////////////////////////////////

}


// Calculate Legendre polynomials at mu
void MGXS_Estimator::calculatePl( const double mu )
{
	Pl[1] = mu;
	for ( int n = 2 ; n <=N ; n++ )
	{ Pl[n] = ( ( 2.0*n - 1.0 ) * mu * Pl[n-1] - ( n - 1.0 ) * Pl[n-2] ) / n; }
}


// Score at events
void MGXS_Estimator::score( const Particle_t& P, const double told, const double track  )
{
	// Flux, Capture, Fission, nuFission, Total, Scater Tallies
	bin->score( P, grid, told, track );
	
	// Scattering matrix Tallies
	// First, simulate scattering reaction to get post-scattering particle information
	Particle_t P_final = P;
	P.cell()->simulate_scatter( P_final );

	// Compute scattering cosine and considered legendre polynomials
	const double mu = P.dir().x*P_final.dir().x + P.dir().y*P_final.dir().y + P.dir().z*P_final.dir().z;
	calculatePl( mu );
	
	// Then, find the final energy bin location
	const int loc = Binary_Search( P_final.energy(), grid );
	
	// Final energy bin location is set, now score into the appropriate bins
	// Iterate over considered Legendre components
	for ( int i = 0 ; i <= N ; i++ )
	{
		tensor_bin[i][loc]->score( P, grid, told, track*Pl[i] );
	}
}


// Closeout history
// Update the sum and sum of squared, and restart history sum of all tallies
void MGXS_Estimator::endHistory()
{
    nhist++;
    // Flux, Capture, Fission, nuFission, Total, Scater Tallies
    for ( int i = 0 ; i < Nscore ; i++ ){
        for ( int j = 0 ; j < Nbin ; j++ ){
	    tally_endHistory( bin->tally[j][i] );
	}
    }
	
    // Scattering matrix Tallies
    for ( int i = 0 ; i < Nbin ; i++ ){
        for ( int j = 0 ; j < Nbin ; j++ ){
	    for ( int n = 0 ; n <= N ; n++ ){
		tally_endHistory( tensor_bin[n][i]->tally[j][0] );
	    }
	}
    }
}

void MGXS_Estimator::endCycle( const double tracks )
{
    ncycle++;
    // Flux, Capture, Fission, nuFission, Total, Scater Tallies
    for ( int i = 0 ; i < Nscore ; i++ ){
        for ( int j = 0 ; j < Nbin ; j++ ){
	    tally_endCycle( bin->tally[j][i], tracks );
	}
    }
	
    // Scattering matrix Tallies
    for ( int i = 0 ; i < Nbin ; i++ ){
        for ( int j = 0 ; j < Nbin ; j++ ){
	    for ( int n = 0 ; n <= N ; n++ ){
		tally_endCycle( tensor_bin[n][i]->tally[j][0], tracks );
	    }
	}
    }
    nhist = 0.0;
}

// Report results
void MGXS_Estimator::report( std::ostringstream& output )
{
	// Compute mean, variance, uncertainty, and F.O.M of all tallies
	// Flux, Capture, Fission, nuFission, Total, Scater Tallies
	for ( int i = 0 ; i < Nscore ; i++ )
	{
		for ( int j = 0 ; j < Nbin ; j++ )
		{
			tally_average( bin->tally[j][i] );
		}	
	}	
	// Scattering matrix Tallies
	for ( int i = 0 ; i < Nbin ; i++ )
	{
		for ( int j = 0 ; j < Nbin ; j++ )
		{
			for ( int n = 0 ; n <= N ; n++ )
			{
				tally_average( tensor_bin[n][i]->tally[j][0] );
			}	
		}	
	}	
	

	// Convert tally (except flux) into homogenized multigroup constant [mean and uncertainty]
	for ( int j = 0 ; j < Nbin ; j++ )
	{
		const double B = bin->tally[j][0].uncer / bin->tally[j][0].mean;
		for ( int i = 1 ; i < Nscore ; i++ )
		{
			const double A = bin->tally[j][i].uncer / bin->tally[j][i].mean;
			bin->tally[j][i].uncer = std::sqrt( A*A + B*B );

			bin->tally[j][i].mean = bin->tally[j][i].mean / bin->tally[j][0].mean;
			bin->tally[j][i].uncer = bin->tally[j][i].uncer * bin->tally[j][i].mean;
		}	
		for ( int i = 0 ; i < Nbin ; i++ )
		{
			const double C = tensor_bin[0][i]->tally[j][0].uncer / tensor_bin[0][i]->tally[j][0].mean;
			tensor_bin[0][i]->tally[j][0].uncer = std::sqrt( C*C + B*B );

			tensor_bin[0][i]->tally[j][0].mean = tensor_bin[0][i]->tally[j][0].mean / bin->tally[j][0].mean;
			tensor_bin[0][i]->tally[j][0].uncer = tensor_bin[0][i]->tally[j][0].uncer * tensor_bin[0][i]->tally[j][0].mean;
		}	
	}	
	

	output << "\n\n";
	output << "Estimator report: " + e_name + "\n";
	for ( int i = 0 ; i < e_name.length()+18 ; i++ ) { output << "="; }
	output << "\n";
	
	// Note energy group number is reversed for convenience
	if ( Nbin != 0 )
	{
		output << "Group structure," << std::endl;

		output << "g\t";
		output << "upper(" + bin->unit + ")\tlower(" + bin->unit + ")\t";
		output << std::setw(12) << std::left << scores[0]->name() << "\t"; 
		output << std::setw(12) << std::left << "uncertainty\t"; 
		output << "\n";
	
		output << "----\t" << "------------\t" << "------------\t" << "------------\t" << "------------\t";
		output << "\n";
	
		for ( int j = 0 ; j < Nbin ; j++ )
		{
			output << j+1;
			output << "\t" << std::setw(12) << std::left << grid[Nbin-j];
		        output << "\t" << std::setw(12) << std::left << grid[Nbin-j-1]; 
		        output << "\t" << std::setw(12) << bin->tally[Nbin-j-1][0].mean; 
		        output << "\t" << std::setw(12) << bin->tally[Nbin-j-1][0].uncer; 
			output << "\n";
		}
		output << "\n\nHomogenized MG Constants," << std::endl;

		output << "g\t";
		for ( int i = 1 ; i < Nscore ; i++ ) 
		{ 
			output << std::setw(12) << std::left << scores[i]->name() << "\t"; 
			output << std::setw(12) << std::left << "uncertainty\t"; 
		}
		output << "\n";
	
		output << "----\t" << "------------\t";
		for ( int i = 1 ; i < Nscore ; i++ ) 
		{ 
			output << "------------\t" << "------------\t"; 
		}
		output << "\n";
	
		for ( int j = 0 ; j < Nbin ; j++ )
		{
			output << j+1;
			output << "\t";

			for ( int i = 1 ; i < Nscore ; i++ )
			{
				output << std::setw(12) << bin->tally[Nbin-j-1][i].mean << "\t";
			        output << std::setw(12) << bin->tally[Nbin-j-1][i].uncer << "\t";
			}
			output << "\n";
		}
		
		output << "\n\nLegendre Scattering Matrix Components," << std::endl;
		
		for ( int n = 0 ; n <= N ; n++ )
		{
			output << "Sigma_s" << n << "\n";

			output << "g\t";
			for ( int j = 0 ; j < Nbin ; j++ )
			{
				const std::string g_prime = std::to_string(j+1);
				output << std::setw(12) << std::left << "g -> " + g_prime << "\t"; 
				output << std::setw(12) << std::left << "uncertainty\t"; 
			}
			output << "\n";
	
			output << "----\t";
			for ( int j = 0 ; j < Nbin ; j++ )
			{ 
				output << "------------\t" << "------------\t"; 
			}
			output << "\n";
	
			for ( int j = 0 ; j < Nbin ; j++ )
			{
				output << j+1;
				output << "\t";
				for ( int i = 0 ; i < Nbin ; i++ )
				{
					output << std::setw(12) << tensor_bin[n][Nbin-i-1]->tally[Nbin-j-1][0].mean << "\t"; 
					output << std::setw(12) << tensor_bin[n][Nbin-i-1]->tally[Nbin-j-1][0].uncer << "\t"; 
				}
				output << "\n";
			}
			output << "\n";
		}
	}
}
*/
