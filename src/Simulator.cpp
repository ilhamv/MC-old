#include <vector>  // vector
#include <iostream>// cout
#include <cstring> // string, strcmp
#include <memory>  // shared_ptr, make_shared
#include <stack>   // stack
#include <cmath>   // exp
#include <sstream> // ostringstream
#include <fstream> // ifstream

#include "Simulator.h"
#include "pugixml.hpp"
#include "Algorithm.h"
#include "Random.h"
#include "H5Cpp.h"


//=============================================================================
// Search Cell
//=============================================================================

std::shared_ptr<Cell> Simulator::search_cell( const Point& p )
{
    for( const auto& C : Cells ){
        if ( C->test_point( p ) ){ return C; }
    }
    std::cout<< "[WARNING] A particle is lost:\n( x, y, z )  (" << p.x << ", " 
             << p.y << ", " << p.z << " )\n";
    std::exit(EXIT_FAILURE);
}



//=============================================================================
//=============================================================================
// Start of Simulator Set Up
//=============================================================================
//=============================================================================

template< typename T >
std::shared_ptr<T> Simulator::find_by_name( 
        const std::vector<std::shared_ptr<T>>& vec,
        const std::string name )
{
    for ( auto& v : vec ){
	if ( v->name() == name ) { return v; }
    }
    return nullptr;
}

void Simulator::set_nuclide( const std::string name, const std::string label, 
                             std::shared_ptr<Nuclide_t>& Nuc )
{
    std::string dirName = "./xs_library/" + name + ".txt";
    std::ifstream xs_file (dirName);
    if ( !xs_file ){
	std::cout << "unknown nuclide " << name << std::endl;
	throw;
    }

    // Data loading
    double A;    // Nuclide mass
    std::vector<double> a; // Chi spectrum parameter a
    std::vector<double> b; // Chi spectrum parameter b
    auto E = std::make_shared< std::vector<double> >();
    std::vector<double> sigmaS;
    std::vector<double> sigmaC;
    std::vector<double> sigmaF;
    std::vector<double> nu;
    double c1, c2, c3, c4, c5;

    // Nuclide mass, 1st line
    if ( xs_file >> c1 ) { A = c1; }
    else { std::cout<< "Failed to read A in library file " << name << "\n"; throw; }
    Nuc = std::make_shared<Nuclide_t> ( label, A );

    // Chi sectrum parameters, 2nd to 4th lines
    for ( int i = 0 ; i < 3 ; i++ )
    {
	if ( xs_file >> c1 >> c2 )
	{
	    a.push_back( c1 );	
	    b.push_back( c2 );	
	}
	else { std::cout<< "Faled to read ab in library file " << name << "\n"; throw; }
    }

    // Cross sections
    // 	Column --> what?
    // 	1 --> energy (eV)
    // 	2 --> sigmaS
    // 	3 --> sigmaC
    // 	4 --> sigmaF
    // 	5 --> nu
    while ( xs_file >> c1 >> c2 >> c3 >> c4 >> c5 )
    {
	E->push_back(c1);
        sigmaS.push_back(c2);
	sigmaC.push_back(c3);
	sigmaF.push_back(c4);
	nu.push_back(c5);
    }
        
    Nuc->setTable( E );
    auto XS_S = std::make_shared<Table_XSec> ( E, sigmaS );
    auto XS_C = std::make_shared<Table_XSec> ( E, sigmaC );
    auto XS_F = std::make_shared<Table_XSec> ( E, sigmaF );
    auto XS_nu = std::make_shared<Table_XSec> ( E, nu );
        
    // Set reactions
    Nuc->addReaction( std::make_shared<Capture_Reaction> ( XS_C ) );
    Nuc->addReaction( std::make_shared< Scatter_Reaction > ( XS_S, std::make_shared< DistributionIsotropicScatter > (), A ) );
    // Fissionable check
    if ( nu[ nu.size() / 2 ] > 0.0 )
    {
	Nuc->addReaction( std::make_shared< Fission_Reaction > ( XS_F, XS_nu, std::make_shared< DistributionWatt > ( a, b ) ) );
    }
}

Simulator::Simulator( const std::string input_dir )
{
    io_dir = input_dir+"/";


    // XML input file
    std::string input_name = io_dir + "input.xml";
    pugi::xml_document input_file;
    pugi::xml_parse_result load_result = input_file.load_file( input_name.c_str() );

    // Able to load file?
    if ( ! load_result ) 
    {
        std::cout<< "Unable to load input file " << input_name << ":\n";
	std::cout<< load_result.description() << std::endl;
        std::exit(EXIT_FAILURE);
    }

    
//=============================================================================
// Set simulation
//=============================================================================

pugi::xml_node input_simulation  = input_file.child("simulation");    
pugi::xml_node input_description = input_simulation.child("description");
pugi::xml_node input_ksearch     = input_simulation.child("ksearch");
pugi::xml_node input_tdmc        = input_simulation.child("tdmc");

// Description
simulation_name = input_description.attribute("name").value();
Nsample         = input_description.attribute("samples").as_double();

// Cut-off
if( input_file.child("cut-off").attribute("energy") ){
    Ecut_off = input_file.child("cut-off").attribute("energy").as_double();
}
if( input_file.child("cut-off").attribute("time") ){
    tcut_off = input_file.child("cut-off").attribute("time").as_double();
}

// ksearch
if( input_ksearch ){
    ksearch  = true;
    unsigned long long active = input_ksearch.attribute("active_cycles")
                                .as_double();
    Npassive = input_ksearch.attribute("passive_cycles").as_double();
    Ncycle   = active + Npassive;
    mode = "k-eigenvalue";
    k_estimator = std::make_shared<EstimatorK>(Ncycle, Ncycle-Npassive, 
                                               Nsample);
}

//=============================================================================
// TDMC
//=============================================================================

if( input_tdmc ){
    tdmc = true;
    const std::string grid_string = input_tdmc.attribute("time").value();
    std::istringstream  iss( grid_string );
    for( double s; iss >> s; ) { tdmc_time.push_back(s); }
    if( input_ksearch ){
        std::cout<<"ksearch and tdmc could not coexist\n";
        std::exit(EXIT_FAILURE);
    }
    tcut_off = tdmc_time.back();
    if( input_tdmc.attribute("split") ){
        tdmc_split = input_tdmc.attribute("split").as_double();
    }
    mode = "time-dependent";
}
        
	// Set user distributuions
  	pugi::xml_node input_distributions = input_file.child("distributions");

  	// Find total number of distributions
  	unsigned int num_distributions = 0;
  	for ( const auto& d : input_distributions ) { num_distributions++; }

  	// Since distributions may depend on other distributions, need to iterate
  	unsigned int set_distributions = 0;
  	while ( set_distributions < num_distributions ) 
	{
    		int previous_set_distributions = set_distributions;

    		for ( const auto& d : input_distributions ) 
		{
      			const std::string type = d.name();
      			const std::string name = d.attribute("name").value();
      			const std::string data = d.attribute("datatype").value();

      			// Double
			if ( data == "double" ) 
			{
        			// Skip rest of loop if distribution already done
        			if ( find_by_name( Distribution_Double, name ) ) { continue; }
	        		std::shared_ptr<Distribution<double>> Dist;
        			
				// Delta-double
				if ( type == "delta" ) 
				{
          				const double val = d.attribute("val").as_double();
          				Dist = std::make_shared< DistributionDelta< double > > ( val, name );
	        		}

        			// Uniform-double
				else if ( type == "uniform" ) 
				{
          				const double a = d.attribute("a").as_double();
          				const double b = d.attribute("b").as_double();
	          			Dist = std::make_shared< DistributionUniform > ( a, b, name );
        			}
        			
        		        // Watt spectrum-double
                		else if (type == "watt" )
		                {
					// Still U-235 by default
					// Next: <watt name="" dtype="double" fissile="U-235"
					// 	Build a function("nuclide name") returning a pair of vectors, a and b
					std::vector<double> a;
					std::vector<double> b;
					a.push_back( 0.988 );
					a.push_back( 0.988 );
					a.push_back( 1.028 );
					b.push_back( 2.249 );
					b.push_back( 2.249 );
					b.push_back( 2.084 );
					Dist = std::make_shared< DistributionWatt > ( a, b, name );
		                }
        	        
		                // Unknown
        		        else 
	                	{
	        	  		std::cout << "unsupported distribution with data type " << data << std::endl;
        	  			throw;
	        		}
				
				// Push new double-distribution
		        	Distribution_Double.push_back( Dist );
		      	}

			// 3D point
			else if ( data == "point" ) 
			{
        			// Skip rest of loop if distribution already done
        			if ( find_by_name( Distribution_Point, name ) ) { continue; }
        			std::shared_ptr< Distribution< Point > > Dist;
		        	
				// Delta-point
				if ( type == "delta" ) 
				{
	          			const double x = d.attribute("x").as_double(); 
	          			const double y = d.attribute("y").as_double(); 
        	  			const double z = d.attribute("z").as_double();         
	          			Dist = std::make_shared< DistributionDelta< Point > > ( Point( x, y, z ), name );
	        		}
        			
				// Isotropic-point
				else if ( type == "isotropic" ) 
				{
          				Dist = std::make_shared< DistributionIsotropicDirection > ( name );
        			}
        			
				// XYZ-point
        			else if ( type == "independentXYZ" ) 
				{
		          		std::shared_ptr< Distribution<double> > distX = find_by_name( Distribution_Double, d.attribute("x").value() ); 
        		  		std::shared_ptr< Distribution<double> > distY = find_by_name( Distribution_Double, d.attribute("y").value() ); 
          				std::shared_ptr< Distribution<double> > distZ = find_by_name( Distribution_Double, d.attribute("z").value() ); 

          				// if any of these distributions have not yet been resolved, skip to the end of the loop
	          			if ( !distX || !distY || !distZ ) { continue; }

		          		Dist = std::make_shared< DistributionIndepndentXYZ > ( distX, distY, distZ, name );
		        	}
		        	
				// Unknown
				else 
				{
          				std::cout << "unsupported " << data << " distribution of type " << type << std::endl;
		          		throw;
		        	}
				
				// Push new point-distribution
        			Distribution_Point.push_back( Dist );
      			}
      			
			// Unknown datatype
			else 
			{
        			std::cout << "unsupported distribution with data type " << data << std::endl;
        			throw;
      			}
      			
			// if we reach here, assume distribution has been set
      			set_distributions++;
		}
    		
		// check to see if number of distributions has increased, if not, caught in an infinite loop
    		if ( previous_set_distributions == set_distributions ) 
		{ 
      			std::cout << "distributions could not be resolved. " << std::endl;
      			throw;
    		}
  	}	

    // Set nuclides
    pugi::xml_node input_nuclides = input_file.child("nuclides");
    for ( const auto& n : input_nuclides )
    {
	std::shared_ptr<Nuclide_t> Nuc;
        const  std::string name   = n.attribute("name").value(); // Nuclide name (or label)
	double             Amass = 1e9;                          // Default nuclide mass
	const  std::string n_tag = n.name();                    // Nuclide tag
		
        if ( n_tag != "nuclide" )
        {
            std::cout<< "[ERROR-INPUT] Unsupported tag under nuclides\n";
            std::exit(EXIT_FAILURE);
        }

        // Standard ZAID nuclide
        if ( n.attribute("ZAID") )
        {
	    set_nuclide( n.attribute("ZAID").value(), name, Nuc );
        }

        // User defined nuclide
        else
        {
    	    // Provided nuclide mass input
	    if ( n.attribute("A") ) 
    	    {
	        Amass = n.attribute("A").as_double();
    	    }

	    Nuc   = std::make_shared<Nuclide_t> ( name, Amass );

    	    // Add nuclide reactions
    	    for ( const auto& r : n.children() ) 
	    {
            	const std::string       rxn_type = r.name();
			
		// Set XSec
		std::shared_ptr<XSec_t> XS;
		
		// Constant XSec
		if ( r.attribute("xs") ) 
    		{
      		    const double xs = r.attribute("xs").as_double();
		    XS = std::make_shared<Constant_XSec> ( xs );
    		}
		// Table look-up XSec
		else if ( r.attribute("xs_file") ) 
		{
                    std::string filename;
	            std::string dirName = r.attribute("xs_file").value();
                
	            // cross section loading
	            std::ifstream xs_file (dirName);
		    auto E_vec = std::make_shared<std::vector<double>>();
        	    std::vector<double> XS_vec;
        	    double c1,c2;
                
		    if (xs_file.is_open())
                    {
                	while(xs_file >> c1 >> c2)
	        	{
        	    	    E_vec->push_back(c1);
                    	    XS_vec.push_back(c2); //4th column is scatter
                	}
                	    xs_file.close();
                	    XS = std::make_shared<Table_XSec> ( E_vec, XS_vec );
	                Nuc->setTable( E_vec );
                    }
                    else
		    { 
			std::cout << "[ERROR-INPUT] Unable to open xs table file for reaction " << rxn_type << " in nuclide " << name << std::endl;
                        std::exit(EXIT_FAILURE);
		    }
		}
		// Unknown XSec type
		else
		{
		    std::cout << "[ERROR-INPUT] Appropriate cross section type for reaction " << rxn_type << " is required" << std::endl;
                    std::exit(EXIT_FAILURE);
		}
		
		// Capture
		if ( rxn_type == "capture" )
		{
        	    Nuc->addReaction( std::make_shared<Capture_Reaction> ( XS ) );
	      	}      

		// Scatter
		else if ( rxn_type == "scatter" )
		{
		    // Set scattering cosine distribution
		    if ( !r.attribute("distribution") ) 
		    { 
			std::cout << "[ERROR-INPUT] Scattering cosine distribution is required for scattering reaction " << std::endl;
                        std::exit(EXIT_FAILURE);
		    }
		
		    const std::string dist_name = r.attribute("distribution").value();
  		    std::shared_ptr< Distribution<double> > f_mu = find_by_name( Distribution_Double, dist_name );

                    const std::string mod = r.attribute("model").value();
                    if ( mod == "zero" )
                    {
		        Nuc->addReaction( std::make_shared< Scatter_Zero_Reaction > ( XS, f_mu, Amass ) );
                    }
                    else
                    {
		        Nuc->addReaction( std::make_shared< Scatter_Reaction > ( XS, f_mu, Amass ) );
                    }
      		}
		
		// Fission
		else if ( rxn_type == "fission" )
		{
		    auto nubar = std::make_shared<Constant_XSec> ( r.attribute("nubar").as_double() );
		    
                    const std::string dist_name = r.attribute("chi").value();
  		    std::shared_ptr< Distribution<double> > watt = find_by_name( Distribution_Double, dist_name );
                
		    Nuc->addReaction( std::make_shared< Fission_Reaction > ( XS, nubar, watt ) );
		} 
	    } // End reactions
    	} // End user defined nuclide
		
	// Push new nuclide
	Nuclides.push_back( Nuc );
    }

  	// Set materials
  	pugi::xml_node input_materials = input_file.child("materials");
  	for ( const auto& m : input_materials )
	{
    		const           std::string name = m.attribute("name").value();
    		std::shared_ptr<Material>  Mat = std::make_shared<Material> ( name );

    		// Add material nuclides
    		for ( const auto& n : m.children() )
		{
			if ( (std::string) n.name() == "nuclide" ) 
			{
        			const std::string                   nuclide_name = n.attribute("name").value();
        			const double                        density      = n.attribute("density").as_double();
				const std::shared_ptr<Nuclide_t>    nucPtr       = find_by_name( Nuclides, nuclide_name );
				
      				if ( nucPtr ) 
				{
        				Mat->addNuclide( nucPtr, density );
      				}
      				else
		       		{
        				std::cout << "unknown nuclide " << nuclide_name << " in material " << name << std::endl;
        				throw;
      				}
      			}
    		}
    		
		// Push new material
		Materials.push_back( Mat );
  	}
  
	// Set surfaces
  	pugi::xml_node input_surfaces = input_file.child("surfaces");
  	for ( const auto& s : input_surfaces )
	{
    		std::shared_ptr<Surface> S;
    		const std::string type = s.name();
      		const std::string name = s.attribute("name").value();
		std::string bc   = "transmission"; // Default boundary condition
		
		// Provided B.C. input
		if ( s.attribute("bc") ) 
		{ 
			bc = s.attribute("bc").value(); 
			if ( bc != "transmission" && bc != "reflective" )
			{
				std::cout<< "unknown boundary condition " << bc << " for surface " << name << std::endl;
				throw;
			}

			if ( bc == "reflective" && type != "plane_x" && type != "plane_y" && type != "plane_z" && type != "plane" )
			{
				std::cout<< "reflective boundary condition is only supported by plane surfaces";
				throw;
			}
		}

                int bc_type;
                if( bc == "transmission" ) { bc_type = 0; }
                if( bc == "reflective" ) { bc_type = 1; }
    		
		// Plane-x
		if ( type == "plane_x" ) 
		{
      			const double x = s.attribute("x").as_double();
			S = std::make_shared< SurfacePlaneX > ( name, Surfaces.size(), bc_type, x );
    		}

		// Plane-y
		else if ( type == "plane_y" ) 
		{
      			const double y = s.attribute("y").as_double();
			S = std::make_shared< SurfacePlaneY > ( name, Surfaces.size(), bc_type, y );
    		}
		
		// Plane-z
		else if ( type == "plane_z" ) 
		{
      			const double z = s.attribute("z").as_double();
			S = std::make_shared< SurfacePlaneZ > ( name, Surfaces.size(), bc_type, z );
    		}

		// Generic plane
		else if ( type == "plane" ) 
		{
      			const double a = s.attribute("a").as_double();
      			const double b = s.attribute("b").as_double();
      			const double c = s.attribute("c").as_double();
      			const double d = s.attribute("d").as_double();
			S = std::make_shared< SurfacePlane > ( name, Surfaces.size(), bc_type, a, b, c, d );
    		}
    		
		// Sphere
		else if ( type == "sphere" )
		{
      			const double x = s.attribute("x").as_double();
      			const double y = s.attribute("y").as_double();
      			const double z = s.attribute("z").as_double();
      			const double r = s.attribute("r").as_double();
      			S = std::make_shared< SurfaceSphere > ( name, Surfaces.size(), bc_type, x, y, z, r );
		}
		
		// Cylinder-x
		else if ( type == "cylinder_x" )
		{
      			const double y = s.attribute("y").as_double();
      			const double z = s.attribute("z").as_double();
      			const double r = s.attribute("r").as_double();
      			S = std::make_shared< SurfaceCylinderX > ( name, Surfaces.size(), bc_type, y, z, r );
		}
		
		// Cylinder-z
		else if ( type == "cylinder_z" )
		{
      			const double x = s.attribute("x").as_double();
      			const double y = s.attribute("y").as_double();
      			const double r = s.attribute("r").as_double();
      			S = std::make_shared< SurfaceCylinderZ > ( name, Surfaces.size(), bc_type, x, y, r );
		}

		// Unknown surface type
		else 
		{
      			std::cout << " unkown surface type " << type << std::endl;
      			throw;
    		}
    	
		// Push new surface
		Surfaces.push_back( S );
  	}
  
	// Set cells
  	pugi::xml_node input_cells = input_file.child("cells");
  	for ( const auto& r : input_cells ) 
	{
    		std::shared_ptr<Cell> Reg;
    		const std::string name       = r.attribute("name").value();
		    double            importance = 1.0;  // default

    		// Modify cell importance
    		if ( r.attribute("importance") ) {
                importance = r.attribute("importance").as_double();
            }
        
    		Reg  = std::make_shared<Cell> ( name, Cells.size(), importance );

    		// Set cell material
    		if ( r.attribute("material") ) 
		{
			const std::shared_ptr<Material> matPtr = find_by_name( Materials, r.attribute("material").value() );
      			if ( matPtr ) 
			{
        			Reg->set_material( matPtr );
      			}
			else
		       	{
        			std::cout << "unknown material " << r.attribute("material").value() << " in cell " << name << std::endl;
        			throw;
      			} 
   		}
   
    		// Set cell bounding surfaces
    		for ( const auto& s : r.children() ) 
		{
      			if ( (std::string) s.name() == "surface" ) 
			{
        			std::string name  = s.attribute("name").value();
        			const int   sense = s.attribute("sense").as_int();

        			std::shared_ptr<Surface> SurfPtr = find_by_name( Surfaces, name );
        			
				if ( SurfPtr ) 
				{
          				Reg->add_surface( find_by_name( Surfaces, name ), sense );
        			}
				else 
				{
          				std::cout << "unknown surface with name " << name << std::endl;
          				throw;
        			}
      			}
			else 
			{
        			std::cout << "unknown data type " << s.name() << " in cell " << name << std::endl;
        			throw;
      			}
    		} 
    		
		// Push new cell
		Cells.push_back( Reg );
  	}
    
//==========================================================================
// Set estimators
//==========================================================================
    
// Loop over estimators
for ( auto& e : input_file.child("estimators").children("estimator") ){
    std::shared_ptr<Estimator> set_estimator;

    // Estimator name
    std::string e_name = "Estimator #" + std::to_string(Estimators.size()+1);
    if ( e.attribute("name") ) { e_name = e.attribute("name").value(); }
    set_estimator = std::make_shared<Estimator> ( e_name, Nsample, 
                                                  Ncycle-Npassive );

    // Estimator type
    std::string e_type = "TL";
    std::shared_ptr<ScoreKernel> e_sk;
    if( e.attribute("type") ) { e_type = e.attribute("type").value(); }
    if( e_type == "TL" ){ 
        e_sk = std::make_shared<ScoreKernelTrackLength>();
    } else if( e_type == "C" ){
        e_sk = std::make_shared<ScoreKernelCollision>();
    } else{
        std::cout<< "[ERROR] Unsupported score type in estimator " 
                 << e_name << "\n";
        std::exit(EXIT_FAILURE);
    }    
    if(tdmc){
        e_sk = std::make_shared<ScoreKernelVelocity>();
    }
    
    // Estimator scores
    std::vector<std::string> e_scores;
    if( !e.attribute("scores") ){
        std::cout<< "[ERROR] There is no score in estimator " << e_name 
                 << "\n";
        std::exit(EXIT_FAILURE);
    }
    std::istringstream iss( e.attribute("scores").value() );
    for( std::string s; iss >> s; ) { e_scores.push_back(s); }
    for( auto& s : e_scores ){
        std::shared_ptr<Score> e_score;
        if( s == "flux" ){
            e_score = std::make_shared<ScoreFlux>(s,e_sk);
        } else if ( s == "absorption" ){ 
            e_score = std::make_shared<ScoreAbsorption>(s,e_sk);
        } else if ( s == "scatter" ){ 
            e_score = std::make_shared<ScoreScatter>(s,e_sk); 
        } else if ( s == "capture" ){ 
            e_score = std::make_shared<ScoreCapture>(s,e_sk);
        } else if ( s == "fission" ){ 
            e_score = std::make_shared<ScoreFission>(s,e_sk); 
        } else if ( s == "nu-fission" ){ 
            e_score = std::make_shared<ScoreNuFission>(s,e_sk);
        } else if ( s == "total" ){ 
            e_score = std::make_shared<ScoreTotal>(s,e_sk);
        } else if ( s == "cross" ){ 
            e_sk = std::make_shared<ScoreKernelNeutron>();
            e_score = std::make_shared<ScoreFlux>(s,e_sk);
        } else{
            std::cout<< "[ERROR] Unsuported score type " << s 
                     << " in estimator " << e_name << "\n";
            std::exit(EXIT_FAILURE);
        }
        set_estimator->add_score( e_score );
    }

    // Estimator filters
    std::vector<double>     f_grid;
    std::shared_ptr<Filter> e_filter;
    // TDMC filter
    if(tdmc){
        f_grid = tdmc_time;
        e_filter = std::make_shared<FilterTDMC> (f_grid);
        set_estimator->add_filter(e_filter);
    }
    // Attach estimator on geometries (and build the corresponding filter)
    f_grid.clear();
    for( auto& surface : e.children("surface") ){
        const std::string s_name = surface.attribute("name").value();
        const std::shared_ptr<Surface> s_ptr = find_by_name( Surfaces, s_name );
        if ( !s_ptr ){
    	std::cout << "[ERROR] Unknown surface label " << s_name 
                      << " in estimator " << e_name << "\n";
            std::exit(EXIT_FAILURE);
   	}
   	s_ptr->attach_estimator( set_estimator );
        f_grid.push_back(s_ptr->ID());
    }
    for( auto& c : e.children("cell") ){
        const std::string c_name = c.attribute("name").value();
        const std::shared_ptr<Cell> c_ptr = find_by_name( Cells, c_name );
        if ( !c_ptr ){
    	std::cout << "[ERROR] Unknown cell label " << c_name 
                      << " in estimator " << e_name << "\n";
            std::exit(EXIT_FAILURE);
   	}
        if( e_type == "TL" ) { c_ptr->attach_estimator_TL( set_estimator ); }
        if( e_type == "C" ) { c_ptr->attach_estimator_C( set_estimator ); }
        f_grid.push_back(c_ptr->ID());
    }
    if( e.child("surface") ){ 
        set_estimator->add_filter( std::make_shared<FilterSurface>(f_grid) );
    } else if( e.child("cell") ){
        set_estimator->add_filter( std::make_shared<FilterCell>(f_grid) );
    } else{
        std::cout<< "[ERROR] Estimator " << e_name 
                 << " needs to be attached somewhere\n";
        std::exit(EXIT_FAILURE);
    }

    // Other filters
    for( auto& f : e.children("filter") ){
        // Filter grid
        f_grid.clear();
	if( f.attribute("grid") ){
	    const std::string   grid_string = f.attribute("grid").value();
	    std::istringstream  iss( grid_string );
	    for( double s; iss >> s; ) { f_grid.push_back(s); }
	} else if( f.attribute("grid_linear") ){
	    const std::string grid_string = f.attribute("grid_linear").value();
	    double a, b, step;
	    std::istringstream  iss( grid_string );
	    iss >> a >> step >> b;
	    f_grid.push_back(a);
	    while( f_grid.back() < b )
	    {
		f_grid.push_back( f_grid.back() + step );
	    }
	    f_grid.pop_back();
	    f_grid.push_back(b);
	} else if( f.attribute("grid_lethargy") ){
	    const std::string grid_string = f.attribute("grid_lethargy").value();
	    double a, b, N, step;
	    std::istringstream  iss( grid_string );
	    iss >> a >> b >> N;
            step = std::log(b/a)/N;
	    f_grid.push_back(0.0);
	    while( f_grid.size()!=N+1 )
	    {
		f_grid.push_back( f_grid.back() + step );
	    }
            std::reverse(f_grid.begin(), f_grid.end());
            for( int i = 0; i < f_grid.size(); i++ ){
                f_grid[i] = b * std::exp(-f_grid[i]);
            }
	} else{
            std::cout<< "[ERROR] Need filter grid for estimator " << e_name 
                     << "\n";
            std::exit(EXIT_FAILURE);
        }
        // Filter type
        if ( !f.attribute("type") ){
            std::cout<< "[ERROR] Need filter type for estimator " << e_name 
                     << "\n";
            std::exit(EXIT_FAILURE);
        }
        std::string f_name = f.attribute("type").value();
        if( f_name == "energy" ){
            e_filter = std::make_shared<FilterEnergy> (f_grid);
        } else if( f_name == "time" ){
            e_filter = std::make_shared<FilterTime> (f_grid);
        } else{
            std::cout<< "[ERROR] Unknown filter type for estimator " << e_name
                     << "\n";
            std::exit(EXIT_FAILURE);
        }
        set_estimator->add_filter(e_filter);
    }
    // Push new estimator
    set_estimator->initialize_tallies();
    Estimators.push_back( set_estimator );
}

//==========================================================================
// TRMM Estimator
//==========================================================================

pugi::xml_node input_trmm = input_file.child("trmm");

if(input_trmm){
    trmm = true;
    if(!ksearch){
        std::cout<< "[ERROR] TRMM should be run in ksearch mode\n" ;
        std::exit(EXIT_FAILURE);
    }
    std::shared_ptr<ScoreKernel> trmm_sk;
    std::shared_ptr<Score>       trmm_score;
    std::vector<double>          trmm_grid;
    std::shared_ptr<Filter>      trmm_filter;

    // Estimator
    trmm_estimator_collision = std::make_shared<Estimator>
        ( "MG", Nsample, Ncycle-Npassive );
    trmm_estimator_scatter = std::make_shared<EstimatorScatter>
        ( "MG_scatter", Nsample, Ncycle-Npassive );
    trmm_estimator_fission = std::make_shared<EstimatorFission>
        ( "MG_fission", Nsample, Ncycle-Npassive );

    // Score
    trmm_sk = std::make_shared<ScoreKernelTrackLengthVelocity>();
    trmm_score = std::make_shared<ScoreTotal>("collision",trmm_sk);
    trmm_estimator_collision->add_score( trmm_score );
    trmm_sk = std::make_shared<ScoreKernelTrackLength>();
    trmm_score = std::make_shared<ScoreFlux>("flux",trmm_sk);
    trmm_estimator_collision->add_score( trmm_score );
    trmm_sk = std::make_shared<ScoreKernelTrackLengthVelocity>();
    trmm_score = std::make_shared<ScoreScatter>("InScatter",trmm_sk);
    trmm_estimator_scatter->add_score( trmm_score );
    trmm_score = std::make_shared<ScoreNuFission>("NuFission",trmm_sk);
    trmm_estimator_fission->add_score( trmm_score );

    // Filters
    for( auto& c : input_trmm.children("cell") ){
        const std::string c_name = c.attribute("name").value();
        const std::shared_ptr<Cell> c_ptr = find_by_name( Cells, c_name );
        if ( !c_ptr ){
            std::cout << "[ERROR] Unknown cell label " << c_name 
                      << " in trmm\n";
            std::exit(EXIT_FAILURE);
   	}
        c_ptr->attach_estimator_TL( trmm_estimator_collision );
        c_ptr->attach_estimator_TL( trmm_estimator_scatter );
        c_ptr->attach_estimator_TL( trmm_estimator_fission );
        trmm_grid.push_back(c_ptr->ID());
    }
    trmm_estimator_collision->add_filter( std::make_shared<FilterCell>(trmm_grid));
    trmm_estimator_scatter->
        add_filter( std::make_shared<FilterCell>(trmm_grid) );
    trmm_estimator_fission->
        add_filter( std::make_shared<FilterCell>(trmm_grid) );
    for( auto& f : input_trmm.children("filter") ){
        // Filter grid
        trmm_grid.clear();
	if( f.attribute("grid") ){
	    const std::string   grid_string = f.attribute("grid").value();
	    std::istringstream  iss( grid_string );
	    for( double s; iss >> s; ){ 
                if( s < trmm_grid.back() ){
                    std::cout << "[ERROR] filter grid value should be ascending\n";
                    std::exit(EXIT_FAILURE);
                }
                trmm_grid.push_back(s); 
            }
	} else if( f.attribute("grid_linear") ){
	    const std::string grid_string = f.attribute("grid_linear").value();
	    double a, b, step;
	    std::istringstream  iss( grid_string );
	    iss >> a >> step >> b;
	    trmm_grid.push_back(a);
	    while( trmm_grid.back() < b )
	    {
		trmm_grid.push_back( trmm_grid.back() + step );
	    }
	    trmm_grid.pop_back();
	    trmm_grid.push_back(b);
	} else if( f.attribute("grid_lethargy") ){
	    const std::string grid_string = f.attribute("grid_lethargy").value();
	    double a, b, N, step;
	    std::istringstream  iss( grid_string );
	    iss >> a >> b >> N;
            step = std::log(b/a)/N;
	    trmm_grid.push_back(0.0);
	    while( trmm_grid.size()!=N+1 )
	    {
		trmm_grid.push_back( trmm_grid.back() + step );
	    }
            std::reverse(trmm_grid.begin(), trmm_grid.end());
            for( int i = 0; i < trmm_grid.size(); i++ ){
                trmm_grid[i] = b * std::exp(-trmm_grid[i]);
            }
	} else{
            std::cout<< "[ERROR] Need filter grid for trmm\n";
            std::exit(EXIT_FAILURE);
        }
        // Filter type
        if ( !f.attribute("type") ){
            std::cout<< "[ERROR] Need filter type for trmm\n";
            std::exit(EXIT_FAILURE);
        }
        std::string f_name = f.attribute("type").value();
        if( f_name == "energy" ){
            trmm_filter = std::make_shared<FilterEnergyOld> (trmm_grid);
            trmm_estimator_scatter->add_filter(trmm_filter);
            trmm_estimator_fission->add_filter(trmm_filter);
            trmm_filter = std::make_shared<FilterEnergy> (trmm_grid);
        } else if( f_name == "time" ){
            trmm_filter = std::make_shared<FilterTime> (trmm_grid);
        } else{
            std::cout<< "[ERROR] Unknown filter type for trmm\n";
            std::exit(EXIT_FAILURE);
        }
        trmm_estimator_collision->add_filter(trmm_filter);
        trmm_estimator_scatter->add_filter(trmm_filter);
        trmm_estimator_fission->add_filter(trmm_filter);
    }
    // Push new estimator
    trmm_estimator_collision->initialize_tallies();
    trmm_estimator_scatter->initialize_tallies();
    trmm_estimator_fission->initialize_tallies();
    Estimators.push_back( trmm_estimator_collision );
    Estimators.push_back( trmm_estimator_scatter );
    Estimators.push_back( trmm_estimator_fission );
}


//==========================================================================
// Source Bank
//==========================================================================

pugi::xml_node input_sources = input_file.child("sources");
for( const auto& s : input_sources.children() ){
    std::shared_ptr<Source> S;
    std::string s_type = s.name();

    // defaults
    double prob = 1.0;
    std::shared_ptr<Distribution<Point> > s_dir = 
        std::make_shared<DistributionIsotropicDirection> ();
    std::shared_ptr<Distribution<double>> s_energy = 
        std::make_shared<DistributionDelta<double>>(2e6);

    // supplied distributions
    std::string s_dir_name;
    std::string s_energy_name;
    if ( s.attribute("probability") ){
        prob = s.attribute("probability").as_double();
    }
    if ( s.attribute("direction") ){
  	s_dir_name = s.attribute("direction").value();
  	s_dir = find_by_name( Distribution_Point, s_dir_name );
    }
    if ( s.attribute("energy") ){
  	s_energy_name = s.attribute("energy").value();
  	s_energy      = find_by_name( Distribution_Double, s_energy_name );
    }
    if ( !s_dir || !s_energy ){
    	std::cout<< "[ERROR] unknown direction distribution in source.\n";
        std::exit(EXIT_FAILURE);
    }
 
    if ( s_type == "point" ){
	const double x = s.attribute("x").as_double();
	const double y = s.attribute("y").as_double();
	const double z = s.attribute("z").as_double();
        Point p(x,y,z);
	S = std::make_shared<SourcePoint>( p, search_cell(p), s_dir, 
                                            s_energy );
    }
    else{
        std::cout << "[INPUT ERROR] Unknown source type...\n";
        std::exit(EXIT_FAILURE);
    }
		
    Fbank.add_source( S, prob );
}   

    
}


//=============================================================================
//=============================================================================
// End of Simulator Set Up
//=============================================================================
//=============================================================================



//=============================================================================
// Move particle
//=============================================================================

void Simulator::move_particle( Particle& P, const double l )
{
    P.move( l );
    if(ksearch) { k_estimator->estimate_TL(P,l); }
    if(tally){ 
        for( auto& e : P.cell()->estimators_TL ) { e->score( P, l ); }
    }
    Ntrack++;
}


//=============================================================================
// Collision
//=============================================================================

void Simulator::collision( Particle& P )
{
    if (ksearch) { k_estimator->estimate_C(P); }
    if(tally){ 
        for( auto& e : P.cell()->estimators_C ) { e->score( P, 0 ); }
    }
    P.cell()->collision( P, Pbank, ksearch, Fbank, k );			
}


//=============================================================================
// Surface Hit
//=============================================================================
        
void Simulator::surface_hit( Particle& P, const std::shared_ptr<Surface>& S )
{
    P.set_surface_old(S);
    if ( S->bc() == 0 ){
	P.move( EPSILON_float );
	P.set_cell( search_cell( P.pos() ) );
    } else{
	S->reflect( P );
        P.set_cell( P.cell() );
	P.move( EPSILON_float );
    }
    if (tally){
	for ( auto& e : S->estimators ){ 
            e->score( P, 0.0 ); 
        }
    }
}


//=============================================================================
// Cut-off and weight rouletting
//=============================================================================

void Simulator::cut_off( Particle& P )
{
    if ( P.energy() <= Ecut_off || P.time() >= tcut_off || P.weight() == 0.0 ){
        P.kill();
    }
    else{
        // Weight rouletting
        if( P.weight() < wr ){
            if( Urand() < P.weight() / ws ) { P.set_weight(ws); }
            else { P.kill(); }
        }
    }
}


//=============================================================================
// Cut-off and weight rouletting
//=============================================================================

bool Simulator::test_point( const Point& p, const std::shared_ptr<Cell>& C )
{
    // Loop over surfaces in cell, if not on correct side return false
    for ( const auto& S : C->surfaces() ) {
    	if ( S.first->eval( p ) * S.second < 0 ) { return false; }
    }
    return true;
}


//=============================================================================
// THE Simulation
//=============================================================================

void Simulator::start()
{
    // Simulation loop
    for ( icycle = 0; icycle < Ncycle ; icycle++ ){
        if ( icycle == Npassive ) { tally = true; }
        Sbank = Fbank; Fbank.reset();

        // Cycle loop
        for ( isample = 0 ; isample < Nsample ; isample++ ){
            Pbank.push( Sbank.get_source() );
                    
            // History loop
            while ( !Pbank.empty() ){
                Particle P = Pbank.top(); // Working particle
                Pbank.pop();

                // Particle loop
                while ( P.alive() ){
                     // To hold nearest surface and its distance
                    std::pair< std::shared_ptr< Surface >, double > SnD;
                                    
                    // Determine nearest surface and its distance
                    SnD = P.cell()->surface_intersect( P );

                    // Determine collision distance
                    double dcol = P.cell()->collision_distance( P.energy() );

                    // Exceeding TDMC time boundary?
                    if(tdmc){
                        double dbound = (tdmc_time[P.tdmc()] - P.time()) 
                                        * P.speed();
                        if( std::min(SnD.second,dcol) > dbound ){
                            move_particle( P, dbound );
                            P.increment_tdmc();
                            cut_off( P );
                            // Time splitting
                            if(P.alive()){
                                P.set_weight(P.weight()/tdmc_split);
                                for( int i = 0; i < tdmc_split - 1; i++ ){
                                    Pbank.push(P);
                                }
                            }
                            continue;
                        }
                    }
                                    
                    // Hit surface?
                    if ( dcol > SnD.second ){	
                        move_particle( P, SnD.second );
                        surface_hit( P, SnD.first );
                        split_roulette( P, Pbank );
                    }
                    // Collide!!
                    else{
                        move_particle( P, dcol );
                        collision( P );
                    }        
                    cut_off( P );
                } 
            }

            // Estimator history closeout
            if (tally)
            { for ( auto& E : Estimators ) { E->end_history(); } }
            if (ksearch) { k_estimator->end_history(); }
        }

        // Estimator cycle closeout
        if (tally) { for ( auto& E : Estimators ) { E->end_cycle(); } }
        if (ksearch){ 
            k_estimator->report_cycle(tally);
            k = k_estimator->k;
        }
    
    } // All cycles are done, end of simulation loop
    for ( auto& E : Estimators ) { E->end_simulation(); }
}


//=============================================================================
// Report results
//=============================================================================

void Simulator::report()
{
    // H5 output
    io_dir += "output.h5";
    H5std_string FILE_NAME(io_dir);
    H5::H5File output(FILE_NAME, H5F_ACC_TRUNC);
    H5::DataSet dataset;
    H5::Group group;
    H5::DataSpace space_scalar(H5S_SCALAR);
    H5::DataType type_ull    = H5::PredType::NATIVE_ULLONG;
    H5::DataType type_double = H5::PredType::NATIVE_DOUBLE;
    H5::StrType type_str(0, H5T_VARIABLE);

    // Summary
    group = output.createGroup("/summary");
    dataset = group.createDataSet( "Ncycle", type_ull, space_scalar );
    dataset.write(&Ncycle, type_ull);
    dataset = group.createDataSet( "Nsample", type_ull, space_scalar );
    dataset.write(&Nsample, type_ull);
    dataset = group.createDataSet( "Npassive",type_ull, space_scalar );
    dataset.write(&Npassive, type_ull);
    dataset = group.createDataSet( "Ntrack",type_ull, space_scalar );
    dataset.write(&Ntrack, type_ull);
    dataset = group.createDataSet( "mode", type_str, space_scalar );
    dataset.write(mode, type_str);
    dataset = group.createDataSet( "cut_off-E", type_double, space_scalar);
    dataset.write(&Ecut_off, type_double);
    dataset = group.createDataSet( "cut_off-t", type_double, space_scalar);
    dataset.write(&tcut_off, type_double);
    group = output.createGroup("/summary/survival_roulette");
    dataset = group.createDataSet( "wr", type_double, space_scalar);
    dataset.write(&wr, type_double);
    dataset = group.createDataSet( "ws", type_double, space_scalar);
    dataset.write(&ws, type_double);
    if(tdmc){
        group = output.createGroup("/summary/tdmc");
        dataset = group.createDataSet( "split", type_ull, space_scalar);
        dataset.write(&tdmc_split, type_ull);
        hsize_t dimsv[1]; dimsv[0] = tdmc_time.size();
        H5::DataSpace data_spacev(1,dimsv);
        dataset = group.createDataSet( "time", type_double, data_spacev);
        dataset.write(tdmc_time.data(), type_double);
    }

    // Report estimators
    for ( auto& E : Estimators ) { E->report( output ); }
    if(ksearch){k_estimator->report(output);}

    if(!trmm){return;}

    // Set TRM
    unsigned long long trm_N = trmm_estimator_collision->tally_size()/2;
    Eigen::MatrixXd TRM,TRM_real;
    TRM = Eigen::MatrixXd(trm_N,trm_N);

    for( int f = 0; f < trm_N; f++ ){
        for( int i = 0; i < trm_N; i++ ){
            if( i == f ){
                TRM(i,i)  = -trmm_estimator_collision->tally(i).mean;
                TRM(i,i) +=  trmm_estimator_scatter->tally(i+i*trm_N).mean;
                TRM(i,i) +=  trmm_estimator_fission->tally(i+i*trm_N).mean;
                TRM(i,i) /=  trmm_estimator_collision->tally(i+trm_N).mean;
            } else{
                TRM(f,i)  = trmm_estimator_scatter->tally(f+i*trm_N).mean;
                TRM(f,i) += trmm_estimator_fission->tally(f+i*trm_N).mean;
                TRM(f,i) /= trmm_estimator_collision->tally(i+trm_N).mean;
            }
        }
    }
    TRM_real = TRM.real();

    // Solve eigenvalue of TRM
    Eigen::MatrixXcd phi_mode;
    Eigen::VectorXcd alpha;
    Eigen::EigenSolver<Eigen::MatrixXd> eSolve(TRM);
    phi_mode = eSolve.eigenvectors();
    alpha    = eSolve.eigenvalues();
    std::vector<double> alpha_real(trm_N);
    std::vector<double> alpha_imag(trm_N);
    for( int i = 0; i < trm_N; i++ ){
        alpha_real[i] = alpha[i].real();
        alpha_imag[i] = alpha[i].imag();
    }

    // Solve coefficients via initial condition
    Eigen::VectorXcd phi0;
    Eigen::VectorXcd A;
    phi0 = Eigen::VectorXcd::Zero(trm_N);
    phi0(trm_N-1) = 1.0 * std::sqrt( 14.1E6 * 191312955.067 ) * 100.0;
    Eigen::ColPivHouseholderQR<Eigen::MatrixXcd> dec(phi_mode);
    A = dec.solve(phi0);

    // Construct solution in time
    std::vector<double> t = {0.0, 3E-8, 15E-8, 4E-6, 1E-4};
    Eigen::MatrixXcd phi = Eigen::MatrixXcd::Zero(t.size(),trm_N);
    std::vector<double> phi_real(trm_N*t.size());

    unsigned long long idx = 0;
    for( int i = 0; i < t.size(); i++){
        for(int g = 0; g < trm_N; g++){
            for(int n = 0; n < trm_N; n++){
                phi(i,g) += A(n) * phi_mode(g,n) * std::exp( alpha[n] * t[i] );
            }
            phi_real[idx] = phi(i,g).real();
            idx++;
        }
    }
    
    // TRMM results
    group = output.createGroup("/TRMM");
    hsize_t dimsM[2]; dimsM[0] = trm_N; dimsM[1] = trm_N;
    H5::DataSpace data_spaceM(2,dimsM);
    dataset = group.createDataSet( "TRM", type_double, data_spaceM);
    dataset.write(TRM_real.data(), type_double);
    hsize_t dims[2]; dims[0] = t.size(); dims[1] = trm_N;
    H5::DataSpace data_spacev(2,dims);
    dataset = group.createDataSet( "flux", type_double, data_spacev);
    dataset.write(phi_real.data(), type_double);
    
    hsize_t dims_alpha[1]; dims_alpha[0] = trm_N;
    H5::DataSpace space_alpha(1,dims_alpha);
    group = group.createGroup("alpha");
    dataset = group.createDataSet( "real", type_double, space_alpha);
    dataset.write(alpha_real.data(), type_double);
    dataset = group.createDataSet( "imag", type_double, space_alpha);
    dataset.write(alpha_imag.data(), type_double);
}
