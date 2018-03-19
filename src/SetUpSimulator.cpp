#include <vector>  
#include <iostream>
#include <cstring> 
#include <memory>  
#include <stack>   
#include <cmath>   
#include <fstream> 

#include "pugixml.hpp"
#include "Simulator.h"


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


//=============================================================================
// Constructor: Simulator Set Up Start
//=============================================================================

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
// Basic parametres
//=============================================================================

pugi::xml_node input_simulation  = input_file.child("simulation");    
pugi::xml_node input_description = input_simulation.child("description");
pugi::xml_node input_ksearch     = input_simulation.child("ksearch");
pugi::xml_node input_tdmc        = input_simulation.child("tdmc");

// Description
simulation_name = input_description.attribute("name").value();
Nsample         = input_description.attribute("samples").as_double();

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
    tdmc_interval.push_back(0.0);
    for( double s; iss >> s; ){
        tdmc_time.push_back(s); 
        tdmc_interval.push_back(s-tdmc_interval.back()); 
    }
    tdmc_interval.erase(tdmc_interval.begin());
    if( input_ksearch ){
        std::cout<<"ksearch and tdmc could not coexist\n";
        std::exit(EXIT_FAILURE);
    }
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

//=============================================================================
// Nuclides
//=============================================================================
pugi::xml_node input_nuclides = input_file.child("nuclides");
for( const auto& n : input_nuclides.children("nuclide") ){
    std::shared_ptr<Nuclide> N;
    const std::string n_name = n.attribute("name").value();
    double            n_A = MAX_float;
    std::shared_ptr<Reaction> n_capture = nullptr;
    std::shared_ptr<Reaction> n_absorb = nullptr;
    std::shared_ptr<Reaction> n_total = nullptr;
    std::shared_ptr<ReactionScatter> n_scatter = nullptr;
    std::shared_ptr<ReactionFission> n_fission = nullptr;
    std::vector<double> n_E;
    std::vector<std::shared_ptr<Distribution<double>>> n_ChiD; n_ChiD.resize(6,nullptr);
	
    // Standard ZAID nuclide
    if ( n.attribute("ZAID") ){
        const std::string n_ZAID = n.attribute("ZAID").value();
        std::string dirName = "./xs_library/" + n_ZAID + ".txt";
        std::ifstream xs_file (dirName);

        // Data loading
        double A;    // Nuclide mass
        std::vector<double> a; // Chi spectrum parameter a
        std::vector<double> b; // Chi spectrum parameter b
        std::vector<double> sigmaS;
        std::vector<double> sigmaC;
        std::vector<double> sigmaF;
        std::vector<double> sigmaA;
        std::vector<double> sigmaT;
        std::vector<double> nu;
        std::vector<double> beta;
        std::vector<double> lambda;
        std::vector<double> fraction;
        std::vector<double> f_lambda;
        double c1, c2, c3, c4, c5, c6;

        // Nuclide mass, 1st line
        if ( xs_file >> c1 ) { n_A = c1; }
        else { std::cout<< "Failed to read A in library file " << dirName << "\n"; throw; }

        // Chi sectrum parameters, 2nd to 4th lines
        for ( int i = 0 ; i < 3 ; i++ ){
            if ( xs_file >> c1 >> c2 ){
                a.push_back( c1 );	
                b.push_back( c2 );	
            }
            else { std::cout<< "Faled to read ab in library file " << dirName  << "\n"; throw; }
        }

        // Cross sections
        // 	Column --> what?
        // 	1 --> energy (eV)
        // 	2 --> sigmaS
        // 	3 --> sigmaC
        // 	4 --> sigmaF
        // 	5 --> nu
        // 	6 --> nu delayed
        while ( xs_file >> c1 >> c2 >> c3 >> c4 >> c5 >> c6 ){
            n_E.push_back(c1);
            sigmaS.push_back(c2);
            sigmaC.push_back(c3);
            sigmaF.push_back(c4);
            sigmaA.push_back(c3+c4);
            sigmaT.push_back(c2+c3+c4);
            nu.push_back(c5);
            if(c6!=0) { beta.push_back(c6/c5); }
            else { beta.push_back(c6); }
        }
        if(sigmaF[0] != 0){
            std::vector<double> c; c.resize(7);
            std::vector<double> d_E;
            std::vector<std::vector<double>> cdf;
            std::vector<double> dummy = {0.0};
            cdf.resize(6,dummy);
            dirName = "./xs_library/" + n_ZAID + "D.txt";
            std::ifstream d_file (dirName);
            d_file >> c[0] >> c[1] >> c[2] >> c[3] >> c[4] >> c[5];
            for( int i = 0; i < 6; i++ ){
                lambda.push_back(c[i]);
            }
            d_file >> c[0] >> c[1] >> c[2] >> c[3] >> c[4] >> c[5];
            for( int i = 0; i < 6; i++ ){
                fraction.push_back(c[i]);
                f_lambda.push_back(lambda[i]*c[i]);
            }
            while( d_file >> c[0] >> c[1] >> c[2] >> c[3] >> c[4] >> c[5] >> c[6] ){
                d_E.push_back(c[0]);
                for( int i = 0; i < 6; i++ ){
                    cdf[i].push_back(c[i+1]);
                }
            }
            for ( int j = 0; j < 6; j++ ){
                for( int i = 1; i < d_E.size() + 1; i++ ){
                    cdf[j][i] = cdf[j][i-1] + cdf[j][i] * (d_E[i]-d_E[i-1]);
                }
            }
            for( int i = 0; i < 6; i++ ){
                n_ChiD[i] = std::make_shared<DistributionDelayedNeutron>(cdf[i],d_E);
            }
        }
            
        auto XS_S = std::make_shared<XSTable> ( sigmaS );
        auto XS_C = std::make_shared<XSTable> ( sigmaC );
        auto XS_F = std::make_shared<XSTable> ( sigmaF );
        auto XS_A = std::make_shared<XSTable> ( sigmaA );
        auto XS_T = std::make_shared<XSTable> ( sigmaT );
        auto XS_nu = std::make_shared<XSTable> ( nu );
        auto XS_beta = std::make_shared<XSTable> ( beta );
            
        // Set reactions
        n_capture = std::make_shared<Reaction>(XS_C);
        n_absorb = std::make_shared<Reaction>(XS_A);
        n_total = std::make_shared<Reaction>(XS_T);
        n_scatter = std::make_shared<ReactionScatter>
            ( XS_S, std::make_shared<DistributionIsotropicScatter>(), n_A );
        n_fission = std::make_shared<ReactionFission>( XS_F, XS_nu,
                std::make_shared<DistributionWatt>( a, b ),
                n_ChiD,
                XS_beta, lambda, fraction, f_lambda );
    }
    // User defined nuclide
    else{
	if ( n.attribute("A") ){
	    n_A = n.attribute("A").as_double();
    	}

    	// Add nuclide reactions
    	for ( const auto& r : n.children() ){
            const std::string rxn_type = r.name();
			
	    // Set XSec
	    std::shared_ptr<XS> XS;
	    if ( r.attribute("xs") ){
      		const double xs = r.attribute("xs").as_double();
		XS = std::make_shared<XSConstant> ( xs );
    	    }
	    else{
	        std::cout << "[ERROR-INPUT] Unknown XS type...\n";
                std::exit(EXIT_FAILURE);
	    }	
	    
            // The Reaction
	    if( rxn_type == "capture" ){
        	n_capture = std::make_shared<Reaction> ( XS );
	    } else{
                std::cout<<"User defined nuclide only support capture now\n";
                std::exit(EXIT_FAILURE);
            }
    	}
		
    }
    N = std::make_shared<Nuclide> ( n_name, n_A, n_capture, n_scatter, 
                                    n_fission, n_absorb, n_total, n_E );
    Nuclides.push_back( N );
}


//=============================================================================
// Materials
//=============================================================================

pugi::xml_node input_materials = input_file.child("materials");
for( const auto& m : input_materials.children("material") ){
    const std::string m_name = m.attribute("name").value();
    std::vector<std::pair<std::shared_ptr<Nuclide>,double>> m_nuclides;

    for( const auto& n : m.children("nuclide") ){
	const std::string n_name = n.attribute("name").value();
	const double n_density = n.attribute("density").as_double();
	const std::shared_ptr<Nuclide> N = find_by_name( Nuclides, n_name );	
	if(N){
	    m_nuclides.push_back( std::make_pair( N, n_density ) );
	}else{
	    std::cout << "[INPUT_ERROR] Unknown nuclide found...\n";
            std::exit(EXIT_FAILURE);
	}
    }
    std::shared_ptr<Material> M =std::make_shared<Material>( m_name, m_nuclides );
    Materials.push_back( M );
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
  

//=============================================================================
// Cells
//=============================================================================

pugi::xml_node input_cells = input_file.child("cells");
for( const auto& c : input_cells.children("cell") ) 
{
    std::shared_ptr<Cell> C;
    std::string c_name = "Cell " + std::to_string(Cells.size() + 1 );
    double c_importance = 1.0;
    std::shared_ptr<Material> c_material = nullptr;
    std::vector<std::pair<std::shared_ptr<Surface>,int>> c_surfaces;

    if( c.attribute("name") ){
        c_name = c.attribute("name").value();
    }
    if( c.attribute("importance") ){
        c_importance = c.attribute("importance").as_double();
    }
    if( c.attribute("material") ){
        const std::string m_name = c.attribute("material").value();
	c_material = find_by_name( Materials,m_name );
        if( !c_material ){
            std::cout << "[INPUT_ERROR] Unknown material in cell\n";
            std::exit(EXIT_FAILURE);
        }
    }
    for ( const auto& s : c.children("surface") ){
        const std::string s_name  = s.attribute("name").value();
        const int s_sense = s.attribute("sense").as_int();
        std::shared_ptr<Surface> S = find_by_name( Surfaces, s_name );        
        if (S){
            c_surfaces.push_back( std::make_pair( S, s_sense ) );
        }
        else {
            std::cout << "[INPUT_ERROR] Unknown surface\n";
            std::exit(EXIT_FAILURE);
        }
    } 
    C = std::make_shared<Cell>( c_name, Cells.size(), c_importance, c_material,
                                c_surfaces );
    Cells.push_back(C);
}
    

//==========================================================================
// Estimators
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
    trmm_score = std::make_shared<ScoreScatterOld>("InScatter",trmm_sk);
    trmm_estimator_scatter->add_score( trmm_score );
    trmm_score = std::make_shared<ScoreNuFissionOld>("NuFission",trmm_sk);
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


