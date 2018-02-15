#include <vector>       // vector
#include <iostream>     // cout
#include <cstring>      // strcmp
#include <memory>       // shared_ptr, make_shared
#include <stack>        // stack
#include <cmath>        // exp
#include <sstream>      // istringstream
#include <fstream>      // file stream
#include <dirent.h>     // open a folder

#include "VReduction.h" // Split_Roulette
#include "Const.h"      // MAX
#include "pugixml.hpp"
#include "Geometry.h"
#include "Particle.h"
#include "Distribution.h"
#include "Source.h"
#include "Nuclide.h"
#include "Material.h"
#include "Reaction.h"
#include "Estimator.h"
#include "XSec.h"
#include "XMLparser.h"


// Function that returns an item from a vector of objects of type T by name provided
// the object class has a string and a method called name() allowing for it to be returned
template< typename T >
std::shared_ptr< T > findByName( const std::vector< std::shared_ptr< T > >& vec, const std::string name ) 
{
	for ( auto& v : vec ) 
	{
		if ( v->name() == name ) { return v; }
	}
	return nullptr;
}


// Function that set nuclide based on library availability
// Current approximation: isotropic in C.O.M.
void setNuclide( const std::string name, const std::string label, std::shared_ptr<Nuclide_t>& Nuc )
{
	// Check availability
	std::string dirName = "./xs_library/" + name + ".txt"; // Library file location
	std::ifstream xs_file (dirName);
	if ( !xs_file )
	{
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
	Nuc->addReaction( std::make_shared< Scatter_Reaction > ( XS_S, std::make_shared< IsotropicScatter_Distribution > (), A ) );
	// Fissionable check
	if ( nu[ nu.size() / 2 ] > 0.0 )
	{
		Nuc->addReaction( std::make_shared< Fission_Reaction > ( XS_F, XS_nu, std::make_shared< Watt_Distribution > ( a, b ) ) );
	}
}


// XML input pasrese
void XML_input
(
    std::string                                              file_dir,
    std::string&                                             simulation_name,
    unsigned long long&                                      Nsample,          
    bool&                                                    ksearch,
    bool&                                                    tdmc,
    unsigned long long&                                      Ncycle,          
    unsigned long long&                                      Npassive,          
    double&                                                  Ecut_off,
    double&                                                  tcut_off,
    Source_Bank&                                             Sbank,
    std::vector < std::shared_ptr<Surface_t>   >&            Surface,     
    std::vector < std::shared_ptr<Cell>      >&              cell,    
    std::vector < std::shared_ptr<Nuclide_t>   >&            Nuclide,   
    std::vector < std::shared_ptr<Material_t>  >&            Material, 
    std::vector < std::shared_ptr<Estimator> >&            estimator,
    std::vector < std::shared_ptr<Distribution_t<double>> >& Distribution_Double,
    std::vector < std::shared_ptr<Distribution_t<Point_t>>>& Distribution_Point,
    std::vector<double>& tdmc_time,
    int& tdmc_split
)
{
    // XML input file
    file_dir += "input.xml";
    pugi::xml_document input_file;
    pugi::xml_parse_result load_result = input_file.load_file( file_dir.c_str() );

    // Able to load file?
    if ( ! load_result ) 
    {
        std::cout<< "Unable to load input file " << file_dir << ":\n";
	std::cout<< load_result.description() << std::endl;
        std::exit(EXIT_FAILURE);
    }

    
//==============================================================================
// Set simulation
//==============================================================================

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
}

// TDMC
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
        tdmc_split = input_tdmc.attribute("split").as_int();
    }
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
        			if ( findByName( Distribution_Double, name ) ) { continue; }
	        		std::shared_ptr<Distribution_t<double>> Dist;
        			
				// Delta-double
				if ( type == "delta" ) 
				{
          				const double val = d.attribute("val").as_double();
          				Dist = std::make_shared< Delta_Distribution< double > > ( val, name );
	        		}

        			// Uniform-double
				else if ( type == "uniform" ) 
				{
          				const double a = d.attribute("a").as_double();
          				const double b = d.attribute("b").as_double();
	          			Dist = std::make_shared< Uniform_Distribution > ( a, b, name );
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
					Dist = std::make_shared< Watt_Distribution > ( a, b, name );
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
        			if ( findByName( Distribution_Point, name ) ) { continue; }
        			std::shared_ptr< Distribution_t< Point_t > > Dist;
		        	
				// Delta-point
				if ( type == "delta" ) 
				{
	          			const double x = d.attribute("x").as_double(); 
	          			const double y = d.attribute("y").as_double(); 
        	  			const double z = d.attribute("z").as_double();         
	          			Dist = std::make_shared< Delta_Distribution< Point_t > > ( Point_t( x, y, z ), name );
	        		}
        			
				// Isotropic-point
				else if ( type == "isotropic" ) 
				{
          				Dist = std::make_shared< IsotropicDirection_Distribution > ( name );
        			}
        			
				// XYZ-point
        			else if ( type == "independentXYZ" ) 
				{
		          		std::shared_ptr< Distribution_t<double> > distX = findByName( Distribution_Double, d.attribute("x").value() ); 
        		  		std::shared_ptr< Distribution_t<double> > distY = findByName( Distribution_Double, d.attribute("y").value() ); 
          				std::shared_ptr< Distribution_t<double> > distZ = findByName( Distribution_Double, d.attribute("z").value() ); 

          				// if any of these distributions have not yet been resolved, skip to the end of the loop
	          			if ( !distX || !distY || !distZ ) { continue; }

		          		Dist = std::make_shared< IndependentXYZ_Distribution > ( distX, distY, distZ, name );
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
	    setNuclide( n.attribute("ZAID").value(), name, Nuc );
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
  		    std::shared_ptr< Distribution_t<double> > f_mu = findByName( Distribution_Double, dist_name );

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
  		    std::shared_ptr< Distribution_t<double> > watt = findByName( Distribution_Double, dist_name );
                
		    Nuc->addReaction( std::make_shared< Fission_Reaction > ( XS, nubar, watt ) );
		} 
	    } // End reactions
    	} // End user defined nuclide
		
	// Push new nuclide
	Nuclide.push_back( Nuc );
    }

  	// Set materials
  	pugi::xml_node input_materials = input_file.child("materials");
  	for ( const auto& m : input_materials )
	{
    		const           std::string name = m.attribute("name").value();
    		std::shared_ptr<Material_t>  Mat = std::make_shared<Material_t> ( name );

    		// Add material nuclides
    		for ( const auto& n : m.children() )
		{
			if ( (std::string) n.name() == "nuclide" ) 
			{
        			const std::string                   nuclide_name = n.attribute("name").value();
        			const double                        density      = n.attribute("density").as_double();
				const std::shared_ptr<Nuclide_t>    nucPtr       = findByName( Nuclide, nuclide_name );
				
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
		Material.push_back( Mat );
  	}
  
	// Set surfaces
  	pugi::xml_node input_surfaces = input_file.child("surfaces");
  	for ( const auto& s : input_surfaces )
	{
    		std::shared_ptr< Surface_t > S;
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
    		
		// Plane-x
		if ( type == "plane_x" ) 
		{
      			const double x = s.attribute("x").as_double();
			S = std::make_shared< PlaneX_Surface > ( name, Surface.size(), bc, x );
    		}

		// Plane-y
		else if ( type == "plane_y" ) 
		{
      			const double y = s.attribute("y").as_double();
			S = std::make_shared< PlaneY_Surface > ( name, Surface.size(), bc, y );
    		}
		
		// Plane-z
		else if ( type == "plane_z" ) 
		{
      			const double z = s.attribute("z").as_double();
			S = std::make_shared< PlaneZ_Surface > ( name, Surface.size(), bc, z );
    		}

		// Generic plane
		else if ( type == "plane" ) 
		{
      			const double a = s.attribute("a").as_double();
      			const double b = s.attribute("b").as_double();
      			const double c = s.attribute("c").as_double();
      			const double d = s.attribute("d").as_double();
			S = std::make_shared< Plane_Surface > ( name, Surface.size(), bc, a, b, c, d );
    		}
    		
		// Sphere
		else if ( type == "sphere" )
		{
      			const double x = s.attribute("x").as_double();
      			const double y = s.attribute("y").as_double();
      			const double z = s.attribute("z").as_double();
      			const double r = s.attribute("r").as_double();
      			S = std::make_shared< Sphere_Surface > ( name, Surface.size(), bc, x, y, z, r );
		}
		
		// Cylinder-x
		else if ( type == "cylinder_x" )
		{
      			const double y = s.attribute("y").as_double();
      			const double z = s.attribute("z").as_double();
      			const double r = s.attribute("r").as_double();
      			S = std::make_shared< CylinderX_Surface > ( name, Surface.size(), bc, y, z, r );
		}
		
		// Cylinder-z
		else if ( type == "cylinder_z" )
		{
      			const double x = s.attribute("x").as_double();
      			const double y = s.attribute("y").as_double();
      			const double r = s.attribute("r").as_double();
      			S = std::make_shared< CylinderZ_Surface > ( name, Surface.size(), bc, x, y, r );
		}

		// Cone-X
		else if ( type == "cone_x" )
		{
      			const double x = s.attribute("x").as_double();
      			const double y = s.attribute("y").as_double();
      			const double z = s.attribute("z").as_double();
      			const double r = s.attribute("r").as_double();
      			S = std::make_shared< ConeX_Surface > ( name, Surface.size(), bc, x, y, z, r );
		}
		
		// Cone-Y
		else if ( type == "cone_y" )
		{
      			const double x = s.attribute("x").as_double();
      			const double y = s.attribute("y").as_double();
      			const double z = s.attribute("z").as_double();
      			const double r = s.attribute("r").as_double();
      			S = std::make_shared< ConeX_Surface > ( name, Surface.size(), bc, x, y, z, r );
		}

		// Cone-Z
		else if ( type == "cone_z" )
		{
      			const double x = s.attribute("x").as_double();
      			const double y = s.attribute("y").as_double();
      			const double z = s.attribute("z").as_double();
      			const double r = s.attribute("r").as_double();
      			S = std::make_shared< ConeX_Surface > ( name, Surface.size(), bc, x, y, z, r );
		}
		
		// Unknown surface type
		else 
		{
      			std::cout << " unkown surface type " << type << std::endl;
      			throw;
    		}
    	
		// Push new surface
		Surface.push_back( S );
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
        
    		Reg  = std::make_shared<Cell> ( name, cell.size(), importance );

    		// Set cell material
    		if ( r.attribute("material") ) 
		{
			const std::shared_ptr<Material_t> matPtr = findByName( Material, r.attribute("material").value() );
      			if ( matPtr ) 
			{
        			Reg->setMaterial( matPtr );
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

        			std::shared_ptr<Surface_t> SurfPtr = findByName( Surface, name );
        			
				if ( SurfPtr ) 
				{
          				Reg->addSurface( findByName( Surface, name ), sense );
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
		cell.push_back( Reg );
  	}
    
//==========================================================================
// Set estimators
//==========================================================================
    
// Loop over estimators
for ( auto& e : input_file.child("estimators").children("estimator") ){
    std::shared_ptr<Estimator> set_estimator;

    // Estimator name
    std::string e_name = "Estimator #" + std::to_string(estimator.size()+1);
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
        const std::shared_ptr<Surface_t> s_ptr = findByName( Surface, s_name );
        if ( !s_ptr ){
    	std::cout << "[ERROR] Unknown surface label " << s_name 
                      << " in estimator " << e_name << "\n";
            std::exit(EXIT_FAILURE);
   	}
   	s_ptr->attach_estimator_C( set_estimator );
        f_grid.push_back(s_ptr->ID());
    }
    for( auto& c : e.children("cell") ){
        const std::string c_name = c.attribute("name").value();
        const std::shared_ptr<Cell> c_ptr = findByName( cell, c_name );
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
    estimator.push_back( set_estimator );
}

//==========================================================================
// TRMM
//==========================================================================

pugi::xml_node input_trmm = input_file.child("trmm");

if(input_trmm){
    if(!ksearch){
        std::cout<< "[ERROR] TRMM should be run in ksearch mode\n" ;
        std::exit(EXIT_FAILURE);
    }
    std::shared_ptr<Estimator>   trmm_estimator;
    std::shared_ptr<ScoreKernel> trmm_sk;
    std::shared_ptr<Score>       trmm_score;
    std::vector<double>          trmm_grid;
    std::shared_ptr<Filter>      trmm_filter;

    // Estimator
    trmm_estimator = std::make_shared<Estimator>
        ( "TRMM[g]", Nsample, Ncycle-Npassive );

    // Score
    trmm_sk = std::make_shared<ScoreKernelTrackLengthVelocity>();
    trmm_score = std::make_shared<ScoreTotal>("Collision",trmm_sk);
    trmm_estimator->add_score( trmm_score );
    trmm_sk = std::make_shared<ScoreKernelTrackLength>();
    trmm_score = std::make_shared<ScoreFlux>("Flux",trmm_sk);
    trmm_estimator->add_score( trmm_score );

    // Filters
    for( auto& c : input_trmm.children("cell") ){
        const std::string c_name = c.attribute("name").value();
        const std::shared_ptr<Cell> c_ptr = findByName( cell, c_name );
        if ( !c_ptr ){
            std::cout << "[ERROR] Unknown cell label " << c_name 
                      << " in trmm\n";
            std::exit(EXIT_FAILURE);
   	}
        c_ptr->attach_estimator_TL( trmm_estimator );
        trmm_grid.push_back(c_ptr->ID());
    }
    trmm_estimator->add_filter( std::make_shared<FilterCell>(trmm_grid) );
    for( auto& f : input_trmm.children("filter") ){
        // Filter grid
        trmm_grid.clear();
	if( f.attribute("grid") ){
	    const std::string   grid_string = f.attribute("grid").value();
	    std::istringstream  iss( grid_string );
	    for( double s; iss >> s; ) { trmm_grid.push_back(s); }
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
            trmm_filter = std::make_shared<FilterEnergy> (trmm_grid);
        } else if( f_name == "time" ){
            trmm_filter = std::make_shared<FilterTime> (trmm_grid);
        } else{
            std::cout<< "[ERROR] Unknown filter type for trmm\n";
            std::exit(EXIT_FAILURE);
        }
        trmm_estimator->add_filter(trmm_filter);
    }
    // Push new estimator
    trmm_estimator->initialize_tallies();
    estimator.push_back( trmm_estimator );
}

        // Set source bank
  	pugi::xml_node input_sources = input_file.child("sources");
  	for ( const auto& s : input_sources )
	{
		std::shared_ptr<Source_t> Src;
		
		// Default parameters
		double prob = 1.0;                                                                                            // probability or ratio
  		std::shared_ptr< Distribution_t<Point_t> > dirDist  = std::make_shared< IsotropicDirection_Distribution > (); // direction distribution
  		std::string dir_dist_name;
  		std::shared_ptr< Distribution_t<double> >  enrgDist = std::make_shared< Delta_Distribution<double> > ( 2e6 ); // energy distribution
  		std::string enrg_dist_name;
  		std::shared_ptr< Distribution_t<double> >  timeDist = std::make_shared< Delta_Distribution<double> > ( 0.0 ); // time distribution
  		std::string time_dist_name;
    		
		// Input provided probability
    		if ( s.attribute("probability") )
    		{
	    		prob = s.attribute("probability").as_double();
    		}

		// Input provided direction distribution
    		if ( s.attribute("direction") )
    		{
  			dir_dist_name = s.attribute("direction").value();
  			dirDist       = findByName( Distribution_Point, dir_dist_name );
    		}
    		
		// Input provided energy distribution
		if ( s.attribute("energy") )
    		{
  			enrg_dist_name = s.attribute("energy").value();
  			enrgDist       = findByName( Distribution_Double, enrg_dist_name );
    		}

		// Input provided time distribution
    		if ( s.attribute("time") )
    		{
  			time_dist_name = s.attribute("time").value();
  			timeDist       = findByName( Distribution_Double, time_dist_name );
    		}
		
		// Check distribution availability
            	if ( ! dirDist || ! enrgDist || ! timeDist ) 
		{
    			if ( ! dirDist )  { std::cout << " unknown direction distribution " << dir_dist_name  << " in source " << std::endl; }
    			if ( ! enrgDist ) { std::cout << " unknown energy distribution "    << enrg_dist_name << " in source " << std::endl; }
    			if ( ! timeDist ) { std::cout << " unknown time distribution "      << time_dist_name << " in source " << std::endl; }
    			throw;
 		}

		// Point source
		if ( (std::string) s.name() == "point" )
		{
			const double x = s.attribute("x").as_double();
			const double y = s.attribute("y").as_double();
			const double z = s.attribute("z").as_double();
			Src = std::make_shared<Point_Source> ( x, y, z, dirDist, enrgDist, timeDist );
    		}
		
		// Spherical shell source
		else if ( (std::string) s.name() == "sphere_shell" )
		{
			const double x  = s.attribute("x").as_double();
			const double y  = s.attribute("y").as_double();
			const double z  = s.attribute("z").as_double();
			const double ri = s.attribute("ri").as_double();
			const double ro = s.attribute("ro").as_double();
			Src = std::make_shared<Sphere_Shell_Source> ( x, y, z, ri, ro, dirDist, enrgDist, timeDist );
    		}

		// Disk-x
		else if ( (std::string) s.name() == "disk_x" )
		{
			const double x     = s.attribute("x").as_double();
			const double y     = s.attribute("y").as_double();
			const double z     = s.attribute("z").as_double();
			const double r     = s.attribute("r").as_double();
			Src = std::make_shared<DiskX_Source> ( x, y, z, r, dirDist, enrgDist, timeDist );
    		}

		// Disk-z
		else if ( (std::string) s.name() == "disk_z" )
		{
			const double x     = s.attribute("x").as_double();
			const double y     = s.attribute("y").as_double();
			const double z     = s.attribute("z").as_double();
			const double r     = s.attribute("r").as_double();
			Src = std::make_shared<DiskZ_Source> ( x, y, z, r, dirDist, enrgDist, timeDist );
    		}
		
		// Generic source
		else if ( (std::string) s.name() == "source" )
		{
		
  			std::string pos_dist_name  = s.attribute("position").value();

  			std::shared_ptr< Distribution_t<Point_t> > posDist  = findByName( Distribution_Point , pos_dist_name );

  			if ( ! posDist )
			{ 
				std::cout << " unknown position distribution "  << pos_dist_name  << " in source " << std::endl;
    				throw;
  			}
			
			Src = std::make_shared<Generic_Source> ( posDist, dirDist, enrgDist, timeDist );
		}

		// Unknown source type
		else 
		{
            		std::cout << "unknown source type " << (std::string) s.name() << std::endl;
			throw;
          	}

		// Push new source
		Sbank.addSource( Src, prob );
  	}
    
}


