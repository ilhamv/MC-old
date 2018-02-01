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
    const std::string                                        file_name,
    std::string&                                             simName,
    unsigned long long&                                      nSample,          
    bool&                                                    ksearch,
    unsigned long long&                                      nCycle,          
    unsigned long long&                                      nPassive,          
    double&                                                  Ecut_off,
    double&                                                  tcut_off,
    Source_Bank&                                             Sbank,
    std::vector < std::shared_ptr<Surface_t>   >&            Surface,     
    std::vector < std::shared_ptr<Cell_t>      >&            Cell,    
    std::vector < std::shared_ptr<Nuclide_t>   >&            Nuclide,   
    std::vector < std::shared_ptr<Material_t>  >&            Material, 
    std::vector < std::shared_ptr<Estimator_t> >&            Estimator,
    std::vector < std::shared_ptr<Distribution_t<double>> >& Distribution_Double,
    std::vector < std::shared_ptr<Distribution_t<Point_t>>>& Distribution_Point
)
{
    // XML input file
    pugi::xml_document input_file;
    pugi::xml_parse_result load_result = input_file.load_file( file_name.c_str() );

    // Able to load file?
    if ( ! load_result ) 
    {
        std::cout<< "Unable to load input file " << file_name << ":\n";
	std::cout<< load_result.description() << std::endl;
        std::exit(EXIT_FAILURE);
    }
	
    // Set simulation descriptions
    pugi::xml_node input_simulation = input_file.child("simulation");
    	
    for ( const auto& s : input_simulation )
    {
	if( (std::string) s.name() == "description" )
	{
	    simName = s.attribute("name").value();         
	    nSample   = s.attribute("samples").as_double();
	}
	else if ( (std::string) s.name() == "cut-off" )
	{
	    if ( s.attribute("energy") ) { Ecut_off = s.attribute("energy").as_double(); }
	    if ( s.attribute("time") )   { tcut_off = s.attribute("time").as_double(); }
	}
	else if ( (std::string) s.name() == "ksearch" )
	{
            unsigned long long active;
            ksearch  = true;
	    active   = s.attribute("active_cycles").as_double();
	    nPassive = s.attribute("passive_cycles").as_double();
            nCycle   = active + nPassive;
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
			S = std::make_shared< PlaneX_Surface > ( name, bc, x );
    		}

		// Plane-y
		else if ( type == "plane_y" ) 
		{
      			const double y = s.attribute("y").as_double();
			S = std::make_shared< PlaneY_Surface > ( name, bc, y );
    		}
		
		// Plane-z
		else if ( type == "plane_z" ) 
		{
      			const double z = s.attribute("z").as_double();
			S = std::make_shared< PlaneZ_Surface > ( name, bc, z );
    		}

		// Generic plane
		else if ( type == "plane" ) 
		{
      			const double a = s.attribute("a").as_double();
      			const double b = s.attribute("b").as_double();
      			const double c = s.attribute("c").as_double();
      			const double d = s.attribute("d").as_double();
			S = std::make_shared< Plane_Surface > ( name, bc, a, b, c, d );
    		}
    		
		// Sphere
		else if ( type == "sphere" )
		{
      			const double x = s.attribute("x").as_double();
      			const double y = s.attribute("y").as_double();
      			const double z = s.attribute("z").as_double();
      			const double r = s.attribute("r").as_double();
      			S = std::make_shared< Sphere_Surface > ( name, bc, x, y, z, r );
		}
		
		// Cylinder-x
		else if ( type == "cylinder_x" )
		{
      			const double y = s.attribute("y").as_double();
      			const double z = s.attribute("z").as_double();
      			const double r = s.attribute("r").as_double();
      			S = std::make_shared< CylinderX_Surface > ( name, bc, y, z, r );
		}
		
		// Cylinder-z
		else if ( type == "cylinder_z" )
		{
      			const double x = s.attribute("x").as_double();
      			const double y = s.attribute("y").as_double();
      			const double r = s.attribute("r").as_double();
      			S = std::make_shared< CylinderZ_Surface > ( name, bc, x, y, r );
		}

		// Cone-X
		else if ( type == "cone_x" )
		{
      			const double x = s.attribute("x").as_double();
      			const double y = s.attribute("y").as_double();
      			const double z = s.attribute("z").as_double();
      			const double r = s.attribute("r").as_double();
      			S = std::make_shared< ConeX_Surface > ( name, bc, x, y, z, r );
		}
		
		// Cone-Y
		else if ( type == "cone_y" )
		{
      			const double x = s.attribute("x").as_double();
      			const double y = s.attribute("y").as_double();
      			const double z = s.attribute("z").as_double();
      			const double r = s.attribute("r").as_double();
      			S = std::make_shared< ConeX_Surface > ( name, bc, x, y, z, r );
		}

		// Cone-Z
		else if ( type == "cone_z" )
		{
      			const double x = s.attribute("x").as_double();
      			const double y = s.attribute("y").as_double();
      			const double z = s.attribute("z").as_double();
      			const double r = s.attribute("r").as_double();
      			S = std::make_shared< ConeX_Surface > ( name, bc, x, y, z, r );
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
    		std::shared_ptr<Cell_t> Reg;
    		const std::string name       = r.attribute("name").value();
		    double            importance = 1.0;  // default

    		// Modify cell importance
    		if ( r.attribute("importance") ) {
                importance = r.attribute("importance").as_double();
            }
        
    		Reg  = std::make_shared<Cell_t> ( name, importance );

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
		Cell.push_back( Reg );
  	}
  	
	// Set estimators
  	pugi::xml_node input_estimators = input_file.child("estimators");
  	for ( const auto& e : input_estimators )
	{
    		std::shared_ptr<Estimator_t> Est;

    		const std::string e_type     = e.name();
		const std::string name       = e.attribute("name").value();
		
		if ( e_type == "estimator" )
		{
			const std::string        score_string = e.attribute("score").value();
			std::vector<std::string> scores;
			std::istringstream iss( score_string );
			for(std::string s; iss >> s; )
			{ scores.push_back(s); }
		
			// Generic estimator
				Est = std::make_shared<Generic_Estimator> ( name );
				for ( auto& s : scores )
				{
					if      ( s == "current" )    { Est->addScore( std::make_shared<Current_Score>()    ); }
					else if ( s == "flux" )       { Est->addScore( std::make_shared<Flux_Score>()       ); }
					else if ( s == "absorption" ) { Est->addScore( std::make_shared<Absorption_Score>() ); }
					else if ( s == "scatter" )    { Est->addScore( std::make_shared<Scatter_Score>()    ); }
					else if ( s == "capture" )    { Est->addScore( std::make_shared<Capture_Score>()    ); }
					else if ( s == "fission" )    { Est->addScore( std::make_shared<Fission_Score>()    ); }
					else if ( s == "total" )      { Est->addScore( std::make_shared<Total_Score>()      ); }
					else if ( s == "count" )    
					{
        					std::cout << "score 'count' in estimator " << name << " has to be the only score" << std::endl;
						throw;
	       				}
       					else 
					{
        					std::cout << "unsuported score type " << s << " in estimator " << name << std::endl;
						throw;
	       				}
				}
		}

		else if ( e_type == "mgxs" ) 
		{ 
			unsigned int N = 1; // Default Legendre scattering components considered
		       	if ( e.attribute("N") ) { N = e.attribute("N").as_uint(); }
			Est = std::make_shared<MGXS_Estimator> ( name, N ); 
		}
		else { std::cout << "unknown estimator type " << name << std::endl; throw; }

      		for ( const auto& eChild : e.children() )
		{
			// Add estimator to surface
        		if ( (std::string) eChild.name() == "surface" )
			{
          			const std::string                s_name  = eChild.attribute("name").value();
          			const std::shared_ptr<Surface_t> SurfPtr = findByName( Surface, s_name );
          			
				if ( SurfPtr ) 					
				{
       					SurfPtr->addEstimator( Est );
       				}
       				else 
				{
        				std::cout << "unknown surface label " << s_name << " in estimator " << name << std::endl;
					throw;
       				}
       			}
			
			// Add estimator to cell
        		else if ( (std::string) eChild.name() == "cell" )
			{
          			const std::string          r_name  = eChild.attribute("name").value();
          			std::shared_ptr<Cell_t>  RegPtr  = findByName( Cell, r_name );
          				
				if ( RegPtr ) 
				{
            				RegPtr->addEstimator( Est );
          			}
          			else 
				{
            				std::cout << "unknown cell label " << r_name << " in estimator " << name << std::endl;
					throw;
          			}
        		}
            
			
			// Set bin (for generic estimator) or group (for MGXS)
        		else if ( (std::string) eChild.name() == "bin" || (std::string) eChild.name() == "group" )
			{
				// Construct bin grid
				std::vector<double> bin_grid;
				// Just grids
				if( eChild.attribute("grid") )
				{
					const std::string   bin_string = eChild.attribute("grid").value();
					std::istringstream  iss( bin_string );
					for( double s; iss >> s; )
					{ bin_grid.push_back(s); }
				}
				// Linear spaced grid
				if( eChild.attribute("grid_linear") )
				{
					const std::string   bin_string = eChild.attribute("grid_linear").value();
					double a, b, range; // Begin, end, step size
					std::istringstream  iss( bin_string );

					iss >> a >> range >> b;

					bin_grid.push_back(a);
					while( bin_grid.back() < b )
					{
						bin_grid.push_back( bin_grid.back() + range );
					}
					bin_grid.pop_back();
					bin_grid.push_back(b);
				}
				
        			// Generic estimator bin type
				if ( (std::string) eChild.name() == "bin" )
				{
					// Type
					std::string type = eChild.attribute("type").value();
        				if ( type != "energy" && type != "time" ) 
					{
            					std::cout << "unsuported bin type " << type << " in estimator " << name << std::endl;
						throw;
					}
					Est->setBin( type, bin_grid ); 
				}
        			// MGXS group
				else
				{
					if ( bin_grid[0] > 0.0 ) { bin_grid.insert( bin_grid.begin(), 0.0 ); }
					if ( bin_grid.back() > 3.1e7 ) { bin_grid[bin_grid.size()-1] = 3.1e7; }
					Est->setBin( "energy", bin_grid ); 
				}
			}
    		}
    		
		// Push new estimator
    		Estimator.push_back( Est );
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

