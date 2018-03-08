#include "XSec.h"
#include "Algorithm.h" 


double XSConstant::xs( const unsigned long long idx, const double E, 
                       const std::vector<double>& E_vec )
{ return val; }

double XSTable::xs( const unsigned long long idx, const double E, 
                    const std::vector<double>& E_vec )
{
    // If it is at the same Energy
    if ( E == E_current ) { ; }

    // If incident energy exceed the table
    //   use the last data corresponding to the highest energy provided
    else if ( idx == E_vec.size() - 1 ){
	E_current = E;
	XS_current = XSdata.back();
    }
    	
    // Similarly for energy below the lowest energy in the table
    else if ( idx == -1 ){
	E_current = E;
	XS_current = XSdata[0];
    }
	
    // Interpolate with the new given bin index
    else{
    	double E1,E2,XS1,XS2,XS;
    	E1  = E_vec.at(idx); E2  = E_vec.at(idx+1);
	XS1 = XSdata[idx];    XS2 = XSdata[idx+1];
    	XS  = interpolate(E,E1,E2,XS1,XS2);
		
	E_current = E;
	XS_current = XS;
    }

    return XS_current;
}
