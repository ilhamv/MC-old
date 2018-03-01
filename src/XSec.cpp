#include "XSec.h"
#include "Algorithm.h" // interpolate
#include <iostream>

double Table_XSec:: xs( const double E, const unsigned long long bin /*= 0*/ )
{
    // If it is at the same Energy
    if ( E == E_current ) { ; }

    // If incident energy exceed the table
    //   use the last data corresponding to the highest energy provided
    else if ( bin == Edata->size() - 1 )
    {
	E_current = E;
	XS_current = XSdata.back();
    }
    	
    // Similarly for energy below the lowest energy in the table
    else if ( bin == -1 )
    {
	E_current = E;
	XS_current = XSdata[0];
    }
	
    // Interpolate with the new given bin index
    else
    {
    	double E1,E2,XS1,XS2,XS;
    	E1  = Edata->at(bin); E2  = Edata->at(bin+1);
	XS1 = XSdata[bin];    XS2 = XSdata[bin+1];
    	XS  = interpolate(E,E1,E2,XS1,XS2);
		
	E_current = E;
	XS_current = XS;
    }

    return XS_current;
}
