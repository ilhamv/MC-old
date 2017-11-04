#ifndef _XSEC_HEADER_
#define _XSEC_HEADER_


#include <cmath> // sqrt
#include <vector>
#include <memory>
#include "Solver.h"


class XSec_t
{
	public:
     		 XSec_t() {};
    		~XSec_t() {};
		
		// Get cross section in energy
		virtual double xs( const double E, const unsigned long long idx = 0 ) = 0;
};


// Constant cross section
class Constant_XSec : public XSec_t
{
	const double val;

	public:
		 Constant_XSec( const double x ) : val(x) {};
		~Constant_XSec() {};

		double xs( const double E, const unsigned long long idx = 0 ) { return val; }
};


// 1/v cross section
class OverV_XSec : public XSec_t
{
	const double a, b; // a + b / sqrt(E)
	                   // E is in eV

	public:
		 OverV_XSec( const double p1, const double p2 ) : a(p1), b(p2) {};
		~OverV_XSec() {};

		double xs( const double E, const unsigned long long idx = 0 ) { return a + b / std::sqrt(E); }
};


// Table look-up cross section
class Table_XSec : public XSec_t
{
	private:
		std::shared_ptr< std::vector<double> > Edata;
		std::vector<double>                    XSdata;
		double                                 E_current = 0.0;  // Current energy and XS store (to save time in repeating energy binary search)
		double                                 XS_current = 0.0;

	public:
		 Table_XSec( const std::shared_ptr< std::vector<double> >& p1, const std::vector<double> p2 ) : Edata(p1), XSdata(p2) {};
		~Table_XSec() {};

		double xs( const double E, const unsigned long long idx = 0 );
        
};

#endif
