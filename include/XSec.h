#ifndef _XS_H
#define _XS_H

#include <vector>


//=============================================================================
// XS
//=============================================================================

class XS
{
    public:
     	XS() {};
    	~XS() {};
		
	virtual double xs( const unsigned long long idx, const double E, 
                           const std::vector<double>& E_vec ) = 0;
};

//=============================================================================
// Constant
//=============================================================================
class XSConstant : public XS
{
    const double val;

    public:
	 XSConstant( const double x ) : val(x) {};
	~XSConstant() {};

	double xs( const unsigned long long idx, const double E, 
                   const std::vector<double>& E_vec );
};

//=============================================================================
// Table
//=============================================================================
class XSTable : public XS
{
    private:
	const std::vector<double> XSdata;
	double E_current = 0.0; 
	double XS_current = 0.0;

    public:
	XSTable( const std::vector<double> XS ): XSdata(XS) {};
	~XSTable() {};

	double xs( const unsigned long long idx, const double E, 
                   const std::vector<double>& E_vec );
};

#endif // XS_H
