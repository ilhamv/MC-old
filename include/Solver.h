#ifndef _SOLVER_HEADER_
#define _SOLVER_HEADER_

#include <vector>

#include "Const.h" // MAX
#include "Random.h"
#include "Point.h"


// Quad solver for geometry-point evaluation
// return smallest positive real root if it exists; if it does not, return very big number
double solve_quad( const double a, const double b, const double c );


// Binary search a double location in bin grids
// Note:
// 	value < lowest  grid --> -1
// 	value > highest grid --> vector.size - 1 (or number of bins)
// 	value = grid points  --> location of bin whose upper bound is the value
// 	                         (-1 if value = lowest grid)
int Binary_Search( const double x, const std::vector<double>& vec );


// Scatter direction
// return final direction dir_f after scatter initial direction dir_i with scattering polar cosine mu
Point_t scatter_direction( const Point_t dir_i, const double mu0 );


// Linear interpolation
double Linterpolate( const double x, const double x1, const double x2, const double y1, const double y2 );


// Shannon entropy
class Shannon_Entropy_Mesh
{
    private:
	double xmin, xmax, ymin, ymax, zmin, zmax;
	int    x_nmesh, y_nmesh, z_nmesh, total;

	std::vector<std::vector<std::vector<double>>> mesh;

    public:
	Shannon_Entropy_Mesh( double x, double X, int xn, double y, double Y, int yn, double z, double Z, int zn ) :
            xmin(x), xmax(X), x_nmesh(xn), ymin(y), ymax(Y), y_nmesh(yn), zmin(z), zmax(Z), z_nmesh(zn)
	    {
		std::vector <double> v;
		for ( int k = 0; k < z_nmesh; k++ ) { v.push_back(0.0); }
		std::vector < std::vector <double> > vv;
		for (int j = 0; j < y_nmesh; j++ ) { vv.push_back(v); }
		for (int i = 0; i < x_nmesh; i++ ) { mesh.push_back(vv); }
                total = 0;
	    }
	~Shannon_Entropy_Mesh() {};

	void clear();
	void update( const Point_t& p, const double N );
	double entropy();
};


#endif
