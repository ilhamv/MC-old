#ifndef _POINT_HEADER_
#define _POINT_HEADER_


// Point in 3d space
class Point_t
{
    public:
        double x, y, z;

        // Constructor: pass the point xyz value
        Point_t( const double a = 0.0, const double b = 0.0, const double c = 0.0 ) : x(a), y(b), z(c) {};
        ~Point_t() {};		

        void normalize();
};


#endif
