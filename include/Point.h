#ifndef POINT_H
#define POINT_H


class Point
{
    public:
        double x, y, z;

        // Constructor: pass the point xyz value
        Point( const double a = 0.0, const double b = 0.0, 
                 const double c = 0.0 ) : x(a), y(b), z(c) {};
        ~Point() {};		
};


#endif // POINT_H
