#ifndef ENTROPY_H
#define ENTROPY_H

#include <vector>

#include "Point.h"


class Entropy
{
    public:
        Entropy() {};
        ~Entropy() {};

        virtual void score( const Point& p, const int N ) = 0;
        virtual double calculate_H() = 0;
};

class ShannonEntropy : public Entropy
{
    private:
        const std::vector<double> x;
        const std::vector<double> y;
        const std::vector<double> z;
        std::vector<double> p;
        int I, Ix, Iy, Iz, Izy;

    public:
        ShannonEntropy( const std::vector<double>& vx,
                        const std::vector<double>& vy,
                        const std::vector<double>& vz );
        ~ShannonEntropy() {};

        void score( const Point& p, const int N );
        double calculate_H();
};

class EntropyNone : public Entropy
{
    public:
        EntropyNone() {};
        ~EntropyNone() {};

        void score( const Point& p, const int N ) {;}
        double calculate_H() {return 0.0;}
};

#endif // ENTROPY_H
