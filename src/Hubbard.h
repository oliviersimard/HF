#ifndef Hubbard_H_
#define Hubbard_H_

#include "param.h"

#define VERBOSE = 0

namespace HubbardM{

    struct Integral1D{
        Integral1D(double, std::complex<double>); // Constructor
        ~Integral1D(){}; // Destructor
        Integral1D(const Integral1D& rhs); // Copy constructor
        Integral1D& operator=(const Integral1D& source); // assignment operator
        Integral1D operator+(const Integral1D& rhs) const;
        Integral1D operator-(const Integral1D& rhs) const;

        protected:
            double qx;
            std::complex<double> iwn;
    };

    struct Integral2D : Integral1D{
        Integral2D(double, double, std::complex<double>); // Constructor
        ~Integral2D(){}; // Destructor
        Integral2D(const Integral2D& rhs); // Copy constructor
        Integral2D& operator=(const Integral2D& source); // assignment operator
        Integral2D operator+(const Integral2D& rhs) const;
        Integral2D operator-(const Integral2D& rhs) const;

        protected:
            double qx, qy;
            std::complex<double> iwn;
    };

    class HubbardC : public Param{
        public:
            std::vector< std::complex<double> > _matsubara_grid, _matsubara_grid_boson;

            HubbardC(Json_utils json_utilsObj, std::string filename);
            ~HubbardC(){};
    };

}


#endif /* Hubbard_H_ */
