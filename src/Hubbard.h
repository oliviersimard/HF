#ifndef Hubbard_H_
#define Hubbard_H_

#include "param.h"
#include <exception>

namespace HubbardM{

    struct Integral1D{
        Integral1D(double, std::complex<double>); // Constructor
        ~Integral1D(){}; // Destructor
        Integral1D(const Integral1D& rhs); // Copy constructor
        Integral1D& operator=(const Integral1D& source); // assignment operator
        Integral1D operator+(const Integral1D& rhs) const;
        Integral1D operator-(const Integral1D& rhs) const;

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

        double qx, qy;
        std::complex<double> iwn;
    };

    struct Integrals {
        Integral1D* _int1D; 
        Integral2D* _int2D;
        Integrals(Integral1D* int1D){
            this->_int1D = int1D;
            this->_int2D = nullptr;
        }
        Integrals(Integral2D* int2D){
            this->_int1D = nullptr;
            this->_int2D = int2D;
        }
        ~Integrals()=default;
        Integrals(const Integrals& rhs);
        Integrals& operator=(const Integrals& source);
        Integrals operator+(const Integrals& source) const throw();
        Integrals operator-(const Integrals& source) const throw();
    };

    class HubbardC : public Param{
        public:
            using functor_t = std::complex<double> (*)(int, double, int, int);
            std::vector< std::complex<double> > _matsubara_grid, _matsubara_grid_boson;
            static const arma::Mat< std::complex<double> > II_;
            static const arma::Mat< std::complex<double> > ZEROS_;
            HubbardC(Json_utils json_utilsObj, std::string filename, functor_t w);
            ~HubbardC(){};
            std::complex<double> callFunctor(int N_tau, double beta, int n, int l){
                return this->_w(N_tau, beta, n, l);
            }
        private:
            functor_t _w;
    };

    class Hubbard{
        public:
            Hubbard() = default;
            ~Hubbard() = default;

            double epsilonk(Integrals kk) throw();
            arma::Mat< std::complex<double> >& swap(arma::Mat< std::complex<double> >& M);
            arma::Mat< std::complex<double> > initGk(HubbardC model, Integrals kk, Integrals qq, int n, int l) throw();
            arma::Mat< std::complex<double> > Gk(HubbardC model, Integrals kk, Integrals qq, int n, int l, arma::Mat< std::complex<double> > SE) throw();
            arma::Mat< std::complex<double> > frec(arma::Mat< std::complex<double> > (*funct)(int,int), int n, int l) throw();
        private:
            int _stop = 0;
    };

    template<typename T, typename C, typename Q>
    class HubbardSelfCon{
        public:
            arma::Mat< std::complex<double> > tmpSelf(functorStruct<T,C,Q> functObj, HubbardC model, int ii, int ll, arma::Mat< std::complex<double> > SE) throw();
    };

}


#endif /* Hubbard_H_ */
