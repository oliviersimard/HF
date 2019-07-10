#ifndef Param_H_
#define Param_H_

#include <iostream>
#include <complex>
#include <vector>
#include "json_utils.h"

const std::complex<double> im(0.0,1.0);
const double t_=1.0; const double tp_=-0.3; const double tpp_=0.2;

class Param{
    friend std::ostream& operator<<(std::ostream& os, const Param& param_);
    friend class FFT;
    private:
        double* _db_ptr; int* _integ_ptr; bool* _bool_ptr;
    public:
        double _U, _V, _beta, _temperature, _d_tau, _q_1D;
        int _N_tau, _N_it, _dims, _gridK; // N_tau has same value (meaning) has Niwn in julia program.
        bool _precompK, _twoDMap;
        double* _q2D_ptr;

        Param(Json_utils, std::string);
        ~Param();
};

#endif /* Param_H_ */
