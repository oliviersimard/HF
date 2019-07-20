#include "integral_utils.h"

template<typename T, typename C, typename Q>
functorStruct<T,C,Q>::functorStruct(funct_init_t initFunct, funct_con_t conFunct) : _initFunct(initFunct), _conFunct(conFunct){};

template<typename T, typename C, typename Q>
arma::Mat< std::complex<double> > functorStruct<T,C,Q>::callInitFunct(C& obj, Q model, T kk, T qq, int n, int l){
    return (obj.*_initFunct)(model, kk, qq, n, l);
}

template<typename T, typename C, typename Q>
arma::Mat< std::complex<double> > functorStruct<T,C,Q>::callConFunct(C& obj, Q model, T kk, T qq, int n, int l, arma::Mat< std::complex<double> > SE){
    return (obj.*_conFunct)(model, kk, qq, n, l, SE);
}