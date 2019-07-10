#include "param.h"

Param::Param(Json_utils json_utilsObj, std::string filename){
    MembCarrier contentFile = json_utilsObj.JSONLoading(filename);
    double* dptr = contentFile.db_ptr;
    int* iptr = contentFile.int_prt;
    bool* bptr = contentFile.bool_ptr;
    double* d2ptr = contentFile.db_ptr2;
    //private members
    this->_db_ptr = dptr;
    this->_integ_ptr = iptr;
    this->_bool_ptr = bptr;
    //public members
    this->_q2D_ptr = d2ptr;
    this->_U = _db_ptr[0];
    this->_V = _db_ptr[1];
    this->_beta = _db_ptr[2];
    this->_temperature = 1.0/_beta;
    this->_N_tau = _integ_ptr[0];
    this->_d_tau = _beta/_N_tau;
    this->_q_1D = _db_ptr[3];
    this->_N_it = _integ_ptr[1];
    this->_dims = _integ_ptr[2];
    this->_gridK = _integ_ptr[3];
    this->_precompK = _bool_ptr[0];
    this->_twoDMap = _bool_ptr[1];
}

Param::~Param(){
    std::cout << "param instance destructed" << std::endl;
}

std::ostream& operator<<(std::ostream& os, const Param& param_){
    return os << "The content is: \n" << "U: " << param_._U << std::endl << "V: " << param_._V << std::endl <<
    "beta: " << param_._beta << std::endl << "N_tau: " << param_._N_tau << std::endl <<
    "dims: " << param_._dims << std::endl << "gridK: " << param_._gridK << std::endl;
}