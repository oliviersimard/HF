#include "Hubbard.h"

using namespace HubbardM;

Integral1D::Integral1D(double qxx, std::complex<double> iwnn) : qx(qxx), iwn(iwnn){};
    
Integral1D::Integral1D(const Integral1D& source) : qx(source.qx), iwn(source.iwn){};

Integral1D& Integral1D::operator=(const Integral1D& source){

    if (&source == this){
        return *this;
    }
    this->iwn = source.iwn;
    this->qx = source.qx;

    return *this;
}

Integral1D Integral1D::operator+(const Integral1D& rhs) const{

    Integral1D Obj(qx,iwn);
    Obj.qx = this->qx + rhs.qx;
    Obj.iwn = this->iwn + rhs.iwn;

    return Obj;    
}

Integral1D Integral1D::operator-(const Integral1D& rhs) const{
    
    Integral1D Obj(qx,iwn);
    Obj.qx = this->qx - rhs.qx;
    Obj.iwn = this->iwn - rhs.iwn;

    return Obj;
}

Integral2D::Integral2D(double qxx, double qyy, std::complex<double> iwnn) : Integral1D(qxx,iwnn){
    this->qy = qyy;
}

Integral2D::Integral2D(const Integral2D& source) : Integral1D(source.qx,source.iwn), qy(source.qy){};

Integral2D& Integral2D::operator=(const Integral2D& source){

    if (&source == this){
        return *this;
    }
    this->qx = source.qx;
    this->qy = source.qy;
    this->iwn = source.iwn;

    return *this;
}

Integral2D Integral2D::operator+(const Integral2D& rhs) const{

    Integral2D Obj(qx,qy,iwn);
    Obj.qx = this->qx + rhs.qx;
    Obj.qy = this->qy + rhs.qy;
    Obj.iwn = this->iwn + rhs.iwn;

    return Obj;    
}

Integral2D Integral2D::operator-(const Integral2D& rhs) const{
    
    Integral2D Obj(qx,qy,iwn);
    Obj.qx = this->qx - rhs.qx;
    Obj.qy = this->qy - rhs.qy;
    Obj.iwn = this->iwn - rhs.iwn;

    return Obj;
}

HubbardC::HubbardC(Json_utils json_utilsObj, std::string filename) : Param(json_utilsObj, filename) {
    
    std::vector< std::complex<double> > matsubara_grid_tmp, matsubara_grid_bosons_tmp;

    for (size_t n=0; n<this->_N_tau; n++){
        double valF = (2.0*(double)n + 1.0)*M_PI/this->_beta;
        double valB = (2.0*(double)n)*M_PI/this->_beta;
        matsubara_grid_tmp.push_back(im*valF);
        matsubara_grid_bosons_tmp.push_back(im*valB);
    }

    this->_matsubara_grid = matsubara_grid_tmp;
    this->_matsubara_grid_boson = matsubara_grid_bosons_tmp;
}