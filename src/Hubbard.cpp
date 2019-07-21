#include "Hubbard.h"

using namespace HubbardM;

/* Integral 1D */

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


/* Integral 2D */


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

/* Integrals */

Integrals::Integrals(const Integrals& rhs){
    this->_int1D = rhs._int1D;
    this->_int2D = rhs._int2D;
}

Integrals& Integrals::operator=(const Integrals& source){
    if (&source == this) return *this;
    this->_int1D = source._int1D;
    this->_int2D = source._int2D;

    return *this;
}

Integrals Integrals::operator+(const Integrals& source) const throw(){
    if (source._int1D != nullptr && this->_int1D != nullptr){
        Integrals Obj1D(this->_int1D);
        Obj1D._int1D->iwn = this->_int1D->iwn + source._int1D->iwn;
        Obj1D._int1D->qx = this->_int1D->qx + source._int1D->qx;
        return Obj1D; 
    }
    else if (source._int2D != nullptr && this->_int2D != nullptr){
        Integrals Obj2D(this->_int2D);
        Obj2D._int2D->iwn = this->_int2D->iwn + source._int2D->iwn;
        Obj2D._int2D->qx = this->_int2D->qx + source._int2D->qx;
        Obj2D._int2D->qy = this->_int2D->qy + source._int2D->qy;
        return Obj2D;
    }
    else{
        throw std::invalid_argument("Look operator overloading in Integrals struct.");
    }
}

Integrals Integrals::operator-(const Integrals& source) const throw(){
    if (source._int1D != nullptr && this->_int1D != nullptr){
        Integrals Obj1D(this->_int1D);
        Obj1D._int1D->iwn = this->_int1D->iwn - source._int1D->iwn;
        Obj1D._int1D->qx = this->_int1D->qx - source._int1D->qx;
        return Obj1D; 
    }
    else if (source._int2D != nullptr && this->_int2D != nullptr){
        Integrals Obj2D(this->_int2D);
        Obj2D._int2D->iwn = this->_int2D->iwn - source._int2D->iwn;
        Obj2D._int2D->qx = this->_int2D->qx - source._int2D->qx;
        Obj2D._int2D->qy = this->_int2D->qy - source._int2D->qy;
        return Obj2D;
    }
    else{
        throw std::invalid_argument("Look operator overloading in Integrals struct.");
    }
}

/* HubbardC */

HubbardC::HubbardC(Json_utils json_utilsObj, std::string filename, functor_t w) : Param(json_utilsObj, filename){
    
    std::vector< std::complex<double> > matsubara_grid_tmp, matsubara_grid_bosons_tmp;

    for (size_t n=0; n<this->_N_tau; n++){
        double valF = (2.0*(double)n + 1.0)*M_PI/this->_beta;
        double valB = (2.0*(double)n)*M_PI/this->_beta;
        matsubara_grid_tmp.push_back(im*valF);
        matsubara_grid_bosons_tmp.push_back(im*valB);
    }

    this->_matsubara_grid = matsubara_grid_tmp;
    this->_matsubara_grid_boson = matsubara_grid_bosons_tmp;
    this->_w = w;
}

const arma::Mat< std::complex<double> > HubbardC::II_(2, 2, arma::fill::eye);
const arma::Mat< std::complex<double> > HubbardC::ZEROS_(2, 2, arma::fill::zeros);

/* Hubbard */

arma::Mat< std::complex<double> >& Hubbard::swap(arma::Mat< std::complex<double> >& M){
    std::complex<double> buffer = M(0,0);
    M(0,0) = M(1,1);
    M(1,1) = buffer;
    return M;
}

double Hubbard::epsilonk(Integrals kk) throw(){
    if (kk._int1D != nullptr && kk._int2D == nullptr){
        return -2.0*t_*cos(kk._int1D->qx);
    }
    else if (kk._int1D == nullptr && kk._int2D != nullptr){
        return -2.0*t_*(cos(kk._int2D->qx)+cos(kk._int2D->qy));
    }
    else{
        throw std::invalid_argument("Check epsilonk function in Hubbard.cpp");
    }
    
}

arma::Mat< std::complex<double> > Hubbard::initGk(HubbardC model, Integrals kk, Integrals qq, int n, int l) throw(){
    arma::Mat< std::complex<double> > tmpMat = model.II_;
    double mu = model._U/2.0; // Half-filling
    tmpMat(0,0) = 1.0/(model.callFunctor(model._N_tau,model._beta,n,l) - epsilonk(kk+qq) - mu) - 1.0/(model.callFunctor(model._N_tau,model._beta,n,l)); // Diagonal up part
    tmpMat(1,1) = 1.0/(model.callFunctor(model._N_tau,model._beta,n,l) - epsilonk(kk+qq) + mu) - 1.0/(model.callFunctor(model._N_tau,model._beta,n,l)); // Diagonal down part
    tmpMat(0,1) = tmpMat(1,0) = 0.0+im*0.0;

    return tmpMat;
}

arma::Mat< std::complex<double> > Hubbard::Gk(HubbardC model, Integrals kk, Integrals qq, int n, int l, arma::Mat< std::complex<double> > SE) throw(){
    arma::Mat< std::complex<double> > tmpMat = model.ZEROS_;
    double mu = model._U/2.0; // Half-filling
    tmpMat += arma::inv(1.0/model.callFunctor(model._N_tau,model._beta,n,l)*model.II_ + mu*model.II_ - epsilonk(kk+qq)*model.II_ - SE) - 1.0/(model.callFunctor(model._N_tau,model._beta,n,l))*model.II_;
    
    return tmpMat;
}

arma::Mat< std::complex<double> > Hubbard::frec(std::function< arma::Mat< std::complex<double> >(int,int) > funct, int n, int l) throw(){
    this->_stop = l;
    if (this->_stop == 0){
        return funct(n,0);
    }
    else{
        this->_stop -= 1;
        return funct(n,l) + frec(funct,n,l-1);
    }
}
