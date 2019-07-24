#include<iostream>
#include<fstream>
#include<complex>
#include<vector>
#include<armadillo>

using namespace std;

class FunctorBuildGk;

const complex<double> im(0.0,1.0);
const arma::Mat< complex<double> > II_(2, 2, arma::fill::eye);
const arma::Mat< complex<double> > ZEROS_(2, 2, arma::fill::zeros);
static arma::Mat< complex<double> > statMat(2,2);

#ifdef ONED
double epsk1D(double);
arma::Mat< complex<double> > buildGkAA_1D(int,double,int,double,double,double);
#else
double epsk2D(double,double);
arma::Mat< complex<double> > buildGkAA_2D(int,double,int,double,double,double,double);
#endif

complex<double> w(int,double,int);
arma::Mat< complex<double> > someFunction(FunctorBuildGk& Obj, vector<double> kArr, int Nomega);

class FunctorBuildGk{
    friend std::ostream& operator<<(std::ostream& os, const FunctorBuildGk& obj);
    public:
        FunctorBuildGk(double,int,double,double,vector<double>,int,int,vector< complex<double> >&);
        ~FunctorBuildGk()=default;
        
        arma::Mat< complex<double> > operator()(int, double, double);
        void get_ndo_2D();
        arma::Mat< complex<double> > operator()(int, double);
        void get_ndo_1D();

    private:
        double _mu, _u, _ndo;
        int _beta, _Nit, _Nk;
        vector<double> _kArr;
        complex<double>* _Gup_k;
        size_t _size;
};

int main(int argc, char ** argv){
 
	int Nomega=400;
    int Nk=50;
	double beta=5;
    double ndo_initial=0.6;
    int Niterations=100;
    double mu=0.;
    double ndo=ndo_initial;
    vector<complex<double> > Gup_k(Nomega);

    vector<double> kArr(Nk+1);
    for (int k=0; k<=Nk; k++){
        kArr[k] = -1.0*M_PI/2.0 + k*2.0*M_PI/(2.0*Nk);
    }

  
    for (double u=1; u<8; u+=0.2) {
        mu=u/2.;
        ndo=ndo_initial;

        FunctorBuildGk u_ndo_c(mu,beta,u,ndo,kArr,Niterations,Nk,Gup_k);
        //arma::Mat< complex<double> > test;
        cout << "First: " << u_ndo_c << endl;
        u_ndo_c.get_ndo_2D();
        cout << "After: " << u_ndo_c << endl;
        //test = someFunction(u_ndo_c,kArr,Nomega);
    
    for (int i=0; i<Niterations; i++) {
        
        double ndo_av=0;
        for (int kx=0; kx<=Nk; kx++) {

            for (int ky=0; ky<=Nk; ky++){
    
                // calculate Gup_k in Matsubara space (AFM)
                for (int j=0; j<Gup_k.size(); j++)
                    Gup_k[j] = buildGkAA_2D(j,mu,beta,u,ndo,kArr[kx],kArr[ky])(0,0);
                // calculate ndo_k
                double ndo_k=0;
                for (int j=0; j<Gup_k.size(); j++)
                    ndo_k += (2./beta)*(Gup_k[j]-1./w(j,0.0,beta)).real();
                ndo_k -= 0.5;
                ndo_k *= (-1);
            
                if ( (ky==0) || (ky==Nk) || (kx==0) || (kx==Nk) ){
                    if ( ((kx==0) || (kx==Nk)) && ((ky==0) || (ky==Nk)) ){
                        ndo_av += 0.25*ndo_k;
                    }
                    ndo_av += 0.5*ndo_k;
                }
                else
                    ndo_av += ndo_k;
            }
        }
        ndo_av /= (Nk*Nk);
        ndo = ndo_av;
    }
        
    cout << "u= " << u << " ndo: " << ndo << "\n";
        
    }
  
  return 0;
}


complex<double> w(int n, double mu, int beta){
    return im*(2.0*(double)n+1.0)*M_PI/(double)beta + mu;
}

arma::Mat< complex<double> > someFunction(FunctorBuildGk& Obj, vector<double> kArr, int Nomega){
    arma::Mat< complex<double> > objMat = ZEROS_;
    for (size_t i=0; i<kArr.size(); i+=4){
        for (size_t j=0; j<kArr.size(); j+=4){
            for (int n=0; n<Nomega; n+=4){
                objMat += Obj(n,kArr[i],kArr[j]);
            }
        }
    }
    cout << "Big Mat: " << objMat << endl;
    return objMat;
}

FunctorBuildGk::FunctorBuildGk(double mu,int beta,double u,double ndo,vector<double> kArr,int Nit,int Nk,vector< complex<double> >& Gup_k) : 
_mu(mu), _beta(beta), _u(u), _ndo(ndo), _kArr(kArr), _Nit(Nit), _Nk(Nk){
    this->_Gup_k = &Gup_k.front(); 
    this->_size = Gup_k.size();
}


#ifdef ONED

double epsk1D(double kx){
    return -2.0*cos(kx);
}

arma::Mat< complex<double> > buildGkAA_1D(int j, double mu, int beta, double u, double ndo, double kx){
    statMat(0,0) = 1.0/( w(j,mu,beta) - u*ndo - epsk1D(kx)*epsk1D(kx)/( w(j,mu,beta) - u*(1.0-ndo) ) ); // G^{AA}_{up}
    statMat(1,1) = 1.0/( w(j,mu,beta) - u*(1.0-ndo) - epsk1D(kx)*epsk1D(kx)/( w(j,mu,beta) - u*(ndo) ) ); // G^{AA}_{down}
    statMat(0,1) = 0.0+0.0*im; statMat(1,0) = 0.0+0.0*im;
    return statMat;
}

void FunctorBuildGk::get_ndo_1D(){
    for (int i=0; i<_Nit; i++) {
        
        double ndo_av=0.0;
        for (int kkx=0; kkx<=_Nk; kkx++) {
    
            // calculate Gup_k in Matsubara space (AFM)
            for (int jj=0; jj<_size; jj++)
                *(_Gup_k+jj) = buildGkAA_1D(jj,_mu,_beta,_u,_ndo,_kArr[kkx])(0,0);
            // calculate ndo_k
            double ndo_k=0;
            for (int jj=0; jj<_size; jj++)
                ndo_k += (2./_beta)*( *(_Gup_k+jj)-1./w(jj,0.0,_beta) ).real();
            ndo_k -= 0.5;
            ndo_k *= (-1);
            
            if ((kkx==0) || (kkx==_Nk)){
                ndo_av += 0.5*ndo_k;
            }
            else
                ndo_av += ndo_k;
        }
        ndo_av /= (_Nk);
        _ndo = ndo_av;
    }
}

arma::Mat< complex<double> > FunctorBuildGk::operator()(int j, double kx){
    cout << "Tabarnak osti: " << "u= " << _u << " ndo: " << _ndo << "\n";
    return buildGkAA_1D(j,_mu,_beta,_u,_ndo,kx);
}

#else

double epsk2D(double kx, double ky){
    return -2.0*(cos(kx)+cos(ky));
}

arma::Mat< complex<double> > buildGkAA_2D(int j, double mu, int beta, double u, double ndo, double kx, double ky){
    statMat(0,0) = 1.0/( w(j,mu,beta) - u*ndo - epsk2D(kx,ky)*epsk2D(kx,ky)/( w(j,mu,beta) - u*(1.0-ndo) ) ); // G^{AA}_{up}
    statMat(1,1) = 1.0/( w(j,mu,beta) - u*(1.0-ndo) - epsk2D(kx,ky)*epsk2D(kx,ky)/( w(j,mu,beta) - u*(ndo) ) ); // G^{AA}_{down}
    statMat(0,1) = 0.0+0.0*im; statMat(1,0) = 0.0+0.0*im;
    return statMat;
}

void FunctorBuildGk::get_ndo_2D(){
    for (int i=0; i<_Nit; i++) {
        
        double ndo_av=0.0;
        for (int kkx=0; kkx<=_Nk; kkx++) {

            for (int kky=0; kky<=_Nk; kky++){
    
                // calculate Gup_k in Matsubara space (AFM)
                for (int jj=0; jj<_size; jj++)
                    //_Gup_k[jj] = buildGkAA_2D(jj,_mu,_beta,_u,_ndo,_kArr[kkx],_kArr[kky])(0,0);
                    *(_Gup_k+jj) = buildGkAA_2D(jj,_mu,_beta,_u,_ndo,_kArr[kkx],_kArr[kky])(0,0);
                // calculate ndo_k
                double ndo_k=0;
                for (int jj=0; jj<_size; jj++)
                    //ndo_k += (2./_beta)*( _Gup_k[jj]-1./w(jj,0.0,_beta) ).real();
                    ndo_k += (2./_beta)*( *(_Gup_k+jj)-1./w(jj,0.0,_beta) ).real();
                ndo_k -= 0.5;
                ndo_k *= (-1);
            
                if ( (kky==0) || (kky==_Nk) || (kkx==0) || (kkx==_Nk) ){
                    if ( ((kkx==0) || (kkx==_Nk)) && ((kky==0) || (kky==_Nk)) ){
                        ndo_av += 0.25*ndo_k;
                    }
                    ndo_av += 0.5*ndo_k;
                }
                else
                    ndo_av += ndo_k;
            
            }
        }
        ndo_av /= (_Nk*_Nk);
        _ndo = ndo_av;
    }
}

arma::Mat< complex<double> > FunctorBuildGk::operator()(int j, double kx, double ky){
    cout << "Tabarnak osti: " << "u= " << _u << " ndo: " << _ndo << "\n";
    return buildGkAA_2D(j,_mu,_beta,_u,_ndo,kx,ky);
}

#endif

ostream& operator<<(ostream& os, const FunctorBuildGk& obj){
    return os << "The content is: \n" << "U: " << obj._u << endl << "mu: " << obj._mu << endl <<
    "beta: " << obj._beta << endl << "n_do: " << obj._ndo << endl <<
    "size Matsubara arr: " << obj._size << endl << "gridK: " << obj._Nk << endl;
}
