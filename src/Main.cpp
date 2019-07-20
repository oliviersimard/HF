#include "fft.h"

using namespace std;

std::complex<double> z(int N_tau, double beta, int n, int l); // Function used to compute Matsubara frequencies

int main(){
    
    string filename("params.json");
    Json_utils Json_utilsObj;
    HubbardM::HubbardC param_(Json_utilsObj, filename, &z); // Instantiating

    //cout << typeid(U).name() << "  " << N_it << "  " << V << endl;

    if (VERBOSE > 0){
        cout << param_ << endl;
    }

    FFT FFTObj;
    int N = param_._N_tau;

    HubbardM::Integral1D testIntegral1D(0.4,0.1+im*0.3); HubbardM::Integral2D testIntegral2D(0.3,0.4,0.1+im*0.3);
    HubbardM::Integrals testIntegrals1D(&testIntegral1D); HubbardM::Integrals testIntegrals2D(&testIntegral2D);

    HubbardM::Hubbard HubbardObj;
    arma::Mat< complex<double> > matVal = HubbardObj.initGk(param_,testIntegrals1D,testIntegrals1D,4,4);
    cout << "matVal: " << matVal << endl;
    arma::Mat< complex<double> > SE = arma::randu< arma::Mat< complex<double> > >(2,2);
    arma::Mat< complex<double> > matVal2 = HubbardObj.Gk(param_,testIntegrals1D,testIntegrals1D,4,4,SE);
    cout << "matVal2: " << matVal2 << endl;

    return 0;
}

std::complex<double> z(int N_tau, double beta, int n, int l){
    return im*(2.0*((double)n + 1.0 + (double)l*(double)N_tau))*M_PI/beta;
}