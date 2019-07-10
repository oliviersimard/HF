#include "fft.h"

using namespace std;

template<typename T> vector<T> GenerateVecT(unsigned int numOfEls, double min, double max);
template<> vector< complex<double> > GenerateVecT(unsigned int numOfEls, double min, double max);

int main(){
    
    string filename("params.json");
    Json_utils Json_utilsObj;
    HubbardM::HubbardC param_(Json_utilsObj, filename); // Instantiating

    //cout << typeid(U).name() << "  " << N_it << "  " << V << endl;

    cout << param_ << endl;

    FFT FFTObj;
    vector<double> test_out(30); vector< complex<double> > test_out_inv(30);
    vector< complex<double> > test_in = GenerateVecT< complex<double> >(30,0.5,3.0);
    for (size_t i=0; i<test_in.size(); i++){
        cout << test_in[i] << endl;
    }

    FFTObj.fft_w2t(param_,test_in,test_out);

    for (size_t i=0; i<test_out.size(); i++){
        cout << test_out[i] << endl;
    }

    FFTObj.fft_t2w(param_,test_out,test_out_inv); // Should give back test_in.

    for (size_t i=0; i<test_out_inv.size(); i++){
        cout << test_out_inv[i] << endl;
    }

    return 0;
}

template<typename T> 
vector<T> GenerateVecT(unsigned int numOfEls, double min, double max){
    vector<T> vecValues;
    srand(time(NULL)); // Setting random number generator.
    unsigned int i = 0;
    T randVal = 0;
    while(i < numOfEls){
        double randVal_temp = (double)rand() / RAND_MAX;
        randVal = (T)( min + randVal_temp*(max - min) );
        vecValues.push_back(randVal);
        i++;
    }
    return vecValues;
}

template<>
vector< complex<double> > GenerateVecT(unsigned int numOfEls, double min, double max){
    vector< complex<double> > vecValues;
    srand(time(NULL)); // Setting random number generator.
    unsigned int i = 0;
    complex<double> randVal(0.0,0.0);
    while(i < numOfEls){
        double randVal_temp_im = (double)rand() / RAND_MAX;
        double randVal_temp_re = (double)rand() / RAND_MAX;
        randVal = complex<double>( min + randVal_temp_re*(max - min), min + randVal_temp_im*(max - min) );
        vecValues.push_back(randVal);
        i++;
    }
    return vecValues;
}