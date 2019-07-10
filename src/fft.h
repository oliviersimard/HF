#ifndef fft_H_
#define fft_H_

#include "Hubbard.h"

//extern const std::complex<double> im; 

class FFT{
    public:
        void fft_t2w(HubbardM::HubbardC parm_, std::vector<double> &y, std::vector<std::complex<double> > &z);
        void fft_w2t(HubbardM::HubbardC parm_, std::vector<std::complex<double> > &z, std::vector<double> &y);
        template<typename T> std::vector<T> GenerateVecT(unsigned int numOfEls, int min, int max);
};

#ifdef GENERAL
std::ostream& output = std::cout << "\n\n" << "Chose the general Fourier transform option." << std::endl;

#endif /* General */

#endif /* fft_H_ */