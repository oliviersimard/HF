#ifndef HubbardSelfCon_H_
#define HubbardSelfCon_H_

#include "fft.h"
#include <unordered_map>

using namespace HubbardM;

class HubbardSelfCon{
    public:
        arma::Mat< std::complex<double> > tmpSelf(HubbardC model, int ii, int ll, arma::Mat< std::complex<double> > SE, std::vector<double> boundArr) throw();
        std::unordered_map<int, arma::Mat< std::complex<double> > > iterationProcess(HubbardC model, std::vector<double> boundArr, int ll) throw();
};

#endif /* HubbardSelfCon_H_ */

