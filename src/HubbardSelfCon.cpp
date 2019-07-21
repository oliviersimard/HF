#include "HubbardSelfCon.h"

/* HubbardSelfCon */

//using namespace HubbardM;

arma::Mat< std::complex<double> > HubbardSelfCon::tmpSelf(HubbardC model, int ii, int ll, arma::Mat< std::complex<double> > SE, std::vector<double> boundArr) throw(){
    std::vector<double> tauVecUp(model._N_tau + 1, 0.0); std::vector<double> tauVecDown(model._N_tau + 1, 0.0);
    std::vector< std::complex<double> > iwnVecUp(model._N_tau, 0.0+im*0.0); std::vector< std::complex<double> > iwnVecDown(model._N_tau, 0.0+im*0.0);
    arma::Mat< std::complex<double> > tmpMatSelf = model.ZEROS_;

    Hubbard HubbardObj; // Declaring default constructor of Hubbard class.
    FFT FFTObj; // Declaring default constructor of FFT class.
    
    functorStruct<Integrals,Hubbard,HubbardC> functObj(&Hubbard::initGk,&Hubbard::Gk); // The functions needed in the iteration procedure.
    
    std::vector<double> kArr;
    for (size_t i=0; i<=model._gridK; i++){
        kArr.push_back(boundArr[0]+(double)i*2.0*boundArr[1]/model._gridK);
    }

    if (ii <= 1){
        std::cout << ii << " tmpSelf" << std::endl;
        
        if (model._dims == 1){

            std::function< arma::Mat< std::complex<double> >(int,int) > functKInit;
            functKInit = [&](int n, int l){
                Integral1D dummy(0.0,0.0+0.0*im);
                for (size_t i = 0; i<kArr.size(); i++){
                    Integral1D int1D(kArr[i],model.callFunctor(model._N_tau,model._beta,n,l));
                    Integrals ints1D(&int1D); Integrals ints1Dnull(&dummy);
                    tmpMatSelf += functObj.callInitFunct(HubbardObj,model,ints1D,ints1Dnull,n,l);
                }
                return tmpMatSelf *= 1.0/model._gridK;
            };
            for (int n = 0; n<model._N_tau; n++){
                iwnVecUp[n] = HubbardObj.frec(functKInit,n,ll)(0,0);
                iwnVecDown[n] = HubbardObj.frec(functKInit,n,ll)(1,1);
            }

            FFTObj.fft_w2t_notc(model, iwnVecUp, tauVecUp);
            FFTObj.fft_w2t_notc(model, iwnVecDown, tauVecDown);

            for (size_t i=0; i<=model._N_tau; i++){
                tauVecUp[i] += -0.5;
                tauVecDown[i] += -0.5;
            }

            if (VERBOSE > 0){
                std::cout << "Down HF density" << -1.0*(tauVecDown.back()) << std::endl;
                std::cout << "Up HF density" << -1.0*(tauVecUp.back()) << std::endl;
            }

            tmpMatSelf(0,0) = -1.0*(tauVecDown.back());
            tmpMatSelf(1,1) = -1.0*(tauVecUp.back());

            return tmpMatSelf;
        }
        else if (model._dims == 2){

            std::cout << "to translate from Julia" << std::endl;
            return tmpMatSelf;
        }
        else{
            throw std::invalid_argument("Check functor in tmpSelf function.");
        }
    }
    else if (ii > 1){
        std::cout << ii << " tmpSelf" << std::endl;

        if (model._dims == 1){

            std::function< arma::Mat< std::complex<double> >(int,int) > functKCon;
            functKCon = [&](int n, int l){
                Integral1D dummy(0.0,0.0+0.0*im); // For extended Hubbard model, but 0 in the Hubbard case, because of local interaction. 
                for (size_t i = 0; i<kArr.size(); i++){
                    Integral1D int1D(kArr[i],model.callFunctor(model._N_tau,model._beta,n,l));
                    Integrals ints1D(&int1D); Integrals ints1Dnull(&dummy);
                    tmpMatSelf += functObj.callConFunct(HubbardObj,model,ints1D,ints1Dnull,n,l,SE);
                }
                return tmpMatSelf *= 1.0/model._gridK;
            };
            for (int n = 0; n<model._N_tau; n++){
                iwnVecUp[n] = HubbardObj.frec(functKCon,n,ll)(0,0);
                iwnVecDown[n] = HubbardObj.frec(functKCon,n,ll)(1,1);
            }
            
            FFTObj.fft_w2t_notc(model, iwnVecUp, tauVecUp);
            FFTObj.fft_w2t_notc(model, iwnVecDown, tauVecDown);

            for (size_t i=0; i<=model._N_tau; i++){
                tauVecUp[i] += -0.5;
                tauVecDown[i] += -0.5;
            }

            if (VERBOSE > 0){
                std::cout << "Down HF density" << -1.0*(tauVecDown.back()) << std::endl;
                std::cout << "Up HF density" << -1.0*(tauVecUp.back()) << std::endl;
            }

            tmpMatSelf(0,0) = -1.0*(tauVecDown.back());
            tmpMatSelf(1,1) = -1.0*(tauVecUp.back());

            return tmpMatSelf;
        }
        else if (model._dims == 2){

            std::cout << "to translate from Julia" << std::endl;
            return tmpMatSelf;
        }
        else{
            throw std::invalid_argument("Check functor in tmpSelf function.");
        }
    }
}

std::unordered_map<int, arma::Mat< std::complex<double> > > HubbardSelfCon::iterationProcess(HubbardC model, std::vector<double> boundArr, int ll) throw(){
    std::unordered_map<int, arma::Mat< std::complex<double> > > umap;
    for (size_t it=1; it<=model._N_it; it++){
        if (it == 1){
            arma::Mat< std::complex<double> > selfCplxMat;
            selfCplxMat = tmpSelf(model, it, ll, model.ZEROS_, boundArr);
            umap[it] = selfCplxMat;
        }
        else if (it > 1){
            arma::Mat< std::complex<double> > selfCplxMat;
            selfCplxMat = tmpSelf(model, it, ll, umap[it-1], boundArr);
            umap[it] = selfCplxMat;
        }
    }
    return umap;
}