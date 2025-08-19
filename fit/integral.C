#include "RooAbsPdf.h"
#include "RooExtendPdf.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooRealVar.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooFit.h"  
#include <iostream>

using namespace RooFit;

void integral(){
    RooRealVar B_mass("B_mass", "B_mass", 3.6, 3.8);
    RooRealVar c1("c1", "coefficient #1", 0.131057);
    RooRealVar c2("c2", "coefficient #2", -0.0613971);
    RooChebychev bkg("bkg", "Chebychev background", B_mass, RooArgList(c1, c2));
    
    // Assume B_mass is an externally defined RooRealVar
    B_mass.setRange("fullregion", 3.6, 3.8);
    B_mass.setRange("signalregion", 3.65, 3.72);

    // Create integral objects specifying variable and range
    RooAbsReal* integral_full = bkg.createIntegral(B_mass, NormSet(B_mass), Range("fullregion"));
    RooAbsReal* integral_signal = bkg.createIntegral(B_mass, NormSet(B_mass), Range("signalregion"));

    // Get integral values
    double value_full = integral_full->getVal();
    double value_signal = integral_signal->getVal();

    std::cout << "Integral over full region is: " << value_full << std::endl;
    std::cout << "Integral over signal region is: " << value_signal << std::endl;
    std::cout << "fb=(R1+R3)/R2 is: " << value_signal/(value_full - value_signal) << std::endl;
}
