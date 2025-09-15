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
    RooRealVar B_mass("B_mass", "B_mass", 3.75, 4.0);
    RooRealVar c1("c1", "coefficient #1", -0.22408);
    RooRealVar c2("c2", "coefficient #2", -0.0398719);
    RooChebychev bkg("bkg", "Chebychev background", B_mass, RooArgList(c1, c2));
    
    // Assume B_mass is an externally defined RooRealVar
    B_mass.setRange("fullregion", 3.75, 4.0);
    //B_mass.setRange("PSIsignalregion", 3.66, 3.72);
    B_mass.setRange("Xsignalregion", 3.83, 3.91);

    // Create integral objects specifying variable and range
    RooAbsReal* integral_full = bkg.createIntegral(B_mass, NormSet(B_mass), Range("fullregion"));
    //RooAbsReal* integral_PSI = bkg.createIntegral(B_mass, NormSet(B_mass), Range("PSIsignalregion"));
    RooAbsReal* integral_X = bkg.createIntegral(B_mass, NormSet(B_mass), Range("Xsignalregion"));

    // Get integral values
    double value_full = integral_full->getVal();
    //double value_PSI = integral_PSI->getVal();
    double value_X = integral_X->getVal();

    std::cout << "Integral over full region is: " << value_full << std::endl;
    //std::cout << "Integral over PSI(2S) signal region is: " << value_PSI << std::endl;
    std::cout << "Integral over X(3872) signal region is: " << value_X << std::endl;
    //std::cout << "fb for PSI(2S) is: " << value_PSI/(value_full - value_PSI - value_X) << std::endl;
    std::cout << "fb for X(3872) is: " << value_X/(value_full - value_X) << std::endl;
}
