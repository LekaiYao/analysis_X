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
    RooRealVar c1("c1", "coefficient #1", 0.00734315);
	RooRealVar c2("c2", "coefficient #2", -0.318798);
    RooChebychev bkg("bkg", "Chebychev background", B_mass, RooArgList(c1, c2));
    // 假设 B_mass 是在外部定义的 RooRealVar
    B_mass.setRange("fullregion", 3.6, 3.8);
    B_mass.setRange("signalregion", 3.656,3.718);

    // 创建积分对象，指定变量和区间
    RooAbsReal* integral_full = bkg.createIntegral(B_mass, NormSet(B_mass), Range("fullregion"));
    RooAbsReal* integral_signal = bkg.createIntegral(B_mass, NormSet(B_mass), Range("signalregion"));

    // 获取积分值
    double value_full = integral_full->getVal();
    double value_signal = integral_signal->getVal();

    std::cout << "Integral over full region is: " << value_full << std::endl;
    std::cout << "Integral over signal region is: " << value_signal << std::endl;
    std::cout << "fb=(R1+R3)/R2 is: " << value_signal/(value_full-value_signal) << std::endl;
}