#include <cmath>
#include <iostream>

using namespace std;

void width(){
    double N_gaus1 = 4867.25 , N_gaus2 = 3473.39 ;
    double sigma1 =  0.00999993  , sigma2 =  0.00461782  ;
    double mean = 3.87254;
    double sigma;
    sigma = sqrt((N_gaus1*pow(sigma1,2)+N_gaus2*pow(sigma2,2))/(N_gaus1+N_gaus2));
    cout << "mean = " << mean << endl;
    cout << "sigma = " << sigma << endl;
    cout << "signal region is [" << mean -4*sigma << "," << mean +4*sigma << "]" << endl;
}