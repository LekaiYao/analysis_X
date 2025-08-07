#include <cmath>
#include <iostream>

using namespace std;

void Ratio(){
    double N_psi2S = 44910.0 , N_X = 3296.0 ;
    double sigma_psi2S = 465.0 , sigma_X = 303.0 ;
    double r , sigma_ratio;
    r = N_X/N_psi2S;
    sigma_ratio = r * sqrt(pow((sigma_psi2S/N_psi2S),2)+pow((sigma_X/N_X),2));
    cout << "ratio = N_X/N_psi2S = " << r << endl;
    cout << "uncertainty of ratio is " << sigma_ratio << endl;
}