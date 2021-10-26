#ifndef SWAPTION_G2_H_INCLUDED
#define SWAPTION_G2_H_INCLUDED

#include "product.h"
#include "G2++.h"
#include <vector>
using namespace std;

class Swaption : public product {
    public:

 double strike_;

 double swap_maturity_;

 double tenor_;

 double payOff(const vector<double> &term_structure, const vector<double> &brownian_motion_1, const vector<double> &brownian_motion_2, const vector<double> &short_rate_path,const vector<double> & init_zc);

 Swaption(g2_model *g2mod, double swap_maturity, double maturity, double strike, double tenor);

 double sigma_1_tilde ();
 double sigma_2_tilde ();
 double rho_tilde();
 double mu_1_tilde();
 double mu_2_tilde();
 double A_T_Ti(double T_i,const vector<double> & init_zc,const vector<double> &term_structure);
 double c_i(double T_i);
 double kappa(double x, double T_i);
 double lambda_i(double x,double T_i,const vector<double> & init_zc,const vector<double> &term_structure);
 double x_bar(double x, const vector<double> & init_zc,const vector<double> &term_structure);
 double h_1 (double x, const vector<double> & init_zc,const vector<double> &term_structure);
 double h_2 (double x,double T_i,const vector<double> & init_zc,const vector<double> &term_structure);

 double inte_f(double x,const vector<double> & init_zc,const vector<double> &term_structure);

 double price_theory_Swaption(const vector<double> &term_structure,const vector<double> & init_zc);

 virtual ~Swaption();
};



#endif // SWAPTION_G2_H_INCLUDED
