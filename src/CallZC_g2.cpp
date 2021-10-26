#include "CallZC_g2.h"
#include "Math.h"
#include <math.h>
#include <iostream>
#include <numeric>

CallZC::CallZC(g2_model *g2mod, double zc_maturity, double maturity, double strike) {
	mod_ = g2mod;
	zc_maturity_ = zc_maturity;
	T_ = maturity;
	strike_ = strike;
}

double CallZC::payOff(const vector<double> &term_structure, const vector<double> &brownian_motion_1, const vector<double> &brownian_motion_2, const vector<double> &short_rate_path,const vector<double> & init_zc) {
	double zero_coupon = mod_->zero_coupon_bonds_g2(T_, zc_maturity_, brownian_motion_1, brownian_motion_2, term_structure, init_zc);
//    vector<double> vect(10);
//    double discount_value=mod_->zero_coupon_bonds_g2(0,T_,vect,vect,term_structure,init_zc);
    int N = term_structure.size();
    double tenor = double(term_structure[N - 1]) / double(N - 1);
    int index_T = int (T_/tenor);
    double discount_value = 0.;
//    for (int i=0; i<index_T-1; i++){
//            discount_value  = discount_value + short_rate_path[i]*tenor;
//    }
//    discount_value = exp(-discount_value);
    vector<double> vect(10);
    discount_value=mod_->zero_coupon_bonds_g2(0,T_,vect,vect,term_structure,init_zc);
//    cout<<"zero_coupon_bond" << zero_coupon <<endl;
//    cout<<"discount_value" << discount_value <<endl;

	return max(zero_coupon-strike_, 0.)*discount_value;
}


double CallZC::price_theory_CallZC(const vector<double> & init_zc,const vector<double> & term_s){
   double T=T_;
   double S=zc_maturity_;
   double strike=strike_;

   vector<double> vect(10);

   double P_t_T=mod_->zero_coupon_bonds_g2(0,T,vect,vect,term_s,init_zc);
   cout<<"P_t_T:"<<P_t_T<<endl;
   double P_t_S=mod_->zero_coupon_bonds_g2(0,S,vect,vect,term_s,init_zc);
   cout<<"P_t_S:"<<P_t_S<<endl;
   double sigma_hat_sq = pow(mod_->sigma_1,2)*(1-exp(-2*(mod_->k_1)*T))*pow(mod_->int_1(T,S),2)/2/mod_->k_1
                         + pow(mod_->sigma_2,2)*(1-exp(-2*(mod_->k_2)*T))*pow(mod_->int_2(T,S),2)/2/mod_->k_2
                         + 2*mod_->sigma_1*mod_->sigma_2*mod_->rho*mod_->int_1(T,S)*mod_->int_2(T,S)*mod_->int_1_2(0,T);
   double sigma_hat = sqrt(sigma_hat_sq);
//   cout<<"sigma_1:"<<(pow(mod_->sigma_1,2)*(1-exp(-2*(mod_->k_1)*T))*pow(mod_->int_1(T,S),2)/2/mod_->k_1)<<endl;
//   cout<<"sigma_2:"<<(pow(mod_->sigma_2,2)*(1-exp(-2*(mod_->k_2)*T))*pow(mod_->int_2(T,S),2)/2/mod_->k_2)<<endl;
//   cout<<"sigma_3:"<<(2*mod_->sigma_1*mod_->sigma_2*mod_->rho*mod_->int_1(T,S)*mod_->int_2(T,S)*mod_->int_1_2(0,T))<<endl;
//   cout<<"sigma_hat_sq :"<<sigma_hat_sq <<endl;
   cout<<"sigma_hat:"<<sigma_hat<<endl;
   cout<<"sigma_1:"<<mod_->sigma_1<<endl;
   cout<<"sigma_2:"<<mod_->sigma_2<<endl;
   cout<<"k_1:"<<mod_->k_1<<endl;
   cout<<"k_2:"<<mod_->k_2<<endl;
   cout<<"mod_->rho:"<<mod_->rho<<endl;
   double h =  log(P_t_S/P_t_T/strike)/sigma_hat+sigma_hat/2;
   cout<<"h:"<<h<<endl;
   cout<<"cdf h:"<<cdfGaussian(h)<<endl;
   double price=P_t_S*cdfGaussian(h)-strike*P_t_T*cdfGaussian(h-sigma_hat);
   return price;
}

CallZC::~CallZC() {

}
