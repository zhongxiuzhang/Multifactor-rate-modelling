#include "CallZC.h"
#include "Math.h"
#include <math.h>
#include <iostream>

CallZC::CallZC(HW1 *hjm, double zc_maturity, int maturity, double strike) {
	mod_ = hjm;
	zc_maturity_ = zc_maturity;
	T_ = maturity;
	strike_ = strike;
}

double CallZC::payOff(const vector<double> &term_structure, const vector<double> &brownian_motion, const vector<double> &short_rate_path) {
	double zero_coupon = mod_->zero_coupon_bonds_hw1(T_, zc_maturity_, brownian_motion, term_structure);
	return max(mod_->zero_coupon_bonds_hw1(T_, zc_maturity_, brownian_motion, term_structure), 0.);
}

double CallZC::price_theory_CallZC(const vector<double> &term_structure, const vector<double> & initial_zero_coupon_, double zc_maturity, double Maturity, double strike){
   double a=mod_->a_;
   double Y=(1-exp(-a*(zc_maturity-Maturity)))/a;
   double sigma=mod_->vol_;
  // vector<double> initial_forward_rate=copy_term_structure(mod_->ini_fwd_rate_);
  // vector<double> initial_zero_coupon_=mod_->initial_zero_coupon(initial_forward_rate, term_structure);
   double sigmaP2=sigma*sigma*(1-exp(-2*a*Maturity))*Y*Y/(2*a);
   double sigmaP=sqrt(sigmaP2);
   // find t, T position
	vector<double> term_structure_S = copy_term_structure(term_structure);
	vector<double> term_structure_T = copy_term_structure(term_structure);

	for (size_t i = 0; i<term_structure_S.size(); i++) {
		term_structure_S[i] = abs(term_structure_S[i] - zc_maturity);
	}
	vector<double>::iterator Smallest_S = min_element(term_structure_S.begin(), term_structure_S.end());
	int index_S = distance(term_structure_S.begin(), Smallest_S);

	for (size_t i = 0; i<term_structure_T.size(); i++) {
		term_structure_T[i] = abs(term_structure_T[i] - Maturity);
	}
	vector<double>::iterator Smallest_T = min_element(term_structure_T.begin(), term_structure_T.end());
	int index_T = distance(term_structure_T.begin(), Smallest_T);

	cout<<"index_S:"<<initial_zero_coupon_[index_S]<<endl<<"index_T:"<<initial_zero_coupon_[index_T]<<endl;
		for(int j = 0;j<60;j++){
		cout<<"index:"<<j<<"  zc:"<<initial_zero_coupon_[j]<<endl;
	}

   double h=log(initial_zero_coupon_[index_S]/(initial_zero_coupon_[index_T]*strike))/sigmaP+sigmaP/2;
   double price=initial_zero_coupon_[index_S]*cdfGaussian(h)-strike*initial_zero_coupon_[index_T]*cdfGaussian(h-sigmaP);
   return price;
}

CallZC::~CallZC() {

}
