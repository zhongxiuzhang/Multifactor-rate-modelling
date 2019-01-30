#include "CallZC.h"
#include <math.h>

CallZC::CallZC(double zc_maturity, int maturity, double strike) {
	zc_maturity_ = zc_maturity;
	T_ = maturity; 
	strike_ = strike;
}

double CallZC::payOff(const vector<double> &term_structure, const vector<double> &brownian_motion, const vector<double> &short_rate_path) {
	double zero_coupon = mod_->zero_coupon_bonds_hw1(T_, zc_maturity_, brownian_motion, term_structure);
	return max(mod_->zero_coupon_bonds_hw1(T_, zc_maturity_, brownian_motion, term_structure), 0.);
}

CallZC::~CallZC() {

}