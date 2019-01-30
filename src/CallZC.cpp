#include "CallZC.h"
#include <math.h>

CallZC::CallZC(HW1 *hjm, double zc_maturity, int maturity, double strike) {
	mod_ = hjm;
	zc_maturity_ = zc_maturity;
	T_ = maturity; 
	strike_ = strike;
}

double CallZC::payOff(const vector<double> &term_structure, const vector<double> &brownian_motion, const vector<double> &short_rate_path) {
	double zero_coupon = mod_->zero_coupon_bonds_hw1(T_, zc_maturity_, brownian_motion, term_structure);
	return max(mod_->zero_coupon_bonds_hw1(T_, zc_maturity_, brownian_motion, term_structure) - strike_, 0.);
}

CallZC::~CallZC() {

}