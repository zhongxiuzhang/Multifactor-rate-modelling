#include "MonteCarlo.h"


MonteCarlo::MonteCarlo(HW1 *mod, Option *opt, int nbSamples, double maturity, int step) {
	mod_ = mod; 
	opt_ = opt; 
	nbSamples_ = nbSamples;
	step_ = step;
	maturity_ = maturity;
}

double MonteCarlo::price() {
	vector<double> term_structure = mod_->term_structure_hw1(maturity_, step_);
	vector<double> brownian_motion_path = mod_->vector_movement_brownien(term_structure);
	vector<double> short_rate_path = mod_->short_rate_hw1(brownian_motion_path, term_structure);
	double payOff = opt_->payOff(term_structure, brownian_motion_path, short_rate_path);
	double price = payOff;

	for (int i = 1; i < nbSamples_; i++) {
		cout << i << endl;
		brownian_motion_path = mod_->vector_movement_brownien(term_structure);
		short_rate_path = mod_->short_rate_hw1(brownian_motion_path, term_structure);
		price += opt_->payOff(term_structure, brownian_motion_path, short_rate_path);
	}

	price = price / nbSamples_;
	vector<double> initial_forward_rate = mod_->ini_fwd_rate_;
	vector<double> initial_zero_coupon = mod_->initial_zero_coupon(initial_forward_rate, term_structure);
	vector<double> term_structure_T = copy_term_structure(term_structure);
	for (size_t i = 0; i<term_structure_T.size(); i++) {
		term_structure_T[i] = abs(term_structure_T[i] - opt_->T_);
	}
	vector<double>::iterator Smallest_T = min_element(term_structure_T.begin(), term_structure_T.end());
	int index_T = distance(term_structure_T.begin(), Smallest_T);
	
	price = price * initial_zero_coupon[index_T];
	return price;
}

MonteCarlo::~MonteCarlo() {

}