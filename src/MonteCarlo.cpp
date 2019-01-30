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
	return price;
}

MonteCarlo::~MonteCarlo() {

}