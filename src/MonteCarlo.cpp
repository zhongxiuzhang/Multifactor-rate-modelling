#include "MonteCarlo.h"

MonteCarlo::MonteCarlo(HW1 *mod, Option *opt, int nbSamples, int step) {
	mod_ = mod; 
	opt_ = opt; 
	nbSamples_ = nbSamples;
	step_ = step;
}

double MonteCarlo::price() {
	vector<double> term_structure = mod_->term_structure_hw1(opt_->T_, step_);
	vector<double> brownian_motion_path = mod_->vector_movement_brownien(term_structure);
	vector<double> short_rate_path = mod_->short_rate_hw1(brownian_motion_path, term_structure);
	double payOff = opt_->payOff(short_rate_path);
	double price = payOff;

	for (int i = 1; i < nbSamples_; i++) {
		brownian_motion_path = mod_->vector_movement_brownien(term_structure);
		short_rate_path = mod_->short_rate_hw1(brownian_motion_path, term_structure);
		price += opt_->payOff(short_rate_path);
	}

	price = price / nbSamples_;
	return price;
}

MonteCarlo::~MonteCarlo() {

}