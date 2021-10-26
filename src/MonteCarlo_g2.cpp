#include "MonteCarlo_g2.h"

MonteCarlo::MonteCarlo(g2_model *mod, product *prod, int nbSamples, double maturity, int step) {
	mod_ = mod;
	prod_ = prod;
	nbSamples_ = nbSamples;
	step_ = step;
	maturity_ = maturity;
}

double MonteCarlo::price(const vector<double> & init_zc) {
	vector<double> term_structure = mod_->term_structure_g2(maturity_, step_);
//    int N = term_structure.size();
//	double tenor = double(term_structure[N - 1]) / double(N - 1);
//	cout<<"tenor:"<<tenor<<endl;
//    //double tenor = 0.;
	vector<double> brownian_motion_path_1 = mod_->vector_movement_brownien_1(term_structure);
    vector<double> brownian_motion_path_2 = mod_->vector_movement_brownien_2(brownian_motion_path_1,term_structure);
	vector<double> short_rate_path = mod_->short_rate_g2(brownian_motion_path_1,brownian_motion_path_2, term_structure,init_zc);
	double payOff = prod_->payOff(term_structure, brownian_motion_path_1, brownian_motion_path_2,short_rate_path,init_zc);
	//cout<<"payoff:"<<payOff<<endl;
	double price = payOff;

	for (int i = 1; i < nbSamples_; i++) {
		//cout << i << endl;
		brownian_motion_path_1 = mod_->vector_movement_brownien_1(term_structure);
        brownian_motion_path_2 = mod_->vector_movement_brownien_2(brownian_motion_path_1,term_structure);
		short_rate_path = mod_->short_rate_g2(brownian_motion_path_1,brownian_motion_path_2, term_structure,init_zc);
		price += prod_->payOff(term_structure, brownian_motion_path_1, brownian_motion_path_2,short_rate_path,init_zc);
	}

	price = price / nbSamples_;
	return price;
}

MonteCarlo::~MonteCarlo() {

}
