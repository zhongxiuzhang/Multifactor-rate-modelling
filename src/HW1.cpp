#include "HW1.h"
#include <Eigen/Dense>

#include <algorithm> // max_element, min_element
#include <vector>

using namespace std;
using Eigen::MatrixXd;


HW1::HW1(double a, double vol) {
	a_ = a;
	vol_ = vol;
}

double HW1::gaussian_box_muller() {
	double x = 0.0;
	double y = 0.0;
	double euclid_sq = 0.0;

	do {
		x = 2.0 * rand() / static_cast<double>(RAND_MAX) - 1;
		y = 2.0 * rand() / static_cast<double>(RAND_MAX) - 1;
		euclid_sq = x * x + y * y;
	} while (euclid_sq >= 1.0);

	return x * sqrt(-2 * log(euclid_sq) / euclid_sq);
}

vector<double> HW1::term_structure_hw1(const int& last_year, const int& step) {
	// To generate the term structure of Hull White model.
	// Input: last year of interest rate, the tenor_num denote the number of points in one year.
	// OUT: term structure of interest rate.
	int num = (last_year)*step + 1;
	double tenor = 1 / (double)step;

	vector<double> term_structure;
	for (int i = 0; i<num; i++) {
		term_structure.push_back(tenor*(double)i);
	}
	return term_structure;
}

vector<double>  HW1::vector_movement_brownien(const vector<double> & term_structure)
//function to generate a matrix to save num_trj paths of brownian motion
//IN: term structure, number of paths num_trj
//OUT: a matrix of num_trj brownian motion
{
	int N = term_structure.size();
	double pas = double(term_structure[N - 1]) / double(N - 1);
	vector<double> mov_brow;
	mov_brow.push_back(0.);
	for (int j = 1; j<N; j++) {
		mov_brow.push_back(mov_brow[j - 1] + sqrt(pas)*gaussian_box_muller());
	}
	return mov_brow;
}

vector<double> HW1::initial_forward_rate(const vector<double>& term_structure) {
	//This is the f(0,t) in the short rate simulation which come from the market data.
	//Input: term structure
	//Output: f(0,t)
	int num = term_structure.size();
	vector<double> ini_fwd_rate_;
	for (int i = 0; i<num; i++) {
		ini_fwd_rate_.push_back(0.005*i);
	}
	return ini_fwd_rate_;
}

vector<double> HW1::initial_zero_coupon(const vector<double>& ini_fwd_rate_,const vector<double>& term_structure) {

    int N = term_structure.size();
	double tenor = double(term_structure[N - 1]) / double(N - 1);
	int num = ini_fwd_rate_.size();
	vector<double> initial_zero_coupon_;
	initial_zero_coupon_.push_back(1);
	for (int i = 1; i<num; i++) {
		initial_zero_coupon_.push_back(initial_zero_coupon_[i-1]*exp(-ini_fwd_rate_[i]*tenor));
	}
	return initial_zero_coupon_;
}

vector<double> HW1::short_rate_hw1(const vector<double> & mov_brow, const vector<double> & term_structure)
//function for generating the short rate
//IN: brownian motion matrix, initial forward rate f(0,t), volatility sigma, mean reverting rate a
//OUT: short rate r(T)
{
	int N = term_structure.size();
	double tenor = double(term_structure[N - 1]) / double(N - 1);
	double inte_sto = 0.;
	vector<double> rate;
	rate.push_back(ini_fwd_rate_[0]);
	for (int j = 1; j<N; j++) {
		inte_sto = inte_sto + exp(a_*tenor*j)*(mov_brow[j] - mov_brow[j - 1]);
		rate.push_back(ini_fwd_rate_[j] + pow(vol_, 2) / 2 * pow(a_, 2)*pow(1 - exp(-a_ * tenor*j), 2) + vol_ * exp(-a_ * tenor*j)*inte_sto);
	}
	return rate;
}

vector<double> HW1::beta(const vector<double> & short_rate, const vector<double> & term_structure)
//function for generating beta
//IN: short rate matrix, term structure
//OUT: beta matrix
{
	int N = term_structure.size();
	double tenor = double(term_structure[N - 1]) / double(N - 1);
	double inte = 0.;
	vector<double> beta;
	beta.push_back(1.);
	for (int j = 1; j<N; j++) {
		inte = inte + tenor * short_rate[j];
		beta.push_back(exp(inte));
	}
	return beta;
}


//MatrixXd forward_rate_hw1(double t, double T, const MatrixXd & mov_brow, const vector<double> & ini_fwd_rate,const vector<double>& term_structure, double vol, double a)
double HW1::forward_rate_hw1(double t, double T, const vector<double> & mov_brow, const vector<double>& term_structure)
//function to generating the forward rate f(t,T)
//IN: t, T, brownian motion matrix, initial forward rate f(0,t), volatility sigma, mean reverting rate a
//OUT: forward rate f(t,T)
{
	int N = term_structure.size();
	double tenor = double(term_structure[N - 1]) / double(N - 1);
	double inte_sto = 0.;

	// find t, T position
	vector<double> term_structure_t = copy_term_structure(term_structure);
	vector<double> term_structure_T = copy_term_structure(term_structure);

	for (size_t i = 0; i<term_structure_t.size(); i++) {
		term_structure_t[i] = abs(term_structure_t[i] - t);
	}
	vector<double>::iterator Smallest_t = min_element(term_structure_t.begin(), term_structure_t.end());
	int index_t = distance(term_structure_t.begin(), Smallest_t);

	for (size_t i = 0; i<term_structure_T.size(); i++) {
		term_structure_T[i] = abs(term_structure_T[i] - T);
	}
	vector<double>::iterator Smallest_T = min_element(term_structure_T.begin(), term_structure_T.end());
	int index_T = distance(term_structure_T.begin(), Smallest_T);

	double forward_rate;

	for (int j = 1; j <= index_t; j++) {
		inte_sto = inte_sto + exp(a_*tenor*j)*(mov_brow[j] - mov_brow[j - 1]);
	}
	forward_rate = ini_fwd_rate_[index_T] + (vol_*vol_)*exp(-a_ * term_structure[index_t])*(exp(a_*term_structure[index_t]) - 1) / (a_*a_) + vol_ * exp(-a_ * tenor*index_T)*inte_sto;
	return forward_rate;
}

double HW1::zero_coupon_bonds_hw1(double t, double T, const vector<double> & mov_brow, const vector<double>& term_structure)
//function for generating the zero coupon.
//IN: t, T, brownian motion matrix, initial forward rate f(0,t), volatility sigma, mean reverting rate a
//OUT: zero coupon
{

	int N = term_structure.size();
	double tenor = double(term_structure[N - 1]) / double(N - 1);

	double zero_coupon_bonds_hw1 = 1.;

	// find t, T position
	vector<double> term_structure_t = copy_term_structure(term_structure);
	vector<double> term_structure_T = copy_term_structure(term_structure);

	for (size_t i = 0; i<term_structure_t.size(); i++) {
		term_structure_t[i] = abs(term_structure_t[i] - t);
	}
	vector<double>::iterator Smallest_t = min_element(term_structure_t.begin(), term_structure_t.end());
	int index_t = distance(term_structure_t.begin(), Smallest_t);

	for (size_t i = 0; i<term_structure_T.size(); i++) {
		term_structure_T[i] = abs(term_structure_T[i] - T);
	}
	vector<double>::iterator Smallest_T = min_element(term_structure_T.begin(), term_structure_T.end());
	int index_T = distance(term_structure_T.begin(), Smallest_T);

	double t_fwd_rate = term_structure[index_t];
	double s_fwd_rate;

	double forward_rate_hw1_;

	for (int i = index_t; i<index_T; i++)
	{
		s_fwd_rate = term_structure[i];
		forward_rate_hw1_ = forward_rate_hw1(t_fwd_rate, s_fwd_rate, mov_brow, term_structure);
		zero_coupon_bonds_hw1 = zero_coupon_bonds_hw1 * exp(-forward_rate_hw1_ * tenor);
	}
	return zero_coupon_bonds_hw1;
}




vector<double> copy_term_structure(vector<double> const &term_structure) {
	vector<double> copy_term_strcuture;
	for (int i = 0; i < term_structure.size(); i++) {
		copy_term_strcuture.push_back(term_structure[i]);
	}
	return copy_term_strcuture;
}