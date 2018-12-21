#ifndef HW1_H_INCLUDED
#define HW1_H_INCLUDED

#include <Eigen/Dense>
#include <vector>

using namespace std;
using Eigen::MatrixXd;

class HW1 {
public:
	// For random Gaussian generation
	// Note that there are many ways of doing this, but we will
	// be using the Box-Muller method for demonstration purposes
	double vol_; 
	double a_;
	vector<double> ini_fwd_rate_;
	double gaussian_box_muller();

	HW1(double a, double vol, vector<double> & ini_fwd_rate);
	

	vector<double>  vector_movement_brownien(const vector<double> & term_structure);

	// To generate the term structure of Hull White model. The tenor_num denote the number of points in one year.
	vector<double> term_structure_hw1(const int& last_year, const int& step);

	//This is the f(0,t) in the short rate simulation which come from the market data.
	vector<double> initial_forward_rate(const vector<double>& term_structure);

	vector<double> short_rate_hw1(const vector<double> & mov_brow, const vector<double> & term_structure);

	vector<double> beta(const vector<double> & short_rate, const vector<double> & term_structure);

	double forward_rate_hw1(double t, double T, const vector<double> & mov_brow, const vector<double>& term_structure);

	double zero_coupon_bonds_hw1(double t, double T, const vector<double> & mov_brow, const vector<double>& term_structure);

};
#endif // HW1_H_INCLUDED
