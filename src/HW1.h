#include <vector>
#include <cmath>
#ifndef HW1_H_INCLUDED
#define HW1_H_INCLUDED
#include <Eigen/Dense>

using namespace std;
using Eigen::MatrixXd;
// For random Gaussian generation
// Note that there are many ways of doing this, but we will
// be using the Box-Muller method for demonstration purposes
double gaussian_box_muller();

MatrixXd matrix_movement_brownien(const vector<double> & term_structure, int num_trj);

vector<double> initial_forward_rate(const vector<double>& term_structure);

MatrixXd short_rate_hw1(const MatrixXd & mov_brow, const vector<double> & term_structure, const vector<double> & ini_fwd_rate, double vol, double a);

double forward_rate_hw1();

MatrixXd beta(const MatrixXd & short_rate, const vector<double> & term_structure);

double zero_coupon_bonds_hw1();

vector<double> term_structure_hw1(const int& last_year,const int& tenor_num);


#endif // HW1_H_INCLUDED
