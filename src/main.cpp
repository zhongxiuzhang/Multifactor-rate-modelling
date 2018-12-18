
#include <iostream>
#include <Eigen/Dense>

#include <vector>
#include <cmath>
#include <cstdlib>
#include "HW1.h"

using namespace std;
using Eigen::MatrixXd;



/*int main()
{
    cout << "Loi N "<<gaussian_box_muller() << endl;
    cout << gaussian_box_muller() << endl;
    return 0;
}*/



int main()
{
  int T=10, N=10, num_trj=3;
  /*MatrixXd m(num_trj,N+1);
  m = matrix_movement_brownien(T, N, num_trj);
  std::cout << m << std::endl;*/
    int last_year=4;
    int tenor_num=4;

    vector<double> term_structure_=term_structure_hw1(last_year,tenor_num);
    MatrixXd m;
    m = matrix_movement_brownien(term_structure_, num_trj);
    cout << m << std::endl;

    cout<< term_structure_.size() << endl;
    cout<< term_structure_[116] << endl;
    vector<double> ini_fwd_rate = initial_forward_rate(term_structure_);
    MatrixXd sr = short_rate_hw1(m, term_structure_, ini_fwd_rate, 2., 1.);
    cout << sr << std::endl;
    MatrixXd beta_= beta(sr, term_structure_);
    cout << beta_ << std::endl;
}