#include "HW1.h"
#include <iostream>
#include <cstdio>      /* printf, NULL */
#include <cstdlib>     /* srand, rand */
#include <ctime>       /* time */
#include <list>
#include <cmath>
#include <vector>
#include <Eigen/Dense>

using namespace std;
using Eigen::MatrixXd;

double gaussian_box_muller() {
  double x = 0.0;
  double y = 0.0;
  double euclid_sq = 0.0;

  do {
    x = 2.0 * rand() / static_cast<double>(RAND_MAX)-1;
    y = 2.0 * rand() / static_cast<double>(RAND_MAX)-1;
    euclid_sq = x*x + y*y;
  } while (euclid_sq >= 1.0);

  return x*sqrt(-2*log(euclid_sq)/euclid_sq);
}

vector<double> term_structure_hw1(const int& last_year,const int& step){
// To generate the term structure of Hull White model.
// Input: last year of interest rate, the tenor_num denote the number of points in one year.
// OUT: term structure of interest rate.
    int num=(last_year)*step+1;
    double tenor=1/(double)step;

    vector<double> term_structure;
    for(int i=0;i<num;i++){
          term_structure.push_back(tenor*(double)i);
    }
    return term_structure;
}

MatrixXd matrix_movement_brownien(const vector<double> & term_structure, int num_trj) 
//function to generate a matrix to save num_trj paths of brownian motion
//IN: term structure, number of paths num_trj
//OUT: a matrix of num_trj brownian motion 
{
  int N = term_structure.size();
  double pas = double(term_structure[N-1])/double(N-1);
  MatrixXd mov_brow(num_trj,N); 
  for(int i=0;i<num_trj;i++){
    mov_brow(i,0)=0.;
    for(int j=1;j<N;j++){
      mov_brow(i,j)=mov_brow(i,j-1)+sqrt(pas)*gaussian_box_muller();
    }
  }
  return mov_brow;
  
}

vector<double> initial_forward_rate(const vector<double>& term_structure){
//This is the f(0,t) in the short rate simulation which come from the market data.
//Input: term structure
//Output: f(0,t)
    int num=term_structure.size();
    vector<double> initial_forward_rate;
    for(int i=0;i<num;i++){
        initial_forward_rate.push_back(0.005*i);
    }
    return initial_forward_rate;
}

MatrixXd short_rate_hw1(const MatrixXd & mov_brow, const vector<double> & term_structure, const vector<double> & ini_fwd_rate, double vol, double a)
//function for generating the short rate
//IN: brownian motion matrix, initial forward rate f(0,t), volatility sigma, mean reverting rate a
//OUT: short rate r(T)
{
  int r = mov_brow.rows();
  int c = mov_brow.cols();
  int N = term_structure.size();
  double tenor = double(term_structure[N-1])/double(N-1);
  double last_date = double(term_structure[N-1]);
  double inte_sto =0.;
  MatrixXd rate(r,c);
  for(int i=0;i<r;i++){
    rate(i,0)=ini_fwd_rate[0];
    for(int j=1;j<N;j++){
      inte_sto = inte_sto + exp(a*tenor*j)*(mov_brow(i,j)-mov_brow(i,j-1));
      rate(i,j)=ini_fwd_rate[j]+pow(vol,2)/2*pow(a,2)*pow(1-exp(-a*tenor*j),2)+vol*exp(-a*tenor*j)*inte_sto;
    }
  }
  return rate;
}

MatrixXd beta(const MatrixXd & short_rate, const vector<double> & term_structure)
//function for generating beta
//IN: short rate matrix, term structure
//OUT: beta matrix
{
  int r = short_rate.rows();
  int c = short_rate.cols();
  int N = term_structure.size();
  double tenor = double(term_structure[N-1])/double(N-1);
  double inte = 0.;
  MatrixXd beta(r,c);
  for(int i=0;i<r;i++){
    beta(i,0)=1.;
    for(int j=1;j<N;j++){
      inte = inte + tenor*short_rate(i,j);
      beta(i,j)=exp(inte);
    }
  }
  return beta;
}
