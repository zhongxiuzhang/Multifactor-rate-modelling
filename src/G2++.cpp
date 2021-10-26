#include "G2++.h"
#include <math.h>
#include <stdlib.h>     /* srand, rand */
#include <iostream>

g2_model::g2_model(double k_1_, double k_2_, double sigma_1_, double sigma_2_,double rho_) {
	k_1 = k_1_;
	k_2 = k_2_;
    sigma_1 = sigma_1_;
    sigma_2 = sigma_2_;
    rho = rho_;
}

vector<double> g2_model::term_structure_g2(const int& last_year, const int& step) {
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

double g2_model::gaussian_box_muller() {
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

vector<double>  g2_model::vector_movement_brownien_1(const vector<double> & term_structure)
//function to generate a matrix to save num_trj paths of brownian motion
//IN: term structure, number of paths num_trj
//OUT: a matrix of num_trj brownian motion
{
	int N = term_structure.size();
	double pas = double(term_structure[N - 1]) / double(N - 1);
	vector<double> mov_brow_1;
	mov_brow_1.push_back(0.);
	for (int j = 1; j<N; j++) {
		mov_brow_1.push_back(mov_brow_1[j - 1] + sqrt(pas)*gaussian_box_muller());
	}
	return mov_brow_1;
}

vector<double>  g2_model::vector_movement_brownien_2(const vector<double> & mov_brow_1,const vector<double> & term_structure)
//function to generate a matrix to save num_trj paths of brownian motion
//IN: term structure, number of paths num_trj
//OUT: a matrix of num_trj brownian motion
{
	int N = term_structure.size();
	double pas = double(term_structure[N - 1]) / double(N - 1);
	vector<double> mov_brow_2;
	mov_brow_2.push_back(0.);
	for (int j = 1; j<N; j++) {
		mov_brow_2.push_back(mov_brow_2[j - 1] + rho*(mov_brow_1[j]-mov_brow_1[j-1])+sqrt(1-rho*rho)*sqrt(pas)*gaussian_box_muller());
	}
	return mov_brow_2;
}

double g2_model::int_1(double t,double T){
    double int_1=(1-exp(-k_1*(T-t)))/k_1;
    return int_1;
}

double g2_model::int_2(double t,double T){
    double int_2=(1-exp(-k_2*(T-t)))/k_2;
    return int_2;
}

double g2_model::int_1_2(double t,double T){
    double int_1_2=(1-exp(-(k_1+k_2)*(T-t)))/(k_1+k_2);
    return int_1_2;
}


double g2_model::x_1(double t, const vector<double> & mov_brow_1, const vector<double> & term_structure){
	int N = term_structure.size();
	double tenor = double(term_structure[N - 1]) / double(N - 1);
	// find time interval of t in the term structure
	int t_left_index;
	int t_right_index;
    for (int i = 0; i<int(term_structure.size()); i++) {
            if ( t<0 || t>term_structure.back()){
                cout<< "t is not in the right interval"<<endl;
            }else if ( t>i*tenor && t<=(i+1)*tenor ){
                    t_left_index=i;
                    t_right_index=i+1;
            }else if (t==0){
                    t_left_index=0;
                    t_right_index=1;
            }
	}
	double x_left=0.;
	double x_right=0.;
    for (int i=0; i<t_left_index-1; i++){
            x_left = x_left + exp(k_1*tenor*i)*(mov_brow_1[i+1] - mov_brow_1[i]);
    }
    x_left=sigma_1*exp(-k_1*tenor*t_left_index)*x_left;
    for (int i=0; i<t_right_index-1; i++){
            x_right = x_right + exp(k_1*tenor*i)*(mov_brow_1[i+1] - mov_brow_1[i]);
    }
    x_right=sigma_1*exp(-k_1*tenor*t_right_index)*x_right;

    double x_=0.;
    x_=(t-t_left_index*tenor)/tenor * x_right + (t_right_index*tenor-t)/tenor * x_left;
    return x_;
}

double g2_model::x_2(double t, const vector<double> & mov_brow_2, const vector<double> & term_structure){
	int N = term_structure.size();
	double tenor = double(term_structure[N - 1]) / double(N - 1);
	// find time interval of t in the term structure
	int t_left_index;
	int t_right_index;
    for (int i = 0; i<int(term_structure.size()); i++) {
            if ( t<0 || t>term_structure.back()){
                cout<< "t is not in the right interval"<<endl;
            }else if ( t>i*tenor && t<=(i+1)*tenor ){
                    t_left_index=i;
                    t_right_index=i+1;
            }else if (t==0){
                    t_left_index=0;
                    t_right_index=1;
            }
	}
	double x_left=0.;
	double x_right=0.;
    for (int i=0; i<t_left_index-1; i++){
            x_left = x_left + exp(k_2*tenor*i)*(mov_brow_2[i+1] - mov_brow_2[i]);
    }
    x_left=sigma_2*exp(-k_2*tenor*t_left_index)*x_left;
    for (int i=0; i<t_right_index-1; i++){
            x_right = x_right + exp(k_2*tenor*i)*(mov_brow_2[i+1] - mov_brow_2[i]);
    }
    x_right=sigma_2*exp(-k_2*tenor*t_right_index)*x_right;

    double x_=0.;
    x_=(t-t_left_index*tenor)/tenor * x_right + (t_right_index*tenor-t)/tenor * x_left;
    return x_;
}

//double g2_model::zero_coupon_bonds_g2(double t, double T, const vector<double> & mov_brow_1,
//       const vector<double> & mov_brow_2, const vector<double> & term_structure){
//           double x_1_;
//           double x_2_;
//           double b_1;
//           double b_2;
//           double b_1_2;
//           double m_;
//           double v_;
//           x_1_=x_1(t, mov_brow_1,term_structure);
//           x_2_=x_2(t, mov_brow_2,term_structure);
//           b_1=int_1(t,T);
//           b_2=int_2(t,T);
//           b_1_2=int_1_2(t,T);
//           m_=x_1_*b_1+x_2_*b_2;
//           v_=(T-t-b_1-k_1*b_1*b_1/2)*sigma_1*sigma_1/k_1/k_1;
//           v_=v_+(T-t-b_2-k_2*b_2*b_2/2)*sigma_2*sigma_2/k_2/k_2;
//           v_=v_+(T-t-b_1-b_2+b_1_2)*2*rho*sigma_1*sigma_2/k_1/k_2;
//           double into_phi=0.;
//           double dt=(T-t)/1000;
//           for (int i=0; i<1000; i++){
//            into_phi=into_phi+phi(t+dt*double(i))*dt;
//           }
//           double zero_coupon_=0.;
//           zero_coupon_=exp(-into_phi-m_+v_/2.);
//           return zero_coupon_;
//}



//double g2_model::phi(double t)
//{
//    return 0.00005*t+0.001;
//}

double g2_model::phi_calib(const vector<double> &init_zc,const vector<double> &term_structure,double t)
{
    int N = term_structure.size();
    double tenor = double(term_structure[N - 1]) / double(N - 1);
    int t_left_index;
    int t_right_index;
    for (int i = 0; i<int(term_structure.size()); i++) {
            if ( t<0 || t>term_structure.back()){
                cout<< "t is not in the right interval"<<endl;
            }else if ( t>i*tenor && t<=(i+1)*tenor ){
                    t_left_index=i;
                    t_right_index=i+1;
            }else if (t==0){
                    t_left_index=0;
                    t_right_index=1;
            }
    }
    double phi_left;
    double phi_right;
    if (t_left_index==0){phi_left=0.;
    }else {
        phi_left=(log(init_zc[t_left_index-1])-log(init_zc[t_left_index]))/tenor+pow(sigma_1*int_1(0,t_left_index*tenor),2)/2.+pow(sigma_2*int_2(0,t_left_index*tenor),2)/2.+sigma_1*sigma_2*rho*int_1(0,t_left_index*tenor)*int_2(0,t_left_index*tenor);
    }
    phi_right=(log(init_zc[t_right_index-1])-log(init_zc[t_right_index]))/tenor+pow(sigma_1*int_1(0,t_right_index*tenor),2)/2.+pow(sigma_2*int_2(0,t_right_index*tenor),2)/2.+sigma_1*sigma_2*rho*int_1(0,t_right_index*tenor)*int_2(0,t_right_index*tenor);
    double phi=(t-t_left_index*tenor)/tenor * phi_right + (t_right_index*tenor-t)/tenor * phi_left;

    return phi;
}

vector<double> g2_model::short_rate_g2(const vector<double> & mov_brow_1,const vector<double> & mov_brow_2,const vector<double> & term_structure,const vector<double> & init_zc)
{
    int N = term_structure.size();
    vector<double> short_rate;
    for (int i = 0; i<N; i++){
        short_rate.push_back(phi_calib(init_zc,term_structure,term_structure[i])+x_1(term_structure[i],mov_brow_1,term_structure)+x_2(term_structure[i],mov_brow_2,term_structure));
    }
    return short_rate;
}

double g2_model::zero_coupon_bonds_g2(double t, double T, const vector<double> & mov_brow_1,
       const vector<double> & mov_brow_2, const vector<double> & term_structure,const vector<double> & init_zc){
           double x_1_; double x_2_;
           double b_1; double b_2;
           double b_1_2; double m_; double v_;
           x_1_=x_1(t, mov_brow_1,term_structure);
           //cout<<"x_1_"<<x_1_<<endl;
           x_2_=x_2(t, mov_brow_2,term_structure);
           //cout<<"x_2_"<<x_2_<<endl;
           b_1=int_1(t,T);
           b_2=int_2(t,T);
           b_1_2=int_1_2(t,T);
           m_=x_1_*b_1+x_2_*b_2;
           v_=(T-t-b_1-k_1*b_1*b_1/2)*sigma_1*sigma_1/k_1/k_1;
           v_=v_+(T-t-b_2-k_2*b_2*b_2/2)*sigma_2*sigma_2/k_2/k_2;
           v_=v_+(T-t-b_1-b_2+b_1_2)*2*rho*sigma_1*sigma_2/k_1/k_2;
           double into_phi=0.;
           double dt=(T-t)/1000.;
           for (int i=0; i<1000; i++){
            into_phi=into_phi+phi_calib(init_zc,term_structure,t+dt*double(i))*dt;
           }
           double zero_coupon_=0.;
           zero_coupon_=exp(-into_phi-m_+v_/2.);
           return zero_coupon_;
}

double g2_model::forward_rate_g2(double t, double T, const vector<double> & mov_brow_1, const vector<double> & mov_brow_2, const vector<double> & term_structure,const vector<double> & init_zc)
{
    double x_1_=x_1(t, mov_brow_1, term_structure);
    double x_2_=x_2(t, mov_brow_2, term_structure);
    double fwd_rate=phi_calib(init_zc,term_structure,T) + x_1_*exp(-k_1*(T-t)) + x_2_*exp(-k_2*(T-t)) - pow(sigma_1,2)*pow(int_1(t,T),2)/2 - pow(sigma_2,2)*pow(int_2(t,T),2)/2 \
                    - sigma_1*sigma_2*rho*int_1(t,T)*int_2(t,T);
	return fwd_rate;

}

//vector<double> g2_model::short_rate_g2(const vector<double> & mov_brow_1,const vector<double> & mov_brow_2,const vector<double> & term_structure)
//{
//    int N = term_structure.size();
//    vector<double> short_rate;
//    for (int i = 0; i<N; i++){
//        short_rate.push_back(phi(term_structure[i])+x_1(term_structure[i],mov_brow_1,term_structure)+x_2(term_structure[i],mov_brow_2,term_structure));
//    }
//    return short_rate;
//}
//
//double g2_model::forward_rate_g2(double t, double T, const vector<double> & mov_brow_1, const vector<double> & mov_brow_2, const vector<double> & term_structure)
//{
//    double x_1_=x_1(t, mov_brow_1, term_structure);
//    double x_2_=x_2(t, mov_brow_2, term_structure);
//    double fwd_rate=phi(T) + x_1_*exp(-k_1*(T-t)) + x_2_*exp(-k_2*(T-t)) - pow(sigma_1,2)*pow(int_1(t,T),2)/2 - pow(sigma_2,2)*pow(int_2(t,T),2)/2 \
//                    - sigma_1*sigma_2*rho*int_1(t,T)*int_2(t,T);
//	return fwd_rate;
//
//}

//vector<double> g2_model::x_1(const vector<double> & mv_brow_1, const vector<double> & term_structure){
//    int N = term_structure.size();
// double tenor = double(term_structure[N - 1]) / double(N - 1);
// double inte_sto = 0.;
// vector<double> x_1;
// x_1.push_back(inte_sto);
// for (int j = 1; j<N; j++) {
//        inte_sto = 0.;
//        for (int i=1; i<j; i++){
//                inte_sto = inte_sto + exp(k_1*tenor*i)*(mv_brow_1[i] - mv_brow_1[i - 1]);
//            }
//        inte_sto = sigma_1*exp(-k_1*tenor*j)*inte_sto;
//  x_1.push_back(inte_sto);
// }
// return x_1;
//}
//vector<double> g2_model::x_2(const vector<double> & mv_brow_2, const vector<double> & term_structure){
//    int N = term_structure.size();
// double tenor = double(term_structure[N - 1]) / double(N - 1);
// double inte_sto = 0.;
// vector<double> x_2;
// x_2.push_back(inte_sto);
// for (int j = 1; j<N; j++) {
//        inte_sto = 0.;
//        for (int i=1; i<j; i++){
//                inte_sto = inte_sto + exp(k_2*tenor*i)*(mv_brow_2[i] - mv_brow_2[i - 1]);
//            }
//        inte_sto = sigma_1*exp(-k_2*tenor*j)*inte_sto;
//  x_2.push_back(inte_sto);
// }
// return x_2;
//}

//vector<double> copy_term_structure(vector<double> const &term_structure) {
//	vector<double> copy_term_strcuture;
//	for (int i = 0; i < term_structure.size(); i++) {
//		copy_term_strcuture.push_back(term_structure[i]);
//	}
//	return copy_term_strcuture;
//}

