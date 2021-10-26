#include "Swaption_g2.h"
#include "Math.h"
#include <math.h>
#include <iostream>
#include <numeric>

using namespace std;


Swaption::Swaption(g2_model *g2mod, double swap_maturity, double maturity, double strike, double tenor) {
 mod_ = g2mod;
 swap_maturity_ = swap_maturity;
 T_ = maturity;
 strike_ = strike;
 tenor_ = tenor;
}

double Swaption::payOff(const vector<double> &term_structure, const vector<double> &brownian_motion_1, const vector<double> &brownian_motion_2, const vector<double> &short_rate_path,const vector<double> & init_zc) {
 double N=(swap_maturity_-T_)/tenor_;
    double A_0=0.; double A_t0=0.;
    for(int i=0;i<N;i++){
        A_0=A_0+tenor_*mod_->zero_coupon_bonds_g2(0, T_+(i+1)*tenor_, brownian_motion_1, brownian_motion_2, term_structure, init_zc);
        A_t0=A_t0+tenor_*mod_->zero_coupon_bonds_g2(T_, T_+(i+1)*tenor_, brownian_motion_1, brownian_motion_2, term_structure, init_zc);
    }
    double P_t0_Tn=mod_->zero_coupon_bonds_g2(T_, swap_maturity_, brownian_motion_1, brownian_motion_2, term_structure, init_zc);
    return A_0*max((1-P_t0_Tn)/A_t0-strike_,0.);
}


double Swaption::sigma_1_tilde (){
    double T = T_;
    double sigma_1 = mod_->sigma_1;
    double k_1 = mod_->k_1;
    double simga_1_tilde_=sigma_1*sqrt((1.-exp(-2.*k_1*T))/2./k_1);
	return simga_1_tilde_;
}

double Swaption::sigma_2_tilde (){
    double T = T_;
    double sigma_2 = mod_->sigma_2;
    double k_2 = mod_->k_2;
    double simga_2_tilde_=sigma_2*sqrt((1.-exp(-2.*k_2*T))/2./k_2);
	return simga_2_tilde_;
}


double Swaption::rho_tilde(){
    double T=T_;
    double sigma_1 = mod_->sigma_1;
    double sigma_2 = mod_->sigma_2;
    double rho = mod_->rho;
    double rho_tilde_ = sigma_1*sigma_2*rho*mod_->int_1_2(0.,T)/sigma_1_tilde ()/sigma_2_tilde ();
    return rho_tilde_;
}


double Swaption::mu_1_tilde(){
    double T=T_;
    double sigma_1 = mod_->sigma_1;
    double sigma_2 = mod_->sigma_2;
    double k_1 = mod_->k_1;
    double k_2 = mod_->k_2;
    double rho = mod_->rho;
    double mu_1_tilde_ = sigma_1*sigma_1*(1.-exp(-2.*k_1*T))/2./k_1/k_1 + sigma_1*sigma_2*rho*mod_->int_1_2(0.,T)/k_2
                         -(sigma_1*sigma_1/k_1+sigma_1*sigma_2*rho/k_2)*mod_->int_1(0.,T);
    return mu_1_tilde_;
}

double Swaption::mu_2_tilde(){
    double T=T_;
    double sigma_1 = mod_->sigma_1;
    double sigma_2 = mod_->sigma_2;
    double k_1 = mod_->k_1;
    double k_2 = mod_->k_2;
    double rho = mod_->rho;
    double mu_2_tilde_ = sigma_2*sigma_2*(1.-exp(-2.*k_2*T))/2./k_2/k_2 + sigma_1*sigma_2*rho*mod_->int_1_2(0.,T)/k_1
                         -(sigma_2*sigma_2/k_2+sigma_1*sigma_2*rho/k_1)*mod_->int_2(0.,T);
    return mu_2_tilde_;
}


double Swaption::c_i(double T_i){
    if(abs(swap_maturity_-T_i)<tenor_){
            return 1+strike_*tenor_;

    }else{
            return strike_*tenor_;
    }
}

double Swaption::A_T_Ti(double T_i,const vector<double> & init_zc,const vector<double> &term_structure){
    double T=T_;
    double sigma_1 = mod_->sigma_1;
    double sigma_2 = mod_->sigma_2;
    double k_1 = mod_->k_1;
    double k_2 = mod_->k_2;
    double rho = mod_->rho;

    double b_1=mod_->int_1(T,T_i);
    double b_2=mod_->int_2(T,T_i);
    double b_1_2=mod_->int_1_2(T,T_i);
    double v_=(T_i-T-b_1-k_1*b_1*b_1/2)*sigma_1*sigma_1/k_1/k_1;
    v_=v_+(T_i-T-b_2-k_2*b_2*b_2/2)*sigma_2*sigma_2/k_2/k_2;
    v_=v_+(T_i-T-b_1-b_2+b_1_2)*2*rho*sigma_1*sigma_2/k_1/k_2;

    double into_phi=0.;
    double dt=(T_i-T)/1000.;
    for (int i=0; i<1000; i++){
        into_phi=into_phi+mod_->phi_calib(init_zc,term_structure,T+dt*double(i))*dt;
    }

    double A_T_Ti_=exp(-into_phi+v_/2.);
    return A_T_Ti_;
}

double Swaption::lambda_i(double x,double T_i,const vector<double> & init_zc,const vector<double> &term_structure){
    double T = T_;
    double lambda_i_=c_i(T_i)*A_T_Ti(T_i,init_zc,term_structure)*exp(-x*mod_->int_1(T,T_i));
    //cout<<"A_T_T_i "<<A_T_Ti(T_i,init_zc,term_structure)<<endl;
//    cout<<"x "<<x<<endl;
//    cout<<"exp "<<exp(200*mod_->int_1(T,T_i))<<endl;
    return lambda_i_;
}

double Swaption::x_bar(double x, const vector<double> & init_zc,const vector<double> &term_structure){
    double tol = 0.001;
    double tol_ = 100;
    double x_up_bound = -1000.;
    double x_down_bound = 1000.;
    double x_change = 0.;

    double x_bar_= 0.;
    double N=(swap_maturity_-T_)/tenor_;

    double sum_x_bar = 0.;
    double sum_x_up = 0.;
    double sum_x_down = 0.;
    double T_i_ = 0.;
    double lambda_i_ = 0.;

    while (abs(tol_)>tol){
        x_bar_ =(x_up_bound+x_down_bound)/2.;
        sum_x_bar = 0.;
        sum_x_down = 0.;
        sum_x_up = 0.;
        for(int i=1;i<N+1;i++){
 //           cout<<"i "<<i<<endl;
            T_i_ = T_ + double(tenor_*i);
 //           cout<<"T_i "<<T_i_<<endl;
            lambda_i_= lambda_i(x,T_i_,init_zc,term_structure);
            sum_x_bar += lambda_i_*exp(-x_bar_*mod_->int_2(T_,T_i_));
            sum_x_up += lambda_i_*exp(-x_up_bound*mod_->int_2(T_,T_i_));
            sum_x_down += lambda_i_*exp(-x_down_bound*mod_->int_2(T_,T_i_));
        }
        tol_ = sum_x_bar-1;

//        cout<<"x_bar "<< x_bar_ << " tol xbar "<< tol_ <<endl;
//
//        cout<<"x up "<< x_up_bound <<" x up tol "<<sum_x_up-1<<endl;
//        cout<<"x down "<<x_down_bound <<" x down tol "<<sum_x_down-1<<endl;

        if((sum_x_bar-1)>0 && (sum_x_up-1) >0 &&(sum_x_down-1)<0 ){
            x_up_bound = x_bar_;
        }else if((sum_x_bar-1)<0 && (sum_x_up-1) >0 &&(sum_x_down-1)<0){
            x_down_bound = x_bar_;
        }else if((sum_x_up-1) <0 &&(sum_x_down-1)>0){
            x_change = x_up_bound;
            x_up_bound = x_down_bound;
            x_down_bound = x_change;
        }
    }
   return x_bar_;
}

double Swaption::h_1 (double x, const vector<double> & init_zc,const vector<double> &term_structure){
    double x_bar_= x_bar(x, init_zc,term_structure);
    double mu_2_tilde_ = mu_2_tilde();
    double sigma_2_tilde_ = sigma_2_tilde();
    double rho_tilde_ = rho_tilde();
    double mu_1_tilde_ = mu_1_tilde();
    double sigma_1_tilde_ = sigma_1_tilde();
    double h_1_ = 0.;
    h_1_ = (x_bar_-mu_2_tilde_)/sqrt(1-rho_tilde_*rho_tilde_)/sigma_2_tilde_;
    h_1_ = h_1_ - rho_tilde_*(x-mu_1_tilde_)/sqrt(1-rho_tilde_*rho_tilde_)/sigma_1_tilde_;
	return h_1_;
}

double Swaption::h_2 (double x,double T_i,const vector<double> & init_zc,const vector<double> &term_structure){
    double h_1_ = h_1(x, init_zc,term_structure);
    double b_2_T_T_i_ = mod_->int_2(T_,T_i);
    double sigma_2_tilde_ = sigma_2_tilde();
    double rho_tilde_ = rho_tilde();
    double h_2_ = 0.;
    h_2_ = h_1_ + b_2_T_T_i_*sigma_2_tilde_*sqrt(1-rho_tilde_*rho_tilde_);
	return h_2_;
}
double Swaption::kappa(double x, double T_i){
    double kappa_=0.;
    kappa_ = mu_2_tilde()-sigma_2_tilde()*sigma_2_tilde()*(1-rho_tilde()*rho_tilde())*mod_->int_2(T_,T_i)/2.+rho_tilde()*sigma_2_tilde()*(x-mu_1_tilde())/sigma_1_tilde();
    kappa_ = -kappa_*mod_->int_2(T_,T_i);
    return kappa_;
}

double Swaption::inte_f(double x,const vector<double> & init_zc,const vector<double> &term_structure){
    double h_1_ = h_1(x, init_zc,term_structure);
//    cout<<"h_1"<<endl;
//    cout<<h_1_<<endl;
    double term_1_ = cdfGaussian (-h_1_);
    double N = (swap_maturity_-T_)/tenor_;
    double term_2_ = 0.;
    double T_i = 0.;
    double h_2_ = 0.;
    for(int i=1;i<N+1;i++){
        T_i = T_ + double(tenor_*i);
        h_2_ = h_2(x,T_i,init_zc,term_structure);
//        cout<<"h_2"<<endl;
//        cout<<h_2_<<endl;
        term_2_ += lambda_i(x,T_i,init_zc,term_structure)*exp(kappa(x,T_i))*cdfGaussian (-h_2_);
//        cout<<"lambda "<<lambda_i(x,T_i,init_zc,term_structure)<<endl;
        //cout<<exp(-kappa(x,T_i))<<endl;
    }
    double inte_f_ = 0.;
//    cout<<"term_1"<<endl;
//    cout<<term_1_<<endl;
//    cout<<"term_2"<<endl;
//    cout<<term_2_<<endl;
    inte_f_ = (term_1_-term_2_)*exp(-(x-mu_1_tilde())*(x-mu_1_tilde())/2./sigma_1_tilde()/sigma_1_tilde())/sigma_1_tilde()/sqrt(2*M_PI);
//    cout<<"term_1_-term_2_"<<term_1_-term_2_<<endl;
    return inte_f_;
}


double Swaption::price_theory_Swaption(const vector<double> &term_structure,const vector<double> & init_zc){
   vector<double> vect(10);
   cout<<"T_"<<endl;
   cout<<T_<<endl;
   cout<<"swap_maturity_"<<endl;
   cout<<swap_maturity_<<endl;
   cout<<"tenor_"<<endl;
   cout<<tenor_<<endl;

   double P_t0_Tn=mod_->zero_coupon_bonds_g2(T_, swap_maturity_, vect, vect, term_structure, init_zc);
//   cout<<"P_t0_Tn"<<endl;
//   cout<<P_t0_Tn<<endl;

//   double x = -200.;
//   double integral_= inte_f(x,init_zc,term_structure);
//   cout<<"integral"<<endl;
//   cout<<integral_<<endl;

   double integral_ =0.;
   double x_down = -10.;
   double x_up = 10.;
   double dx_ = 0.2;
   double N = (x_up - x_down)/dx_;
   double x = 0.;
   for(int i=0;i<N+1;i++){
//        cout<<i<<endl;
        x = x_down + double(i*dx_);
        integral_ +=  inte_f(x,init_zc,term_structure)*dx_;
//        cout<<"integral"<<endl;
//        cout<<"integral_f"<<inte_f(x,init_zc,term_structure)<<endl;
    }
//   cout<<"integral"<<endl;
//   cout<<integral_<<endl;

   double price_=P_t0_Tn*integral_;
   return price_;
}

Swaption::~Swaption() {

}
