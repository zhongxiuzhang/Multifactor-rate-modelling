#include <iostream>
#include <math.h>
#include <vector>
#include "Math.h"

#include "G2++.h"
#include "MonteCarlo_g2.h"
#include "CallZC_g2.h"
#include "Swaption_g2.h"

using namespace std;

int main()
{
    g2_model *g2test=new g2_model(0.1,0.2,0.2,0.2,-0.5);

    vector<double> term_structure_ = g2test->term_structure_g2(5, 4);

    double a[21] = {1.0,1.0005,1.001,1.0017,1.0025,1.00248,1.00245,1.00242,1.0024,1.0018,1.001,1.0006,0.9998,0.998,0.9965,0.995,0.9945,0.9925,0.9905,0.988,0.9867};
    vector<double> init_zc_;

    int N = 21;
	for (int j = 0; j<N; j++) {
		init_zc_.push_back(a[j]);
	}


    cout<<term_structure_.back()<<endl;
    cout<<term_structure_[0]<<endl;
    cout<<term_structure_[1]<<endl;

    vector<double> mov_brow_1 = g2test->vector_movement_brownien_1(term_structure_);
    vector<double> mov_brow_2 = g2test->vector_movement_brownien_2(mov_brow_1,term_structure_);
    cout<<mov_brow_1.back()<<endl;
    cout<<mov_brow_2.back()<<endl;

//    double x_1_=g2test->x_1(4,mov_brow_1,term_structure_);
//    double x_2_=g2test->x_2(3,mov_brow_2,term_structure_);
//    cout<<x_1_<<endl;
//    cout<<x_2_<<endl;

    cout<<"cdfGaussian"<<cdfGaussian(2)<<endl;


//    double zero_coupon_=g2test->zero_coupon_bonds_g2(1,5,mov_brow_1,mov_brow_2,term_structure_);
//    cout<<" ZC "<<zero_coupon_<<endl;

    vector<double> short_rate_=g2test->short_rate_g2(mov_brow_1,mov_brow_2,term_structure_,init_zc_);

//    cout<<" short rate " <<endl;
//    for (int i=0; i<short_rate_.size();i++){
//            cout<< short_rate_[i]<<endl;
//    }

//    double forward_rate_=g2test->forward_rate_g2(0,8,mov_brow_1,mov_brow_2,term_structure_,init_zc_);
//    cout<<" forward rate " <<endl;
//    cout<<forward_rate_<<endl;

    cout<<"end"<<endl;

    double zero_coupon_=g2test->zero_coupon_bonds_g2(2.,5.,mov_brow_1,mov_brow_2,term_structure_,init_zc_);
    cout<<"Discount ZC "<<zero_coupon_<<endl;

    CallZC *CallZC_=new CallZC(g2test,5.,1.,0.5);

    cout<<"CallZC"<<endl;
    cout<< CallZC_->price_theory_CallZC(init_zc_,term_structure_)<<endl;
    MonteCarlo *MonteCarlo_=new MonteCarlo(g2test,CallZC_,2000,5,4);
    cout<<"MonteCarlo"<<endl;
    cout<< MonteCarlo_->price(init_zc_)<<endl;

//    Swaption *Swaption_=new Swaption(g2test, 2, 1,0.01,1);
//    cout<<"Swaption"<<endl;
//    MonteCarlo *MonteCarlo_swaption=new MonteCarlo(g2test,Swaption_,2000,5,4);
//    cout<<"MonteCarlo"<<endl;
//    cout<< MonteCarlo_swaption->price(init_zc_)<<endl;
//    cout<< Swaption_->x_bar(1,init_zc_,term_structure_)<<endl;
//    cout<<"Swaption theory price"<<endl;
//    cout<< Swaption_->price_theory_Swaption(term_structure_,init_zc_)<<endl;

    return 0;
}
