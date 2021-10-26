#ifndef G2++_H_INCLUDED
#define G2++_H_INCLUDED
#include <vector>

using namespace std;
class g2_model{
    public:
	// For random Gaussian generation
	// Note that there are many ways of doing this, but we will
	// be using the Box-Muller method for demonstration purposes
	double k_1; double k_2;
	double sigma_1; double sigma_2;
    double rho;
    double gaussian_box_muller();

    g2_model(double k_1_, double k_2_, double sigma_1_, double sigma_2_,double rho_);

	vector<double>  vector_movement_brownien_1(const vector<double> & term_structure);
    vector<double>  vector_movement_brownien_2(const vector<double> & mov_brow_1,const vector<double> & term_structure);
    double x_1(double t, const vector<double> & mov_brow_1, const vector<double> & term_structure);
    double x_2(double t, const vector<double> & mov_brow_2, const vector<double> & term_structure);
    //double phi(double t);
    double phi_calib(const vector<double> &init_zc,const vector<double> &term_structure,double t);
    double int_1(double t,double T);
    double int_2(double t,double T);
    double int_1_2(double t,double T);
    vector<double> short_rate_g2(const vector<double> & mov_brow_1,const vector<double> & mov_brow_2,const vector<double> & term_structure,const vector<double> & init_zc);
    double zero_coupon_bonds_g2(double t, double T, const vector<double> & mov_brow_1, const vector<double> & mov_brow_2, const vector<double> & term_structure,const vector<double> & init_zc);
    double forward_rate_g2(double t, double T,const vector<double> & mov_brow_1, const vector<double> & mov_brow_2, const vector<double> & term_structure,const vector<double> & init_zc);

    //vector<double> short_rate_g2(const vector<double> & mov_brow_1,const vector<double> & mov_brow_2,const vector<double> & term_structure);
    //double zero_coupon_bonds_g2(double t, double T, const vector<double> & mov_brow_1, const vector<double> & mov_brow_2, const vector<double> & term_structure);
    //double forward_rate_g2(double t, double T,const vector<double> & mov_brow_1, const vector<double> & mov_brow_2, const vector<double> & term_structure);

	// To generate the term structure of Hull White model. The tenor_num denote the number of points in one year.
	vector<double> term_structure_g2(const int& last_year, const int& step);

};

vector<double> copy_term_structure(vector<double> const &term_structure);

#endif // G2++_H_INCLUDED
