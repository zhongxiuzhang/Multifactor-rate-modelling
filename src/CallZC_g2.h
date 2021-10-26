#ifndef CALLZC_G2_H_INCLUDED
#define CALLZC_G2_H_INCLUDED

#include "product.h"
#include "G2++.h"
#include <vector>
using namespace std;

class CallZC : public product {
    public:

	double strike_;

	double zc_maturity_;

	double payOff(const vector<double> &term_structure, const vector<double> &brownian_motion_1, const vector<double> &brownian_motion_2, const vector<double> &short_rate_path,const vector<double> & init_zc);

	CallZC(g2_model *g2mod, double zc_maturity, double maturity, double strike);

    double price_theory_CallZC(const vector<double> & init_zc,const vector<double> & term_s);

	virtual ~CallZC();
};

#endif // CALLZC_G2_H_INCLUDED
