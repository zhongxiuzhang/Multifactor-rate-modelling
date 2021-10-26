#ifndef PRODUCT_H_INCLUDED
#define PRODUCT_H_INCLUDED

#include <vector>
#include "G2++.h"
using namespace std;

class product
{
public:
	int T_; //Maturity
	g2_model *mod_;
	virtual ~product() {}

	virtual double payOff(const vector<double> &term_structure, const vector<double> &brownian_motion_1,const vector<double> &brownian_motion_2, const vector<double> &short_rate_path,const vector<double> & init_zc) = 0;
};


#endif // PRODUCT_H_INCLUDED
