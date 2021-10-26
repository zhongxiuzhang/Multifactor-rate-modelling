#ifndef MONTECARLO_G2_H_INCLUDED
#define MONTECARLO_G2_H_INCLUDED

#include <iostream>
#include "product.h"
#include "G2++.h"
using namespace std;

/**
* @brief MonteCarlo class
*/

class MonteCarlo {

public:
	g2_model* mod_; /*! pointer to the model */
	product *prod_; /*! pointer to the product */
	int nbSamples_; /*! number of draws Monte Carlo */
	int step_;
	double maturity_;

	double price(const vector<double> & init_zc);

	virtual ~MonteCarlo();

	MonteCarlo(g2_model *mod, product *prod, int nbSamples, double maturity, int step);

};

#endif // MONTECARLO_G2_H_INCLUDED
