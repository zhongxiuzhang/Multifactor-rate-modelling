#ifndef MONTECARLO_H
#define	MONTECARLO_H

#include <iostream>
#include "Option.h"
#include "HW1.h"
using namespace std;

/**
* @brief MonteCarlo class
*/

class MonteCarlo {

public:
	HW1* mod_; /*! pointer to the model */
	Option *opt_; /*! pointer to the option */
	int nbSamples_; /*! number of draws Monte Carlo */
	int step_;
	double maturity_;

	double price();

	virtual ~MonteCarlo();

	MonteCarlo(HW1 *mod, Option *opt, int nbSamples, double maturity, int step);

					   
};

#endif /* MONTECARLO_H */
