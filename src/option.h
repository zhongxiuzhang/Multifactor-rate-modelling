#ifndef OPTION_H
#define	OPTION_H
#include <vector>
#include "HW1.h" 
using namespace std;

class Option
{
public:
	int T_; //Maturity
	HW1 *mod_; 
	virtual ~Option() {}
	
	virtual double payOff(const vector<double> &term_structure, const vector<double> &brownian_motion, const vector<double> &short_rate_path) = 0;
};


#endif /* OPTION_H */

