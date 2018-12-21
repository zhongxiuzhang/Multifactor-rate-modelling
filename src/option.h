#ifndef OPTION_H
#define	OPTION_H
#include <vector>
using namespace std;

class Option
{
public:
	int T_; //Maturity
	virtual ~Option() {}
	
	virtual double payOff(const vector<double> &short_rate_path) = 0;
};


#endif /* OPTION_H */

