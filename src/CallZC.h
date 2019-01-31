#include "Option.h"
#include "HW1.h"
#include <vector>
using namespace std;

class CallZC : public Option {
public:

	double strike_;

	double zc_maturity_;

	double payOff(const vector<double> &term_structure, const vector<double> &brownian_motion, const vector<double> &short_rate_path);

	CallZC(HW1 *hjm, double zc_maturity, int Maturity, double strike);

	double price_theory_CallZC(const vector<double> &term_structure, const vector<double> & initial_zero_coupon_);

	virtual ~CallZC();
};
