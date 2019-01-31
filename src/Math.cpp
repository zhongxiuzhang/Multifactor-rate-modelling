#include "Math.h"
#include <cmath>
#include <stdlib.h>

using namespace std;

double cdfGaussian (double quantile){
	
	return erfc(-quantile / sqrt(2)) / 2.;
}
