#include "Math.h"
#include <cmath>
#include <stdlib.h>

using namespace std;

double cdfGaussian (double quantile){
	double proba;
	double pi = 3.141592657;
	double b0 = 0.2316419;
    double b1 = 0.319381530;
    double b2 = -0.356563782;
    double b3 = 1.781477937;
    double b4 = -1.821255978;
    double b5 = 1.330274429;
	double t1 = 1./(1.+b0*quantile);

	proba = 1. - (1./sqrt(2.*pi))*exp(-0.5*pow(quantile,2))*(b1*t1+b2*pow(t1,2)+b3*pow(t1,3)+b4*pow(t1,4)+b5*pow(t1,5));
    return proba;

}
