#include "LegendrePoly.h"
#include <cmath>

double P_l(int l, double x) {
	if(l == 0) {
		return 1.0;
	} else if (l == 1) {
		return x;
	} else {
		return (2.0*l - 1.0)/l*x*P_l(l-1, x) - (l-1.0)/l*P_l(l-2, x);
	}
}

double Normed_P_l_sq(int l, double x) {
	return (2.0*l+1.0)/2.0*std::pow(P_l(l, x), 2.0);
}

double P_0(double x) {
	return 1.0;
}

double P_1(double x) {
	return x;
}

double P_2(double x) {
	return 0.5*(3.0*x*x -1.0);
}

double P_l_ROOT(double* x, double* pars) {
	return P_l(pars[0], x[0]);
}