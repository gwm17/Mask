#include "LegendrePoly.h"

double P_l(int l, double x) {
	if(l == 0) {
		return 1.0;
	} else if (l == 1) {
		return x;
	} else {
		return ((2.0*l + 1.0)*x*P_l(l-1, x) - (l-1.0)*P_l(l-2, x))/(double(l));
	}
}