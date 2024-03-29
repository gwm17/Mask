#ifndef LEGENDREPOLY_H
#define LEGENDREPOLY_H

namespace Mask {

	double P_l(int l, double x);
	double P_l_ROOT(double* x, double* pars); //Only for use with ROOT cern in plotting
	double Normed_P_l_sq(int l, double x);
	
	double P_0(double x);
	double P_1(double x);
	double P_2(double x);

}

#endif