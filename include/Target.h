/*

Target.h
A basic target unit for use in the Mask environment. A target
is defined as a single compound with elements Z,A of a given stoichiometry
Holds an energy loss class

Based on code by D.W. Visser written at Yale for the original SPANC

Written by G.W. McCann Aug. 2020

*/
#ifndef TARGET_H
#define TARGET_H

#include <string>
#include <vector>
#include <cmath>
#include "EnergyLoss.h"

namespace Mask {

	class Target {
	
	public:
	 	Target(double thick);
	 	~Target();
	 	void SetElements(std::vector<int>& z, std::vector<int>& a, std::vector<int>& stoich);
	 	bool ContainsElement(int z, int a);
	 	double GetEnergyLossTotal(int zp, int ap, double startEnergy, double angle);
	 	double GetReverseEnergyLossTotal(int zp, int ap, double finalEnergy, double angle);
	 	double GetEnergyLossFractionalDepth(int zp, int ap, double startEnergy, double angle, double percent_depth);
	 	double GetReverseEnergyLossFractionalDepth(int zp, int ap, double finalEnergy, double angle, double percent_depth);
	 	inline const double& GetThickness() { return thickness; }
	 	inline int GetNumberOfElements() { return Z.size(); }
	 	inline int GetElementZ(int index) { return Z[index]; }
	 	inline int GetElementA(int index) { return A[index]; }
	 	inline int GetElementStoich(int index) { return Stoich[index]; }
	
	private:
		EnergyLoss eloss;
		double thickness;
		std::vector<int> Z, A, Stoich;
	
	};

}

#endif 
