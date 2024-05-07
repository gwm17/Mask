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
#include "catima/gwm_integrators.h"
#include "MassLookup.h"

namespace Mask {

	class Target {
	
	public:
	 	Target(const std::vector<int>& z, const std::vector<int>& a, const std::vector<int>& stoich, double thick);
	 	~Target();
	 	bool ContainsElement(int z, int a);
	 	double GetEnergyLossTotal(int zp, int ap, double startEnergy, double angle);
	 	double GetReverseEnergyLossTotal(int zp, int ap, double finalEnergy, double angle);
	 	double GetEnergyLossFractionalDepth(int zp, int ap, double startEnergy, double angle, double percent_depth);
	 	double GetReverseEnergyLossFractionalDepth(int zp, int ap, double finalEnergy, double angle, double percent_depth);
	 	const double& GetThickness() const { return m_thickness; }
	 	int GetNumberOfElements() const { return m_Z.size(); }
	 	int GetElementZ(int index) const { return m_Z[index]; }
	 	int GetElementA(int index) const { return m_A[index]; }
	 	int GetElementStoich(int index) const { return m_stoich[index]; }

	private:
	 	void Init(const std::vector<int>& z, const std::vector<int>& a, const std::vector<int>& stoich);
		
		catima::Material m_material;
		double m_thickness;
		double m_thickness_gcm2;
		std::vector<int> m_Z, m_A, m_stoich;
	
	};

}

#endif 
