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
	 	inline const double& GetThickness() { return m_thickness; }
	 	inline int GetNumberOfElements() { return m_Z.size(); }
	 	inline int GetElementZ(int index) { return m_Z[index]; }
	 	inline int GetElementA(int index) { return m_A[index]; }
	 	inline int GetElementStoich(int index) { return m_stoich[index]; }
	
	private:
	 	void Init(const std::vector<int>& z, const std::vector<int>& a, const std::vector<int>& stoich);
		
		catima::Material m_material;
		double m_thickness;
		double m_thickness_gcm2;
		std::vector<int> m_Z, m_A, m_stoich;
	
	};

}

#endif 
