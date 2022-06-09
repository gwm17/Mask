/*

Target.cpp
A basic target unit for use in the SPANCRedux environment. A target
is defined as a single compound with elements Z,A of a given stoichiometry 
Holds an energy loss class

Based on code by D.W. Visser written at Yale for the original SPANC

Written by G.W. McCann Aug. 2020

*/
#include "Target.h"
#include "catima/nucdata.h"

namespace Mask {

	/*Targets must be of known thickness*/
	Target::Target(double thick) {
		thickness = thick;
		thickness_gcm2 = thickness*1.0e-6;
	}
	
	Target::~Target() {}
	
	/*Set target elements of given Z, A, S*/
	void Target::SetElements(std::vector<int>& z, std::vector<int>& a, std::vector<int>& stoich) {
		Z = z;
		A = a;
		Stoich = stoich;
		MassLookup& masses = MassLookup::GetInstance();
		eloss.SetTargetComponents(Z, A, Stoich);
		for(size_t i=0; i<Z.size(); i++)
		{
			target_material.add_element(masses.FindMassU(Z[i], A[i]), Z[i], Stoich[i]);
		}
	}
	
	/*Element verification*/
	bool Target::ContainsElement(int z, int a) {
		for(unsigned int i=0; i<Z.size(); i++)
			if( z == Z[i] && a == A[i]) 
				return true;
		return false;
	}
	
	/*Calculates energy loss for travelling all the way through the target*/
	double Target::GetEnergyLossTotal(int zp, int ap, double startEnergy, double theta) {
		if(theta == M_PI/2.) 
			return startEnergy;
		else if (theta > M_PI/2.) 
			theta = M_PI - theta;
		catima::Projectile proj(MassLookup::GetInstance().FindMassU(zp, ap), zp, 0.0, 0.0);
		proj.T = startEnergy/proj.A;
		target_material.thickness(thickness_gcm2/fabs(cos(theta)));
		return catima::integrate_energyloss(proj, target_material);
		//return eloss.GetEnergyLoss(zp, ap, startEnergy, thickness/fabs(cos(theta)));
	}

	/*Calculates the energy loss for traveling some fraction through the target*/
	double Target::GetEnergyLossFractionalDepth(int zp, int ap, double finalEnergy, double theta, double percent_depth)
	{
		if(theta == M_PI/2.)
			return finalEnergy;
		else if (theta > M_PI/2.)
			theta = M_PI-theta;

		catima::Projectile proj(MassLookup::GetInstance().FindMassU(zp, ap), zp, 0.0, 0.0);
		proj.T = finalEnergy/proj.A;
		target_material.thickness(thickness_gcm2*percent_depth/fabs(cos(theta)));
		return catima::integrate_energyloss(proj, target_material);
		//return eloss.GetEnergyLoss(zp, ap, finalEnergy, thickness*percent_depth/(std::fabs(std::cos(theta))));
	}
	
	/*Calculates reverse energy loss for travelling all the way through the target*/
	double Target::GetReverseEnergyLossTotal(int zp, int ap, double finalEnergy, double theta) {
		if(theta == M_PI/2.) 
			return finalEnergy;
		else if (theta > M_PI/2.) 
			theta = M_PI - theta;

		catima::Projectile proj(MassLookup::GetInstance().FindMassU(zp, ap), zp, 0.0, 0.0);
		proj.T = finalEnergy/proj.A;
		target_material.thickness(thickness_gcm2/fabs(cos(theta)));
		return catima::reverse_integrate_energyloss(proj, target_material);
		//return eloss.GetReverseEnergyLoss(zp, ap, finalEnergy, thickness/fabs(cos(theta)));
	}

	/*Calculates the reverse energy loss for traveling some fraction through the target*/
	double Target::GetReverseEnergyLossFractionalDepth(int zp, int ap, double finalEnergy, double theta, double percent_depth)
	{
		if(theta == M_PI/2.)
			return finalEnergy;
		else if (theta > M_PI/2.)
			theta = M_PI-theta;
		catima::Projectile proj(MassLookup::GetInstance().FindMassU(zp, ap), zp, 0.0, 0.0);
		proj.T = finalEnergy/proj.A;
		target_material.thickness(thickness_gcm2*percent_depth/fabs(cos(theta)));
		return catima::reverse_integrate_energyloss(proj, target_material);
		//return eloss.GetReverseEnergyLoss(zp, ap, finalEnergy, thickness*percent_depth/(std::fabs(std::cos(theta))));
	}

}
