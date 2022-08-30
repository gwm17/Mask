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
	Target::Target(const std::vector<int>& z, const std::vector<int>& a, const std::vector<int>& stoich, double thick) :
		m_thickness(thick), m_thickness_gcm2(thick*1.0e-6)
	{
		Init(z, a, stoich);
	}
	
	Target::~Target() {}
	
	/*Set target elements of given Z, A, S*/
	void Target::Init(const std::vector<int>& z, const std::vector<int>& a, const std::vector<int>& stoich) 
	{
		m_Z = z;
		m_A = a;
		m_stoich = stoich;
		MassLookup& masses = MassLookup::GetInstance();
		for(size_t i=0; i<m_Z.size(); i++)
		{
			m_material.add_element(masses.FindMassU(m_Z[i], m_A[i]), m_Z[i], m_stoich[i]);
		}
	}
	
	/*Element verification*/
	bool Target::ContainsElement(int z, int a)
	{
		for(std::size_t i=0; i<m_Z.size(); i++)
		{
			if( z == m_Z[i] && a == m_A[i]) 
				return true;
		}
		return false;
	}
	
	/*Calculates energy loss for travelling all the way through the target*/
	double Target::GetEnergyLossTotal(int zp, int ap, double startEnergy, double theta)
	{
		if(theta == M_PI/2.) 
			return startEnergy;
		else if (theta > M_PI/2.) 
			theta = M_PI - theta;

		catima::Projectile proj(MassLookup::GetInstance().FindMassU(zp, ap), zp, 0.0, 0.0);
		proj.T = startEnergy/proj.A;
		m_material.thickness(m_thickness_gcm2/fabs(cos(theta)));
		return catima::integrate_energyloss(proj, m_material);
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
		m_material.thickness(m_thickness_gcm2*percent_depth/fabs(cos(theta)));
		return catima::integrate_energyloss(proj, m_material);
	}
	
	/*Calculates reverse energy loss for travelling all the way through the target*/
	double Target::GetReverseEnergyLossTotal(int zp, int ap, double finalEnergy, double theta)
	{
		if(theta == M_PI/2.) 
			return finalEnergy;
		else if (theta > M_PI/2.) 
			theta = M_PI - theta;

		catima::Projectile proj(MassLookup::GetInstance().FindMassU(zp, ap), zp, 0.0, 0.0);
		proj.T = finalEnergy/proj.A;
		m_material.thickness(m_thickness_gcm2/fabs(cos(theta)));
		return catima::reverse_integrate_energyloss(proj, m_material);
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
		m_material.thickness(m_thickness_gcm2*percent_depth/fabs(cos(theta)));
		return catima::reverse_integrate_energyloss(proj, m_material);
	}

}
