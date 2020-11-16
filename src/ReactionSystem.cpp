#include "ReactionSystem.h"
#include "KinematicsExceptions.h"
#include <iostream>

ReactionSystem::ReactionSystem() :
	m_beamDist(0,0), m_theta1Range(0,0), m_exDist(0,0), generator(nullptr), target_set_flag(false), gen_set_flag(false), rxnLayer(0), m_sys_equation("")
{
}

ReactionSystem::ReactionSystem(std::vector<int>& z, std::vector<int>& a) :
	m_beamDist(0,0), m_theta1Range(0,0), m_exDist(0,0), generator(nullptr), target_set_flag(false), gen_set_flag(false), rxnLayer(0), m_sys_equation("")
{
	SetNuclei(z, a);
}

ReactionSystem::~ReactionSystem() {

}

bool ReactionSystem::SetNuclei(std::vector<int>&z, std::vector<int>& a) {
	if(z.size() != a.size() || z.size() < 3) {
		return false;
	}

	step1.SetNuclei(z[0], a[0], z[1], a[1], z[2], a[2]);
	SetSystemEquation();
	return true;
}

void ReactionSystem::AddTargetLayer(std::vector<int>& zt, std::vector<int>& at, std::vector<int>& stoich, double thickness) {
	target.AddLayer(zt, at, stoich, thickness);
}

void ReactionSystem::LinkTarget() {
	step1.SetLayeredTarget(&target);

	rxnLayer = target.FindLayerContaining(step1.GetTarget().GetZ(), step1.GetTarget().GetA());
	if(rxnLayer != -1) {
		step1.SetRxnLayer(rxnLayer);
		target_set_flag = true;
	} else {
		throw ReactionLayerException();
	}
}

void ReactionSystem::SetSystemEquation()  {
	m_sys_equation = step1.GetTarget().GetIsotopicSymbol();
	m_sys_equation += "(";
	m_sys_equation += step1.GetProjectile().GetIsotopicSymbol();
	m_sys_equation += ", ";
	m_sys_equation += step1.GetEjectile().GetIsotopicSymbol();
	m_sys_equation += ")";
	m_sys_equation += step1.GetResidual().GetIsotopicSymbol();
}

void ReactionSystem::RunSystem() {
	if(!gen_set_flag) return;
	
	//Link up the target if it hasn't been done yet
	if(!target_set_flag) {
		LinkTarget();
	}

	if(step1.IsDecay()) {
		double rxnTheta = acos(generator->Uniform(-1, 1));
		double rxnPhi = generator->Uniform(0, 2.0*M_PI);
		step1.SetPolarRxnAngle(rxnTheta);
		step1.SetAzimRxnAngle(rxnPhi);

		step1.TurnOnResidualEloss();
		step1.Calculate();
	} else {
		//Sample parameters
		double bke = generator->Gaus(m_beamDist.first, m_beamDist.second);
		double rxnTheta = acos(generator->Uniform(cos(m_theta1Range.first), cos(m_theta1Range.second)));
		double rxnPhi = 0;
		double residEx = generator->Gaus(m_exDist.first, m_exDist.second);
	
		step1.SetBeamKE(bke);
		step1.SetPolarRxnAngle(rxnTheta);
		step1.SetAzimRxnAngle(rxnPhi);
		step1.SetExcitation(residEx);
	
		step1.TurnOnResidualEloss();
		step1.Calculate();
	}

}