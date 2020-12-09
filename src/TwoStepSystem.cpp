#include "TwoStepSystem.h"
#include "KinematicsExceptions.h"

TwoStepSystem::TwoStepSystem() :
	ReactionSystem()
{
}

TwoStepSystem::TwoStepSystem(std::vector<int>& z, std::vector<int>& a) :
	ReactionSystem()
{
	SetNuclei(z, a);
}

TwoStepSystem::~TwoStepSystem() {

}

bool TwoStepSystem::SetNuclei(std::vector<int>&z, std::vector<int>& a) {
	if(z.size() != a.size() || z.size() < 4) {
		return false;
	}
	int zr = z[0] + z[1] - z[2];
	int ar = a[0] + a[1] - a[2];

	step1.SetNuclei(z[0], a[0], z[1], a[1], z[2], a[2]);
	step2.SetNuclei(zr, ar, 0, 0, z[3], a[3]);
	SetSystemEquation();
	return true;
}

void TwoStepSystem::LinkTarget() {
	step1.SetLayeredTarget(&target);
	step2.SetLayeredTarget(&target);

	rxnLayer = target.FindLayerContaining(step1.GetTarget().GetZ(), step1.GetTarget().GetA());
	if(rxnLayer != -1) {
		step1.SetRxnLayer(rxnLayer);
		step2.SetRxnLayer(rxnLayer);
		target_set_flag = true;
	} else {
		throw ReactionLayerException();
	}
}

void TwoStepSystem::SetSystemEquation() {
	m_sys_equation = step1.GetTarget().GetIsotopicSymbol();
	m_sys_equation += "(";
	m_sys_equation += step1.GetProjectile().GetIsotopicSymbol();
	m_sys_equation += ", ";
	m_sys_equation += step1.GetEjectile().GetIsotopicSymbol();
	m_sys_equation += ")";
	m_sys_equation += step1.GetResidual().GetIsotopicSymbol();
	m_sys_equation += "-> ";
	m_sys_equation += step2.GetEjectile().GetIsotopicSymbol();
	m_sys_equation += "+";
	m_sys_equation += step2.GetResidual().GetIsotopicSymbol();
}

void TwoStepSystem::RunSystem() {
	if(!gen_set_flag) return;
	
	//Link up the target if it hasn't been done yet
	if(!target_set_flag) {
		LinkTarget();
	}

	//Sample parameters
	double bke = generator->Gaus(m_beamDist.first, m_beamDist.second);
	double rxnTheta = acos(generator->Uniform(cos(m_theta1Range.first), cos(m_theta1Range.second)));
	double rxnPhi = 0;
	double decay1Theta = acos(generator->Uniform(-1, 1));
	double decay1Phi = generator->Uniform(0, M_PI*2.0);
	double residEx = generator->Gaus(m_exDist.first, m_exDist.second);

	step1.SetBeamKE(bke);
	step1.SetPolarRxnAngle(rxnTheta);
	step1.SetAzimRxnAngle(rxnPhi);
	step1.SetExcitation(residEx);

	step2.SetPolarRxnAngle(decay1Theta);
	step2.SetAzimRxnAngle(decay1Phi);
	
	step1.Calculate();

	step2.SetTarget(step1.GetResidual());
	step2.TurnOnResidualEloss();
	step2.Calculate();


}