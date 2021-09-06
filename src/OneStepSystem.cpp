#include "OneStepSystem.h"

namespace Mask {

	OneStepSystem::OneStepSystem() :
		ReactionSystem()
	{
	}
	
	OneStepSystem::OneStepSystem(std::vector<int>& z, std::vector<int>& a) :
		ReactionSystem()
	{
		SetNuclei(z, a);
	}
	
	OneStepSystem::~OneStepSystem() {}
	
	bool OneStepSystem::SetNuclei(std::vector<int>& z, std::vector<int>& a) {
		if(z.size() != a.size() || z.size() != 3) {
			return false;
		}
	
		step1.SetNuclei(z[0], a[0], z[1], a[1], z[2], a[2]);
		SetSystemEquation();
		return true;
	}
	
	void OneStepSystem::LinkTarget() {
		step1.SetLayeredTarget(&target);
	
		rxnLayer = target.FindLayerContaining(step1.GetTarget().GetZ(), step1.GetTarget().GetA());
		if(rxnLayer != -1) {
			step1.SetRxnLayer(rxnLayer);
			target_set_flag = true;
		} else {
			throw ReactionLayerException();
		}
	}
	
	void OneStepSystem::SetSystemEquation() {
		m_sys_equation = step1.GetTarget().GetIsotopicSymbol();
		m_sys_equation += "(";
		m_sys_equation += step1.GetProjectile().GetIsotopicSymbol();
		m_sys_equation += ", ";
		m_sys_equation += step1.GetEjectile().GetIsotopicSymbol();
		m_sys_equation += ")";
		m_sys_equation += step1.GetResidual().GetIsotopicSymbol();
	}
	
	void OneStepSystem::RunSystem() {
		if(!gen_set_flag) return;
		
		//Link up the target if it hasn't been done yet
		if(!target_set_flag) {
			LinkTarget();
		}
	
		//Sample parameters
		double bke = (*m_beamDist)(*generator);
		double rxnTheta = std::acos((*m_theta1Range)(*generator));
		double rxnPhi = (*m_phi1Range)(*generator);
		double residEx = (*m_exDist)(*generator);
		
		step1.SetBeamKE(bke);
		step1.SetPolarRxnAngle(rxnTheta);
		step1.SetAzimRxnAngle(rxnPhi);
		step1.SetExcitation(residEx);
		
		step1.TurnOnResidualEloss();
		step1.Calculate();
	}

}