#include "TwoStepSystem.h"
#include "RandomGenerator.h"
#include "KinematicsExceptions.h"

namespace Mask {

	TwoStepSystem::TwoStepSystem() :
		ReactionSystem(), m_phi2Range(0, 2.0*M_PI)
	{
		nuclei.resize(6);
	}
	
	TwoStepSystem::TwoStepSystem(std::vector<int>& z, std::vector<int>& a) :
		ReactionSystem(), m_phi2Range(0, 2.0*M_PI)
	{
		nuclei.resize(6);
		SetNuclei(z, a);
	}
	
	TwoStepSystem::~TwoStepSystem() {
	
	}
	
	bool TwoStepSystem::SetNuclei(std::vector<int>&z, std::vector<int>& a) {
		if(z.size() != a.size() || z.size() != 4) {
			return false;
		}
		int zr = z[0] + z[1] - z[2];
		int ar = a[0] + a[1] - a[2];
	
		step1.SetNuclei(z[0], a[0], z[1], a[1], z[2], a[2]);
		step2.SetNuclei(zr, ar, 0, 0, z[3], a[3]);
		SetSystemEquation();
		return true;
	}

	const std::vector<Nucleus>& TwoStepSystem::GetNuclei()
	{
		nuclei[0] = step1.GetTarget();
		nuclei[1] = step1.GetProjectile();
		nuclei[2] = step1.GetEjectile();
		nuclei[3] = step1.GetResidual();
		nuclei[4] = step2.GetEjectile();
		nuclei[5] = step2.GetResidual();

		return nuclei;
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
		//Link up the target if it hasn't been done yet
		if(!target_set_flag) {
			LinkTarget();
		}
	
		//Sample parameters
		double bke = (*m_beamDist)(RandomGenerator::GetInstance().GetGenerator());
		double rxnTheta = acos((*m_theta1Range)(RandomGenerator::GetInstance().GetGenerator()));
		double rxnPhi = (*m_phi1Range)(RandomGenerator::GetInstance().GetGenerator());
		double decay1costheta = decay1dist.GetRandomCosTheta();
		double decay1Theta = std::acos(decay1costheta);
		double decay1Phi = m_phi2Range(RandomGenerator::GetInstance().GetGenerator());
		double residEx = (*m_exDist)(RandomGenerator::GetInstance().GetGenerator());
	
		step1.SetBeamKE(bke);
		step1.SetPolarRxnAngle(rxnTheta);
		step1.SetAzimRxnAngle(rxnPhi);
		step1.SetExcitation(residEx);
	
		step2.SetPolarRxnAngle(decay1Theta);
		step2.SetAzimRxnAngle(decay1Phi);
		
		step1.Calculate();
	
		step2.SetTarget(step1.GetResidual());
		if(decay1costheta == -10) {
			step2.ResetEjectile();
			step2.ResetResidual();
			return;
		}
		step2.TurnOnResidualEloss();
		step2.Calculate();
	
	
	}

}
