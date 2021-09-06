#include "DecaySystem.h"

namespace Mask {

	DecaySystem::DecaySystem() :
		ReactionSystem()
	{
	}
	
	DecaySystem::DecaySystem(std::vector<int>& z, std::vector<int>& a) :
		ReactionSystem()
	{
		SetNuclei(z, a);
	}
	
	DecaySystem::~DecaySystem() {}
	
	void DecaySystem::SetRandomGenerator(std::mt19937* gen) {
		generator = gen;
		decay1dist.AttachRandomNumberGenerator(gen);
		gen_set_flag = true;
	}
	
	bool DecaySystem::SetNuclei(std::vector<int>& z, std::vector<int>& a) {
		if(z.size() != a.size() || z.size() != 2) {
			return false;
		}
	
		step1.SetNuclei(z[0], a[0], 0, 0, z[1], a[1]);
		SetSystemEquation();
		return true;
	}
	
	void DecaySystem::LinkTarget() {
		step1.SetLayeredTarget(&target);
	
		rxnLayer = target.FindLayerContaining(step1.GetTarget().GetZ(), step1.GetTarget().GetA());
		if(rxnLayer != -1) {
			step1.SetRxnLayer(rxnLayer);
			target_set_flag = true;
		} else {
			throw ReactionLayerException();
		}
	}
	
	void DecaySystem::SetSystemEquation() {
		m_sys_equation = step1.GetTarget().GetIsotopicSymbol();
		m_sys_equation += "-> ";
		m_sys_equation += step1.GetEjectile().GetIsotopicSymbol();
		m_sys_equation += "+";
		m_sys_equation += step1.GetResidual().GetIsotopicSymbol();
	}
	
	void DecaySystem::RunSystem() {
		if(!gen_set_flag) return;
		
		//Link up the target if it hasn't been done yet
		if(!target_set_flag) {
			LinkTarget();
		}
	
		double rxnTheta = std::acos(decay1dist.GetRandomCosTheta());
		double rxnPhi = (*m_phi1Range)(*generator);
		step1.SetPolarRxnAngle(rxnTheta);
		step1.SetAzimRxnAngle(rxnPhi);
	
		step1.TurnOnResidualEloss();
		step1.Calculate();

	}

}