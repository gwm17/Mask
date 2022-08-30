#include "TwoStepSystem.h"
#include "RandomGenerator.h"
#include "KinematicsExceptions.h"

#include <sstream>

namespace Mask {

	TwoStepSystem::TwoStepSystem() :
		ReactionSystem(), m_phi2Range(0, 2.0*M_PI)
	{
		m_nuclei.resize(6);
	}
	
	TwoStepSystem::TwoStepSystem(const std::vector<int>& z, const std::vector<int>& a) :
		ReactionSystem(), m_phi2Range(0, 2.0*M_PI)
	{
		m_nuclei.resize(6);
		SetNuclei(z, a);
	}
	
	TwoStepSystem::~TwoStepSystem() {}
	
	bool TwoStepSystem::SetNuclei(const std::vector<int>&z, const std::vector<int>& a)
	{
		if(z.size() != a.size() || z.size() != 4)
			return false;

		int zr = z[0] + z[1] - z[2];
		int ar = a[0] + a[1] - a[2];
		int zb = zr - z[3];
		int ab = ar - a[3];

		m_nuclei[0] = CreateNucleus(z[0], a[0]); //target
		m_nuclei[1] = CreateNucleus(z[1], a[1]); //projectile
		m_nuclei[2] = CreateNucleus(z[2], a[2]); //ejectile
		m_nuclei[3] = CreateNucleus(zr, ar); //residual
		m_nuclei[4] = CreateNucleus(z[3], a[3]); //breakup1
		m_nuclei[5] = CreateNucleus(zb, ab); //breakup2
	
		m_step1.BindNuclei(&(m_nuclei[0]), &(m_nuclei[1]), &(m_nuclei[2]), &(m_nuclei[3]));
		m_step2.BindNuclei(&(m_nuclei[3]), nullptr, &(m_nuclei[4]), &(m_nuclei[5]));
		SetSystemEquation();
		return true;
	}

	std::vector<Nucleus>* TwoStepSystem::GetNuclei()
	{
		return &m_nuclei;
	}
	
	void TwoStepSystem::LinkTarget()
	{
		m_step1.SetLayeredTarget(&m_target);
		m_step2.SetLayeredTarget(&m_target);
	
		m_rxnLayer = m_target.FindLayerContaining(m_nuclei[0].Z, m_nuclei[0].A);
		if(m_rxnLayer != m_target.GetNumberOfLayers())
		{
			m_step1.SetRxnLayer(m_rxnLayer);
			m_step2.SetRxnLayer(m_rxnLayer);
			m_isTargetSet = true;
		}
		else
			throw ReactionLayerException();
	}
	
	void TwoStepSystem::SetSystemEquation()
	{
		std::stringstream stream;
		stream << m_nuclei[0].isotopicSymbol << "("
			   << m_nuclei[1].isotopicSymbol << ", "
			   << m_nuclei[2].isotopicSymbol << ")"
			   << m_nuclei[3].isotopicSymbol << "->"
			   << m_nuclei[4].isotopicSymbol << "+"
			   << m_nuclei[5].isotopicSymbol;
		m_sysEquation = stream.str();
	}
	
	void TwoStepSystem::RunSystem()
	{
		//Link up the target if it hasn't been done yet
		if(!m_isTargetSet)
			LinkTarget();
	
		//Sample parameters
		double bke = (*m_beamDist)(RandomGenerator::GetInstance().GetGenerator());
		double rxnTheta = acos((*m_theta1Range)(RandomGenerator::GetInstance().GetGenerator()));
		double rxnPhi = (*m_phi1Range)(RandomGenerator::GetInstance().GetGenerator());
		double decay1costheta = m_step2Distribution.GetRandomCosTheta();
		double decay1Theta = std::acos(decay1costheta);
		double decay1Phi = m_phi2Range(RandomGenerator::GetInstance().GetGenerator());
		double residEx = (*m_exDist)(RandomGenerator::GetInstance().GetGenerator());
	
		m_step1.SetBeamKE(bke);
		m_step1.SetPolarRxnAngle(rxnTheta);
		m_step1.SetAzimRxnAngle(rxnPhi);
		m_step1.SetExcitation(residEx);
	
		m_step2.SetPolarRxnAngle(decay1Theta);
		m_step2.SetAzimRxnAngle(decay1Phi);
		
		m_step1.Calculate();
	
		if(decay1costheta == -10)
		{
			m_nuclei[4].vec4.SetPxPyPzE(0., 0., 0., m_nuclei[4].groundStateMass);
			m_nuclei[5].vec4.SetPxPyPzE(0., 0., 0., m_nuclei[5].groundStateMass);
			return;
		}
		m_step2.SetResidualEnergyLoss(true);
		m_step2.Calculate();
	
	
	}

}
