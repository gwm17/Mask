#include "OneStepSystem.h"
#include "RandomGenerator.h"

#include <sstream>

namespace Mask {

	OneStepSystem::OneStepSystem() :
		ReactionSystem()
	{
		m_nuclei.resize(4);
	}
	
	OneStepSystem::OneStepSystem(const std::vector<int>& z, const std::vector<int>& a) :
		ReactionSystem()
	{
		m_nuclei.resize(4);
		SetNuclei(z, a);
	}
	
	OneStepSystem::~OneStepSystem() {}
	
	bool OneStepSystem::SetNuclei(const std::vector<int>& z, const std::vector<int>& a)
	{
		if(z.size() != a.size() || z.size() != 3)
			return false;

		int zr = z[0] + z[1] - z[2];
		int ar = a[0] + a[1] - a[2];

		m_nuclei[0] = CreateNucleus(z[0], a[0]); //target
		m_nuclei[1] = CreateNucleus(z[1], a[1]); //projectile
		m_nuclei[2] = CreateNucleus(z[2], a[2]); //ejectile
		m_nuclei[3] = CreateNucleus(zr, ar); //residual

		m_step1.BindNuclei(&(m_nuclei[0]), &(m_nuclei[1]), &(m_nuclei[2]), &(m_nuclei[3]));
		SetSystemEquation();
		return true;
	}

	std::vector<Nucleus>* OneStepSystem::GetNuclei()
	{
		return &m_nuclei;
	}
	
	void OneStepSystem::LinkTarget()
	{
		m_step1.SetLayeredTarget(&m_target);
	
		m_rxnLayer = m_target.FindLayerContaining(m_nuclei[0].Z, m_nuclei[0].A);
		if(m_rxnLayer != m_target.GetNumberOfLayers())
		{
			m_step1.SetRxnLayer(m_rxnLayer);
			m_isTargetSet = true;
		}
		else
			throw ReactionLayerException();
	}
	
	void OneStepSystem::SetSystemEquation()
	{
		std::stringstream stream;
		stream << m_nuclei[0].isotopicSymbol << "("
			   << m_nuclei[1].isotopicSymbol << ", "
			   << m_nuclei[2].isotopicSymbol << ")"
			   << m_nuclei[3].isotopicSymbol;
		m_sysEquation = stream.str();
	}
	
	void OneStepSystem::RunSystem() {
		//Link up the target if it hasn't been done yet
		if(!m_isTargetSet)
		{
			LinkTarget();
		}
	
		//Sample parameters
		double bke = (*m_beamDist)(RandomGenerator::GetInstance().GetGenerator());
		double rxnTheta = std::acos((*m_theta1Range)(RandomGenerator::GetInstance().GetGenerator()));
		double rxnPhi = (*m_phi1Range)(RandomGenerator::GetInstance().GetGenerator());
		double residEx = (*m_exDist)(RandomGenerator::GetInstance().GetGenerator());
		
		m_step1.SetBeamKE(bke);
		m_step1.SetPolarRxnAngle(rxnTheta);
		m_step1.SetAzimRxnAngle(rxnPhi);
		m_step1.SetExcitation(residEx);
		
		m_step1.SetResidualEnergyLoss(true);
		m_step1.Calculate();
	}

}