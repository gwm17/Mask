#include "DecaySystem.h"
#include "RandomGenerator.h"

#include <sstream>

namespace Mask {

	DecaySystem::DecaySystem() :
		ReactionSystem()
	{
		m_nuclei.resize(3);
	}
	
	DecaySystem::DecaySystem(const std::vector<int>& z, const std::vector<int>& a) :
		ReactionSystem()
	{
		m_nuclei.resize(3);
		SetNuclei(z, a);
	}
	
	DecaySystem::~DecaySystem() {}
	
	bool DecaySystem::SetNuclei(const std::vector<int>& z, const std::vector<int>& a)
	{
		if(z.size() != a.size() || z.size() != 2)
			return false;

		int zr = z[0] - z[1];
		int ar = a[0] - a[1];

		m_nuclei[0] = CreateNucleus(z[0], a[0]); //target
		m_nuclei[1] = CreateNucleus(z[1], a[1]); //breakup1
		m_nuclei[2] = CreateNucleus(zr, ar); //breakup2

		m_step1.BindNuclei(&(m_nuclei[0]), nullptr, &(m_nuclei[1]), &(m_nuclei[2]));
		SetSystemEquation();
		return true;
	}

	std::vector<Nucleus>* DecaySystem::GetNuclei()
	{
		return &m_nuclei;
	}
	
	void DecaySystem::LinkTarget() {
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
	
	void DecaySystem::SetSystemEquation()
	{
		std::stringstream stream;
		stream << m_nuclei[0].isotopicSymbol << "->"
			   << m_nuclei[1].isotopicSymbol << "+"
			   << m_nuclei[2].isotopicSymbol;
		m_sysEquation = stream.str();
	}
	
	void DecaySystem::RunSystem()
	{
		//Link up the target if it hasn't been done yet
		if(!m_isTargetSet)
			LinkTarget();
	
		double rxnTheta = std::acos(m_step1Distribution.GetRandomCosTheta());
		double rxnPhi = (*m_phi1Range)(RandomGenerator::GetInstance().GetGenerator());
		m_step1.SetPolarRxnAngle(rxnTheta);
		m_step1.SetAzimRxnAngle(rxnPhi);
		m_step1.SetResidualEnergyLoss(true);
		m_step1.Calculate();
	}

}