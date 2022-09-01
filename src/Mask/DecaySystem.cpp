#include "DecaySystem.h"
#include "RandomGenerator.h"

#include <sstream>

namespace Mask {

	DecaySystem::DecaySystem(const std::vector<StepParameters>& params) :
		ReactionSystem()
	{
		m_nuclei.resize(3);
		Init(params);
	}
	
	DecaySystem::~DecaySystem() {}
	
	void DecaySystem::Init(const std::vector<StepParameters>& params)
	{
		if(params.size() != 1 || params[0].rxnType != RxnType::Decay ||
		   params[0].Z.size() != 2 || params[0].A.size() != 2)
		{
			m_isValid = false;
			std::cerr << "Invalid parameters at DecaySystem::Init(), does not match Decay signature!" << std::endl;
			return;
		}

		const StepParameters& step1Params = params[0];

		int zr = step1Params.Z[0] - step1Params.Z[1];
		int ar = step1Params.A[0] - step1Params.A[1];

		m_nuclei[0] = CreateNucleus(step1Params.Z[0], step1Params.A[0]); //target
		m_nuclei[1] = CreateNucleus(step1Params.Z[1], step1Params.A[1]); //breakup1
		m_nuclei[2] = CreateNucleus(zr, ar); //breakup2

		m_step1.BindNuclei(&(m_nuclei[0]), nullptr, &(m_nuclei[1]), &(m_nuclei[2]));
		SetSystemEquation();

		
		AddPhiRange(step1Params.phiMin, step1Params.phiMax);
		AddDecayAngularDistribution(step1Params.angularDistFile);
		AddExcitationDistribution(step1Params.meanResidualEx, step1Params.sigmaResidualEx);

		return;
	}
	
	void DecaySystem::SetLayeredTarget(const LayeredTarget& target)
	{
		m_target = target;
		m_rxnLayer = m_target.FindLayerContaining(m_nuclei[0].Z, m_nuclei[0].A);
		if(m_rxnLayer != m_target.GetNumberOfLayers())
		{
			m_step1.SetLayeredTarget(&m_target);
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
		std::mt19937& gen = RandomGenerator::GetInstance().GetGenerator();
		double rxnTheta = m_decayAngularDistributions[0].GetRandomCosTheta();
		double rxnPhi = m_phiRanges[0](gen);
		double ex = m_exDistributions[0](gen);

		m_step1.SetPolarRxnAngle(rxnTheta);
		m_step1.SetAzimRxnAngle(rxnPhi);
		m_step1.SetExcitation(ex);
		m_step1.SetResidualEnergyLoss(true);
		m_step1.Calculate();
	}

}