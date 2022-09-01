#include "OneStepSystem.h"
#include "RandomGenerator.h"

#include <sstream>

namespace Mask {

	OneStepSystem::OneStepSystem(const std::vector<StepParameters>& params) :
		ReactionSystem()
	{
		m_nuclei.resize(4);
		Init(params);
	}
	
	OneStepSystem::~OneStepSystem() {}
	
	void OneStepSystem::Init(const std::vector<StepParameters>& params)
	{
		if(params.size() != 1 || params[0].rxnType != RxnType::Reaction ||
		   params[0].Z.size() != 3 || params[0].A.size() != 3)
		{
			m_isValid = false;
			std::cerr << "Invalid parameters at OneStepSystem::Init(), does not match OneStep signature!" << std::endl;
			return;
		}

		const StepParameters& step1Params = params[0];

		//Set nuclei

		int zr = step1Params.Z[0] + step1Params.Z[1] - step1Params.Z[2];
		int ar = step1Params.A[0] + step1Params.A[1] - step1Params.A[2];

		m_nuclei[0] = CreateNucleus(step1Params.Z[0], step1Params.A[0]); //target
		m_nuclei[1] = CreateNucleus(step1Params.Z[1], step1Params.A[1]); //projectile
		m_nuclei[2] = CreateNucleus(step1Params.Z[2], step1Params.A[2]); //ejectile
		m_nuclei[3] = CreateNucleus(zr, ar); //residual

		m_step1.BindNuclei(&(m_nuclei[0]), &(m_nuclei[1]), &(m_nuclei[2]), &(m_nuclei[3]));
		SetSystemEquation();

		//Set sampling parameters

		AddBeamDistribution(step1Params.meanBeamEnergy, step1Params.sigmaBeamEnergy);
		m_step1.SetEjectileThetaType(step1Params.thetaType);
		AddThetaRange(step1Params.thetaMin, step1Params.thetaMax);
		AddPhiRange(step1Params.phiMin, step1Params.phiMax);
		AddExcitationDistribution(step1Params.meanResidualEx, step1Params.sigmaResidualEx);
	}
	
	void OneStepSystem::SetLayeredTarget(const LayeredTarget& target)
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
	
	void OneStepSystem::SetSystemEquation()
	{
		std::stringstream stream;
		stream << m_nuclei[0].isotopicSymbol << "("
			   << m_nuclei[1].isotopicSymbol << ", "
			   << m_nuclei[2].isotopicSymbol << ")"
			   << m_nuclei[3].isotopicSymbol;
		m_sysEquation = stream.str();
	}
	
	void OneStepSystem::RunSystem()
	{
		//Sample parameters
		std::mt19937& gen = RandomGenerator::GetInstance().GetGenerator();
		double bke = (m_beamDistributions[0])(gen);
		double rxnTheta = std::acos((m_thetaRanges[0])(gen));
		double rxnPhi = (m_phiRanges[0])(gen);
		double residEx = (m_exDistributions[0])(gen);
		
		m_step1.SetBeamKE(bke);
		m_step1.SetPolarRxnAngle(rxnTheta);
		m_step1.SetAzimRxnAngle(rxnPhi);
		m_step1.SetExcitation(residEx);
		
		m_step1.SetResidualEnergyLoss(true);
		m_step1.Calculate();
	}

}