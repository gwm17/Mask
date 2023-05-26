#include "TwoStepSystem.h"
#include "RandomGenerator.h"
#include "KinematicsExceptions.h"

#include <sstream>

namespace Mask {
	
	TwoStepSystem::TwoStepSystem(const std::vector<StepParameters>& params) :
		ReactionSystem()
	{
		m_nuclei.resize(6);
		Init(params);
	}
	
	TwoStepSystem::~TwoStepSystem() {}
	
	void TwoStepSystem::Init(const std::vector<StepParameters>& params)
	{
		if(params.size() != 2 || params[0].rxnType != RxnType::Reaction || params[1].rxnType != RxnType::Decay ||
		   params[0].Z.size() != 3 || params[0].A.size() != 3 || params[1].Z.size() != 2 || params[1].A.size() != 2)
		{
			m_isValid = false;
			std::cerr << "Invalid parameters at TwoStepSystem::Init(), does not match TwoStep signature!" << std::endl;
			return;
		}

		const StepParameters& step1Params = params[0];
		const StepParameters& step2Params = params[1];

		//Setup nuclei
		int zr = step1Params.Z[0] + step1Params.Z[1] - step1Params.Z[2];
		int ar = step1Params.A[0] + step1Params.A[1] - step1Params.A[2];
		if(zr != step2Params.Z[0] || ar != step2Params.A[0])
		{
			m_isValid = false;
			std::cerr << "Invalid parameters at TwoStepSystem::Init(), step one and step two are not sequential! Step one recoil (Z,A): ("
					  << zr << "," << ar << ") Step two target (Z,A): (" << step2Params.Z[0] << "," << step2Params.A[0] << ")" <<std::endl;
			return;  
		}
		int zb = step2Params.Z[0] - step2Params.Z[1];
		int ab = step2Params.A[0] - step2Params.A[1];

		m_nuclei[0] = CreateNucleus(step1Params.Z[0], step1Params.A[0]); //target
		m_nuclei[1] = CreateNucleus(step1Params.Z[1], step1Params.A[1]); //projectile
		m_nuclei[2] = CreateNucleus(step1Params.Z[2], step1Params.A[2]); //ejectile
		m_nuclei[3] = CreateNucleus(zr, ar); //residual
		m_nuclei[4] = CreateNucleus(step2Params.Z[1], step2Params.A[1]); //breakup1
		m_nuclei[5] = CreateNucleus(zb, ab); //breakup2

		m_step1.BindNuclei(&(m_nuclei[0]), &(m_nuclei[1]), &(m_nuclei[2]), &(m_nuclei[3]));
		m_step2.BindNuclei(&(m_nuclei[3]), nullptr, &(m_nuclei[4]), &(m_nuclei[5]));
		SetSystemEquation();

		//Step one sampling parameters
		AddBeamDistribution(step1Params.meanBeamEnergy, step1Params.sigmaBeamEnergy);
		m_step1.SetEjectileThetaType(step1Params.thetaType);
		AddThetaRange(step1Params.thetaMin, step1Params.thetaMax);
		AddPhiRange(step1Params.phiMin, step1Params.phiMax);
		AddExcitationDistribution(step1Params.meanResidualEx, step1Params.sigmaResidualEx);

		//Step two sampling parameters
		AddPhiRange(step2Params.phiMin, step2Params.phiMax);
		AddDecayAngularDistribution(step2Params.angularDistFile);
		AddExcitationDistribution(step2Params.meanResidualEx, step2Params.sigmaResidualEx);
	}
	
	void TwoStepSystem::SetLayeredTarget(const LayeredTarget& target)
	{
		m_target = target;
		m_rxnLayer = m_target.FindLayerContaining(m_nuclei[0].Z, m_nuclei[0].A);
		if(m_rxnLayer != m_target.GetNumberOfLayers())
		{
			m_step1.SetLayeredTarget(&m_target);
			m_step2.SetLayeredTarget(&m_target);
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

	TwoStepParameters TwoStepSystem::SampleParameters()
	{
		TwoStepParameters params;
		std::mt19937& gen = RandomGenerator::GetInstance().GetGenerator();
		params.beamEnergy = (m_beamDistributions[0])(gen);
		params.rxnTheta = std::acos((m_thetaRanges[0])(gen));
		params.rxnPhi = (m_phiRanges[0])(gen);
		params.cosdecay1Theta = m_decayAngularDistributions[0].GetRandomCosTheta();
		params.decay1Theta = std::acos(params.cosdecay1Theta);
		params.decay1Phi = m_phiRanges[1](gen);
		params.residEx = (m_exDistributions[0])(gen);
		params.decay2Ex = m_exDistributions[1](gen);
		params.rxnDepth = (m_rxnDepthDist(gen));
		return params;
	}
	
	void TwoStepSystem::RunSystem()
	{
		//Sample parameters
		// std::mt19937& gen = RandomGenerator::GetInstance().GetGenerator();
		// double bke = (m_beamDistributions[0])(gen);
		// double rxnTheta = std::acos((m_thetaRanges[0])(gen));
		// double rxnPhi = (m_phiRanges[0])(gen);
		// double decay1costheta = m_decayAngularDistributions[0].GetRandomCosTheta();
		// double decay1Theta = std::acos(decay1costheta);
		// double decay1Phi = m_phiRanges[1](gen);
		// double residEx = (m_exDistributions[0])(gen);
		// double decay2Ex = m_exDistributions[1](gen);
		// double rxnDepth = (m_rxnDepthDist(gen));

		TwoStepParameters params;
		do
		{
			params = SampleParameters();
		}
		while(!(m_step1.CheckReactionThreshold(params.beamEnergy, params.residEx) 
				&& m_step2.CheckDecayThreshold(params.residEx, params.decay2Ex)));
		
		m_step1.SetReactionDepth(params.rxnDepth);
		m_step1.SetBeamKE(params.beamEnergy);
		m_step1.SetPolarRxnAngle(params.rxnTheta);
		m_step1.SetAzimRxnAngle(params.rxnPhi);
		m_step1.SetExcitation(params.residEx);
	
		m_step2.SetReactionDepth(params.rxnDepth);
		m_step2.SetPolarRxnAngle(params.decay1Theta);
		m_step2.SetAzimRxnAngle(params.decay1Phi);
		m_step2.SetExcitation(params.decay2Ex);
		
		m_step1.Calculate();
	
		if(params.cosdecay1Theta == -10)
		{
			m_nuclei[4].vec4.SetPxPyPzE(0., 0., 0., m_nuclei[4].groundStateMass);
			m_nuclei[5].vec4.SetPxPyPzE(0., 0., 0., m_nuclei[5].groundStateMass);
			return;
		}
		m_step2.SetResidualEnergyLoss(true);
		m_step2.Calculate();
	
	
	}

}
