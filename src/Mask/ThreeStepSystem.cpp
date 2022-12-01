#include "ThreeStepSystem.h"
#include "RandomGenerator.h"
#include "KinematicsExceptions.h"

namespace Mask {
	
	ThreeStepSystem::ThreeStepSystem(const std::vector<StepParameters>& params) :
		ReactionSystem()
	{
		m_nuclei.resize(8);
		Init(params);
	}
	
	ThreeStepSystem::~ThreeStepSystem() {}
	
	void ThreeStepSystem::Init(const std::vector<StepParameters>& params)
	{
		if(params.size() != 3 || params[0].rxnType != RxnType::Reaction || params[1].rxnType != RxnType::Decay ||
		  params[2].rxnType != RxnType::Decay || params[0].Z.size() != 3 || params[0].A.size() != 3 || params[1].Z.size() != 2 ||
		  params[1].A.size() != 2 || params[2].Z.size() != 2 || params[2].A.size() != 2)
		{
			m_isValid = false;
			std::cerr << "Invalid parameters at ThreeStepSystem::Init(), does not match ThreeStep signature!" << std::endl;
			return;
		}

		const StepParameters& step1Params = params[0];
		const StepParameters& step2Params = params[1];
		const StepParameters& step3Params = params[2];

		//Setup nuclei
		int zr = step1Params.Z[0] + step1Params.Z[1] - step1Params.Z[2];
		int ar = step1Params.A[0] + step1Params.A[1] - step1Params.A[2];
		if(zr != step2Params.Z[0] || ar != step2Params.A[0])
		{
			m_isValid = false;
			std::cerr << "Invalid parameters at ThreeStepSystem::Init(), step one and step two are not sequential! Step one recoil (Z,A): ("
					  << zr << "," << ar << ") Step two target (Z,A): (" << step2Params.Z[0] << "," << step2Params.A[0] << ")" <<std::endl;
			return;  
		}
		int zb2 = step2Params.Z[0] - step2Params.Z[1];
		int ab2 = step2Params.A[0] - step2Params.A[1];
		if(zb2 != step3Params.Z[0] || ab2 != step3Params.A[0])
		{
			m_isValid = false;
			std::cerr << "Invalid parameters at ThreeStepSystem::Init(), step two and step three are not sequential! Step two heavy (Z,A): ("
					  << zb2 << "," << ab2 << ") Step three target (Z,A): (" << step3Params.Z[0] << "," << step3Params.A[0] << ")" <<std::endl;
			return;
		}
		int zb4 = step3Params.Z[0] - step3Params.Z[1];
		int ab4 = step3Params.A[0] - step3Params.A[1];

		m_nuclei[0] = CreateNucleus(step1Params.Z[0], step1Params.A[0]); //target
		m_nuclei[1] = CreateNucleus(step1Params.Z[1], step1Params.A[1]); //projectile
		m_nuclei[2] = CreateNucleus(step1Params.Z[2], step1Params.A[2]); //ejectile
		m_nuclei[3] = CreateNucleus(zr, ar); //residual
		m_nuclei[4] = CreateNucleus(step2Params.Z[1], step2Params.A[1]); //breakup1
		m_nuclei[5] = CreateNucleus(zb2, ab2); //breakup2
		m_nuclei[5] = CreateNucleus(step3Params.Z[1], step3Params.A[1]); //breakup3
		m_nuclei[5] = CreateNucleus(zb4, ab4); //breakup4

		m_step1.BindNuclei(&(m_nuclei[0]), &(m_nuclei[1]), &(m_nuclei[2]), &(m_nuclei[3]));
		m_step2.BindNuclei(&(m_nuclei[3]), nullptr, &(m_nuclei[4]), &(m_nuclei[5]));
		m_step3.BindNuclei(&(m_nuclei[5]), nullptr, &(m_nuclei[6]), &(m_nuclei[7]));
	
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

		//Step three sampling parameters
		AddPhiRange(step3Params.phiMin, step3Params.phiMax);
		AddDecayAngularDistribution(step3Params.angularDistFile);
		AddExcitationDistribution(step3Params.meanResidualEx, step3Params.sigmaResidualEx);
	}

	void ThreeStepSystem::SetLayeredTarget(const LayeredTarget& target)
	{
		m_target = target;
		m_rxnLayer = m_target.FindLayerContaining(m_nuclei[0].Z, m_nuclei[0].A);
		if(m_rxnLayer != m_target.GetNumberOfLayers())
		{
			m_step1.SetLayeredTarget(&m_target);
			m_step2.SetLayeredTarget(&m_target);
			m_step3.SetLayeredTarget(&m_target);
			m_step1.SetRxnLayer(m_rxnLayer);
			m_step2.SetRxnLayer(m_rxnLayer);
			m_step3.SetRxnLayer(m_rxnLayer);
			m_isTargetSet = true;
		}
		else
			throw ReactionLayerException();
	}
	
	void ThreeStepSystem::SetSystemEquation()
	{
		std::stringstream stream;
		stream << m_nuclei[0].isotopicSymbol << "("
			   << m_nuclei[1].isotopicSymbol << ", "
			   << m_nuclei[2].isotopicSymbol << ")"
			   << m_nuclei[3].isotopicSymbol << "->"
			   << m_nuclei[4].isotopicSymbol << "+"
			   << m_nuclei[5].isotopicSymbol << "->"
			   << m_nuclei[6].isotopicSymbol << "+"
			   << m_nuclei[7].isotopicSymbol;
		m_sysEquation = stream.str();
	}
	
	void ThreeStepSystem::RunSystem()
	{
		//Sample parameters
		std::mt19937& gen = RandomGenerator::GetInstance().GetGenerator();
		double bke = m_beamDistributions[0](gen);
		double rxnTheta = std::acos((m_thetaRanges[0])(gen));
		double rxnPhi = m_phiRanges[0](gen);
		double decay1costheta = m_decayAngularDistributions[0].GetRandomCosTheta();
		double decay1Theta = std::acos(decay1costheta);
		double decay1Phi = m_phiRanges[1](gen);
		double decay2costheta = m_decayAngularDistributions[1].GetRandomCosTheta();
		double decay2Theta = std::acos(decay2costheta);
		double decay2Phi = m_phiRanges[2](gen);
		double residEx = m_exDistributions[0](gen);
		double decay1Ex = m_exDistributions[1](gen);
		double decay2Ex = m_exDistributions[2](gen);
		double rxnDepth = (m_rxnDepthDist(gen));
		
		m_step1.SetReactionDepth(rxnDepth);
		m_step1.SetBeamKE(bke);
		m_step1.SetPolarRxnAngle(rxnTheta);
		m_step1.SetAzimRxnAngle(rxnPhi);
		m_step1.SetExcitation(residEx);
	
		m_step2.SetReactionDepth(rxnDepth);
		m_step2.SetPolarRxnAngle(decay1Theta);
		m_step2.SetAzimRxnAngle(decay1Phi);
		m_step2.SetExcitation(decay1Ex);
	
		m_step3.SetReactionDepth(rxnDepth);
		m_step3.SetPolarRxnAngle(decay2Theta);
		m_step3.SetAzimRxnAngle(decay2Phi);
		m_step3.SetExcitation(decay2Ex);
	
		m_step1.Calculate();
	
		if(decay1costheta == -10)
		{
			m_nuclei[4].vec4.SetPxPyPzE(0., 0., 0., m_nuclei[4].groundStateMass);
			m_nuclei[5].vec4.SetPxPyPzE(0., 0., 0., m_nuclei[5].groundStateMass);
			m_nuclei[6].vec4.SetPxPyPzE(0., 0., 0., m_nuclei[6].groundStateMass);
			m_nuclei[7].vec4.SetPxPyPzE(0., 0., 0., m_nuclei[7].groundStateMass);
			return;
		}
		m_step2.Calculate();
	
		if(decay2costheta == -10)
		{
			m_nuclei[6].vec4.SetPxPyPzE(0., 0., 0., m_nuclei[6].groundStateMass);
			m_nuclei[7].vec4.SetPxPyPzE(0., 0., 0., m_nuclei[7].groundStateMass);
			return;
		}
		m_step3.SetResidualEnergyLoss(true);
		m_step3.Calculate();
	
	}

}