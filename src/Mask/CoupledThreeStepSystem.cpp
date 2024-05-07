#include "CoupledThreeStepSystem.h"

#include "Math/Boost.h"
#include "Math/Vector3D.h"
#include "Math/RotationY.h"
#include "Math/RotationZ.h"

namespace Mask {

    CoupledThreeStepSystem::CoupledThreeStepSystem(const std::vector<StepParameters>& params) :
		ReactionSystem()
	{
		m_nuclei.resize(8);
		Init(params);
	}
	
	CoupledThreeStepSystem::~CoupledThreeStepSystem() {}
	
	void CoupledThreeStepSystem::Init(const std::vector<StepParameters>& params)
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
		m_nuclei[6] = CreateNucleus(step3Params.Z[1], step3Params.A[1]); //breakup3
		m_nuclei[7] = CreateNucleus(zb4, ab4); //breakup4

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

	void CoupledThreeStepSystem::SetLayeredTarget(const LayeredTarget& target)
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
	
	void CoupledThreeStepSystem::SetSystemEquation()
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

	CoupledThreeStepParameters CoupledThreeStepSystem::SampleParameters()
	{
		CoupledThreeStepParameters params;
		std::mt19937& gen = RandomGenerator::GetInstance().GetGenerator();
		params.beamEnergy = m_beamDistributions[0](gen);
		params.rxnTheta = std::acos((m_thetaRanges[0])(gen));
		params.rxnPhi = m_phiRanges[0](gen);
		params.cosdecay1Theta = m_decayAngularDistributions[0].GetRandomCosTheta();
        params.cosRelativeAngle = m_decayAngularDistributions[1].GetRandomCosTheta();
		params.decay1Theta = std::acos(params.cosdecay1Theta);
		params.decay1Phi = m_phiRanges[1](gen);
		params.residEx = m_exDistributions[0](gen);
		params.decay1Ex = m_exDistributions[1](gen);
		params.decay2Ex = m_exDistributions[2](gen);
		params.rxnDepth = (m_rxnDepthDist(gen));
		return params;
	}

    //Called after running step2
    void CoupledThreeStepSystem::SampleCoupling(CoupledThreeStepParameters& params)
    {
        ROOT::Math::Boost boost(m_nuclei[5].vec4.BoostToCM());
		ROOT::Math::PxPyPzEVector a1Vec = (boost * m_nuclei[4].vec4);
		ROOT::Math::XYZVector a1PVec(a1Vec.Px(), a1Vec.Py(), a1Vec.Pz());
        ROOT::Math::PxPyPzEVector a2Vec;
		ROOT::Math::XYZVector a2PVec;
		ROOT::Math::RotationZ zRot;
		ROOT::Math::RotationY yRot;
		a1PVec *= 1.0/a1Vec.P();
        zRot.SetAngle(a1PVec.Phi());
		yRot.SetAngle(a1PVec.Theta());
		double relAngle = std::acos(params.cosRelativeAngle);
		ROOT::Math::XYZVector a1PVecAligned = yRot.Inverse() * (zRot.Inverse() * a1PVec);
		a2PVec.SetXYZ(
			std::sin(relAngle),
			0.0,
			std::cos(relAngle)
		);
		a2PVec = zRot * (yRot * a2PVec);
		params.decay2Theta = a2PVec.Theta();
		params.decay2Phi = a2PVec.Phi();
    }
	
	void CoupledThreeStepSystem::RunSystem()
	{
		CoupledThreeStepParameters params;
		do
		{
			params = SampleParameters();
		}
		while(!(m_step1.CheckReactionThreshold(params.beamEnergy, params.residEx) 
				&& m_step2.CheckDecayThreshold(params.residEx, params.decay1Ex)
				&& m_step3.CheckDecayThreshold(params.decay1Ex, params.decay2Ex)));
		
		m_step1.SetReactionDepth(params.rxnDepth);
		m_step1.SetBeamKE(params.beamEnergy);
		m_step1.SetPolarRxnAngle(params.rxnTheta);
		m_step1.SetAzimRxnAngle(params.rxnPhi);
		m_step1.SetExcitation(params.residEx);
	
		m_step2.SetReactionDepth(params.rxnDepth);
		m_step2.SetPolarRxnAngle(params.decay1Theta);
		m_step2.SetAzimRxnAngle(params.decay1Phi);
		m_step2.SetExcitation(params.decay1Ex);
	
		m_step3.SetReactionDepth(params.rxnDepth);
		
		m_step3.SetExcitation(params.decay2Ex);
	
		m_step1.Calculate();
	
		if(params.cosdecay1Theta == -10)
		{
			m_nuclei[4].vec4.SetPxPyPzE(0., 0., 0., m_nuclei[4].groundStateMass);
			m_nuclei[5].vec4.SetPxPyPzE(0., 0., 0., m_nuclei[5].groundStateMass);
			m_nuclei[6].vec4.SetPxPyPzE(0., 0., 0., m_nuclei[6].groundStateMass);
			m_nuclei[7].vec4.SetPxPyPzE(0., 0., 0., m_nuclei[7].groundStateMass);
			return;
		}
		m_step2.Calculate();
	
		m_step3.SetResidualEnergyLoss(true);

		SampleCoupling(params); //must be called here as frame needs to be set

		m_step3.SetPolarRxnAngle(params.decay2Theta);
		m_step3.SetAzimRxnAngle(params.decay2Phi);
		m_step3.Calculate();
	}
}