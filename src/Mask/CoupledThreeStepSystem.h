#ifndef COUPLEDTHREESTEPSYSTEM_H
#define COUPLEDTHREESTEPSYSTEM_H

#include "ReactionSystem.h"
#include "AngularDistribution.h"

namespace Mask {

	struct CoupledThreeStepParameters
	{
		double beamEnergy = 0.;
		double rxnTheta = 0.;
		double rxnPhi = 0.;
		double cosdecay1Theta = 0.;
		double cosRelativeAngle = 0.;
		double decay1Theta = 0.;
		double decay1Phi = 0.;
		double decay2Theta = 0.;
		double decay2Phi = 0.;
		double residEx = 0.;
		double decay1Ex = 0.;
		double decay2Ex = 0.;
		double rxnDepth = 0.;
	};

	class CoupledThreeStepSystem : public ReactionSystem
	{
	public:
		CoupledThreeStepSystem(const std::vector<StepParameters>& params);
		~CoupledThreeStepSystem();

		virtual void SetLayeredTarget(const LayeredTarget& target) override;
		virtual void RunSystem() override;
		
	protected:
		void Init(const std::vector<StepParameters>& params);
		void SetSystemEquation() override;
		CoupledThreeStepParameters SampleParameters();
        void SampleCoupling(CoupledThreeStepParameters& params);
	
		Reaction m_step1, m_step2, m_step3;
	};

}

#endif