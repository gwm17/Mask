#ifndef THREESTEPSYSTEM_H
#define THREESTEPSYSTEM_H

#include "ReactionSystem.h"
#include "AngularDistribution.h"

namespace Mask {

	struct ThreeStepParameters
	{
		double beamEnergy = 0.;
		double rxnTheta = 0.;
		double rxnPhi = 0.;
		double cosdecay1Theta = 0.;
		double cosdecay2Theta = 0.;
		double decay1Theta = 0.;
		double decay1Phi = 0.;
		double decay2Theta = 0.;
		double decay2Phi = 0.;
		double residEx = 0.;
		double decay1Ex = 0.;
		double decay2Ex = 0.;
		double rxnDepth = 0.;
	};

	class ThreeStepSystem : public ReactionSystem
	{
	public:
		ThreeStepSystem(const std::vector<StepParameters>& params);
		~ThreeStepSystem();

		virtual void SetLayeredTarget(const LayeredTarget& target) override;
		virtual void RunSystem() override;
		
	protected:
		void Init(const std::vector<StepParameters>& params);
		void SetSystemEquation() override;
		ThreeStepParameters SampleParameters();
	
		Reaction m_step1, m_step2, m_step3;
	
	};

}

#endif