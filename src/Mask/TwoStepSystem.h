#ifndef TWOSTEPSYSTEM_H
#define TWOSTEPSYSTEM_H

#include "ReactionSystem.h"
#include "AngularDistribution.h"

namespace Mask {

	struct TwoStepParameters
	{
		double beamEnergy = 0.;
		double rxnTheta = 0.;
		double rxnPhi = 0.;
		double cosdecay1Theta = 0.;
		double decay1Theta = 0.;
		double decay1Phi = 0.;
		double residEx = 0.;
		double decay2Ex = 0.;
		double rxnDepth = 0.;
	};

	class TwoStepSystem : public ReactionSystem
	{
	public:
		TwoStepSystem(const std::vector<StepParameters>& params);
		~TwoStepSystem();

		virtual void SetLayeredTarget(const LayeredTarget& target) override;
		virtual void RunSystem() override;
	
	private:
		void Init(const std::vector<StepParameters>& params);
		void SetSystemEquation() override;
		TwoStepParameters SampleParameters();

		Reaction m_step1, m_step2;
	};

}

#endif