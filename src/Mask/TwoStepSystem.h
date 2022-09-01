#ifndef TWOSTEPSYSTEM_H
#define TWOSTEPSYSTEM_H

#include "ReactionSystem.h"
#include "AngularDistribution.h"

namespace Mask {

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

		Reaction m_step1, m_step2;
	};

}

#endif