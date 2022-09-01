#ifndef THREESTEPSYSTEM_H
#define THREESTEPSYSTEM_H

#include "ReactionSystem.h"
#include "AngularDistribution.h"

namespace Mask {

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
	
		Reaction m_step1, m_step2, m_step3;
	
	};

}

#endif