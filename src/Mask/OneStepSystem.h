#ifndef ONESTEPSYSTEM_H
#define ONESTEPSYSTEM_H

#include "ReactionSystem.h"

namespace Mask {

	class OneStepSystem: public ReactionSystem
	{
	public:
		OneStepSystem(const std::vector<StepParameters>& params);
		~OneStepSystem();
	
		virtual void SetLayeredTarget(const LayeredTarget& target) override;
		void RunSystem() override;
	
	private:
		void Init(const std::vector<StepParameters>& params);
		void SetSystemEquation() override;
	
		Reaction m_step1;
	};

}

#endif