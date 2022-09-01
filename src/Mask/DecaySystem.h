#ifndef DECAYSYSTEM_H
#define DECAYSYSTEM_H

#include "ReactionSystem.h"
#include "AngularDistribution.h"

namespace Mask {

	class DecaySystem: public ReactionSystem
	{
	public:
		DecaySystem(const std::vector<StepParameters>& params);
		~DecaySystem();
	
		virtual void SetLayeredTarget(const LayeredTarget& target) override;
		virtual void RunSystem() override;
	
	private:
		void Init(const std::vector<StepParameters>& params);
		void SetSystemEquation() override;
	
		Reaction m_step1;
	};

}

#endif