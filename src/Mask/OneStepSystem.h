#ifndef ONESTEPSYSTEM_H
#define ONESTEPSYSTEM_H

#include "ReactionSystem.h"

namespace Mask {

	class OneStepSystem: public ReactionSystem {
	public:
		OneStepSystem();
		OneStepSystem(const std::vector<int>& z, const std::vector<int>& a);
		~OneStepSystem();
	
		bool SetNuclei(const std::vector<int>& z, const std::vector<int>& a) override;
		void RunSystem() override;
		std::vector<Nucleus>* GetNuclei() override;
	
		virtual void SetReactionThetaType(RxnThetaType type) override { m_step1.SetEjectileThetaType(type); };
	
	private:
		void LinkTarget() override;
		void SetSystemEquation() override;
	
		Reaction m_step1;
	
	};

}

#endif