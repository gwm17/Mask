#ifndef DECAYSYSTEM_H
#define DECAYSYSTEM_H

#include "ReactionSystem.h"
#include "AngularDistribution.h"

namespace Mask {

	class DecaySystem: public ReactionSystem {
	public:
		DecaySystem();
		DecaySystem(const std::vector<int>& z, const std::vector<int>& a);
		~DecaySystem();
	
		bool SetNuclei(const std::vector<int>& z, const std::vector<int>& a) override;
		void RunSystem() override;
		std::vector<Nucleus>*GetNuclei() override;
	
		virtual void SetDecay1Distribution(const std::string& filename) override { m_step1Distribution.ReadDistributionFile(filename); }
	
		int GetDecay1AngularMomentum() { return m_step1Distribution.GetL(); }
		double GetDecay1BranchingRatio() { return m_step1Distribution.GetBranchingRatio(); }
	
	private:
		void LinkTarget() override;
		void SetSystemEquation() override;
	
		Reaction m_step1;
	
		AngularDistribution m_step1Distribution;
		
	};

}

#endif