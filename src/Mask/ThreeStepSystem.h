#ifndef THREESTEPSYSTEM_H
#define THREESTEPSYSTEM_H

#include "ReactionSystem.h"
#include "AngularDistribution.h"

namespace Mask {

	class ThreeStepSystem : public ReactionSystem {
	public:
		ThreeStepSystem();
		ThreeStepSystem(const std::vector<int>& z, const std::vector<int>& a);
		~ThreeStepSystem();
		bool SetNuclei(const std::vector<int>& z, const std::vector<int>& a) override;
		void RunSystem() override;
		std::vector<Nucleus>* GetNuclei() override;
	
		inline void SetDecay1Distribution(const std::string& filename) { m_step2Distribution.ReadDistributionFile(filename); };
		inline void SetDecay2Distribution(const std::string& filename) { m_step3Distribution.ReadDistributionFile(filename); };
	
		void SetReactionThetaType(int type) { m_step1.SetEjectileThetaType(type); };
	
		int GetDecay1AngularMomentum() { return m_step2Distribution.GetL(); };
		int GetDecay2AngularMomentum(){ return m_step3Distribution.GetL(); };
		double GetDecay1BranchingRatio() { return m_step2Distribution.GetBranchingRatio(); };
		double GetDecay2BranchingRatio(){ return m_step3Distribution.GetBranchingRatio(); };
	
	protected:
		void LinkTarget() override;
		void SetSystemEquation() override;

		std::uniform_real_distribution<double> m_phi2Range;
	
		Reaction m_step1, m_step2, m_step3;
	
		AngularDistribution m_step2Distribution, m_step3Distribution;
	
	};

}

#endif