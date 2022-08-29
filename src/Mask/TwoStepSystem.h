#ifndef TWOSTEPSYSTEM_H
#define TWOSTEPSYSTEM_H

#include "ReactionSystem.h"
#include "AngularDistribution.h"

namespace Mask {

	class TwoStepSystem : public ReactionSystem
	{
	public:
		TwoStepSystem();
		TwoStepSystem(const std::vector<int>& z, const std::vector<int>& a);
		~TwoStepSystem();
		bool SetNuclei(const std::vector<int>& z, const std::vector<int>& a) override;
		void RunSystem() override;
		std::vector<Nucleus>* GetNuclei() override;
	
		virtual void SetDecay1Distribution(const std::string& filename) override { m_step2Distribution.ReadDistributionFile(filename); };
	
		virtual void SetReactionThetaType(RxnThetaType type) override { m_step1.SetEjectileThetaType(type); };

		int GetDecay1AngularMomentum() { return m_step2Distribution.GetL(); };
		double GetDecay1BranchingRatio() { return m_step2Distribution.GetBranchingRatio(); };
	
	private:
		void LinkTarget() override;
		void SetSystemEquation() override;

		std::uniform_real_distribution<double> m_phi2Range;
		
		Reaction m_step1, m_step2;
	
		AngularDistribution m_step2Distribution;
	
	};

}

#endif