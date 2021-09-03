#ifndef THREESTEPSYSTEM_H
#define THREESTEPSYSTEM_H

#include "ReactionSystem.h"
#include "AngularDistribution.h"

namespace Mask {

	class ThreeStepSystem : public ReactionSystem {
	public:
		ThreeStepSystem();
		ThreeStepSystem(std::vector<int>& z, std::vector<int>& a);
		~ThreeStepSystem();
		bool SetNuclei(std::vector<int>& z, std::vector<int>& a) override;
		void RunSystem() override;
		void SetRandomGenerator(std::mt19937* gen) override;
	
		inline void SetDecay1Distribution(const std::string& filename) { decay1dist.ReadDistributionFile(filename); };
		inline void SetDecay2Distribution(const std::string& filename) { decay2dist.ReadDistributionFile(filename); };
	
		inline void SetReactionThetaType(int type) { step1.SetEjectileThetaType(type); };
		inline const Nucleus& GetTarget() const { return step1.GetTarget(); };
		inline const Nucleus& GetProjectile() const { return step1.GetProjectile(); };
		inline const Nucleus& GetEjectile() const { return step1.GetEjectile(); };
		inline const Nucleus& GetResidual() const { return step1.GetResidual(); };
		inline const Nucleus& GetBreakup1() const { return step2.GetEjectile(); };
		inline const Nucleus& GetBreakup2() const { return step2.GetResidual(); };
		inline const Nucleus& GetBreakup3() const { return step3.GetEjectile(); };
		inline const Nucleus& GetBreakup4() const { return step3.GetResidual(); };
	
		inline int GetDecay1AngularMomentum() { return decay1dist.GetL(); };
		inline int GetDecay2AngularMomentum(){ return decay2dist.GetL(); };
		inline double GetDecay1BranchingRatio() { return decay1dist.GetBranchingRatio(); };
		inline double GetDecay2BranchingRatio(){ return decay2dist.GetBranchingRatio(); };
	
	protected:
		void LinkTarget() override;
		void SetSystemEquation() override;

		std::uniform_real_distribution<double> m_phi2Range;
	
		Reaction step1, step2, step3;
	
		AngularDistribution decay1dist, decay2dist;
	
	};

}

#endif