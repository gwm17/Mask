#ifndef TWOSTEPSYSTEM_H
#define TWOSTEPSYSTEM_H

#include "ReactionSystem.h"
#include "AngularDistribution.h"

namespace Mask {

	class TwoStepSystem : public ReactionSystem {
	public:
		TwoStepSystem();
		TwoStepSystem(std::vector<int>& z, std::vector<int>& a);
		~TwoStepSystem();
		bool SetNuclei(std::vector<int>& z, std::vector<int>& a) override;
		void RunSystem() override;
		const std::vector<Nucleus>& GetNuclei() override;
	
		inline void SetDecay1Distribution(const std::string& filename) { decay1dist.ReadDistributionFile(filename); };
	
		inline void SetReactionThetaType(int type) { step1.SetEjectileThetaType(type); };
		inline const Nucleus& GetTarget() const { return step1.GetTarget(); };
		inline const Nucleus& GetProjectile() const { return step1.GetProjectile(); };
		inline const Nucleus& GetEjectile() const { return step1.GetEjectile(); };
		inline const Nucleus& GetResidual() const { return step1.GetResidual(); };
		inline const Nucleus& GetBreakup1() const { return step2.GetEjectile(); };
		inline const Nucleus& GetBreakup2() const { return step2.GetResidual(); };
	
		inline int GetDecay1AngularMomentum() { return decay1dist.GetL(); };
		inline double GetDecay1BranchingRatio() { return decay1dist.GetBranchingRatio(); };
	
	private:
		void LinkTarget() override;
		void SetSystemEquation() override;

		std::uniform_real_distribution<double> m_phi2Range;
		
		Reaction step1, step2;
	
		AngularDistribution decay1dist;
	
	};

}

#endif