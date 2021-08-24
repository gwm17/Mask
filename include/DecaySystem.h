#ifndef DECAYSYSTEM_H
#define DECAYSYSTEM_H

#include "ReactionSystem.h"
#include "AngularDistribution.h"

namespace Mask {

class DecaySystem: public ReactionSystem {
public:
	DecaySystem();
	DecaySystem(std::vector<int>& z, std::vector<int>& a);
	~DecaySystem();

	bool SetNuclei(std::vector<int>& z, std::vector<int>& a) override;
	void RunSystem() override;
	void SetRandomGenerator(TRandom3* gen) override;

	inline void SetDecay1Distribution(const std::string& filename) { decay1dist.ReadDistributionFile(filename); };

	inline const Nucleus& GetTarget() const { return step1.GetTarget(); };
	inline const Nucleus& GetEjectile() const { return step1.GetEjectile(); };
	inline const Nucleus& GetResidual() const { return step1.GetResidual(); };

	inline int GetDecay1AngularMomentum() { return decay1dist.GetL(); };
	inline double GetDecay1BranchingRatio() { return decay1dist.GetBranchingRatio(); };


private:
	void LinkTarget() override;
	void SetSystemEquation() override;

	Reaction step1;

	AngularDistribution decay1dist;
	
};

};

#endif