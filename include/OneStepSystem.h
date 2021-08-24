#ifndef ONESTEPSYSTEM_H
#define ONESTEPSYSTEM_H

#include "ReactionSystem.h"

namespace Mask {

class OneStepSystem: public ReactionSystem {
public:
	OneStepSystem();
	OneStepSystem(std::vector<int>& z, std::vector<int>& a);
	~OneStepSystem();

	bool SetNuclei(std::vector<int>& z, std::vector<int>& a) override;
	void RunSystem() override;

	inline void SetReactionThetaType(int type) { step1.SetEjectileThetaType(type); };
	inline const Nucleus& GetTarget() const { return step1.GetTarget(); };
	inline const Nucleus& GetProjectile() const { return step1.GetProjectile(); };
	inline const Nucleus& GetEjectile() const { return step1.GetEjectile(); };
	inline const Nucleus& GetResidual() const { return step1.GetResidual(); };

private:
	void LinkTarget() override;
	void SetSystemEquation() override;

	Reaction step1;

};

};

#endif