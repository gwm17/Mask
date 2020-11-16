#ifndef THREESTEPSYSTEM_H
#define THREESTEPSYSTEM_H

#include "ReactionSystem.h"

class ThreeStepSystem : public ReactionSystem {
public:
	ThreeStepSystem();
	ThreeStepSystem(std::vector<int>& z, std::vector<int>& a);
	~ThreeStepSystem();
	bool SetNuclei(std::vector<int>& z, std::vector<int>& a);
	void RunSystem();

	inline const Nucleus& GetBreakup1() const { return step2.GetEjectile(); };
	inline const Nucleus& GetBreakup2() const { return step2.GetResidual(); };
	inline const Nucleus& GetBreakup3() const { return step3.GetEjectile(); };
	inline const Nucleus& GetBreakup4() const { return step3.GetResidual(); };

protected:
	void LinkTarget();
	void SetSystemEquation();

	Reaction step2, step3;

};

#endif