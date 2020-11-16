#ifndef TWOSTEPSYSTEM_H
#define TWOSTEPSYSTEM_H

#include "ReactionSystem.h"

class TwoStepSystem : public ReactionSystem {
public:
	TwoStepSystem();
	TwoStepSystem(std::vector<int>& z, std::vector<int>& a);
	~TwoStepSystem();
	bool SetNuclei(std::vector<int>& z, std::vector<int>& a);
	void RunSystem();

	inline const Nucleus& GetBreakup1() const { return step2.GetEjectile(); };
	inline const Nucleus& GetBreakup2() const { return step2.GetResidual(); };

private:
	void LinkTarget();
	void SetSystemEquation();
	
	Reaction step2;

};

#endif