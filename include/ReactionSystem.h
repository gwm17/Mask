#ifndef REACTIONSYSTEM_H
#define REACTIONSYSTEM_H

#include "Reaction.h"
#include <vector>
#include <TRandom3.h>

namespace Mask {

class ReactionSystem {
public:
	ReactionSystem();
	ReactionSystem(std::vector<int>& z, std::vector<int>& a);
	virtual ~ReactionSystem();
	virtual bool SetNuclei(std::vector<int>& z, std::vector<int>& a);
	void AddTargetLayer(std::vector<int>& zt, std::vector<int>& at, std::vector<int>& stoich, double thickness);
	inline void SetRandomGenerator(TRandom3* gen) { generator = gen; gen_set_flag = true; };
	inline void SetBeamDistro(double mean, double sigma) { m_beamDist = std::make_pair(mean, sigma); };
	inline void SetTheta1Range(double min, double max) { m_theta1Range = std::make_pair(min*deg2rad, max*deg2rad); };
	inline void SetExcitationDistro(double mean, double sigma) { m_exDist = std::make_pair(mean, sigma); };
	virtual void RunSystem();

	inline const Nucleus& GetTarget() const { return step1.GetTarget(); };
	inline const Nucleus& GetProjectile() const { return step1.GetProjectile(); };
	inline const Nucleus& GetEjectile() const { return step1.GetEjectile(); };
	inline const Nucleus& GetResidual() const { return step1.GetResidual(); };
	inline const char* GetSystemEquation() const { return m_sys_equation.c_str(); };

protected:
	virtual void LinkTarget();
	virtual void SetSystemEquation();
	
	Reaction step1;
	LayeredTarget target;
	std::pair<double, double> m_beamDist, m_theta1Range, m_exDist;
	TRandom3* generator; //not owned by ReactionSystem

	bool target_set_flag, gen_set_flag;
	int rxnLayer;
	std::string m_sys_equation;
	static constexpr double deg2rad = M_PI/180.0;
};

};

#endif