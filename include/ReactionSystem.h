/*
	ReactionSystem.h
	ReactionSystem is the base class from which all kinematics calculations should inherit. It contains all members and functions
	to perform a single step reaction/decay, which is the basic building block of all subsequent types of reactions in MASK. More
	complicated systems (see TwoStepSystem and ThreeStepSystem) add further data members and override the virtual functions.

	--GWM Jan. 2021
*/
#ifndef REACTIONSYSTEM_H
#define REACTIONSYSTEM_H

#include "Reaction.h"
#include "AngularDistribution.h"
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

	/*Set sampling parameters*/
	void SetRandomGenerator(TRandom3* gen);
	inline void SetBeamDistro(double mean, double sigma) { m_beamDist = std::make_pair(mean, sigma); };
	inline void SetReactionThetaType(int type) { step1.SetEjectileThetaType(type); };
	inline void SetTheta1Range(double min, double max) { m_theta1Range = std::make_pair(min*deg2rad, max*deg2rad); };
	inline void SetPhi1Range(double min, double max) { m_phi1Range = std::make_pair(min*deg2rad, max*deg2rad); };
	inline void SetExcitationDistro(double mean, double sigma) { m_exDist = std::make_pair(mean, sigma); };
	inline void SetDecay1Distribution(const std::string& file) { decay1dist.ReadDistributionFile(file); };
	inline void SetDecay2Distribution (const std::string& file) { decay2dist.ReadDistributionFile(file); };

	virtual void RunSystem();

	inline const Nucleus& GetTarget() const { return step1.GetTarget(); };
	inline const Nucleus& GetProjectile() const { return step1.GetProjectile(); };
	inline const Nucleus& GetEjectile() const { return step1.GetEjectile(); };
	inline const Nucleus& GetResidual() const { return step1.GetResidual(); };
	inline const char* GetSystemEquation() const { return m_sys_equation.c_str(); };
	inline int GetDecay1AngularMomentum() { return decay1dist.GetL(); };
	inline int GetDecay2AngularMomentum(){ return decay2dist.GetL(); };
	inline double GetDecay1BranchingRatio() { return decay1dist.GetBranchingRatio(); };
	inline double GetDecay2BranchingRatio(){ return decay2dist.GetBranchingRatio(); };

protected:
	virtual void LinkTarget();
	virtual void SetSystemEquation();
	
	Reaction step1;
	LayeredTarget target;

	//Sampling information
	std::pair<double, double> m_beamDist, m_theta1Range, m_phi1Range, m_exDist;
	TRandom3* generator; //not owned by ReactionSystem

	AngularDistribution decay1dist, decay2dist;

	bool target_set_flag, gen_set_flag;
	int rxnLayer;
	std::string m_sys_equation;
	static constexpr double deg2rad = M_PI/180.0;
};

};

#endif