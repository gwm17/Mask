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
	inline void SetRandomGenerator(TRandom3* gen) { generator = gen; gen_set_flag = true; };
	inline void SetBeamDistro(double mean, double sigma) { m_beamDist = std::make_pair(mean, sigma); };
	inline void SetTheta1Range(double min, double max) { m_theta1Range = std::make_pair(min*deg2rad, max*deg2rad); };
	inline void SetExcitationDistro(double mean, double sigma) { m_exDist = std::make_pair(mean, sigma); };
	inline void SetDecay1AngularMomentum(double l) { L1 = l; };
	inline void SetDecay2AngularMomentum(double l) { L2 = l; };

	/*Sampling over angular distribution*/
	double GetDecayTheta(int L);

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

	//Sampling information
	std::pair<double, double> m_beamDist, m_theta1Range, m_exDist;
	TRandom3* generator; //not owned by ReactionSystem

	bool target_set_flag, gen_set_flag;
	int rxnLayer;
	int L1, L2;
	std::string m_sys_equation;
	static constexpr double deg2rad = M_PI/180.0;
};

};

#endif