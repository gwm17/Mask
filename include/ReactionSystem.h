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
#include "KinematicsExceptions.h"
#include <vector>
#include <TRandom3.h>

namespace Mask {

class ReactionSystem {
public:
	ReactionSystem();
	virtual ~ReactionSystem();

	virtual bool SetNuclei(std::vector<int>& z, std::vector<int>& a) = 0;
	virtual void RunSystem() = 0;

	void AddTargetLayer(std::vector<int>& zt, std::vector<int>& at, std::vector<int>& stoich, double thickness);

	/*Set sampling parameters*/
	virtual void SetRandomGenerator(TRandom3* gen);
	inline void SetBeamDistro(double mean, double sigma) { m_beamDist = std::make_pair(mean, sigma); };
	inline void SetTheta1Range(double min, double max) { m_theta1Range = std::make_pair(min*deg2rad, max*deg2rad); };
	inline void SetPhi1Range(double min, double max) { m_phi1Range = std::make_pair(min*deg2rad, max*deg2rad); };
	inline void SetExcitationDistro(double mean, double sigma) { m_exDist = std::make_pair(mean, sigma); };

	inline const std::string& GetSystemEquation() const { return m_sys_equation; };

protected:
	virtual void LinkTarget() = 0;
	virtual void SetSystemEquation() = 0;
	
	LayeredTarget target;

	//Sampling information
	std::pair<double, double> m_beamDist, m_theta1Range, m_phi1Range, m_exDist;
	TRandom3* generator; //not owned by ReactionSystem

	bool target_set_flag, gen_set_flag;
	int rxnLayer;
	std::string m_sys_equation;
	static constexpr double deg2rad = M_PI/180.0;
};

};

#endif