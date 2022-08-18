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
#include <random>

namespace Mask {

	class ReactionSystem
	{
	public:
		ReactionSystem();
		virtual ~ReactionSystem();
	
		virtual bool SetNuclei(const std::vector<int>& z, const std::vector<int>& a) = 0;
		virtual void RunSystem() = 0;
		virtual std::vector<Nucleus>* GetNuclei() = 0;
	
		void AddTargetLayer(const std::vector<int>& zt, const std::vector<int>& at, const std::vector<int>& stoich, double thickness);
	
		/*Set sampling parameters*/
		void SetBeamDistro(double mean, double sigma)
		{
			if(m_beamDist)
				delete m_beamDist;
			m_beamDist = new std::normal_distribution<double>(mean, sigma); 
		}
	
		void SetTheta1Range(double min, double max)
		{ 
			if(m_theta1Range)
				delete m_theta1Range;
			m_theta1Range = new std::uniform_real_distribution<double>(std::cos(min*s_deg2rad), std::cos(max*s_deg2rad)); 
		}
	
		void SetPhi1Range(double min, double max)
		{
			if(m_phi1Range)
				delete m_phi1Range;
			m_phi1Range = new std::uniform_real_distribution<double>(min*s_deg2rad, max*s_deg2rad); 
		}
	
		void SetExcitationDistro(double mean, double sigma)
		{
			if(m_exDist)
				delete m_exDist;
			m_exDist = new std::normal_distribution<double>(mean, sigma); 
		}
	
		const std::string& GetSystemEquation() const { return m_sysEquation; }
	
	protected:
		virtual void LinkTarget() = 0;
		virtual void SetSystemEquation() = 0;
		
		LayeredTarget m_target;
	
		//Sampling information
		std::normal_distribution<double> *m_beamDist, *m_exDist;
		std::uniform_real_distribution<double> *m_theta1Range, *m_phi1Range; 
	
		bool m_isTargetSet;
		std::size_t m_rxnLayer;
		std::string m_sysEquation;
		std::vector<Nucleus> m_nuclei;
		static constexpr double s_deg2rad = M_PI/180.0;
	};

}

#endif