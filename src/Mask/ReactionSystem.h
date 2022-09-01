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
#include "RxnType.h"
#include "AngularDistribution.h"
#include <vector>
#include <random>

namespace Mask {

	struct StepParameters
	{
		RxnType rxnType = RxnType::None;
		std::vector<int> Z;
		std::vector<int> A;
		double meanBeamEnergy = -1.0;
		double sigmaBeamEnergy = -1.0;
		RxnThetaType thetaType = RxnThetaType::None;
		double thetaMin = -1.0;
		double thetaMax = -1.0;
		double phiMin = -1.0;
		double phiMax = -1.0;
		double meanResidualEx = -1.0;
		double sigmaResidualEx = -1.0;
		std::string angularDistFile;
	};

	class ReactionSystem
	{
	public:
		ReactionSystem();
		virtual ~ReactionSystem();

		virtual void SetLayeredTarget(const LayeredTarget& target) = 0;
		virtual void RunSystem() = 0;

		std::vector<Nucleus>* GetNuclei() { return &m_nuclei; }
		const std::string& GetSystemEquation() const { return m_sysEquation; }
		bool IsValid() const { return m_isValid; }

	protected:
		virtual void SetSystemEquation() = 0;

		void AddBeamDistribution(double mean, double sigma);
		void AddThetaRange(double min, double max);
		void AddPhiRange(double min, double max);
		void AddExcitationDistribution(double mean, double sigma);
		void AddDecayAngularDistribution(const std::string& filename);
		
		LayeredTarget m_target;
	
		//Sampling information
		std::vector<std::normal_distribution<double>> m_beamDistributions, m_exDistributions;
		std::vector<std::uniform_real_distribution<double>> m_thetaRanges, m_phiRanges;
		std::vector<AngularDistribution> m_decayAngularDistributions;

		bool m_isTargetSet;
		bool m_isValid;

		std::size_t m_rxnLayer;
		std::string m_sysEquation;
		std::vector<Nucleus> m_nuclei;

		static constexpr double s_deg2rad = M_PI/180.0;
	};

	ReactionSystem* CreateSystem(const std::vector<StepParameters>& params);
}

#endif