#include "ReactionSystem.h"
#include "RxnType.h"
#include "KinematicsExceptions.h"
#include "LegendrePoly.h"
#include "DecaySystem.h"
#include "OneStepSystem.h"
#include "TwoStepSystem.h"
#include "ThreeStepSystem.h"

namespace Mask {

	ReactionSystem::ReactionSystem() :
		m_isTargetSet(false), m_isValid(true), m_rxnLayer(0), m_sysEquation(""), m_rxnDepthDist(0.0, 1.0)
	{
	}
	
	ReactionSystem::~ReactionSystem()
	{
	}

	ReactionSystem* CreateSystem(const std::vector<StepParameters>& params)
	{
		switch(params.size())
		{
			case 1:
			{
				if(params[0].rxnType == RxnType::Decay)
					return new DecaySystem(params);
				else if (params[0].rxnType == RxnType::Reaction)
					return new OneStepSystem(params);
			}
			case 2: return new TwoStepSystem(params);
			case 3: return new ThreeStepSystem(params);
		}

		return nullptr;
	}

	/*Set sampling parameters*/
	//Should only ever really be one of these, but they cannot be set a priori
	void ReactionSystem::AddBeamDistribution(double mean, double sigma)
	{
		if(mean == -1.0 || sigma == -1.0)
		{
			m_isValid = false;
			std::cerr << "Invalid parameters at ReactionSystem::AddBeamDistribution(), distribution is invalid -> mean: " << mean
					  << " sigma: " << sigma << std::endl;
			return;
		}
		m_beamDistributions.emplace_back(mean, sigma);
	}

	//Again only really one of these, but same issue as beam. This is the sampling of the primary reaction/decay (step 1)
	void ReactionSystem::AddThetaRange(double min, double max)
	{
		if(min == -1.0 || max == -1.0 || min > max)
		{
			m_isValid = false;
			std::cerr << "Invalid parameters at ReactionSystem::AddThetaRange(), range is invalid -> min: " << min
					  << " max: " << max << std::endl;
			return;
		}
		m_thetaRanges.emplace_back(std::cos(min*s_deg2rad), std::cos(max*s_deg2rad));
	}

	void ReactionSystem::AddPhiRange(double min, double max)
	{
		if(min == -1.0 || max == -1.0 || min > max)
		{
			m_isValid = false;
			std::cerr << "Invalid parameters at ReactionSystem::AddPhiRange(), range is invalid -> min: " << min
					  << " max: " << max << std::endl;
			return;
		}
		m_phiRanges.emplace_back(min*s_deg2rad, max*s_deg2rad);
	}

	//Each reaction step can generate an excited nucleus (for a reaction can make an excited residual, decay can make an excited
	//breakup2 or "heavy")
	void ReactionSystem::AddExcitationDistribution(double mean, double sigma)
	{
		if(mean == -1.0 || sigma == -1.0)
		{
			m_isValid = false;
			std::cerr << "Invalid parameters at ReactionSystem::AddExcitationDistribution(), distribution is invalid -> mean: " << mean
					  << " sigma: " << sigma << std::endl;
			return;
		}
		m_exDistributions.emplace_back(mean, sigma);
	}

	//Apply angular distributions to decay products. Here apply them to the breakup1 or "light" fragment
	void ReactionSystem::AddDecayAngularDistribution(const std::string& filename)
	{
		if(filename.empty())
			m_decayAngularDistributions.emplace_back();
		else
			m_decayAngularDistributions.emplace_back(filename);
	}
}