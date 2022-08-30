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
		m_beamDist(nullptr), m_theta1Range(nullptr), m_phi1Range(nullptr), m_exDist(nullptr), m_isTargetSet(false),
		m_rxnLayer(0), m_sysEquation("")
	{
	}
	
	ReactionSystem::~ReactionSystem()
	{
		delete m_beamDist;
		delete m_theta1Range;
		delete m_phi1Range;
		delete m_exDist;
	}
	
	void ReactionSystem::AddTargetLayer(const std::vector<int>& zt, const std::vector<int>& at, const std::vector<int>& stoich, double thickness)
	{
		m_target.AddLayer(zt, at, stoich, thickness);
	}

	ReactionSystem* CreateSystem(const std::vector<int>& z, const std::vector<int>& a)
	{
		if(z.size() != a.size())
		{
			std::cerr<<"Size of Z list does not equal size of A list!"<<std::endl;
			return nullptr;
		}

		switch(a.size())
		{
			case RxnSize::DecaySize: return new DecaySystem(z, a);
			case RxnSize::OneStepSize: return new OneStepSystem(z, a);
			case RxnSize::TwoStepSize: return new TwoStepSystem(z, a);
			case RxnSize::ThreeStepSize: return new ThreeStepSystem(z, a);
		}

		return nullptr;
	}

}