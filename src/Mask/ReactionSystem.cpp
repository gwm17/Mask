#include "ReactionSystem.h"
#include "KinematicsExceptions.h"
#include "LegendrePoly.h"

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

}