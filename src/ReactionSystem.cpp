#include "ReactionSystem.h"
#include "KinematicsExceptions.h"
#include "LegendrePoly.h"

namespace Mask {

	ReactionSystem::ReactionSystem() :
		m_beamDist(nullptr), m_theta1Range(nullptr), m_phi1Range(nullptr), m_exDist(nullptr), generator(nullptr), target_set_flag(false),
		gen_set_flag(false), rxnLayer(0), m_sys_equation("")
	{
	}
	
	ReactionSystem::~ReactionSystem() {
		delete m_beamDist;
		delete m_theta1Range;
		delete m_phi1Range;
		delete m_exDist;
	}
	
	void ReactionSystem::SetRandomGenerator(std::mt19937* gen) {
		generator = gen;
		gen_set_flag = true;
	}
	
	void ReactionSystem::AddTargetLayer(std::vector<int>& zt, std::vector<int>& at, std::vector<int>& stoich, double thickness) {
		target.AddLayer(zt, at, stoich, thickness);
	}

}