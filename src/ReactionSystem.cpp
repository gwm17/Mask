#include "ReactionSystem.h"
#include "KinematicsExceptions.h"
#include "LegendrePoly.h"

namespace Mask {

ReactionSystem::ReactionSystem() :
	m_beamDist(0,0), m_theta1Range(0,0), m_phi1Range(0,0), m_exDist(0,0), generator(nullptr), target_set_flag(false), gen_set_flag(false), rxnLayer(0), m_sys_equation("")
{
}

ReactionSystem::~ReactionSystem() {}

void ReactionSystem::SetRandomGenerator(TRandom3* gen) {
	generator = gen;
	gen_set_flag = true;
}

void ReactionSystem::AddTargetLayer(std::vector<int>& zt, std::vector<int>& at, std::vector<int>& stoich, double thickness) {
	target.AddLayer(zt, at, stoich, thickness);
}

};