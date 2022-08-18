#include "RandomGenerator.h"

namespace Mask {
	RandomGenerator* RandomGenerator::s_instance = new RandomGenerator();

	RandomGenerator::RandomGenerator()
	{
		std::random_device rd;
		rng.seed(rd());
	}

	RandomGenerator::~RandomGenerator() {}
}