#include "RandomGenerator.h"

namespace Mask {
	
	RandomGenerator::RandomGenerator()
	{
		std::random_device rd;
		rng.seed(rd());
	}

	RandomGenerator::~RandomGenerator() {}
}