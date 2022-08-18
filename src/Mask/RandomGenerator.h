#ifndef RANDOMGENERATOR_H
#define RANDOMGENERATOR_H

#include <random>

namespace Mask {

	class RandomGenerator
	{
	public:
		~RandomGenerator();
		std::mt19937& GetGenerator() { return rng; }
		static RandomGenerator& GetInstance() { return *s_instance; }

	private:
		RandomGenerator();
		
		static RandomGenerator* s_instance;
		std::mt19937 rng;
	};

}

#endif