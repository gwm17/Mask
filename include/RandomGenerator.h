#ifndef RANDOMGENERATOR_H
#define RANDOMGENERATOR_H

#include <random>

namespace Mask {

	class RandomGenerator {
	public:
		~RandomGenerator();
		inline std::mt19937& GetGenerator() { return rng; }

		inline static RandomGenerator& GetInstance() {
			static RandomGenerator s_instance;
			return s_instance;
		}

	private:
		RandomGenerator();
		std::mt19937 rng;
	};

}

#endif