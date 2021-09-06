#ifndef ANGULARDISTRIBUTION_H
#define ANGULARDISTRIBUTION_H

#include <string>
#include <vector>
#include <random>

namespace Mask {

	class AngularDistribution {
	public:
		AngularDistribution();
		AngularDistribution(const std::string& file);
		~AngularDistribution();
		void ReadDistributionFile(const std::string& file);
		inline void AttachRandomNumberGenerator(std::mt19937* random) { generator = random; }
		double GetRandomCosTheta();
		inline int GetL() { return L; }
		inline double GetBranchingRatio() { return branchingRatio; }
	
	private:
		inline bool IsIsotropic() { return isoFlag; }
		inline bool IsGeneratorSet() {
			if(generator) {
				return true;
			} else {
				return false;
			}
		}
	
		std::mt19937* generator; //NOT OWNED BY ANGULAR DISTRIBUTION
		std::uniform_real_distribution<double> uniform_cosine_dist;
		std::uniform_real_distribution<double> uniform_prob_dist;
	
		double branchingRatio;
		int L;
		std::vector<double> constants;
		bool isoFlag;
	};

}

#endif