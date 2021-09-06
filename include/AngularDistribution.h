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
		double GetRandomCosTheta();
		inline int GetL() { return L; }
		inline double GetBranchingRatio() { return branchingRatio; }
	
	private:
		inline bool IsIsotropic() { return isoFlag; }
		
		std::uniform_real_distribution<double> uniform_cosine_dist;
		std::uniform_real_distribution<double> uniform_prob_dist;
	
		double branchingRatio;
		int L;
		std::vector<double> constants;
		bool isoFlag;
	};

}

#endif