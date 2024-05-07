#ifndef ANGULARDISTRIBUTION_H
#define ANGULARDISTRIBUTION_H

#include <string>
#include <vector>
#include <random>

namespace Mask {

	class AngularDistribution
	{
	public:
		AngularDistribution();
		AngularDistribution(const std::string& file);
		~AngularDistribution();
		void ReadDistributionFile(const std::string& file);
		double GetRandomCosTheta();
		int GetL() { return m_L; }
		double GetBranchingRatio() { return m_branchingRatio; }
		double GetProbability(double cosTheta);
	
	private:
		bool IsIsotropic() { return m_isIsotropic; }
		
		std::uniform_real_distribution<double> m_uniformCosineDist;
		std::uniform_real_distribution<double> m_uniformProbDist;
	
		double m_branchingRatio;
		int m_L;
		std::vector<double> m_constants;
		bool m_isIsotropic;
	};

}

#endif