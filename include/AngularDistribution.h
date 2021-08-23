#ifndef ANGULARDISTRIBUTION_H
#define ANGULARDISTRIBUTION_H

#include <TRandom3.h>
#include <string>
#include <vector>

class AngularDistribution {
public:
	AngularDistribution();
	AngularDistribution(const std::string& file);
	~AngularDistribution();
	void ReadDistributionFile(const std::string& file);
	void AttachRandomNumberGenerator(TRandom3* random) { generator = random; };
	double GetRandomCosTheta();
	int GetL() { return L; };
	double GetBranchingRatio() { return branchingRatio; };

private:
	bool IsIsotropic();
	bool IsGeneratorSet() {
		if(generator) {
			return true;
		} else {
			return false;
		}
	}

	TRandom3* generator; //NOT OWNED BY ANGULAR DISTRIBUTION

	double branchingRatio;
	int L;
	std::vector<double> constants;
	bool isoFlag;
};

#endif