#include "AngularDistribution.h"
#include <fstream>
#include <cmath>
#include <iostream>
#include "LegendrePoly.h"

namespace Mask {

	AngularDistribution::AngularDistribution() :
		generator(nullptr), uniform_cosine_dist(-1.0, 1.0), uniform_prob_dist(0.0, 1.0), branchingRatio(1.0), L(0), isoFlag(true)
	{
	}
	
	AngularDistribution::AngularDistribution(const std::string& file) :
		generator(nullptr), branchingRatio(1.0), L(0), isoFlag(true)
	{
		ReadDistributionFile(file);
	}
	
	AngularDistribution::~AngularDistribution() {}
	
	void AngularDistribution::ReadDistributionFile(const std::string& file) {
	
		if(file == "none" || file == "") {
			L=0;
			branchingRatio=1.0;
			constants.clear();
			constants.push_back(0.5);
			isoFlag = true;
			return;
		}
	
		std::ifstream input(file);
		std::string junk;
		int l;
		double par;
	
		if(!input.is_open()) {
			std::cerr<<"Unable to open distribution file. All values reset to default."<<std::endl;
			L=0;
			branchingRatio=1.0;
			constants.clear();
			constants.push_back(0.5);
			isoFlag = true;
			return;
		}
	
		input>>junk>>l;
		while(input>>junk) {
			input>>par;
			constants.push_back(par);
		}
		input.close();
	
		if(constants.size() != ((unsigned int) l+1)) {
			std::cerr<<"Unexpected number of constants for given angular momentum! Expected "<<l+1<<" and given "<<constants.size()<<std::endl;
			std::cerr<<"Setting all values to default."<<std::endl;
			branchingRatio=1.0;
			constants.clear();
			constants.push_back(0.5);
			isoFlag = true;
			return;
		}
	
		//Total branching ratio
		branchingRatio = constants[0]*2.0;
		L = l;
	
		//Renormalize distribution such that total prob is 1.0.
		//Test branching ratio to see if we "make" a decay particle,
		//then use re-normalized distribution to pick an angle. 
		if(constants[0] < 0.5) {
			double norm = 0.5/constants[0];
			for(auto& value : constants)
				value *= norm;
		}
	
		isoFlag = false;
	
	}
	
	double AngularDistribution::GetRandomCosTheta() {
	
		if(!IsGeneratorSet()) {
			std::cerr<<"Random number generator is not set in AngularDistribution! Returning default value of 0"<<std::endl;
			return 0.0;
		}
	
		if(isoFlag) return uniform_cosine_dist(*generator);
	
		double test, probability;
		double costheta;
	
		test = uniform_prob_dist(*generator);
		if(test > branchingRatio) return -10;
	
		do {
			probability = 0.0;
			costheta = uniform_cosine_dist(*generator);
			test = uniform_prob_dist(*generator);
			for(unsigned int i=0; i<constants.size(); i++)
				probability += constants[i]*P_l(i*2, costheta);
		} while(test > probability);
	
		return costheta;
	}

}