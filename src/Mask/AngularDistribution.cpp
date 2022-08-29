#include "AngularDistribution.h"
#include "RandomGenerator.h"
#include <fstream>
#include <cmath>
#include <iostream>
#include "LegendrePoly.h"

namespace Mask {

	AngularDistribution::AngularDistribution() :
		m_uniformCosineDist(-1.0, 1.0), m_uniformProbDist(0.0, 1.0), m_branchingRatio(1.0), m_L(0), m_isIsotropic(true)
	{
	}
	
	AngularDistribution::AngularDistribution(const std::string& file) :
		m_branchingRatio(1.0), m_L(0), m_isIsotropic(true)
	{
		ReadDistributionFile(file);
	}
	
	AngularDistribution::~AngularDistribution() {}
	
	void AngularDistribution::ReadDistributionFile(const std::string& file)
	{
	
		if(file == "none" || file == "")
		{
			m_L=0;
			m_branchingRatio=1.0;
			m_constants.clear();
			m_constants.push_back(0.5);
			m_isIsotropic = true;
			return;
		}
	
		std::ifstream input(file);
		std::string junk;
		int l;
		double par;
	
		if(!input.is_open())
		{
			std::cerr<<"Unable to open distribution file. All values reset to default."<<std::endl;
			m_L=0;
			m_branchingRatio=1.0;
			m_constants.clear();
			m_constants.push_back(0.5);
			m_isIsotropic = true;
			return;
		}
	
		input>>junk>>l;
		while(input>>junk)
		{
			input>>par;
			m_constants.push_back(par);
		}
		input.close();
	
		if(m_constants.size() != ((std::size_t) l+1))
		{
			std::cerr<<"Unexpected number of constants for given angular momentum! Expected "<<l+1<<" and given "<<m_constants.size()<<std::endl;
			std::cerr<<"Setting all values to default."<<std::endl;
			m_L=0;
			m_branchingRatio=1.0;
			m_constants.clear();
			m_constants.push_back(0.5);
			m_isIsotropic = true;
			return;
		}
	
		//Total branching ratio
		m_branchingRatio = m_constants[0]*2.0;
		m_L = l;

		std::cout<<"Angular distribution from "<<file<<" branching ratio: "<<m_branchingRatio<<std::endl;
		std::cout<<"Angular distribution from "<<file<<" L: "<<m_L<<std::endl;
	
		//Renormalize distribution such that total prob is 1.0.
		//Test branching ratio to see if we "make" a decay particle,
		//then use re-normalized distribution to pick an angle. 
		if(m_constants[0] < 0.5)
		{
			double norm = 0.5/m_constants[0];
			for(auto& value : m_constants)
				value *= norm;
		}
	
		m_isIsotropic = false;
	
	}
	
	double AngularDistribution::GetRandomCosTheta()
	{
		if(m_isIsotropic) 
			return m_uniformCosineDist(RandomGenerator::GetInstance().GetGenerator());
	
		double test, probability;
		double costheta;
	
		test = m_uniformProbDist(RandomGenerator::GetInstance().GetGenerator());
		if(test > m_branchingRatio)
			return -10;
	
		do
		{
			probability = 0.0;
			costheta = m_uniformCosineDist(RandomGenerator::GetInstance().GetGenerator());
			test = m_uniformProbDist(RandomGenerator::GetInstance().GetGenerator());
			for(std::size_t i=0; i<m_constants.size(); i++)
				probability += m_constants[i]*P_l(i*2, costheta);
		}
		while(test > probability);
	
		return costheta;
	}

}