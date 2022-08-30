/*

LayeredTarget.h
Functional unit for targets in the Mask environment. Contains a
set (read: vector) of Targets for use in reaction calculations. In this
way handles situations such as carbon backed targets

Based on code by D.W. Visser written at Yale for the original SPANC

Written by G.W. McCann Aug. 2020

*/
#include "LayeredTarget.h"
#include <iostream>

namespace Mask {

	LayeredTarget::LayeredTarget() :
		m_name(""), m_fractionalDepthDistribution(0.0, 1.0)
	{
	}
	
	LayeredTarget::~LayeredTarget() {}
	
	/*Add in a Target made of a compound defined by a set of Zs, As, Ss, and a thickness*/
	void LayeredTarget::AddLayer(const std::vector<int>& Z, const std::vector<int>& A, const std::vector<int>& stoich, double thickness)
	{
		m_layers.emplace_back(Z, A, stoich, thickness);
	}
	
	/*
	  Here projectile refers to the incoming reactant particle (i.e. the beam)
	  Calculates energy loss assuming that the reaction occurs in the middle of the target layer
	  Note that the layer order can matter!
	*/
	double LayeredTarget::GetProjectileEnergyLoss(int zp, int ap, double startEnergy, std::size_t rxnLayer, double angle)
	{
		if(rxnLayer > m_layers.size())
		{
			std::cerr<<"Reaction layer in eloss calculation is not in range! Returning 0"<<std::endl;
			return 0.0;
		}
	
		double eloss = 0.0;
		double newEnergy = startEnergy;
		double frac;
		for(std::size_t i=0; i<=rxnLayer; i++)
		{
			if(i == rxnLayer)
			{
				frac = m_fractionalDepthDistribution(RandomGenerator::GetInstance().GetGenerator());
				eloss += m_layers[i].GetEnergyLossFractionalDepth(zp, ap, newEnergy, angle, frac);
				newEnergy = startEnergy - eloss;
			}
			else
			{
				eloss += m_layers[i].GetEnergyLossTotal(zp, ap, newEnergy, angle);
				newEnergy = startEnergy-eloss;
			}
		}
	
		return eloss;
	}
	
	/*
	  Here ejectile refers to the outgoing reactant particle
	  Calculates energy loss assuming that the reaction occurs in the middle of the target
	  Note that the layer order can matter!
	*/
	double LayeredTarget::GetEjectileEnergyLoss(int ze, int ae, double startEnergy, std::size_t rxnLayer, double angle) {
	
		if(rxnLayer > m_layers.size())
		{
			std::cerr<<"Reaction layer in eloss calculation is not in range! Returning 0"<<std::endl;
			return 0.0;
		}
	
		double eloss = 0.0;
		RandomGenerator& gen = RandomGenerator::GetInstance();
		if(angle < M_PI/2.0)
		{
			double newEnergy = startEnergy;
			eloss += m_layers[rxnLayer].GetEnergyLossFractionalDepth(ze, ae, newEnergy, angle, m_fractionalDepthDistribution(gen.GetGenerator()));
			newEnergy = startEnergy - eloss;
			if(rxnLayer == m_layers.size())
				return eloss;

			for(std::size_t i=rxnLayer; i<m_layers.size(); i++)
			{
				eloss += m_layers[i].GetEnergyLossTotal(ze, ae, newEnergy, angle);
				newEnergy = startEnergy - eloss;
			}
		}
		else
		{ //Travelling backwards through target
			double newEnergy = startEnergy;
			eloss += m_layers[rxnLayer].GetEnergyLossFractionalDepth(ze, ae, newEnergy, angle, m_fractionalDepthDistribution(gen.GetGenerator()));
			newEnergy = startEnergy - eloss;
			if(rxnLayer == 0)
				return eloss;

			for(std::size_t i=rxnLayer-1; i>0; i--) //unsigned ints cant be less than 0
			{
				eloss += m_layers[i].GetEnergyLossTotal(ze, ae, newEnergy, angle);
				newEnergy = startEnergy - eloss;
			}
			eloss += m_layers[0].GetEnergyLossTotal(ze, ae, newEnergy, angle);
			newEnergy = startEnergy - eloss;
		}
	
		return eloss;
	}
	
	/*ReverseEnergyLoss version of GetEjectileEnergyLoss*/
	double LayeredTarget::GetEjectileReverseEnergyLoss(int ze, int ae, double startEnergy, std::size_t rxnLayer, double angle)
	{
		if(rxnLayer > m_layers.size())
		{
			std::cerr<<"Reaction layer in eloss calculation is not in range! Returning 0"<<std::endl;
			return 0.0;
		}	
	
		double eloss = 0.0;
		RandomGenerator& gen = RandomGenerator::GetInstance();
		if(angle < M_PI/2.0)
		{
			double newEnergy = startEnergy;
			for(std::size_t i=(m_layers.size()-1); i>rxnLayer; i--)
			{
				eloss += m_layers[i].GetReverseEnergyLossTotal(ze, ae, newEnergy, angle);
				newEnergy = startEnergy + eloss;
			}
			eloss += m_layers[rxnLayer].GetReverseEnergyLossFractionalDepth(ze, ae, newEnergy, angle, m_fractionalDepthDistribution(gen.GetGenerator()));
			newEnergy = startEnergy + eloss;
		}
		else
		{
			double newEnergy = startEnergy;
			for(std::size_t i=0; i < rxnLayer; i++)
			{
				eloss += m_layers[i].GetReverseEnergyLossTotal(ze, ae, newEnergy, angle);
				newEnergy = startEnergy + eloss;
			}
			eloss += m_layers[rxnLayer].GetReverseEnergyLossFractionalDepth(ze, ae, newEnergy, angle, m_fractionalDepthDistribution(gen.GetGenerator()));
			newEnergy = startEnergy + eloss;
		}
	
		return eloss;
	}
	
	std::size_t LayeredTarget::FindLayerContaining(int Z, int A)
	{
		for(std::size_t i=0; i<m_layers.size(); i++)
		{
			if(m_layers[i].ContainsElement(Z, A)) 
				return i;
		}

		return m_layers.size();
	}
    
    
}