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
		name(""), m_fractional_depth_dist(0.0, 1.0)
	{
	}
	
	LayeredTarget::~LayeredTarget() {}
	
	/*Add in a Target made of a compound defined by a set of Zs, As, Ss, and a thickness*/
	void LayeredTarget::AddLayer(std::vector<int>& Z, std::vector<int>& A, std::vector<int>& stoich, double thickness) {
		Target t(thickness);
		t.SetElements(Z, A, stoich);
		layers.push_back(t);
	}
	
	/*
	  Here projectile refers to the incoming reactant particle (i.e. the beam)
	  Calculates energy loss assuming that the reaction occurs in the middle of the target layer
	  Note that the layer order can matter!
	*/
	double LayeredTarget::GetProjectileEnergyLoss(int zp, int ap, double startEnergy, int rxnLayer, double angle) {
	
		if(rxnLayer < 0 || ((unsigned int) rxnLayer) > layers.size()) {
			std::cerr<<"Reaction layer in eloss calculation is not in range! Returning 0"<<std::endl;
			return 0.0;
		}
	
		double eloss = 0.0;
		double newEnergy = startEnergy;
		double frac;
		for(int i=0; i<=rxnLayer; i++) {
			if(i == rxnLayer) {
				frac = m_fractional_depth_dist(RandomGenerator::GetInstance().GetGenerator());
				eloss += layers[i].GetEnergyLossFractionalDepth(zp, ap, newEnergy, angle, frac);
				newEnergy = startEnergy - eloss;
			} else {
				eloss += layers[i].GetEnergyLossTotal(zp, ap, newEnergy, angle);
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
	double LayeredTarget::GetEjectileEnergyLoss(int ze, int ae, double startEnergy, int rxnLayer, double angle) {
	
		if(rxnLayer < 0 || ((unsigned int) rxnLayer) > layers.size()) {
			std::cerr<<"Reaction layer in eloss calculation is not in range! Returning 0"<<std::endl;
			return 0.0;
		}
	
		double eloss = 0.0;
		if(angle < M_PI/2.0) {
			double newEnergy = startEnergy;
			for(unsigned int i=rxnLayer; i<layers.size(); i++) {
				if(i == ((unsigned int)rxnLayer)) {
					eloss += layers[i].GetEnergyLossFractionalDepth(ze, ae, newEnergy, angle, m_fractional_depth_dist(RandomGenerator::GetInstance().GetGenerator()));
					newEnergy = startEnergy - eloss;
				} else {
					eloss += layers[i].GetEnergyLossTotal(ze, ae, newEnergy, angle);
					newEnergy = startEnergy - eloss;
				}
			}
		} else { //Travelling backwards through target
			double newEnergy = startEnergy;
			for(int i=rxnLayer; i>=0; i--) {
				if(i == ((unsigned int)rxnLayer)) {
					eloss += layers[i].GetEnergyLossFractionalDepth(ze, ae, newEnergy, angle, m_fractional_depth_dist(RandomGenerator::GetInstance().GetGenerator()));
					newEnergy = startEnergy - eloss;
				} else {
					eloss += layers[i].GetEnergyLossTotal(ze, ae, newEnergy, angle);
					newEnergy = startEnergy - eloss;
				}
			}
		}
	
		return eloss;
	}
	
	/*ReverseEnergyLoss version of GetEjectileEnergyLoss*/
	double LayeredTarget::GetEjectileReverseEnergyLoss(int ze, int ae, double startEnergy, int rxnLayer, double angle) {
	
		if(rxnLayer < 0 || ((unsigned int) rxnLayer) > layers.size()) {
			std::cerr<<"Reaction layer in eloss calculation is not in range! Returning 0"<<std::endl;
			return 0.0;
		}	
	
		double eloss = 0.0;
		if(angle < M_PI/2.0) {
			double newEnergy = startEnergy;
			for(int i=(layers.size()-1); i>=rxnLayer; i--) {
				if(i == rxnLayer) {
					eloss += layers[i].GetReverseEnergyLossFractionalDepth(ze, ae, newEnergy, angle, m_fractional_depth_dist(RandomGenerator::GetInstance().GetGenerator()));
					newEnergy = startEnergy + eloss;
				} else {
					eloss += layers[i].GetReverseEnergyLossTotal(ze, ae, newEnergy, angle);
					newEnergy = startEnergy + eloss;
				}
			}
		} else {
			double newEnergy = startEnergy;
			for(int i=0; i <= rxnLayer; i++) {
				if(i == rxnLayer) {
					eloss += layers[i].GetReverseEnergyLossFractionalDepth(ze, ae, newEnergy, angle, m_fractional_depth_dist(RandomGenerator::GetInstance().GetGenerator()));
					newEnergy = startEnergy + eloss;
				} else {
					eloss += layers[i].GetReverseEnergyLossTotal(ze, ae, newEnergy, angle);
					newEnergy = startEnergy + eloss;
				}
			}
		}
	
		return eloss;
	}
	
	int LayeredTarget::FindLayerContaining(int Z, int A) {
		for(unsigned int i=0; i<layers.size(); i++)
			if(layers[i].ContainsElement(Z, A)) 
				return i;

		return -1;
	}
    
    
}